// sim_compartment.cpp

#include "sim_compartment.h"
#include "parameters.h"
#include "reporter.h"
#include "randomizer.h"
#include "user_defined.h"

//
// MODEL DYNAMICS
//

Population::Population(Parameters& P, unsigned int pindex)
 : seed_row(0), p(pindex)
{
    // Set up built-in compartments
    N  = P.pop[p].size;
    S  = N;
    E  = vector<Compartment>(S.size());
    Ea = vector<Compartment>(S.size());
    Ip = vector<Compartment>(S.size());
    Ia = vector<Compartment>(S.size());
    Is = vector<Compartment>(S.size());
    C  = vector<Compartment>(S.size());
    R  = vector<double>(S.size(), 0.);
    V  = vector<double>(S.size(), 0.);
    V2 = vector<double>(S.size(), 0.);

    // Initial immunity
    for (unsigned int a = 0; a < S.size(); ++a) {
        double imm = 0;
        if (P.deterministic) imm = S[a] * P.pop[p].imm0[a];
        else                 imm = (unsigned int)(S[a] * P.pop[p].imm0[a] + 0.5);
        S[a] -= imm;
        R[a] += imm;
    }

    // Set up user-specified processes
    unsigned int n_pc = 0;
    for (auto& p : P.processes)
        n_pc += p.ids.size();
    pc = vector<vector<Compartment>>(n_pc, vector<Compartment>(S.size()));
    pci = vector<double>(pc.size(), 0.);
    pco = vector<double>(pc.size(), 0.);
}

// Do seeding and calculate contagiousness
void Population::Contagiousness(Parameters& P, Randomizer& Rand, double t, vector<double>& contag)
{
    auto add = [&](unsigned int age, double n)
    {
        n = min(n, S[age]);
        S[age] -= n;
        E[age].Add(P, Rand, n, P.pop[p].dE);
    };

    // Do seeding
    while (seed_row < P.pop[p].seed_times.size() && t >= P.pop[p].seed_times[seed_row])
    {
        if (P.deterministic)
        {
            for (unsigned int a = 0; a < S.size(); ++a)
                add(a, P.pop[p].dist_seed_ages.weights[a]);
        }
        else
        {
            Rand.Multinomial(1, P.pop[p].dist_seed_ages.weights, P.pop[p].dist_seed_ages.storage);
            for (unsigned int a = 0; a < S.size(); ++a)
            {
                if (P.pop[p].dist_seed_ages.storage[a] == 1)
                {
                    add(a, 1);
                    break;
                }
            }
        }
        ++seed_row;
    }

    // Calculate contagiousness from this population
    for (unsigned int a = 0; a < contag.size(); ++a)
        contag[a] = (N[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip[a].Size() + P.pop[p].fIa[a] * Ia[a].Size() + P.pop[p].fIs[a] * Is[a].Size()) / N[a];
}

// Execute one time step's events
void Population::Tick(Parameters& P, Randomizer& Rand, double t, vector<double>& infec, Reporter& rep)
{
    // Calculate force of infection in this compartment
    lambda.assign(infec.size(), 0.0);
    for (unsigned int a = 0; a < lambda.size(); ++a)
        for (unsigned int b = 0; b < lambda.size(); ++b)
            lambda[a] += P.pop[p].u[a] * P.pop[p].cm(a,b) * infec[b];

    // Account for seasonality
    if (P.pop[p].season_A[0] != 0)
    {
        double f = 1.0 + P.pop[p].season_A[0] * cos(2. * M_PI * (t - P.pop[p].season_phi[0]) / P.pop[p].season_T[0]);
        for (unsigned int a = 0; a < lambda.size(); ++a)
            lambda[a] = lambda[a] * f;
    }

    // Account for importation
    for (unsigned int a = 0; a < lambda.size(); ++a)
        lambda[a] += P.pop[p].omega[a];

    // Helpers
    auto multinomial = [&](double n, vector<double>& p, vector<double>& nd_out, vector<unsigned int>& ni_out) {
        nd_out.resize(p.size(), 0.);
        if (P.deterministic)
        {
            for (unsigned int i = 0; i < p.size(); ++i)
                nd_out[i] = n * p[i];
        }
        else
        {
            ni_out.resize(p.size(), 0);
            Rand.Multinomial(n, p, ni_out);
            for (unsigned int i = 0; i < p.size(); ++i)
                nd_out[i] = ni_out[i];
        }
    };

    auto poisson = [&](double l) {
        if (P.deterministic)
            return l;
        else
            return (double)Rand.Poisson(l);
    };

    auto binomial = [&](double n, double p) {
        if (P.deterministic)
            return n * p;
        else
            return (double)Rand.Binomial(n, p);
    };

    auto num = [&](double n) {
        if (P.deterministic)
            return n;
        else
            return round(n);
    };

    // Do state transitions and reporting for each age group
    for (unsigned int a = 0; a < lambda.size(); ++a)
    {
        // 0. Report prevalences
        if (t == (int)t)
        {
            // Built-in states
            rep(t, p, a, 0) = S[a];
            rep(t, p, a, 1) = E[a].Size();
            rep(t, p, a, 2) = Ip[a].Size();
            rep(t, p, a, 3) = Is[a].Size();
            rep(t, p, a, 4) = Ia[a].Size();
            rep(t, p, a, 5) = R[a];
            rep(t, p, a, 6) = V[a];
            rep(t, p, a, 7) = V2[a];
            rep(t, p, a, 8) = Ea[a].Size();
            rep(t, p, a, 12) = lambda[a];

            // User-specified processes
            for (auto& process : P.processes)
            {
                for (unsigned int i = 0; i < process.p_cols.size(); ++i)
                    rep(t, p, a, process.p_cols[i]) = pc[process.p_ids[i]][a].Size();
            }

        }

        // 1. Built-in states

        // Vaccination and waning of natural immunity and vaccine protection
        // S -> V, V -> S, R -> S
        double nS_V   = min(S[a],        num(P.pop[p].v[a]   * P.pop[p].ev[a]  * S[a] / N[a] * P.time_step));
        double nS_V2  = min(S[a] - nS_V, num(P.pop[p].v2[a]  * P.pop[p].ev2[a] * S[a] / N[a] * P.time_step));
        double nV_V2  = min(V[a],        num(P.pop[p].v12[a] * P.pop[p].ev2[a] * S[a] / N[a] * P.time_step));
        double nV_S   = binomial(V[a],  1.0 - exp(-P.pop[p].wv[a]  * P.time_step));
        double nV2_S  = binomial(V2[a], 1.0 - exp(-P.pop[p].wv2[a] * P.time_step));
        double nR_S   = binomial(R[a],  1.0 - exp(-P.pop[p].wn[a]  * P.time_step));

        S[a]  -= nS_V;
        V[a]  += nS_V;
        V[a]  -= nV_S;
        S[a]  += nV_S;
        V[a]  -= nV_V2;
        V2[a] += nV_V2;
        S[a]  -= nS_V2;
        V2[a] += nS_V2;
        V2[a] -= nV2_S;
        S[a]  += nV2_S;
        R[a]  -= nR_S;
        S[a]  += nR_S;

        // S -> E
        double nS_E = binomial(S[a], 1.0 - exp(-lambda[a] * P.time_step));
        S[a] -= nS_E;
        E[a].Add(P, Rand, nS_E, P.pop[p].dE);

        // R -> E/Ea
        double nR_EEa = binomial(R[a], 1.0 - exp(-lambda[a] * (1 - P.pop[p].pi_r[a]) * P.time_step));
        double nR_Ea = binomial(nR_EEa, P.pop[p].pd_ri[a]);
        double nR_E = nR_EEa - nR_Ea;
        R[a] -= nR_EEa;
        E[a].Add(P, Rand, nR_E, P.pop[p].dE);
        Ea[a].Add(P, Rand, nR_Ea, P.pop[p].dEa);

        // V -> E/Ea
        double nV_EEa = binomial(V[a], 1.0 - exp(-lambda[a] * (1 - P.pop[p].ei_v[a]) * P.time_step));
        double nV_Ea = binomial(nV_EEa, P.pop[p].ed_vi[a]);
        double nV_E = nV_EEa - nV_Ea;
        V[a] -= nV_EEa;
        E[a].Add(P, Rand, nV_E, P.pop[p].dE);
        Ea[a].Add(P, Rand, nV_Ea, P.pop[p].dEa);

        // V2 -> E/Ea
        double nV2_EEa = binomial(V2[a], 1.0 - exp(-lambda[a] * (1 - P.pop[p].ei_v2[a]) * P.time_step));
        double nV2_Ea = binomial(nV2_EEa, P.pop[p].ed_vi2[a]);
        double nV2_E = nV2_EEa - nV2_Ea;
        V2[a] -= nV2_EEa;
        E[a].Add(P, Rand, nV2_E, P.pop[p].dE);
        Ea[a].Add(P, Rand, nV2_Ea, P.pop[p].dEa);

        // E -> Ip/Ia
        double nE_Ipa = E[a].Mature();
        double nE_Ip = binomial(nE_Ipa, P.pop[p].y[a]);
        double nE_Ia = nE_Ipa - nE_Ip;
        Ip[a].Add(P, Rand, nE_Ip, P.pop[p].dIp);
        Ia[a].Add(P, Rand, nE_Ia, P.pop[p].dIa);

        // Ea -> Ia
        double nEa_Ia = Ea[a].Mature();
        Ia[a].Add(P, Rand, nEa_Ia, P.pop[p].dEa);

        // Ip -> Is -- also, true case onsets
        double nIp_Is = Ip[a].Mature();
        Is[a].Add(P, Rand, nIp_Is, P.pop[p].dIs);

        // Reported cases
        double n_to_report = binomial(nIp_Is, P.pop[p].rho[a]);
        C[a].Add(P, Rand, n_to_report, P.pop[p].dC);
        double n_reported = C[a].Mature();

        // Is -> R
        double nIs_R = Is[a].Mature();
        R[a] += nIs_R;

        // Ia -> R
        double nIa_R = Ia[a].Mature();
        R[a] += nIa_R;

        // 2. User-specified processes
        fill(pco.begin(), pco.end(), -1.);

        for (auto& process : P.processes)
        {
            // Determine number of individuals entering the process
            double n_entering = 0.;
            switch (process.source_id)
            {
                case srcS:
                    n_entering = nS_E; break;
                case srcNewE:
                    n_entering = nS_E + nR_E + nV_E + nV2_E; break;
                case srcNewEa:
                    n_entering = nR_Ea + nV_Ea + nV2_Ea; break;
                case srcNewEEa:
                    n_entering = nS_E + nR_EEa + nV_EEa + nV2_EEa; break;
                case srcE:
                    n_entering = nE_Ipa; break;
                case srcEp:
                    n_entering = nE_Ip; break;
                case srcEa:
                    n_entering = nE_Ia; break;
                case srcIp:
                    n_entering = nIp_Is; break;
                case srcIs:
                    n_entering = nIs_R; break;
                case srcIa:
                    n_entering = nIa_R; break;
                case srcI:
                    n_entering = nIs_R + nIa_R; break;
                case srcCasesReported:
                    n_entering = n_to_report; break;
                default:
                    n_entering = pco[process.source_id];
                    if (n_entering < 0)
                        throw logic_error("Process sourced from unset user process. Have user processes been specified in the right order?");
                    break;
            }

            multinomial(n_entering, process.prob[a], nd_out, ni_out);

            // Seed and mature this process's compartments
            unsigned int c = 0;
            for (unsigned int compartment_id : process.ids)
            {
                if (compartment_id != Null)
                {
                    pc[compartment_id][a].Add(P, Rand, nd_out[c], process.delays[c]);
                    pci[compartment_id] = nd_out[c];
                    pco[compartment_id] = pc[compartment_id][a].Mature();
                }
                ++c;
            }
        }

        // 3. Report incidence / outcidence
        rep(t, p, a, 9)  += nIp_Is;         // cases
        rep(t, p, a, 10) += n_reported;     // reported cases
        rep(t, p, a, 11) += nE_Ia + nEa_Ia; // subclinical infections

        // User-specified processes
        for (auto& process : P.processes)
        {
            for (unsigned int i = 0; i < process.i_cols.size(); ++i)
                rep(t, p, a, process.i_cols[i]) += pci[process.i_ids[i]];
            for (unsigned int i = 0; i < process.o_cols.size(); ++i)
                rep(t, p, a, process.o_cols[i]) += pco[process.o_ids[i]];
        }
    }

    // Births, deaths, aging
    double Ntot = accumulate(N.begin(), N.end(), 0.0);
    for (unsigned int a = N.size() - 1; ; --a)
    {
        // Births
        double B = poisson(Ntot * (exp(P.pop[p].B[a] * P.time_step) - 1.));

        // Deaths
        double death_prob = 1.0 - exp(-P.pop[p].D[a] * P.time_step);
        double DS   = binomial(S[a],  death_prob);
        double DV   = binomial(V[a],  death_prob);
        double DV2  = binomial(V2[a], death_prob);
        double DR   = binomial(R[a],  death_prob);

        // Changes
        N[a] += B;
        S[a] += B;

        S[a]  -= DS;
        V[a]  -= DV;
        V2[a] -= DV2;
        double DE  = E[a] .RemoveProb(P, Rand, death_prob);
        double DEa = Ea[a].RemoveProb(P, Rand, death_prob);
        double DIp = Ip[a].RemoveProb(P, Rand, death_prob);
        double DIa = Ia[a].RemoveProb(P, Rand, death_prob);
        double DIs = Is[a].RemoveProb(P, Rand, death_prob);
        R[a]  -= DR;

        N[a]  -= DS + DV + DV2 + DE + DEa + DIp + DIa + DIs + DR;

        // Agings
        if (a != lambda.size() - 1)
        {
            double age_prob = 1.0 - exp(-P.pop[p].A[a] * P.time_step);
            double AS   = binomial(S[a],  age_prob);
            double AV   = binomial(V[a],  age_prob);
            double AV2  = binomial(V2[a], age_prob);
            double AR   = binomial(R[a],  age_prob);

            S[a]      -= AS;
            S[a + 1]  += AS;
            V[a]      -= AV;
            V[a + 1]  += AV;
            V2[a]     -= AV2;
            V2[a + 1] += AV2;
            double AE  = E[a] .MoveProb(E [a + 1], P, Rand, age_prob);
            double AEa = Ea[a].MoveProb(Ea[a + 1], P, Rand, age_prob);
            double AIp = Ip[a].MoveProb(Ip[a + 1], P, Rand, age_prob);
            double AIa = Ia[a].MoveProb(Ia[a + 1], P, Rand, age_prob);
            double AIs = Is[a].MoveProb(Is[a + 1], P, Rand, age_prob);
            R[a]      -= AR;
            R[a + 1]  += AR;

            N[a]      -= AS + AV + AV2 + AE + AEa + AIp + AIa + AIs + AR;
            N[a + 1]  += AS + AV + AV2 + AE + AEa + AIp + AIa + AIs + AR;
        }

        if (a == 0)
            break;
    }
}

// Print full population details
void Population::DebugPrint() const
{
    auto vecprint = [&](const vector<double>& vec, string name) {
        cout << name;
        for (auto& v : vec)
            cout << " " << v;
        cout << "\n";
    };

    auto comprint = [&](const vector<Compartment>& comp, string name) {
        cout << name;
        for (unsigned int c = 0; c < comp.size(); ++c) {
            cout << "element " << c << "\n";
            comp[c].DebugPrint();
        }
    };

    vecprint(lambda, "lambda");
    vecprint(N, "N");
    vecprint(S, "S");
    vecprint(R, "R");
    vecprint(V, "V");
    vecprint(V2, "V2");
    comprint(E, "E");
    comprint(Ea, "Ea");
    comprint(Ip, "Ip");
    comprint(Ia, "Ia");
    comprint(Is, "Is");
    comprint(C, "C");
    cout << "seed_row " << seed_row << " p " << p << "\n";
    for (auto& c : pc)
        comprint(c, "User");

    cout << "\n\n";
}


Metapopulation::Metapopulation(Parameters& P)
{
    P.changes.Capture(P);

    for (unsigned int i = 0; i < P.pop.size(); ++i)
        pops.push_back(Population(P, i));
}

// Execute one time step's events
bool Metapopulation::Tick(Parameters& P, Randomizer& Rand, double t, unsigned int ts, Reporter& rep)
{
    // Apply any changes to parameters
    P.changes.Apply(P, t);

    unsigned int n_ages = P.pop[0].size.size();

    // Calculate contagiousness from each population
    // NOTE -- 'contag' subscripted first by j, then by a.
    // It's the effective number of infectious individuals FROM subpop j of age a.
    contag.assign(pops.size(), vector<double>(n_ages, 0.0));
    for (unsigned int j = 0; j < pops.size(); ++j)
        pops[j].Contagiousness(P, Rand, t, contag[j]);

    // note -- 'infec' subscripted first by i, then by a
    // It's the effective number of infectious individuals who are CURRENTLY IN subpop i of age a.
    infec.assign(pops.size(), vector<double>(n_ages, 0.0));
    for (unsigned int i = 0; i < pops.size(); ++i)
        for (unsigned int j = 0; j < pops.size(); ++j)
            for (unsigned int a = 0; a < n_ages; ++a)
                infec[i][a] += P.travel(j, i) * contag[j][a] * (j != i ? P.pop[j].tau[a] : 1.0);

    // Update populations
    //#pragma omp parallel for schedule(dynamic) reduction(&&:keep_going)
    for (unsigned int i = 0; i < pops.size(); ++i)
        pops[i].Tick(P, Rand, t, infec[i], rep);

    // Run observer at the last time step of each day.
    if (t + P.time_step == int(t + P.time_step))
        return CppObserver(P, Rand, rep, (int)t, x);

    return true;
}

void Metapopulation::Run(Parameters& P, Randomizer& Rand, Reporter& rep, vector<double> x_fit)
{
    x = x_fit;

    #ifdef _OPENMP
    omp_set_num_threads(6);
    #endif

    // Run simulation
    unsigned int time_steps = (1 + P.time1 - P.time0) / P.time_step;
    for (unsigned int ts = 0; ts < time_steps; ++ts)
    {
        if (!Tick(P, Rand, P.time0 + ts * P.time_step, ts, rep))
            break;
    }
}

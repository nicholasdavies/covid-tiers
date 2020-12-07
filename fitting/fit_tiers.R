library(data.table)
library(ggplot2)
library(lubridate)
library(here)
library(cowplot)
library(readxl)
library(sn)
library(qs)
library(stringr)

N_THREADS = 36
REP_START = 7
N_REPS = 8
BURN_IN = 1500
BURN_IN_FINAL = 2500

uk_covid_data_path = "../uk_covid_data/";
datapath = function(x) paste0(uk_covid_data_path, x)

# Sys.setenv(MP_DUPLICATE_LIB_OK=TRUE) # RCB needs to override this environment variable to run

#
# SETUP
#

# set up covidm
cm_path = "./covidm_for_fitting/";
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R")) # RCB needs to move (backed up) Makevars file into .R repo for this to work with -fopenmp
#source("./locate.R");
popUK = readRDS(datapath("fitting/data/popNHS.rds"));
matricesUK = readRDS(datapath("fitting/data/matricesNHS.rds"));

cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
source("./distribution_fit.R");
source("./spim_output.R");

delay_normal = function(mu, sd, t_max, t_step)
{
    t_points = seq(0, t_max, by = t_step);
    heights = pnorm(t_points + t_step/2, mu, sd) -
        pnorm(pmax(0, t_points - t_step/2), mu, sd);
    return (data.table(t = t_points, p = heights / sum(heights)))
}

#
# DATA
#

nhs_regions = popUK[, unique(name)]
pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100

all_data = qread(datapath("processed-data-2020-10-23-fix.qs"))
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]]
sero = all_data[[4]]

#
# FITTING
#

# NUMBER OF REGIONS TO FIT
N_REG = 12;

# Build parameters for NHS regions
params = cm_parameters_SEI3R(nhs_regions[1:N_REG], deterministic = T, date_start = "2020-01-01", date_end = "2020-10-24",
    dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
    dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
    dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
    dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p);
params = cm_split_matrices_ex_in(params, 15)

# school terms
school_close =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
school_reopen = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");

# Load age-varying symptomatic rate
covid_scenario = qread(datapath("fitting/data/2-linelist_both_fit_fIa0.5-rbzvih.qs"));
covu = unname(rep(colMeans(covid_scenario[,  5:12]), each = 2));
covy = unname(rep(colMeans(covid_scenario[, 13:20]), each = 2));

for (i in seq_along(params$pop)) {
    params$pop[[i]]$u = covu / mean(covu);
    params$pop[[i]]$y = covy;
}

# Health burden processes

# Probability of ICU given hospitalisation (derived from CO-CIN)
picu_cocin_func = function(age)
{
    x = c(-0.1309118, 0, 17.2398874, 65.7016492, 100)
    y = c(-2.1825091, -2.1407043, -1.3993552, -1.2344361, -8.8191062)
    p = splinefun(x, y)(age)
    exp(p) / (1 + exp(p))
}
picu_cocin = picu_cocin_func(0:85)

# Infection fatality rate (derived from Levin et al., preprint)
ifr_levin = 100 * exp(-7.56 + 0.121 * 0:85) / (100 + exp(-7.56 + 0.121 * 0:85)) / 100;

# Infection hospitalisation rate (derived from Salje et al., Science)
ihr_salje = exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85));

# Amalgamate probabilities
probabilities = data.table(age = 0:85, ihr = ihr_salje, ifr = ifr_levin, picu = picu_cocin)
probabilities[, age_group := pmin(15, age %/% 5)]
probabilities = probabilities[, lapply(.SD, mean), by = age_group, .SDcols = 2:4]

# Create model burden processes
P.hosp     = probabilities[, ihr];
P.critical = probabilities[, ihr * picu];
P.severe   = probabilities[, ihr * (1 - picu)];
P.death    = probabilities[, ifr];

burden_processes = list(
    list(source = "E", type = "multinomial", names = c("to_hosp", "null"), report = c("", ""),
        prob = matrix(c(P.hosp, 1 - P.hosp), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(6.0 + 2.5, 0.71, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "E", type = "multinomial", names = c("to_icu", "null"), report = c("", ""),
        prob = matrix(c(P.critical, 1 - P.critical), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(9.6 + 2.5, 1.91, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "to_hosp", type = "multinomial", names = "hosp", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(11.08, 1.202, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "to_icu", type = "multinomial", names = "icu", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(13.33, 1.25, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "E", type = "multinomial", names = c("death", "null"), report = c("o", ""),
        prob = matrix(c(P.death, 1 - P.death), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_lnorm(15, 0.9, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    # note -- delay should actually be gamma (symptom onset) plus normal -- to be fixed.
    # from Borremans et al
    list(source = "E", type = "multinomial", names = "to_lfia_positive", report = c(""),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(delay_normal(11.9 + 2.5, 5.3, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "to_lfia_positive", type = "multinomial", names = "lfia_positive", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(224, 1, 730, 0.25)$p, nrow = 1, byrow = T)),
    
    list(source = "to_hosp", type = "multinomial", names = "hosp_undetected", report = c("po"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(6, 1, 60, 0.25)$p, nrow = 1, byrow = T)),
    
    list(source = "S", type = "multinomial", names = "to_pcr_positive", report = c(""),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(2.76, 4.79, 60, 0.25)$p), nrow = 1, byrow = T)),

    list(source = "to_pcr_positive", type = "multinomial", names = "pcr_positive", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(8.47, 1.96, 730, 0.25)$p, nrow = 1, byrow = T))

)
params$processes = burden_processes

# changes
schedule_all = readRDS(datapath("schedule3-2020-11-29.rds"));
schedule = list();
for (i in seq_along(schedule_all)) {
    if (schedule_all[[i]]$pops < N_REG) {
        schedule[[length(schedule) + 1]] = schedule_all[[i]]
    }
}

# Remove NAs
for (i in seq_along(schedule)) {
    for (j in seq_along(schedule[[i]]$values)) {
        if (any(is.na(schedule[[i]]$values[[j]]))) {
            schedule[[i]]$values[[j]] = ifelse(is.na(schedule[[i]]$values[[j]]), prev, schedule[[i]]$values[[j]])
        }
        prev = schedule[[i]]$values[[j]];
    }
}
params$schedule = schedule


#
# Individual fits
#

# create cpp changes
cpp_chgI = function()
{
    c(
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double x, vector<double>& curve) {',
        '    if (x < 0) return curve[0];',
        '    if (x >= (curve.size() - 1) * 0.01) return curve.back();',
        '    unsigned int i = (unsigned int)(x * 100);',
        '    double f = x * 100 - i;',
        '    return f * curve[i] + (1 - f) * curve[i + 1];',
        '};',
        
        'auto odds = [&](double x, double lo) {',
        '    double a = x / (1 - x);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0] = odds(P.processes[0].prob[g][0], x[7]);',
        '    P.processes[0].prob[g][1] = 1 - P.processes[0].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0] = odds(P.processes[1].prob[g][0], x[7] + x[6]);',
        '    P.processes[1].prob[g][1] = 1 - P.processes[1].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0] = min(1.0, P.processes[4].prob[g][0] * x[5]);',
        '    P.processes[4].prob[g][1] = 1 - P.processes[4].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g] = P.pop[0].u[g] * x[1];',
        '}',

        '// Delays to admission to hospital & ICU',
        'P.processes[0].delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[1].delays[1] = delay_gamma(x[8], 1.91, 60, 0.25);',
        ### '// Lengths of stay',
        ### 'P.processes[2].delays[0] = delay_lnorm(x[9], 1.35, 60, 0.25);',
        ### 'P.processes[3].delays[0] = delay_lnorm(x[8], 1.31, 60, 0.25);',
        '// Death',
        'P.processes[4].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        ### '// Waning of antibodies',
        ### 'P.processes[6].delays[0] = delay_gamma(x[14], 1, 730, 0.25);',
        'P.pop[0].rho = vector<double>(P.pop[0].rho.size(), 1);',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[9], -x[10], x[11]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[15]);',
        ### 'P.changes.ch[5].times[0] = x[24];',
        'P.changes.ch[5].times[0] = x[17];',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 || t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}'
    )
}

# create c++ likelihood components
cpp_likI = function(params, ld, sitreps, sero, virus, popid)
{
    refdate = ld[, max(date)] + 1;

    ret = NULL;

    add = function(ret, new)
    {
        c(ret, new);
    }

    if (any(!is.na(ld$N))) {
        ret = add(ret,
            ld[!is.na(N), paste0(
                'll += nbinom_gammaconf(', N, ', max(0.1, dyn("death_o", ', as.numeric(ymd(date) - ymd(params$date0)),
                ', {', 0, '}, {})), 10, ', as.numeric(refdate - date), ', 4.5, 2.7);')]) # XXX NOTE CHANGE 4.5, 2.7 was x[10], x[11], from deathconf_delay and deathconf_shape
    }

    if (any(!is.na(sitreps$n_in_itu))) {
        ret = add(ret,
            sitreps[!is.na(n_in_itu), paste0(
                'll += nbinom(', n_in_itu, ', max(0.1, dyn("icu_p", ', as.numeric(ymd(date) - ymd(params$date0)),
                ', {', 0, '}, {})), 20);')])
    }

    if (any(!is.na(sitreps$n_in_all_beds))) {
        ret = add(ret,
            sitreps[!is.na(n_in_all_beds), paste0(
                'll += nbinom(', n_in_all_beds, ', max(0.1, dyn("hosp_p", ', as.numeric(ymd(date) - ymd(params$date0)), ', {0}, {})',
                ' - dyn("hosp_undetected_p", ', as.numeric(ymd(date) - ymd(params$date0)), ', {0}, {})), 20);')])
    }

    if (any(!is.na(sitreps$n_admitted_diagnosed))) {
        ret = add(ret,
            sitreps[!is.na(n_admitted_diagnosed), paste0(
                'll += nbinom(', n_admitted_diagnosed, ', max(0.1, dyn("hosp_undetected_o", ', as.numeric(ymd(date) - ymd(params$date0)), ', {0}, {})), 20);')])
    }

    if (nrow(virus) > 0) {
        ret = add(ret,
            virus[, paste0(
                '{ double t0 = ', as.numeric(Start.date - ymd(params$date0)), ', t1 = ', as.numeric(End.date - ymd(params$date0)), '; ',
                'double virus_prev = 0; for (double tt = t0; tt <= t1; ++tt) virus_prev += dyn("pcr_positive_p", tt, {0}, {}) / (t1 - t0 + 1); ',
                'll += 10 * skewnorm(virus_prev / ', sum(params$pop[[1]]$size), ', ', xi, ', ', omega, ', ', alpha, '); }')])
    }

    if (nrow(sero) > 0) {
        ret = add(ret,
            sero[, paste0(
                '{ double t0 = ', as.numeric(Start.date - ymd(params$date0)), ', t1 = ', as.numeric(End.date - ymd(params$date0)), '; ',
                'double sero_prev = 0; for (double tt = t0; tt <= t1; ++tt) sero_prev += dyn("lfia_positive_p", tt, {0}, {}) / (t1 - t0 + 1); ',
                'll += skewnorm(sero_prev / ', sum(params$pop[[1]]$size), ', ', xi, ', ', omega, ', ', alpha, '); }')])
    }
    
    return (ret)
}

# Observer function (basic)
cpp_obsI = function(P.death)
{
    c(
        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        'dyn.Obs(t, 0, 0, 0) = estimate_Rt(P, dyn, t, 0, 50);',
        'dyn.Obs(t, 0, 3, 0) = estimate_R0(P, t, 0, 50);',

        # increase in young person mobility
        'if (t == 182) {',
        '    double mode = 0.2;',
        '    double conc = x[12];',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a)',
        '        P.pop[0].u[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '}',
        'if (t == 213) {',
        '    double mode = 0.2;',
        '    double conc_prev = x[12];',
        '    double conc = x[13];',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        'if (t == 244) {',
        '    double mode = 0.2;',
        '    double conc_prev = x[13];',
        '    double conc = x[14];',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        # changing CFR
        'if ((int)t % 7 == 0) {',
        '    double ifr_table[] = {', paste(P.death, collapse = ", "), '};',
        '    for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '        double ifr = ifr_table[g];',
        '        P.processes[4].prob[g][0] = asc(min(365.0, t/365.0), ifr * x[5], ifr * x[5] * x[16], -2.9, 7.8);',
        '        P.processes[4].prob[g][1] = 1 - P.processes[3].prob[g][0];',
        '    }',
        '}',
        # changing detection rate in patients admitted to hospital
        ### 'double detection = asc(min(t / 365.0, 1.0), x[20], x[21], -x[22], x[23]);',
        'double detection = asc(min(t / 365.0, 1.0), 14, 1, -5.86, 33.4);',
        'P.processes[7].delays[0] = delay_gamma(detection, 0.59, 60, 0.25);'
    )
}

# Fitting
priorsI = list(
    tS = "U 0 60",
    u = "N 0.07 0.01 T 0.04 0.2",
    death_mean = "N 15 2 T 5 30",    # <<< co-cin was N 15 2
    death_shape = "N 1.9 0.2 T 0.1 3", # <<< co-cin was N 1.9 0.2
    admission = "N 8 1 T 4 20", # was N 7.5 1 T 4 20 <<< co-cin
    cfr_rel = "N 1 0.1 T 0.1 4", # <<< co-cin
    icu_rlo = "N 0 0.1 T -2 2", # <<< co-cin
    hosp_rlo = "N 0 0.1 T -2 2", #"L 0 4 T 0.1 1.5",
    ###icu_los = "N 14 0.5 T 7 18",      # <<< co-cin
    ###nonicu_los = "N 11.5 0.5 T 5 15.5",    # <<< co-cin
    icu_admission = "N 12.5 1 T 8 14", # was icu_admission = "N 11.1 1 T 8 14", <<< co-cin
    contact_final = "N 1 0.1 T 0 1",
    contact_s0 = "E 0.1 0.1",
    contact_s1 = "E 0.1 0.1",
    ###waning = "N 224 14 T 0 500",
    concentration1 = "N 2 .3 T 2 10", # was .5
    concentration2 = "N 2 .2 T 2 10", # was .4
    concentration3 = "N 2 .1 T 2 10", # was .2
    sep_boost = "N 1 0.05",
    cfr_rel2 = "N 0.45 0.1 T 0 1", # <<< co-cin was N 0.45 0.1
    ### detect1 = "N 0 100 T 0 200",
    ### detect2 = "N 0 0.25 T 0 5",
    ### detect_s0 = "N 0 20 T 0 40",
    ### detect_s1 = "N 0 20 T 0 40",
    sep_when = "U 224 264"
)

posteriorsI = list()
dynamicsI = list()
parametersI = list()

# Remove problematic virus entries
virus = virus[omega > 1e-9]
#virus = virus[Data.source != "REACT-1 R1"] # very low...

existing_file = paste0("../fit_pp", REP_START - 1, "-NEW.qs");
if (file.exists(existing_file)) {
    saved = qread(existing_file)
    posteriorsI = saved[[1]]
    parametersI = saved[[2]]
    rm(saved)
}

for (replic in REP_START:N_REPS) {

init_previous = TRUE
init_previous_amount = 1

# RCB checking execution time to test multithreading
time1 <- Sys.time()

# Loop through regions of England
for (p in c(1,3,4,5,6,9,10)) {
    paramsI = rlang::duplicate(params);
    paramsI$pop = list(rlang::duplicate(params$pop[[p]]));
    paramsI$travel = matrix(1, nrow = 1, ncol = 1);
    paramsI$schedule = list();
    j = 1;
    for (i in seq_along(params$schedule)) {
        if (p - 1 == params$schedule[[i]]$pops) {
            paramsI$schedule[[j]] = rlang::duplicate(params$schedule[[i]]);
            paramsI$schedule[[j]]$pops = 0;
            j = j + 1;
        }
    }

    # contact placeholder for tier 2
    paramsI$schedule[[3]] = rlang::duplicate(paramsI$schedule[[1]]);
    for (i in seq_along(paramsI$schedule[[3]]$values)) {
        paramsI$schedule[[3]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  0.2497655 / 100;
        paramsI$schedule[[3]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -0.2307939 / 100;
        paramsI$schedule[[3]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -1.5907698 / 100;
        paramsI$schedule[[3]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -3.4866544 / 100;
        paramsI$schedule[[3]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -3.4524518 / 100;
    }
    paramsI$schedule[[3]]$mode = "bypass";

    # contact placeholder for tier 3
    paramsI$schedule[[4]] = rlang::duplicate(paramsI$schedule[[1]]);
    for (i in seq_along(paramsI$schedule[[4]]$values)) {
        paramsI$schedule[[4]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  2.080457 / 100;
        paramsI$schedule[[4]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -8.045226 / 100;
        paramsI$schedule[[4]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -2.476266 / 100;
        paramsI$schedule[[4]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -10.144043 / 100;
        paramsI$schedule[[4]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -7.681244 / 100;
    }
    paramsI$schedule[[4]]$mode = "bypass";

    # contact multiplier for gradual contact change
    paramsI$schedule[[5]] = list(
        parameter = "contact",
        pops = 0,
        mode = "multiply",
        values = rep(list(rep(1, 8)), 366),
        times = 0:365
    )

    # contact multiplier for september boost
    paramsI$schedule[[6]] = list(
        parameter = "contact",
        pops = 0,
        mode = "multiply",
        values = list(rep(1, 8)),
        times = c(244)
    )

    ldI = rlang::duplicate(ld);
    ldI = ldI[pid == p - 1];
    sitrepsI = rlang::duplicate(sitreps);
    sitrepsI = sitrepsI[pid == p - 1];
    seroI = rlang::duplicate(sero);
    #seroI = seroI[pid == p - 1 & Data.source == "REACT-2"]; # first version
    seroI = seroI[pid == p - 1 & Data.source != "NHSBT"];   # second version
    virusI = rlang::duplicate(virus);
    #virusI = virusI[pid == p - 1 & Data.source == "ONS-CIS"]; # ONS fitting
    virusI = virusI[pid == p - 1 & Data.source %like% "REACT"]; # REACT-1 fitting

    # load user defined functions
    cm_source_backend(
        user_defined = list(
            model_v2 = list(
                cpp_changes = cpp_chgI(),
                cpp_loglikelihood = cpp_likI(paramsI, ldI, sitrepsI, seroI, virusI, p),
                cpp_observer = cpp_obsI(P.death)
            )
        )
    )
    
    priorsI2 = rlang::duplicate(priorsI)
    if (init_previous) {
        for (k in seq_along(priorsI2)) {
            pname = names(priorsI2)[k];
            if (pname == "waning" & replic == 42) {
                cat("Skipping waning.");
                next;
            }
            if (length(posteriorsI) >= p && pname %in% names(posteriorsI[[p]])) {
                init_values = quantile(posteriorsI[[p]][[pname]], c(0.25, 0.75));
                cat(paste0("Using IQR ", init_values[1], " - ", init_values[2], " for initial values of parameter ", pname, 
                    " with probability ", init_previous_amount, "\n"));
                priorsI2[[pname]] = paste0(priorsI2[[pname]], " I ", init_values[1], " ", init_values[2], " ", init_previous_amount);
                cat(paste0(priorsI2[[pname]], "\n"));
            } else {
                cat(paste0("Could not find init values for parameter ", pname, "\n"));
                cat(paste0(priorsI2[[pname]], "\n"));
            }
        }
    }

    postI = cm_backend_mcmc_test(cm_translate_parameters(paramsI), priorsI2,
        seed = 0, burn_in = ifelse(replic == N_REPS, BURN_IN_FINAL, BURN_IN), 
        iterations = 500, n_threads = N_THREADS, classic_gamma = T); # default burn_in = 2500, iterations = 1000, n_threads = 6
    setDT(postI)
    posteriorsI[[p]] = postI

    # # Sampling fits
    # paramsI2 = rlang::duplicate(paramsI)
    # paramsI2$time1 = as.character(ymd(paramsI$time1) + 56);
    # test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), postI, 500, seed = 0);
    # 
    # test = rbindlist(test)
    # test = test[, population := p]
    # dynamicsI[[p]] = test

    parametersI[[p]] = rlang::duplicate(paramsI)
    qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("../fit_pp", replic, "-progress-NEW.qs"))

    print(p)
}

# RCB timing check again
time2 <- Sys.time()
print(time2-time1)
# 45 mins for England

qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("../fit_pp", replic, "-NEW.qs"))

}

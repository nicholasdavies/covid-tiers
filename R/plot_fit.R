library(data.table)
library(ggplot2)
library(lubridate)
library(here)
library(cowplot)
library(readxl)
library(sn)
library(qs)

# Sys.setenv(MP_DUPLICATE_LIB_OK=TRUE) # RCB needs to override this environment variable to run

#
# SETUP
#

# set up covidm
cm_path = "~/Documents/covidm_MTPs/covidm_for_fitting/";
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R")) # RCB needs to move (backed up) Makevars file into .R repo for this to work with -fopenmp
source("~/Documents/covidm_MTPs/locate.R");
cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
source("~/Documents/covidm_MTPs/distribution_fit.R");
source("~/Documents/covidm_MTPs/spim_output.R");



#
# DATA
#

nhs_regions = popUK[, unique(name)]

# temp
new_data = fread("~/Dropbox/uk_covid_data/data-2020-10-23.csv")
new_data[, pid := match(location, nhs_regions) - 1]
ld = new_data[indicator == "death_inc_line", .(date, N = value, name = location, pid)];
sitreps = new_data[indicator != "death_inc_line", .(date = ymd(date), value_type = indicator, value, name = location, pid)]
sitreps = dcast(sitreps, date + name + pid ~ value_type, value.var = "value", fill = NA, fun.aggregate = function(x) x[1])
sitreps = sitreps[, .(date, n_in_itu = icu_prev, n_in_all_beds = hospital_prev, n_admitted_diagnosed = hospital_inc, name, pid)];

pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100
sero = fread("~/Dropbox/uk_covid_data/data/seroprev_nhs_regions_20201026.csv")
sero = cbind(sero, sero[, skew_normal_solve(pct(Central.estimate), pct(Lower.bound), pct(Upper.bound))]);
sero[, Start.date := dmy(Start.date)]
sero[, End.date := dmy(End.date)]
virus = fread("~/Dropbox/uk_covid_data/data/virusprev_nhs_regions_20201026.csv")
virus[, Start.date := dmy(Start.date)]
virus[, End.date := dmy(End.date)]
virus = cbind(virus, virus[, skew_normal_solve(pct(Central.estimate), pct(Lower.bound), pct(Upper.bound))]);

sero[, pid := match(NHS.region, nhs_regions) - 1]
virus[, pid := match(NHS.region, nhs_regions) - 1]

ggplot(virus) + 
    geom_linerange(aes(x = Start.date + as.numeric(End.date - Start.date) / 2, 
        ymin = pct(Lower.bound), ymax = pct(Upper.bound), colour = str_remove_all(Data.source, " R[0-9]"))) + 
    facet_wrap(~NHS.region) + 
    theme(legend.position = "bottom") +
    labs(x = "Date", y = "Prevalence", colour = "Study")

# Add England to deaths series
ld = rbind(ld, 
    ld[!name %in% c("Northern Ireland", "Scotland", "Wales"), .(N = sum(N), name = "England", pid = 1), by = date]
)

# Add England to sitrep series
sitreps = rbind(sitreps,
    sitreps[!name %in% c("Northern Ireland", "Scotland", "Wales"), 
        .(n_in_itu = sum(n_in_itu, na.rm = T), n_in_all_beds = sum(n_in_all_beds, na.rm = T), n_admitted_diagnosed = sum(n_admitted_diagnosed, na.rm = T),
            name = "England", pid = 1),
        by = date]
)

all_data = qread("~/Dropbox/uk_covid_data/processed-data-2020-10-23.qs")
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]]
sero = all_data[[4]]


# SUMMARISING FUNCS
SPIM_output_1var = function(summ, varname, outname, cyear, cmonth, cday, ymd_from, ymd_to)
{
    quants = c(0.025, 0.5, 0.975)

    # new code...
    rq = data.table(run = 1:summ[, max(run)], q = runif(summ[, max(run)]))
    summ = merge(summ, rq, by = "run")
    summ[, rv := qnbinom(q, size = 20, mu = get(varname))]
    qsumm = summ[, as.list(quantile(rv, quants)), by = .(t, population)]

    qsumm[, date := t + ymd("2020-01-01")]

    data_out = qsumm[date %between% c(ymd_from, ymd_to),
        .(Group = "LSHTM", Model = "Transmission", ModelType = "Multiple", Version = 2,
        `Creation Day` = cday, `Creation Month` = cmonth, `Creation Year` = cyear,
        `Day of Value` = day(date), `Month of Value` = month(date), `Year of Value` = year(date),
        Geography = population, ValueType = outname, Value = `50%`,
        `2.5%`, `50%`, `97.5%`)]

    names(data_out)[names(data_out) %like% "\\%$"] = paste("Quantile", quants);
    data_out
}

SPIM_output_full = function(test, cyear, cmonth, cday, ymd_from, ymd_to)
{
    cat("Summarizing variables...\n");
    summ = test[, .(death_o = sum(death_o)), by = .(run, t, population)]
    summ3 = test[, .(icu_p = sum(icu_p)), by = .(run, t, population)]
    summ4 = test[, .(bed_p = sum(pmax(0, hosp_p - hosp_undetected_p))), by = .(run, t, population)]
    summ34 = test[, .(admissions = sum(hosp_undetected_o)), by = .(run, t, population)]
    summ_pcr_p = test[, .(pcr_p = sum(pcr_positive_p)), by = .(run, t, population)]
    summ_sero_p = test[, .(sero_p = sum(lfia_positive_p)), by = .(run, t, population)]

    cat("Running quantiles...\n");
    w = rbind(
        SPIM_output_1var(summ, "death_o", "type28_death_inc_line", cyear, cmonth, cday, ymd_from, ymd_to),
        SPIM_output_1var(summ3, "icu_p", "icu_prev", cyear, cmonth, cday, ymd_from, ymd_to),
        SPIM_output_1var(summ4, "bed_p", "hospital_prev", cyear, cmonth, cday, ymd_from, ymd_to),
        SPIM_output_1var(summ34, "admissions", "hospital_inc", cyear, cmonth, cday, ymd_from, ymd_to),
        SPIM_output_1var(summ_pcr_p, "pcr_p", "pcr_prev", cyear, cmonth, cday, ymd_from, ymd_to),
        SPIM_output_1var(summ_sero_p, "sero_p", "sero_prev", cyear, cmonth, cday, ymd_from, ymd_to)
    )
    return (w)
}

pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100

make_data = function(ld, sitreps, virus, sero)
{
    rbind(
        ld[, .(ValueType = "type28_death_inc_line", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = N, ymax = NA)],
        sitreps[, .(ValueType = "icu_prev", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_in_itu, ymax = NA)],
        sitreps[, .(ValueType = "hospital_prev", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_in_all_beds, ymax = NA)],
        sitreps[, .(ValueType = "hospital_inc", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_admitted_diagnosed, ymax = NA)],
        virus[Data.source %like% "REACT", .(ValueType = "pcr_prev", Geography = NHS.region,
            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))],
        sero[!Data.source %like% "NHS BT", .(ValueType = "sero_prev", Geography = NHS.region,
            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))]
    )
}



# Load fit
# saved = qread("~/Desktop/uk-fit-nick-2020-11-02.qs")
# posteriorsI = saved[[1]]
# dynamicsI = saved[[2]]
# parametersI = saved[[3]]
# rm(saved)

saved = qread("./fit_pp3.qs")
posteriorsI = saved[[1]]
parametersI = saved[[2]]
rm(saved)

# Recover dynamics
# TODO this needs to be saved in the qs file above instead of being copied out here.
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
        'P.processes[1].delays[1] = delay_gamma(x[10], 1.91, 60, 0.25);',
        '// Lengths of stay',
        'P.processes[2].delays[0] = delay_lnorm(x[9], 1.35, 60, 0.25);',
        'P.processes[3].delays[0] = delay_lnorm(x[8], 1.31, 60, 0.25);',
        '// Death',
        'P.processes[4].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        '// Waning of antibodies',
        'P.processes[6].delays[0] = delay_gamma(x[14], 1, 730, 0.25);',
        'P.pop[0].rho = vector<double>(P.pop[0].rho.size(), 1);',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[11], -x[12], x[13]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[18]);',
        'P.changes.ch[5].times[0] = x[24];',

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
        
        # # seasonality
        # 'P.pop[0].season_phi[0] = 0;',
        # 'P.pop[0].season_T[0] = 365.25;',
        # 'P.pop[0].season_A[0] = x[24];'
    )
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
        # # pressure month
        # 'if (t == 80 /*(int)x[11]*/) {',
        # '    P.processes[1].delays[0] = delay_gamma(x[8] * x[10], 3, 60, 0.25);',
        # '    P.processes[2].delays[0] = delay_gamma(x[9] * x[10], 3, 60, 0.25);',
        # '} else if (t == 80 + 30 /*(int)x[11] + (int)x[12]*/) {',
        # '    P.processes[1].delays[0] = delay_gamma(x[8], 3, 60, 0.25);',
        # '    P.processes[2].delays[0] = delay_gamma(x[9], 3, 60, 0.25);',
        # '}',
        # increase in young person mobility
        'if (t == 182) {',
        '    double mode = 0.2;',
        '    double conc = x[15];',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a)',
        '        P.pop[0].u[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '}',
        'if (t == 213) {',
        '    double mode = 0.2;',
        '    double conc_prev = x[15];',
        '    double conc = x[16];',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        'if (t == 244) {',
        '    double mode = 0.2;',
        '    double conc_prev = x[16];',
        '    double conc = x[17];',
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
        '        P.processes[4].prob[g][0] = asc(min(365.0, t/365.0), ifr * x[5], ifr * x[5] * x[19], -2.9, 7.8);',
        '        P.processes[4].prob[g][1] = 1 - P.processes[3].prob[g][0];',
        '    }',
        '}',
        # changing detection rate in patients admitted to hospital
        'double detection = asc(min(t / 365.0, 1.0), x[20], x[21], -x[22], x[23]);',
        'P.processes[7].delays[0] = delay_gamma(detection, 0.59, 60, 0.25);'
        
        # # test - susceptibility of children as norm
        # 'if (t == 244) {',
        # '    P.pop[0].u[0] = P.pop[0].u[4];',
        # '    P.pop[0].u[1] = P.pop[0].u[4];',
        # '    P.pop[0].u[2] = P.pop[0].u[4];',
        # '    P.pop[0].u[3] = P.pop[0].u[4];',
        # '}'
    )

}


# Loop through regions of England
dynamicsI = list()
for (p in c(1,3,4,5,6,9,10))  {

    # load user defined functions
    cm_source_backend(
        user_defined = list(
            model_v2 = list(
                cpp_changes = cpp_chgI(),
                cpp_loglikelihood = "",
                cpp_observer = cpp_obsI(P.death)
            )
        )
    )
    # Sampling fits
    paramsI2 = rlang::duplicate(parametersI[[p]])
    paramsI2$time1 = as.character(ymd(paramsI2$time1) + 56);
    test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[p]], 500, seed = 0);

    test = rbindlist(test)
    test = test[, population := p]
    dynamicsI[[p]] = test

    print(p)
}

lapply(posteriorsI, function(x) list(maxlp = max(x$lp), meanlp = mean(x$lp)))




#
# PLOT FITS.
#

# Remove problematic virus entries
virus = virus[omega > 1e-9]
virus = virus[Data.source != "REACT-1 R1"] # very low...

# Calculate total population
popsize = NULL
for (i in seq_along(parametersI)) {
    if (!is.null(parametersI[[i]])) {
        popsize = rbind(popsize,
            data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = sum(parametersI[[i]]$pop[[1]]$size))
        )
    }
}


# Create formatted output
test = rbindlist(dynamicsI, fill = TRUE)
test[, population := nhs_regions[population]]
output = SPIM_output_full(test[!population %in% c("England", "United Kingdom", "Wales", "Scotland", "Northern Ireland")], 2020, 11, 1, "2020-01-01", "2020-10-23")
output[, d := make_date(`Year of Value`, `Month of Value`, `Day of Value`)]

# Make data to output
data = make_data(ld, sitreps, virus, sero)[!Geography %in% c("England", "United Kingdom", "Wales", "Scotland", "Northern Ireland")]
data = merge(data, popsize, by = "Geography")
data[ValueType %in% c("pcr_prev", "sero_prev"), ymin := ymin * population_size]
data[ValueType %in% c("pcr_prev", "sero_prev"), y    := y    * population_size]
data[ValueType %in% c("pcr_prev", "sero_prev"), ymax := ymax * population_size]

# Relabel outcomes
output[ValueType == "hospital_inc", ValueType := "Admissions"]
output[ValueType == "hospital_prev", ValueType := "Hospital beds"]
output[ValueType == "icu_prev", ValueType := "ICU beds"]
output[ValueType == "pcr_prev", ValueType := "PCR positives"]
output[ValueType == "sero_prev", ValueType := "Seropositives"]
output[ValueType == "type28_death_inc_line", ValueType := "Deaths (28 d)"]
output[, ValueType := factor(ValueType, levels = c("Deaths (28 d)", "Admissions", "Hospital beds", "ICU beds", "PCR positives", "Seropositives"))]

data[ValueType == "hospital_inc", ValueType := "Admissions"]
data[ValueType == "hospital_prev", ValueType := "Hospital beds"]
data[ValueType == "icu_prev", ValueType := "ICU beds"]
data[ValueType == "pcr_prev", ValueType := "PCR positives"]
data[ValueType == "sero_prev", ValueType := "Seropositives"]
data[ValueType == "type28_death_inc_line", ValueType := "Deaths (28 d)"]
data[, ValueType := factor(ValueType, levels = c("Deaths (28 d)", "Admissions", "Hospital beds", "ICU beds", "PCR positives", "Seropositives"))]

# Make plot
theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))
# pop1 = c("East of England", "London", "Midlands")
# pop2 = c("North East and Yorkshire", "North West", "South East", "South West")
# 
# ggplot(output[d > "2020-03-01" & Geography %in% pop1]) +
#     geom_ribbon(aes(x = d, ymin = `Quantile 0.025` / 1000, ymax = `Quantile 0.975` / 1000, fill = ValueType), alpha = 0.5) +
#     geom_line(aes(x = d, y = Value / 1000, colour = ValueType)) +
#     geom_point(data = data[Geography %in% pop1], aes(x = d, y = y / 1000), size = 0.01, shape = 20) +
#     geom_linerange(data = data[Geography %in% pop1], aes(x = d, ymin = ymin / 1000, ymax = ymax / 1000), size = 0.2) +
#     geom_linerange(data = data[Geography %in% pop1], aes(xmin = dmin, xmax = dmax, y = y / 1000), size = 0.2) +
#     facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
#     theme(legend.position = "none", strip.placement = "outside") +
#     labs(x = "Date", y = "Value (thousands)")
# 
# ggsave("./figures/fits1.pdf", width = 18, height = 16, units = "cm", useDingbats = F)
# ggsave("./figures/fits1.png", width = 18, height = 16, units = "cm")
# 
# ggplot(output[d > "2020-03-01" & Geography %in% pop2]) +
#     geom_ribbon(aes(x = d, ymin = `Quantile 0.025` / 1000, ymax = `Quantile 0.975` / 1000, fill = ValueType), alpha = 0.5) +
#     geom_line(aes(x = d, y = Value / 1000, colour = ValueType)) +
#     geom_point(data = data[Geography %in% pop2], aes(x = d, y = y / 1000), size = 0.01, shape = 20) +
#     geom_linerange(data = data[Geography %in% pop2], aes(x = d, ymin = ymin / 1000, ymax = ymax / 1000), size = 0.2) +
#     geom_linerange(data = data[Geography %in% pop2], aes(xmin = dmin, xmax = dmax, y = y / 1000), size = 0.2) +
#     facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
#     theme(legend.position = "none", strip.placement = "outside") +
#     labs(x = "Date", y = "Value (thousands)")
# 
# ggsave("./figures/fits2.pdf", width = 18, height = 16, units = "cm", useDingbats = F)
# ggsave("./figures/fits2.png", width = 18, height = 16, units = "cm")

linetypes = c("Deaths (28 d)", "Admissions", "Hospital beds", "ICU beds")

ggplot(output[d > "2020-03-01"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.025` / 1000, ymax = `Quantile 0.975` / 1000, fill = ValueType), alpha = 0.5) +
    geom_line(aes(x = d, y = Value / 1000, colour = ValueType)) +
    geom_line(data = data[ValueType %in% linetypes], aes(x = d, y = y / 1000), size = 0.2) +
    geom_point(data = data[!ValueType %in% linetypes], aes(x = d, y = y / 1000), size = 0.01, shape = 20) +
    geom_linerange(data = data, aes(x = d, ymin = ymin / 1000, ymax = ymax / 1000), size = 0.2) +
    geom_linerange(data = data, aes(xmin = dmin, xmax = dmax, y = y / 1000), size = 0.2) +
    facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
    theme_cowplot(font_size = 6) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
    labs(x = "Date", y = "Value (thousands)")

ggsave("./figures/fits.pdf", width = 18, height = 9, units = "cm", useDingbats = F)
ggsave("./figures/fits.png", width = 18, height = 9, units = "cm")

# Zoomed in version of infection prevalence
ggplot(output[d > "2020-04-15" & ValueType == "PCR positives"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.025` / 1000, ymax = `Quantile 0.975` / 1000, fill = ValueType), alpha = 0.5) +
    geom_line(aes(x = d, y = Value / 1000, colour = ValueType)) +
    geom_line(data = data[ValueType %in% linetypes & ValueType == "PCR positives"], aes(x = d, y = y / 1000), size = 0.2) +
    geom_point(data = data[!ValueType %in% linetypes & ValueType == "PCR positives"], aes(x = d, y = y / 1000), size = 0.01, shape = 20) +
    geom_linerange(data = data[ValueType == "PCR positives"], aes(x = d, ymin = ymin / 1000, ymax = ymax / 1000), size = 0.2) +
    geom_linerange(data = data[ValueType == "PCR positives"], aes(xmin = dmin, xmax = dmax, y = y / 1000), size = 0.2) +
    facet_wrap(~Geography, scales = "free", nrow = 2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    theme_cowplot(font_size = 6) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
    labs(x = "Date", y = "Value (thousands)")

ggsave("./figures/prevalence.pdf", width = 18, height = 8, units = "cm", useDingbats = F)
ggsave("./figures/prevalence.png", width = 18, height = 8, units = "cm")




#
# POSTERIOR DISTRIBUTIONS
#

for (i in seq_along(posteriorsI)) {
    if (!is.null(posteriorsI[[i]])) {
        posteriorsI[[i]][, population := nhs_regions[i]]
    }
}
posterior = rbindlist(posteriorsI, fill = TRUE)
posteriorm = melt(posterior[!population %in% c("England", "United Kingdom", "Wales", "Scotland", "Northern Ireland")], 
    id.vars = c("trial", "lp", "chain", "ll", "population"))
posteriorm = posteriorm[!is.na(value)]
# make sure it is in order...

# Add final CFR to posteriorm -- cfr_rel2 is multiplied by cfr_rel to get final CFR.
cfr_final = posteriorm[variable == "cfr_rel2"]
cfr_final$value = cfr_final$value * posteriorm[variable == "cfr_rel", value]
cfr_final$variable = "cfr_final"
posteriorm = rbind(posteriorm, cfr_final)


# Different types of posterior plot
post_extract = function(p, what)
{
    if (is.character(what)) {
        p[variable == what, value]
    } else {
        what
    }
}

# Return first character value
pickname = function(...) {
    args = list(...)
    for (i in seq_along(args)) {
        if (class(args[[i]]) == "character") {
            return (args[[i]])
        }
    }
    stop("There's been a bad error.")
}

# time on x axis, value on y axis
post_tscatter = function(p, x, y, names)
{
    plots = NULL
    for (i in seq_along(x)) {
        pp = data.table(
                x = ymd("2020-01-01") + post_extract(p, x[[i]]),
                y = post_extract(p, y[[i]]),
                name = names[i],
                population = p[variable == pickname(x[[i]], y[[i]]), population])
        
        plots[[i]] = ggplot(pp) +
            geom_point(aes(x = x, y = y)) +
            facet_wrap(~population, nrow = 1) +
            labs(x = NULL, y = names[i]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    }
    
    return (plots)
}

# gamma distributions
dgamma2 = function(x, mu, shape)
{
    scale = mu / shape;
    dgamma(x, shape = shape, scale = scale)
}

post_gammadist = function(p, x0, x1, mean, shape, names)
{
    plots = NULL
    for (i in seq_along(mean)) {
        pp = data.table(
                x0 = x0[[i]],
                x1 = x1[[i]],
                mean = post_extract(p, mean[[i]]),
                shape = post_extract(p, shape[[i]]),
                name = names[i],
                population = p[variable == pickname(x0[[i]], x1[[i]], mean[[i]], shape[[i]]), population])
        
        pp = pp[, .SD[sample.int(1000)], by = population];
        pops = pp[, unique(population)]

        d = data.table(k = rep(1:1000, length(pops)), population = rep(pops, each = 1000))
        d[, as.character(0:100) := as.list(rep(0, 101))]
        for (j in 1:nrow(pp)) {
            d[j, as.character(0:100) := as.list(dgamma2(pp[j, x0] + ((0:100) / 100) * (pp[j, x1] - pp[j, x0]), pp[j, mean], pp[j, shape]))]
        }
        d = melt(d, id.vars = 1:2)
        d[, x := x0[[i]] + (as.numeric(as.character(variable)) / 100) * (x1[[i]] - x0[[i]])]
        d = d[, as.list(quantile(value, c(0.025, 0.5, 0.975))), by = .(population, variable, x)]
        plots[[i]] = ggplot(d) +
            geom_ribbon(aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "#eeeeee") +
            geom_line(aes(x = x, y = `50%`), size = 0.2) +
            facet_wrap(~population, nrow = 1) +
            labs(x = NULL, y = names[i]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    }
    
    return (plots)
}

# log-normal distributions
dlnorm2 = function(x, mu, cv)
{
    meanlog = log(mu / sqrt(1 + cv^2));
    sdlog = sqrt(log(1 + cv^2));
    return (dlnorm(x, meanlog, sdlog))
}

post_lnormdist = function(p, x0, x1, mu, cv, names)
{
    plots = NULL
    for (i in seq_along(mu)) {
        pp = data.table(
                x0 = x0[[i]],
                x1 = x1[[i]],
                mu = post_extract(p, mu[[i]]),
                cv = post_extract(p, cv[[i]]),
                name = names[i],
                population = p[variable == pickname(x0[[i]], x1[[i]], mu[[i]], cv[[i]]), population])
        
        pp = pp[, .SD[sample.int(1000)], by = population];
        pops = pp[, unique(population)]

        d = data.table(k = rep(1:1000, length(pops)), population = rep(pops, each = 1000))
        d[, as.character(0:100) := as.list(rep(0, 101))]
        for (j in 1:nrow(pp)) {
            d[j, as.character(0:100) := as.list(dlnorm2(pp[j, x0] + ((0:100) / 100) * (pp[j, x1] - pp[j, x0]), pp[j, mu], pp[j, cv]))]
        }
        d = melt(d, id.vars = 1:2)
        d[, x := x0[[i]] + (as.numeric(as.character(variable)) / 100) * (x1[[i]] - x0[[i]])]
        d = d[, as.list(quantile(value, c(0.025, 0.5, 0.975))), by = .(population, variable, x)]
        plots[[i]] = ggplot(d) +
            geom_ribbon(aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "#eeeeee") +
            geom_line(aes(x = x, y = `50%`), size = 0.2) +
            facet_wrap(~population, nrow = 1) +
            labs(x = NULL, y = names[i]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    }
    
    return (plots)
}

# asc curves
asc = function(x, y0, y1, s0, s1)
{
    xx = s0 + x * (s1 - s0)
    h0 = exp(s0) / (1 + exp(s0))
    h1 = exp(s1) / (1 + exp(s1))
    h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0)
    y0 + (y1 - y0) * h
}

post_asc = function(p, x0, x1, y0, y1, s0, s1, names)
{
    plots = NULL
    for (i in seq_along(x0)) {
        pp = data.table(
                x0 = x0[[i]],
                x1 = x1[[i]],
                y0 = post_extract(p, y0[[i]]),
                y1 = post_extract(p, y1[[i]]),
                s0 = post_extract(p, s0[[i]]),
                s1 = post_extract(p, s1[[i]]),
                name = names[i],
                population = p[variable == pickname(x0[[i]], x1[[i]], y0[[i]], y1[[i]], s0[[i]], s1[[i]]), population])
        
        pp = pp[, .SD[sample.int(1000)], by = population];
        pops = pp[, unique(population)]

        d = data.table(k = rep(1:1000, length(pops)), population = rep(pops, each = 1000))
        d[, as.character(0:100) := as.list(rep(0, 101))]
        for (j in 1:nrow(pp)) {
            d[j, as.character(0:100) := as.list(asc((0:100) / 100, pp[j, y0], pp[j, y1], pp[j, -s0], pp[j, s1]))]
        }
        d = melt(d, id.vars = 1:2)
        d[, x := x0[[i]] + (as.numeric(as.character(variable)) / 100) * (x1[[i]] - x0[[i]])]
        d = d[, as.list(quantile(value, c(0.025, 0.5, 0.975))), by = .(population, variable, x)]
        plots[[i]] = ggplot(d) +
            geom_ribbon(aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "#eeeeee") +
            geom_line(aes(x = x, y = `50%`), size = 0.2) +
            facet_wrap(~population, nrow = 1) +
            labs(x = NULL, y = names[i]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    }
    
    return (plots)
}

# Just a simple histogram
post_hist = function(p, x, names)
{
    plots = NULL
    for (i in seq_along(x)) {
        pp = data.table(
                x = post_extract(p, x[[i]]),
                name = names[i],
                population = p[variable == pickname(x[[i]]), population])
        
        plots[[i]] = ggplot(pp) +
            geom_histogram(aes(x)) +
            facet_wrap(~population, nrow = 1) +
            labs(x = NULL, y = names[i])
    }
    
    return (plots)
}

# 2 plots
theme_set(cowplot::theme_cowplot(font_size = 6) + theme(strip.background = element_blank()))
p1 = post_tscatter(posteriorm, 
    list("tS", "sep_when"), 
    list("u", "sep_boost"), 
    c("Start of transmission\nand susceptibility", "Late summer\nboost"))

# 4 plots
p2 = post_gammadist(posteriorm, 
    list(0, 0, 0, 0), 
    list(30, 30, 30, 30), 
    list("death_mean", "admission", "icu_admission", "waning"), 
    list("death_shape", 0.71, 1.91, 1), 
    c("Delay from infection\nto death", "Delay from infection\nto hospital admission", "Delay from infection\nto ICU admission", "Duration of\nseropositivity"))

# 3 plots
p3 = post_asc(posteriorm, 
    list(ymd("2020-01-01"), ymd("2020-01-01"), ymd("2020-01-01")), 
    list(ymd("2020-12-31"), ymd("2020-12-31"), ymd("2020-12-31")), 
    list(1, "cfr_rel", "detect1"), 
    list("contact_final", "cfr_final", "detect2"), 
    list("contact_s0", 2.9, "detect_s0"), 
    list("contact_s1", 7.8, "detect_s1"), 
    c("Contact\nmultiplier", "Relative\nfatality rate", "Time to test\nin hospital"))

# 5 plots
p4 = post_hist(posteriorm, 
    list("hosp_rlo", "icu_rlo", "concentration1", "concentration2", "concentration3"), 
    c("Relative log-odds\nof hospitalisation", "Relative log-odds\nof ICU admission", "Younger-age\ntransmission, July", "Younger-age\ntransmission, August", "Younger-age\ntransmission, September"))


cowplot::plot_grid(plotlist = c(p1, p4), ncol = 1, labels = letters, label_size = 10, align = "hv")
ggsave("./figures/posteriors1.pdf", width = 18, height = 24, units = "cm", useDingbats = FALSE)
ggsave("./figures/posteriors1.png", width = 18, height = 24, units = "cm")

cowplot::plot_grid(plotlist = c(p2, p3), ncol = 1, labels = letters, label_size = 10, align = "hv")
ggsave("./figures/posteriors2.pdf", width = 18, height = 24, units = "cm", useDingbats = FALSE)
ggsave("./figures/posteriors2.png", width = 18, height = 24, units = "cm")



# PLAIN HISTOGRAM STYLE POSTERIORS
theme_set(cowplot::theme_cowplot(font_size = 6) + theme(strip.background = element_blank()))
variables = posteriorm[, unique(variable)];
ggplot(posteriorm[variable %in% variables[1:12]]) +
    geom_histogram(aes(x = value, y = after_stat(density))) +
    facet_wrap(variable ~ population, scales = "free", ncol = 7) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("./figures/posteriors_hist1.pdf", width = 18, height = 24, units = "cm", useDingbats = FALSE)
ggsave("./figures/posteriors_hist1.png", width = 18, height = 24, units = "cm")

ggplot(posteriorm[variable %in% variables[13:26] & variable != "cfr_final"]) +
    geom_histogram(aes(x = value, y = after_stat(density))) +
    facet_wrap(variable ~ population, scales = "free", ncol = 7) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("./figures/posteriors_hist2.pdf", width = 18, height = 24, units = "cm", useDingbats = FALSE)
ggsave("./figures/posteriors_hist2.png", width = 18, height = 24, units = "cm")


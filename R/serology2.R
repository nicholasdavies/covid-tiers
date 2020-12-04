library(data.table)
library(ggplot2)
library(lubridate)
library(here)
library(cowplot)
library(readxl)
library(sn)
library(qs)
library(RCppMCMC)
library(glue)

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
england_regions = nhs_regions[c(1,3,4,5,6,9,10)]
all_data = qread("~/Dropbox/uk_covid_data/processed-data-2020-10-23.qs")
# for more recent data:
# all_data = qread("~/Dropbox/uk_covid_data/processed-data-2020-11-27.qs")
#ld = all_data[[1]]
sitreps = all_data[[2]]
#virus = all_data[[3]]
sero = all_data[[4]]


pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100
ggplot(sero[NHS.region %in% england_regions & Data.source != "NHSBT"]) + 
    geom_pointrange(aes(x = Start.date + (End.date - Start.date) / 2, ymin = pct(Lower.bound), ymax = pct(Upper.bound), y = pct(Central.estimate),
        colour = ifelse(Data.source == "Biobank", "Biobank", "ONS/REACT"))) +
    geom_linerange(aes(xmin = Start.date, xmax = End.date, y = pct(Central.estimate), colour = ifelse(Data.source == "Biobank", "Biobank", "ONS/REACT"))) +
    facet_wrap(~NHS.region) + labs(colour = "Source", x = "Date", y = "Estimate")

dsero = sero[NHS.region %in% england_regions]#### & Data.source != "NHSBT"]
dsitreps = sitreps[name %in% england_regions, .(date, n_admitted_diagnosed, name, pid)]

dsero[, pid := match(NHS.region, england_regions) - 1]
dsitreps[, pid := match(name, england_regions) - 1]

filler = data.table(date = rep(ymd("2020-03-01") + 0:as.numeric(ymd("2020-11-01")-ymd("2020-03-01")), length(england_regions)))
filler[, pid := rep(0:6, each = as.numeric(ymd("2020-11-01")-ymd("2020-03-01")) + 1)]
filler[, name := england_regions[pid + 1]]
dsitreps = merge(dsitreps, filler, by = c("date", "pid", "name"), all = T)
dsitreps[is.na(n_admitted_diagnosed), n_admitted_diagnosed := 0]

dsero[, d0 := as.numeric(Start.date - ymd("2020-03-01"))]
dsero[, d1 := as.numeric(End.date - ymd("2020-03-01"))]
dsitreps[, d := as.numeric(date - ymd("2020-03-01"))]

dsero = dsero[order(pid, d0)]
dsitreps = dsitreps[order(pid, d)]

popsize = fread(
"               Geography, population_size
          East of England,         6426982
                   London,         8908081
                 Midlands,         9497532
 North East and Yorkshire,         8152502
               North West,         6656813
               South East,         8192817
               South West,         5554296")

pids = unique(dsitreps[, .(pid, Geography = name)])
popsize = merge(popsize, pids)
popsize = popsize[order(pid)]

dsero
dsitreps

cpp_vec = function(x) 
{
    paste("{", paste(x, collapse = ", "), "}")
}


sero_model = make_model("serology",
    priors = c(
        IHR = "B 1 10",
        waning_inv = "N 200 100 T 0 1000"
    ),
    likelihood_code = glue(
        "using namespace std;",
        # functions
        "// normal log density",
        "auto norm = [&](double x, double mean, double sd)",
        "{",
        "    return -(x - mean) * (x - mean) / (2 * sd * sd) - log(sd) - 0.918938533204673;",
        "};",

        "// skew normal log density",
        "auto skewnorm = [&](double x, double xi, double omega, double alpha)",
        "{",
        "    double X = (x - xi) / omega;",
        "    return -log(omega) + norm(X, 0., 1.) + gsl_sf_log_erfc(-alpha * X / 1.414213562373095);",
        "};",

        # model
        "vector<vector<double>> x(7, vector<double>());", 
        "x[0] = ${cpp_vec(dsitreps[pid == 0, n_admitted_diagnosed])};",
        "x[1] = ${cpp_vec(dsitreps[pid == 1, n_admitted_diagnosed])};",
        "x[2] = ${cpp_vec(dsitreps[pid == 2, n_admitted_diagnosed])};",
        "x[3] = ${cpp_vec(dsitreps[pid == 3, n_admitted_diagnosed])};",
        "x[4] = ${cpp_vec(dsitreps[pid == 4, n_admitted_diagnosed])};",
        "x[5] = ${cpp_vec(dsitreps[pid == 5, n_admitted_diagnosed])};",
        "x[6] = ${cpp_vec(dsitreps[pid == 6, n_admitted_diagnosed])};",
        
        "for (unsigned int pid = 0; pid < x.size(); ++pid) {",
        "    vector<double>& y = x[pid];",
        "    for (unsigned int t = 0; t < y.size(); ++t) {",
        "        y[t] = (t > 0 ? y[t - 1] : 0) + (y[t] / IHR);",
        "        y[t] *= exp(-1.0 / waning_inv);",
        "    }",
        "}",
        
        # likelihood
        "vector<unsigned int> pid = ${cpp_vec(dsero$pid)};",
        "vector<unsigned int> d0 = ${cpp_vec(dsero$d0)};",
        "vector<unsigned int> d1 = ${cpp_vec(dsero$d1)};",
        "vector<double> xi = ${cpp_vec(dsero$xi)};",
        "vector<double> omega = ${cpp_vec(dsero$omega)};",
        "vector<double> alpha = ${cpp_vec(dsero$alpha)};",
        "vector<double> popsize = ${cpp_vec(popsize$population_size)};",
        
        "for (unsigned int i = 0; i < pid.size(); ++i) {",
        "    double sero_prev = 0;",
        "    for (unsigned int d = d0[i]; d <= d1[i]; ++d) {",
        "        sero_prev += x[pid[i]][d] / (d1[i] - d0[i] + 1.0);",
        "    }",
        "    ll += skewnorm(sero_prev / popsize[pid[i]], xi[i], omega[i], alpha[i]);",
        "}",
        .sep = "\n", .open = "${", .close = "}"
    )
)

results = RCppMCMC(sero_model, 10000, 10000, verbose = TRUE)
setDT(results)
ggplot(results) + geom_histogram(aes(waning_inv))
ggplot(results) + geom_histogram(aes(IHR))
results

sample_fit = function(dsitreps, waning_inv, IHR)
{
    dat = dsitreps
    dat[, sero := 0]
    for (pid in unique(dat$pid))
    {
        for (day in setdiff(unique(dat$d), 0))
        {
            previous = dat[pid == pid & d == day - 1, sero]
            dat[pid == pid & d == day, sero := (..previous + n_admitted_diagnosed / IHR) * exp(-1.0 / waning_inv)];
        }
    }
    dat
}

do_sampling = function(n_samp)
{
    cat("Sampling")
    set.seed(12345)
    dat = NULL
    for (s in 1:10) {
        row = sample(nrow(results), 1)
        dat0 = sample_fit(dsitreps, results[row, waning_inv], results[row, IHR])
        dat = rbind(dat, cbind(dat0, run = s))
        cat(".")
    }
    
    return (dat)
}

dat = do_sampling(10)

dsero2 = dsero[, .(name = NHS.region, d0, d1, y0 = pct(Lower.bound), y1 = pct(Upper.bound), ymid = pct(Central.estimate), source = Data.source)]
dsero2[, dmid := (d0 + d1) / 2]
dsero2[, date0 := d0 + ymd("2020-03-01")]
dsero2[, date1 := d1 + ymd("2020-03-01")]
dsero2[, datemid := dmid + ymd("2020-03-01")]

dat = merge(dat, popsize, by = "pid")
dat[, sero := sero / population_size]
dat = dat[, as.list(quantile(sero, c(0.025, 0.5, 0.975))), keyby = .(name, date)]

ggplot(dat) + 
    geom_ribbon(aes(x = date, ymin = `2.5%`, ymax = `97.5%`)) +
    geom_linerange(data = dsero2, aes(x = datemid, ymin = y0, ymax = y1, colour = source)) +
    geom_linerange(data = dsero2, aes(y = ymid, xmin = date0, xmax = date1, colour = source)) +
    facet_wrap(~name)










# experimental -- IHR for each place.

dsero[, omega := 0.005]
dsero[, alpha := 0]

sero_model = make_model("serology",
    priors = c(
        IHR0 = "B 1 10",
        IHR1 = "B 1 10",
        IHR2 = "B 1 10",
        IHR3 = "B 1 10",
        IHR4 = "B 1 10",
        IHR5 = "B 1 10",
        IHR6 = "B 1 10",
        waning_inv = "N 200 100 T 0 1000"
    ),
    likelihood_code = glue(
        "using namespace std;",
        # functions
        "// normal log density",
        "auto norm = [&](double x, double mean, double sd)",
        "{",
        "    return -(x - mean) * (x - mean) / (2 * sd * sd) - log(sd) - 0.918938533204673;",
        "};",

        "// skew normal log density",
        "auto skewnorm = [&](double x, double xi, double omega, double alpha)",
        "{",
        "    double X = (x - xi) / omega;",
        "    return -log(omega) + norm(X, 0., 1.) + gsl_sf_log_erfc(-alpha * X / 1.414213562373095);",
        "};",

        # model
        "vector<vector<double>> x(7, vector<double>());", 
        "x[0] = ${cpp_vec(dsitreps[pid == 0, n_admitted_diagnosed])};",
        "x[1] = ${cpp_vec(dsitreps[pid == 1, n_admitted_diagnosed])};",
        "x[2] = ${cpp_vec(dsitreps[pid == 2, n_admitted_diagnosed])};",
        "x[3] = ${cpp_vec(dsitreps[pid == 3, n_admitted_diagnosed])};",
        "x[4] = ${cpp_vec(dsitreps[pid == 4, n_admitted_diagnosed])};",
        "x[5] = ${cpp_vec(dsitreps[pid == 5, n_admitted_diagnosed])};",
        "x[6] = ${cpp_vec(dsitreps[pid == 6, n_admitted_diagnosed])};",
        
        "vector<double> IHRs { IHR0, IHR1, IHR2, IHR3, IHR4, IHR5, IHR6 };",
        
        "for (unsigned int pid = 0; pid < x.size(); ++pid) {",
        "    vector<double>& y = x[pid];",
        "    for (unsigned int t = 0; t < y.size(); ++t) {",
        "        y[t] = (t > 0 ? y[t - 1] : 0) + (y[t] / IHRs[pid]);",
        "        y[t] *= exp(-1.0 / waning_inv);",
        "    }",
        "}",
        
        # likelihood
        "vector<unsigned int> pid = ${cpp_vec(dsero$pid)};",
        "vector<unsigned int> d0 = ${cpp_vec(dsero$d0)};",
        "vector<unsigned int> d1 = ${cpp_vec(dsero$d1)};",
        "vector<double> xi = ${cpp_vec(dsero$xi)};",
        "vector<double> omega = ${cpp_vec(dsero$omega)};",
        "vector<double> alpha = ${cpp_vec(dsero$alpha)};",
        "vector<double> popsize = ${cpp_vec(popsize$population_size)};",
        
        "for (unsigned int i = 0; i < pid.size(); ++i) {",
        "    double sero_prev = 0;",
        "    for (unsigned int d = d0[i]; d <= d1[i]; ++d) {",
        "        sero_prev += x[pid[i]][d] / (d1[i] - d0[i] + 1.0);",
        "    }",
        "    ll += skewnorm(sero_prev / popsize[pid[i]], xi[i], omega[i], alpha[i]);",
        "}",
        .sep = "\n", .open = "${", .close = "}"
    )
)

results = RCppMCMC(sero_model, 10000, 10000, verbose = TRUE)
setDT(results)
results2 = melt(results, id.vars = 1:4)
results2[variable %like% "IHR", region := england_regions[as.numeric(str_remove(variable, "IHR")) + 1]]
ggplot(results2[variable %like% "IHR"]) + geom_density(aes(value, colour = region))
#ggplot(results2[variable %like% "waning"]) + geom_density(aes(value))
ggplot(results2) + geom_line(aes(x = trial, y = value, colour = chain, group = chain)) + facet_wrap(~variable, scales = "free")

sample_fit = function(dsitreps, waning_inv, IHRs)
{
    dat = dsitreps
    dat[, sero := 0]
    for (pid in unique(dat$pid))
    {
        for (day in setdiff(unique(dat$d), 0))
        {
            previous = dat[pid == pid & d == day - 1, sero]
            dat[pid == pid & d == day, sero := (..previous + n_admitted_diagnosed / IHRs[pid + 1]) * exp(-1.0 / waning_inv)];
        }
    }
    dat
}

do_sampling = function(n_samp)
{
    cat("Sampling")
    set.seed(12345)
    dat = NULL
    for (s in 1:10) {
        row = sample(nrow(results), 1)
        dat0 = sample_fit(dsitreps, results[row, waning_inv], results[row, c(IHR0, IHR1, IHR2, IHR3, IHR4, IHR5, IHR6)])
        dat = rbind(dat, cbind(dat0, run = s))
        cat(".")
    }
    
    return (dat)
}

dat = do_sampling(10)

dsero2 = dsero[, .(name = NHS.region, d0, d1, y0 = pct(Lower.bound), y1 = pct(Upper.bound), ymid = pct(Central.estimate), source = Data.source)]
dsero2[, dmid := (d0 + d1) / 2]
dsero2[, date0 := d0 + ymd("2020-03-01")]
dsero2[, date1 := d1 + ymd("2020-03-01")]
dsero2[, datemid := dmid + ymd("2020-03-01")]

dat = merge(dat, popsize, by = "pid")
dat[, sero := sero / population_size]
dat = dat[, as.list(quantile(sero, c(0.025, 0.5, 0.975))), keyby = .(name, date)]

ggplot(dat) + 
    geom_ribbon(aes(x = date, ymin = `2.5%`, ymax = `97.5%`)) +
    geom_linerange(data = dsero2, aes(x = datemid, ymin = y0, ymax = y1, colour = source)) +
    geom_linerange(data = dsero2, aes(y = ymid, xmin = date0, xmax = date1, colour = source)) +
    facet_wrap(~name)


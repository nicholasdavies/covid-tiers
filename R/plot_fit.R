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
all_data = qread("~/Dropbox/uk_covid_data/processed-data-2020-10-23-fix.qs")
# for more recent data:
# all_data = qread("~/Dropbox/uk_covid_data/processed-data-2020-11-27-new.qs")
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]]
sero = all_data[[4]]

# SUMMARISING FUNCS
output_1var = function(summ, varname, outname, cyear, cmonth, cday, ymd_from, ymd_to, nbinom = TRUE)
{
    quants = c(0.025, 0.5, 0.975)

    # new code...
    rq = data.table(run = 1:summ[, max(run)], q = runif(summ[, max(run)]))
    summ = merge(summ, rq, by = "run")
    if (nbinom) {
        summ[, rv := qnbinom(q, size = 20, mu = get(varname))]
    } else {
        summ[, rv := qnorm(q, mean = get(varname), sd = 1e-12)]
    }
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

output_full = function(test, cyear, cmonth, cday, ymd_from, ymd_to)
{
    cat("Summarizing variables...\n");
    summ = test[, .(death_o = sum(death_o)), by = .(run, t, population)]
    summ3 = test[, .(icu_p = sum(icu_p)), by = .(run, t, population)]
    summ4 = test[, .(bed_p = sum(pmax(0, hosp_p - hosp_undetected_p))), by = .(run, t, population)]
    summ34 = test[, .(admissions = sum(hosp_undetected_o)), by = .(run, t, population)]
    summ_pcr_p = test[, .(pcr_p = sum(pcr_positive_p)), by = .(run, t, population)]
    summ_sero_p = test[, .(sero_p = sum(lfia_positive_p)), by = .(run, t, population)]
    summ_attackrate = test[, .(attackrate = sum(pcr_positive_i)), by = .(run, t, population)]
    summ_attackrate[, attackrate := cumsum(attackrate), by = .(run, population)]

    cat("Running quantiles...\n");
    w = rbind(
        output_1var(summ, "death_o", "type28_death_inc_line", cyear, cmonth, cday, ymd_from, ymd_to),
        output_1var(summ3, "icu_p", "icu_prev", cyear, cmonth, cday, ymd_from, ymd_to),
        output_1var(summ4, "bed_p", "hospital_prev", cyear, cmonth, cday, ymd_from, ymd_to),
        output_1var(summ34, "admissions", "hospital_inc", cyear, cmonth, cday, ymd_from, ymd_to),
        output_1var(summ_pcr_p, "pcr_p", "pcr_prev", cyear, cmonth, cday, ymd_from, ymd_to, FALSE),
        output_1var(summ_sero_p, "sero_p", "sero_prev", cyear, cmonth, cday, ymd_from, ymd_to, FALSE),
        output_1var(summ_attackrate, "attackrate", "attackrate", cyear, cmonth, cday, ymd_from, ymd_to, FALSE)
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
        sero[!Data.source %like% "NHSBT", .(ValueType = "sero_prev", Geography = NHS.region,
            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))]
    )
}



# Load fit
source("./R/latest_fit.R")

# Ensure no undefined parameters
for (i in seq_along(parametersI)) {
    if (!is.null(parametersI[[i]])) {
        parametersI[[i]]$pop[[1]]$dEa = parametersI$pop[[1]]$dE
        
        parametersI[[i]]$pop[[1]]$ev = rep(1, 16)
        parametersI[[i]]$pop[[1]]$ev2 = rep(1, 16)
        
        parametersI[[i]]$pop[[1]]$ei_v = rep(1, 16)
        parametersI[[i]]$pop[[1]]$ed_vi = rep(1, 16)

        parametersI[[i]]$pop[[1]]$ei_v2 = rep(1, 16)
        parametersI[[i]]$pop[[1]]$ed_vi2 = rep(1, 16)
        
        parametersI[[i]]$pop[[1]]$pi_r = rep(1, 16)
        parametersI[[i]]$pop[[1]]$pd_ri = rep(1, 16)
        
        parametersI[[i]]$pop[[1]]$wv = rep(0, 16)
        parametersI[[i]]$pop[[1]]$wv2 = rep(0, 16)
        
        parametersI[[i]]$pop[[1]]$v = rep(0, 16)
        parametersI[[i]]$pop[[1]]$v12 = rep(0, 16)
        parametersI[[i]]$pop[[1]]$v2 = rep(0, 16)
    }
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

source("./R/cpp_funcs.R")

# Recover dynamics


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
    test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[p]], 50, seed = 0);
    warning("50 fits only.")

    test = rbindlist(test)
    test = test[, population := p]
    dynamicsI[[p]] = test

    print(p)
}

lapply(posteriorsI, function(x) list(maxlp = max(x$lp), meanlp = mean(x$lp)))


#
# PLOT FITS.
#

# # Remove problematic virus entries
# virus = virus[omega > 1e-9]

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
# Use next line for plotting projections against data
# test = proj_ldE4o
output = output_full(test[!population %in% c("England", "United Kingdom", "Wales", "Scotland", "Northern Ireland")], 2020, 11, 1, "2020-01-01", "2020-10-23")
# output = output_full(test[!population %in% c("England", "United Kingdom", "Wales", "Scotland", "Northern Ireland")], 2020, 11, 1, "2020-01-01", "2020-11-30")
output[, d := make_date(`Year of Value`, `Month of Value`, `Day of Value`)]
output = merge(output, popsize, by = "Geography")

adj_output = function(output, val_type, div, pop = 0) {
    output[ValueType == val_type, Value := Value / (div + population_size * pop)]
    output[ValueType == val_type, `Quantile 0.025` := `Quantile 0.025` / (div + population_size * pop)]
    output[ValueType == val_type, `Quantile 0.5`   := `Quantile 0.5`   / (div + population_size * pop)]
    output[ValueType == val_type, `Quantile 0.975` := `Quantile 0.975` / (div + population_size * pop)]
}

adj_output(output, "hospital_inc", 1)
adj_output(output, "hospital_prev", 1)
adj_output(output, "icu_prev", 1)
adj_output(output, "pcr_prev", 0, 0.01)
adj_output(output, "sero_prev", 0, 0.01)
adj_output(output, "type28_death_inc_line", 1)
adj_output(output, "attackrate", 0, 0.01)

# Make data to output
data = make_data(ld, sitreps, virus, sero)[!Geography %in% c("England", "United Kingdom", "Wales", "Scotland", "Northern Ireland")]
data = merge(data, popsize, by = "Geography")

adj_data = function(data, val_type, div, pop = 0) {
    data[ValueType == val_type, ymin := ymin / (div + population_size * pop)]
    data[ValueType == val_type, y    := y    / (div + population_size * pop)]
    data[ValueType == val_type, ymax := ymax / (div + population_size * pop)]
}

adj_data(data, "hospital_inc", 1)
adj_data(data, "hospital_prev", 1)
adj_data(data, "icu_prev", 1)
adj_data(data, "pcr_prev", 0.01)
adj_data(data, "sero_prev", 0.01)
adj_data(data, "type28_death_inc_line", 1)

# Relabel outcomes
output[ValueType == "hospital_inc", ValueType := "Admissions"]
output[ValueType == "hospital_prev", ValueType := "Hosp beds occupied"]
output[ValueType == "icu_prev", ValueType := "ICU beds occupied"]
output[ValueType == "pcr_prev", ValueType := "PCR positivity, %"]
output[ValueType == "sero_prev", ValueType := "Seropositivity, %"]
output[ValueType == "type28_death_inc_line", ValueType := "Deaths"]
output[ValueType == "attackrate", ValueType := "Attack rate, %"]
output[, ValueType := factor(ValueType, levels = c("Deaths", "Admissions", "Hosp beds occupied", 
    "ICU beds occupied", "PCR positivity, %", "Seropositivity, %", "Attack rate, %"))]

data[ValueType == "hospital_inc", ValueType := "Admissions"]
data[ValueType == "hospital_prev", ValueType := "Hosp beds occupied"]
data[ValueType == "icu_prev", ValueType := "ICU beds occupied"]
data[ValueType == "pcr_prev", ValueType := "PCR positivity, %"]
data[ValueType == "sero_prev", ValueType := "Seropositivity, %"]
data[ValueType == "type28_death_inc_line", ValueType := "Deaths"]
data[, ValueType := factor(ValueType, levels = c("Deaths", "Admissions", "Hosp beds occupied", 
    "ICU beds occupied", "PCR positivity, %", "Seropositivity, %", "Attack rate, %"))]

# Make plot
theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

linetypes = c("Deaths", "Admissions", "Hosp beds occupied", "ICU beds occupied")

ggplot(output[d > "2020-09-01"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.025`, ymax = `Quantile 0.975`, fill = ValueType), alpha = 0.5) +
    geom_line(aes(x = d, y = Value, colour = ValueType)) +
    geom_line(data = data[ValueType %in% linetypes & d > "2020-09-01"], aes(x = d, y = y), size = 0.2) +
    geom_point(data = data[!ValueType %in% linetypes & d > "2020-09-01"], aes(x = d, y = y), size = 0.01, shape = 20) +
    geom_linerange(data = data[d > "2020-09-01"], aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
    geom_linerange(data = data[d > "2020-09-01"], aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
    geom_vline(aes(xintercept = ymd("2020-10-14")), size = 0.2, linetype = "42") +
    facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
    theme_cowplot(font_size = 6) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
    labs(x = NULL, y = NULL)

ggsave("./figures/forward_fits_E.pdf", width = 18, height = 14, units = "cm", useDingbats = F)
ggsave("./figures/forward_fits_E.png", width = 18, height = 14, units = "cm")


ggplot(output[d > "2020-03-01"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.025`, ymax = `Quantile 0.975`, fill = ValueType), alpha = 0.5) +
    geom_line(aes(x = d, y = Value, colour = ValueType)) +
    geom_line(data = data[ValueType %in% linetypes], aes(x = d, y = y), size = 0.2) +
    geom_point(data = data[!ValueType %in% linetypes], aes(x = d, y = y), size = 0.01, shape = 20) +
    geom_linerange(data = data, aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
    geom_linerange(data = data, aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
    facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
    theme_cowplot(font_size = 6) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
    labs(x = NULL, y = NULL)

ggsave("./figures/fits.pdf", width = 18, height = 14, units = "cm", useDingbats = F)
ggsave("./figures/fits.png", width = 18, height = 14, units = "cm")

ggplot(output[d > "2020-04-15" & ValueType == "PCR positivity, %"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.025`, ymax = `Quantile 0.975`), fill = "#6F9AF8", alpha = 0.5) +
    geom_line(aes(x = d, y = Value), colour = "#6F9AF8") +
    geom_line(data = data[ValueType %in% linetypes & ValueType == "PCR positivity, %"], aes(x = d, y = y), size = 0.2) +
    geom_point(data = data[!ValueType %in% linetypes & ValueType == "PCR positivity, %"], aes(x = d, y = y), size = 0.01, shape = 20) +
    geom_linerange(data = data[ValueType == "PCR positivity, %"], aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
    geom_linerange(data = data[ValueType == "PCR positivity, %"], aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
    facet_wrap(~Geography, scales = "free", nrow = 2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    theme_cowplot(font_size = 6) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
    labs(x = NULL, y = "PCR positivity, %")

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

ggplot(posteriorm) + geom_density(aes(x = value, colour = population)) + facet_wrap(~variable, scales = "free")

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
            geom_bin2d(aes(x = x, y = y), bins = 60) +
            facet_wrap(~population, nrow = 1) +
            labs(x = NULL, y = names[i]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none") +
            scale_fill_gradient(low = "#ffffff", high = "#000000", na.value = "#ffffff")
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

# 3 plots
p2 = post_gammadist(posteriorm, 
    list(0, 0, 0), 
    list(30, 30, 30), 
    list("death_mean", "admission", "icu_admission"), 
    list("death_shape", 0.71, 1.91), 
    c("Delay from infection\nto death", "Delay from infection\nto hospital admission", "Delay from infection\nto ICU admission"))

# 3 plots
p3 = post_asc(posteriorm, 
    list(ymd("2020-01-01"), ymd("2020-01-01"), ymd("2020-01-01")), 
    list(ymd("2020-12-31"), ymd("2020-12-31"), ymd("2020-12-31")), 
    list(1, "cfr_rel", "detect1"), 
    list("contact_final", "cfr_final", "detect2"), 
    list("contact_s0", 2.9, "detect_s0"), 
    list("contact_s1", 7.8, "detect_s1"), 
    c("Contact\nmultiplier", "Relative\nfatality rate", "Time to test\nin hospital"))

# 8 plots
p4 = post_hist(posteriorm, 
    list("nonicu_los", "icu_los", "hosp_rlo", "icu_rlo", 
        "concentration1", "concentration2", "concentration3", "waning"), 
    c("Length of\nhospital stay", "Length of\nICU stay", "Relative log-odds\nof hospitalisation", "Relative log-odds\nof ICU admission", 
        "Younger-age\ntransmission, July", "Younger-age\ntransmission, August", "Younger-age\ntransmission, September", "Duration of\nseropositivity"))


cowplot::plot_grid(plotlist = c(p1, p2, p3), ncol = 1, labels = letters, label_size = 10, align = "hv")
ggsave("./figures/posteriors1.pdf", width = 18, height = 24, units = "cm", useDingbats = FALSE)
ggsave("./figures/posteriors1.png", width = 18, height = 24, units = "cm")

cowplot::plot_grid(plotlist = p4, ncol = 1, labels = letters, label_size = 10, align = "hv")
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




# PREVIOUS MTPs
theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))
m = rbind(
    fread("~/Dropbox/uk_covid_data/forecast_2020-10-27B.csv"),
    fread("~/Dropbox/uk_covid_data/forecast_2020-10-20.csv"),
    fread("~/Dropbox/uk_covid_data/forecast_2020-10-14.csv"),
    fread("~/Dropbox/uk_covid_data/forecast_2020-10-05b.csv"),
    fread("~/Dropbox/uk_covid_data/forecast_2020-09-29b.csv"),
    fread("~/Dropbox/uk_covid_data/LSHTM_Transmission_forecast_2020-11-02B.csv"),
    use.names = FALSE
)[AgeBand == "All"]
m[, DateStamp := make_date(`Creation Year`, `Creation Month`, `Creation Day`)]
m[, Date := make_date(`Year of Value`, `Month of Value`, `Day of Value`)]
m = m[ValueType != "infections_inc"]

m[ValueType == "hospital_inc", ValueType := "Admissions"]
m[ValueType == "type28_death_inc_line", ValueType := "Deaths"]
m[ValueType == "hospital_prev", ValueType := "Hosp beds occupied"]
m[ValueType == "icu_prev", ValueType := "ICU beds occupied"]

m[, min_date := min(Date), by = DateStamp]
m = m[Date <= min_date + 7 * 6]

data2 = fread("~/Dropbox/uk_covid_data/data-2020-11-27.csv")
data2 = data2[location %in% c("East of England", "London", "Midlands", "North East and Yorkshire", "North West", "South East", "South West")]
names(data2)[2] = "Geography"
names(data2)[3] = "ValueType"
data2[ValueType == "death_inc_line", ValueType := "Deaths"]
data2[ValueType == "hospital_inc", ValueType := "Admissions"]
data2[ValueType == "hospital_prev", ValueType := "Hosp beds occupied"]
data2[ValueType == "icu_prev", ValueType := "ICU beds occupied"]

data2
data_england = data2[, .(value = sum(value)), by = .(ValueType, date)]
data_england

ggplot(m[Geography == "England"]) + 
    geom_ribbon(aes(x = Date, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, group = DateStamp, fill = as.factor(DateStamp)), alpha = 0.5) + 
    geom_line(data = data_england[date > "2020-09-15"], aes(x = date, y = value)) +
    facet_wrap(~ValueType, scales = "free") +
    theme(legend.position = c(0.02, 0.85)) +
    labs(x = "Date", y = NULL, fill = "Date of projection") +
    theme(legend.title = element_text(size = 8)) +
    scale_color_viridis_d(aesthetics = c("colour", "fill"), option = "D")

ggsave("./figures/MTPs.pdf", width = 18, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./figures/MTPs.png", width = 18, height = 12, units = "cm")



#
# Comparison of mobility/tiers to predicted
#

library(data.table)
library(ggplot2)
library(zoo)
library(mgcv)
library(lubridate)
library(stringr)
library(cowplot)
library(qs)
library(ogwrangler)

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

# Load Google Mobility data and interface with ogwrangler
gm = qread("./data/google_mobility_uk.qs");
CreateCache();
gm[, name := ifelse(sub_region_2 != "", paste0(sub_region_2, ", ", sub_region_1), sub_region_1)];
gm = gm[name != ""];
gm_match = data.table(code = ogcode("*", "gmcty"));
gm_match[, name := ogwhat(code, "name")];
gm_match[, nhser := ogwhat(code, "nhser")];
gm_match = gm_match[nhser %like% "^E"];
gm_match[, nhser_name := ogwhat(nhser)];
gm_match[, pop := ogwhat(code, "pop2019")];

gm = merge(gm, gm_match, by = "name")
names(gm) = str_replace_all(names(gm), "_percent_change_from_baseline", "")
gm = melt(gm, id.vars = c("date", "code", "nhser_name", "pop"), measure.vars = c("retail_and_recreation", "grocery_and_pharmacy", "parks", "transit_stations", "workplaces", "residential"))

gm = gm[, .(value = weighted.mean(value, pop, na.rm = T)), keyby = .(date, nhser_name, variable)]
gm[, rollval := rollmean(value, 7, fill = NA), by = .(nhser_name, variable)]

# Compare predicted to actual tiering
tiers = proj_ldW4o[, .(tier = obs0[2]), by = .(run, t, population)]

tiers = tiers[, .(`0` = sum(tier == 0) / .N,
    `1` = sum(tier == 1) / .N,
    `2` = sum(tier == 2) / .N,
    `3` = sum(tier == 3) / .N), by = .(t, population)]

tiers = melt(tiers, id.vars = 1:2, variable.name = "tier", value.name = "frac")
tiers[, tier := as.numeric(as.character(tier))]
tiers[, date := ymd("2020-01-01") + t]

# Empirical
tiers2 = fread("./data/tiers2.txt")
lads = ogcode("*", "lad")
lads = lads[lads %like% "^E"]
measures = data.table(lad = rep(lads, each = 35), date = rep(ymd("2020-10-01") + 0:34, length(lads)))
measures[, pop := ogwhat(lad, "pop2019")]
measures[, population := ogwhat(ogwhat(lad, "nhser"))]
measures = merge(measures, tiers2[, .(lad = ladcd, date, tier)], by = c("lad", "date"), all = TRUE)
measures[date == "2020-10-01", tier := 0]
measures[date == "2020-10-14" & is.na(tier), tier := 1]
measures[, tier := nafill(tier, "locf"), by = .(lad)]

measures = measures[, .(
    `0` = weighted.mean(tier == 0, pop),
    `1` = weighted.mean(tier == 1, pop),
    `2` = weighted.mean(tier == 2, pop),
    `3` = weighted.mean(tier == 3, pop)), by = .(date, population)]

measures = melt(measures, id.vars = 1:2, variable.name = "tier", value.name = "frac")
measures[, tier := as.numeric(as.character(tier))]

ggplot() +
    geom_ribbon(data = measures[date > "2020-10-01" & date < "2020-11-05"],
        aes(x = date, ymin = tier - frac/3, ymax = tier + frac/3, group = tier, fill = "Actual"), colour = NA) +
    geom_line(data = tiers[date > "2020-10-01" & date < "2020-11-05" & frac != 0], 
        aes(x = date, y = tier - frac/3, group = tier, colour = "Projected")) +
    geom_line(data = tiers[date > "2020-10-01" & date < "2020-11-05" & frac != 0], 
        aes(x = date, y = tier + frac/3, group = tier, colour = "Projected")) +
    scale_fill_manual(values = "grey") +
    scale_colour_manual(values = "red") +
    facet_wrap(~population) +
    labs(x = "Date", y = "Tier", fill = NULL, colour = NULL) +
    theme(legend.position = "bottom")



# Comparing mobility indices predicted to actual
get_mobility = function(proj_ldW4o, parametersI)
{
    control = proj_ldW4o[, .(tier = obs0[2], ld = obs0[3]), by = .(run, t, population)]
    control[ld == 1, tier := 0]
    control[, sched := tier + 1]
    control[sched == 2, sched := 1]
    control[, pid := match(population, nhs_regions)]
    
    lookup = NULL
    for (p in unique(control$pid)) {
        times = parametersI[[p]]$schedule[[1]]$times
        
        s1 = matrix(unlist(parametersI[[p]]$schedule[[1]]$values), ncol = 8, byrow = T)
        s3 = matrix(unlist(parametersI[[p]]$schedule[[3]]$values), ncol = 8, byrow = T)
        s4 = matrix(unlist(parametersI[[p]]$schedule[[4]]$values), ncol = 8, byrow = T)
        lookup = rbind(lookup,
            cbind(t = times, data.table(s1), sched = 1, pid = p),
            cbind(t = times, data.table(s3), sched = 3, pid = p),
            cbind(t = times, data.table(s4), sched = 4, pid = p)
        )
    }
    
    control = merge(control, lookup, by = c("pid", "t", "sched"))
    names(control)[9:12] = c("workplaces", "grocery_and_pharmacy", "retail_and_recreation", "transit_stations")
    control[, V1 := NULL]
    control[, V6 := NULL]
    control[, V7 := NULL]
    control[, V8 := NULL]
    control[, date := ymd(parametersI[[1]]$date0) + t]
    control = melt(control, id.vars = c("population", "run", "date"), measure.vars = c("workplaces", "grocery_and_pharmacy", "retail_and_recreation", "transit_stations"))
    names(control)[1] = "nhser_name"
    control[, value := value * 100 - 100]
    
    gm = fread("./data/gm2.txt")
    names(gm) = c("nhser_name", "variable", "date", "t", "rollval")
    gm[, rollval := rollval * 100 - 100]
    gm = gm[!nhser_name %in% c("England", "Northern Ireland", "Scotland", "Wales", "United Kingdom")]
    
    ggplot() + 
        geom_line(data = control[run == 1 & date >= "2020-09-01" & date <= "2020-12-01"], 
            aes(x = date, y = value, colour = variable), linetype = "41") +
        geom_line(data = gm[!variable %in% c("residential", "parks") & date >= "2020-09-01" & date <= "2020-12-01"], 
            aes(x = date, y = rollval, group = variable, colour = variable), size = 0.2) +
        facet_wrap(~nhser_name)
}

parametersI = qread("./figures/params_ldW4o.qs");
get_mobility(proj_ldW4o, parametersI)

parametersI = qread("./figures/params_ldE4o.qs");
get_mobility(proj_ldE4o, parametersI)

nhs_regions
parametersI[[1]]$schedule[[1]]$times



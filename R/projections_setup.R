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
# source("~/Documents/covidm_MTPs/distribution_fit.R");
# source("~/Documents/covidm_MTPs/spim_output.R");


#
# DATA
#

nhs_regions = popUK[, unique(name)]




#
# Individual fits
#

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
    )

}



# LOAD FITS
#saved = qread("~/Desktop/uk-fit-nick-2020-11-02.qs")
# posteriorsI = saved[[1]]
# dynamicsI = saved[[2]]
# parametersI = saved[[3]]
# rm(saved)

saved = qread("./fit_pp15.qs")
posteriorsI = saved[[1]]
parametersI = saved[[2]]
rm(saved)

# Determine population sizes and set new vaccine parameters
popsize = NULL
for (i in seq_along(parametersI)) {
    if (!is.null(parametersI[[i]])) {
        popsize = rbind(popsize,
            data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = sum(parametersI[[i]]$pop[[1]]$size))
        )
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
popsize = rbind(data.table(Geography = "England", population_size = sum(popsize$population_size)), popsize)
setkey(popsize, Geography)



# 2. Bringing in of interventions

# Extract first column from each of tables_list with scenario names
england_only = function(tables_list, names)
{
    t = tables_list[[1]][, 1]
    for (i in seq_along(tables_list)) {
        t = cbind(t, tables_list[[i]][, 2]);
        names(t)[i + 1] = names[i];    
    }
    
    return (t)
}

arrange_projection = function(proj, cumulative_deaths = FALSE, from_date = NULL, england = FALSE)
{
    ft = 0;
    if (!is.null(from_date)) {
        ft = as.numeric(ymd(from_date) - ymd("2020-01-01"))
    }
    
    if (england == FALSE) {
        w = proj[t >= ft, .(deaths = sum(death_o), admissions = sum(hosp_undetected_o),
            beds = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p), 
            Rt = obs0[1], tier = obs0[2], cb = obs0[3]), keyby = .(run, t, population)]
    } else {
        w = proj[t >= ft, .(population = "England", deaths = sum(death_o), admissions = sum(hosp_undetected_o),
            beds = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p), 
            Rt = mean(obs0[seq(1, .N, by = 16)]), tier = mean(obs0[seq(2, .N, by = 16)]), cb = mean(obs0[seq(3, .N, by = 16)])), keyby = .(run, t)]
    }
    
    if (cumulative_deaths) {
        w[, cum_deaths := cumsum(deaths), by = .(run, population)]
    }
    w = melt(w, id.vars = 1:3)
    w = w[, as.list(quantile(value, c(0.025, 0.5, 0.975))), by = .(population, variable, t)]
    w[variable == "deaths", variable := "Deaths"]
    w[variable == "cum_deaths", variable := "Cumulative deaths"]
    w[variable == "admissions", variable := "Admissions"]
    w[variable == "beds", variable := "Hospital beds"]
    w[variable == "icu", variable := "ICU beds"]
}

plot_projection_lt = function(proj, proj_compare, from_date, proj_compare2 = NULL)
{
    p = arrange_projection(proj)
    if (!is.null(proj_compare)) pc = arrange_projection(proj_compare);
    if (!is.null(proj_compare2)) pc2 = arrange_projection(proj_compare2);

    plot = ggplot(p[!variable %in% c("tier", "cb") & ymd("2020-01-01") + t >= from_date]);
    
    cbs = p[variable == "cb", .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population)][, 
        .(tmin = min(t), tmax = max(t)), by = population]
    if (nrow(cbs) > 0) {
        plot = plot +
            geom_rect(data = cbs, aes(xmin = tmin, xmax = tmax, ymin = -Inf, ymax = Inf), fill = "black", alpha = 0.1)
    }
    
    if (!is.null(proj_compare)) {
        plot = plot + geom_ribbon(data = pc[!variable %in% c("tier", "cb") & ymd("2020-01-01") + t >= from_date],
            aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, colour = variable), 
            linetype = "51", fill = NA, size = 0.2, alpha = 0.5)
    }

    if (!is.null(proj_compare2)) {
        plot = plot + geom_ribbon(data = pc2[!variable %in% c("tier", "cb") & ymd("2020-01-01") + t >= from_date],
            aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, colour = variable), 
            linetype = "dotted", fill = NA, size = 0.2, alpha = 0.5)
    }
    
    rline = data.table(variable = factor("Rt", levels(p$variable)), y = 1)
    
    plot +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = variable), alpha = 0.5) +
        geom_hline(data = rline, aes(yintercept = y), size = 0.3) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = variable)) +
        facet_grid(variable ~ population, switch = "y", scales = "free") +
        cowplot::theme_cowplot(font_size = 8) + 
        theme(strip.background = element_blank(), strip.placement = "outside", legend.position = "none",
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.4)) +
        labs(x = NULL, y = NULL) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        ylim(0, NA)
}

plot_projection = function(proj_list, proj_names, from_date, pal = "Accent")
{
    p = NULL
    for (i in seq_along(proj_list)) {
        pp = arrange_projection(proj_list[[i]]);
        pp[, name := proj_names[[i]]];
        p = rbind(p, pp);
    }
    p[, name := factor(name, unique(name))]

    plot = ggplot(p[!variable %in% c("tier", "cb") & ymd("2020-01-01") + t >= from_date]);
    
    cbs = p[variable == "cb", .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population)][, 
        .(tmin = min(t), tmax = max(t)), by = population]
    if (nrow(cbs) > 0) {
        plot = plot +
            geom_rect(data = cbs, aes(xmin = tmin, xmax = tmax, ymin = -Inf, ymax = Inf), fill = "black", alpha = 0.1)
    }
    
    rline = data.table(variable = factor("Rt", levels(p$variable)), y = 1)
    
    plot +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = name), alpha = 0.5) +
        geom_hline(data = rline, aes(yintercept = y), size = 0.3) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = name)) +
        facet_grid(variable ~ population, switch = "y", scales = "free") +
        cowplot::theme_cowplot(font_size = 8) + 
        theme(strip.background = element_blank(), strip.placement = "outside", 
            legend.position = "bottom",
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.4)) +
        labs(x = NULL, y = NULL, colour = NULL, fill = NULL) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal) +
        ylim(0, NA)
}

plot_cum_deaths = function(proj_list, proj_names, from_date, ystart = 26, ydiff = -1.2, titl = "Cumulative deaths")
{
    p = NULL
    for (i in seq_along(proj_list)) {
        w = arrange_projection(proj_list[[i]], cumulative_deaths = TRUE, from_date = from_date, england = TRUE)
        p = rbind(p,
            cbind(w[variable == "Cumulative deaths" | variable == "cb"], name = proj_names[i])
        )
    }
    
    p[, name := factor(name, levels = unique(name))]
    
    cbs = p[variable == "cb", .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population, name)][, 
        .(tmin = min(t), tmax = max(t)), by = .(population, name)]
    
    plot = ggplot(p[variable != "cb"]) +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%` / 1000, ymax = `97.5%` / 1000, fill = name), alpha = 0.5) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%` / 1000, colour = name), size = 0.2) +
        facet_wrap(~population, scales = "free") +
        theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = c(0.02, 0.8),
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2)) +
        labs(x = NULL, y = "Cumulative deaths\n(thousands)", colour = NULL, fill = NULL, title = titl) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")
    
    if (nrow(cbs) > 0) {
        yy = ystart
        for (nm in cbs[, unique(name)]) {
            plot = plot +
                geom_linerange(data = cbs[name == nm], aes(xmin = tmin, xmax = tmax, colour = name), y = yy, show.legend = FALSE, linetype = "32", size = 0.2) +
                geom_linerange(data = cbs[name == nm], aes(x = tmin, colour = name), ymin = yy - (ydiff * 0.2), ymax = yy + (ydiff * 0.2), show.legend = FALSE, size = 0.2) +
                geom_linerange(data = cbs[name == nm], aes(x = tmax, colour = name), ymin = yy - (ydiff * 0.2), ymax = yy + (ydiff * 0.2), show.legend = FALSE, size = 0.2)
            yy = yy + ydiff
        }
    }

    return (plot)
}

plot_icu = function(proj_list, proj_names, from_date)
{
    p = NULL
    for (i in seq_along(proj_list)) {
        w = arrange_projection(proj_list[[i]], cumulative_deaths = TRUE, from_date = from_date, england = TRUE)
        p = rbind(p,
            cbind(w[variable == "ICU beds"], name = proj_names[i])
        )
    }
    
    p[, name := factor(name, levels = unique(name))]
    
    plot = ggplot(p) +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = name), alpha = 0.5) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = name), size = 0.2) +
        facet_wrap(~population, scales = "free") +
        theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "none",
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2)) +
        labs(x = NULL, y = "ICU beds occupied", colour = NULL, fill = NULL) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")

    return (plot)
}


summarize_projection = function(proj, from_date, popsize, scenario_output = NULL, to_date = NULL, wh = "after")
{
    if (!is.null(to_date)) {
        proj = rlang::duplicate(proj)
        proj = proj[ymd("2020-01-01") + t <= to_date]
    }
    
    geo = c("England", proj[, unique(population)])
    
    # Build summaries
    totals = proj[, .(cases = sum(cases), infections = sum(infections_i), 
        deaths = sum(death_o), admissions = sum(hosp_undetected_o)), 
        by = .(population, run, when = ifelse(ymd("2020-01-01") + t < from_date, "before", "after"))]
    totals = rbind(totals[, .(population = "England", cases = sum(cases), infections = sum(infections),
        deaths = sum(deaths), admissions = sum(admissions)), by = .(run, when)], totals, fill = TRUE);
    
    prev_peaks = proj[ymd("2020-01-01") + t < from_date,
        .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
        by = .(population, run, t)][,
        .(all = max(all), icu = max(icu)),
        by = .(population, run)][,
        .(peak_all = mean(all), peak_icu = mean(icu)), keyby = population]
    prev_peaks = rbind(prev_peaks[, .(population = "England", peak_all = sum(peak_all), peak_icu = sum(peak_icu))], prev_peaks, fill = TRUE);
    setkey(prev_peaks, population)
    
    if (wh == "before") {
        post_peaks = proj[ymd("2020-01-01") + t < from_date,
            .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
            by = .(population, run, t)][,
            .(all = max(all), icu = max(icu)),
            keyby = .(population, run)]
    } else {
        post_peaks = proj[ymd("2020-01-01") + t >= from_date,
            .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
            by = .(population, run, t)][,
            .(all = max(all), icu = max(icu)),
            keyby = .(population, run)]
    }
    post_peaks = rbind(post_peaks[, .(population = "England", all = sum(all), icu = sum(icu)), by = run], post_peaks, fill = TRUE);
    setkey(post_peaks, population)
    
    peaks = post_peaks[,
        .(all = all, rel_all = all / prev_peaks[population, peak_all],
          icu = icu, rel_icu = icu / prev_peaks[population, peak_icu]),
        keyby = .(population, run)]

    hi = proj[ymd("2020-01-01") + t >= from_date,
        .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
        keyby = .(population, run, t)][,
        .(all_hi = sum(all > prev_peaks[population, peak_all] * 0.5) / 7,
          icu_hi = sum(icu > prev_peaks[population, peak_icu] * 0.5) / 7), 
        keyby = .(population, run)]
    
    tier = proj[ymd("2020-01-01") + t >= from_date,
        .(tier = obs0[2], cb = obs0[3]), keyby = .(population, run, t)][,
        .(tier2 = sum(tier == 2 & cb == 0) / 7, tier3 = sum(tier == 3 & cb == 0) / 7, cb = sum(cb == 1) / 7),
        keyby = .(population, run)]
    
    hi2 = merge(hi, popsize, by.x = "population", by.y = "Geography")
    hi2 = hi2[, .(population = "England", all_hi = weighted.mean(all_hi, population_size), 
        icu_hi = weighted.mean(icu_hi, population_size)), by = run]
    hi = rbind(hi2, hi)
    
    tier2 = merge(tier, popsize, by.x = "population", by.y = "Geography")
    tier2 = tier2[, .(population = "England", tier2 = weighted.mean(tier2, population_size),
        tier3 = weighted.mean(tier3, population_size), cb = weighted.mean(cb, population_size)), by = run]
    tier = rbind(tier2, tier)

    # Factor
    totals[, population := factor(population, geo)]
    peaks[, population := factor(population, geo)]
    hi[, population := factor(population, geo)]
    tier[, population := factor(population, geo)]
    
    if (is.null(scenario_output)) { # Nice table output
        tab = rbind(
            totals[when == wh, .(indicator = "Deaths", value = niceq(deaths)), keyby = population],
            totals[when == wh, .(indicator = "Admissions", value = niceq(admissions)), keyby = population],
            peaks[, .(indicator = "Peak ICU requirement", value = niceq(icu)), keyby = population],
            peaks[, .(indicator = "Peak ICU (rel. W1)", value = nicepc(rel_icu)), keyby = population],
            hi[, .(indicator = "Weeks of high ICU occupancy", value = niceq(icu_hi)), keyby = population],
            tier[, .(indicator = "Weeks in Tier 2", value = niceq(tier2)), keyby = population],
            tier[, .(indicator = "Weeks in Tier 3", value = niceq(tier3)), keyby = population],
            tier[, .(indicator = "Weeks in lockdown", value = niceq(cb)), keyby = population]
        )
        return (dcast(tab, indicator ~ population))
    } else { # Data output
        tab = rbind(
            totals[when == "after", .(indicator = "Cases",      value = quantile(cases,      c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            totals[when == "after", .(indicator = "Infections", value = quantile(infections, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            totals[when == "after", .(indicator = "Deaths",     value = quantile(deaths,     c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            totals[when == "after", .(indicator = "Admissions", value = quantile(admissions, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            
            prev_peaks[, .(indicator = "Peak hospital beds", value = peak_all, q = "ref"), keyby = population],
            prev_peaks[, .(indicator = "Peak ICU beds", value = peak_icu, q = "ref"), keyby = population],
            peaks[, .(indicator = "Peak hospital beds", value = quantile(all, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            peaks[, .(indicator = "Peak ICU beds", value = quantile(icu, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            
            hi[, .(indicator = "Weeks of high ICU occupancy", value = quantile(icu_hi, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            hi[, .(indicator = "Weeks of high hospital occupancy", value = quantile(all_hi, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            
            tier[, .(indicator = "Weeks in Tier 2", value = quantile(tier2, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            tier[, .(indicator = "Weeks in Tier 3", value = quantile(tier3, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            tier[, .(indicator = "Weeks in lockdown",  value = quantile(cb, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population]
        )
        
        return (cbind(dcast(tab, population + indicator ~ q), scenario = scenario_output))
    }
}




cpp_vec = function(x) 
{
    paste("{", paste(x, collapse = ", "), "}")
}

cpp_obsI_tiers = function(tier_start_date, tier2_threshold, tier3_threshold, pop_size)
{
    tier_start_time = as.numeric(ymd(tier_start_date) - ymd("2020-01-01"));
    
    glue::glue(
        'if (t == 0) { dyn.scratch["tier"] = 0; dyn.scratch["revt"] = 0; dyn.scratch["cb"] = 0; }', # current tier and review time, circuit breaker yes/no
        'if (dyn.scratch["cb"] == 1) {',
        '        P.changes.ch[2].mode = Change::Bypass;',
        '        P.changes.ch[3].mode = Change::Bypass;',
        '}',
        'else if (t >= ${tier_start_time}) {',
        '    double inf_per100k_perweek = dyn("infections_p", t, {}, {}) / ${pop_size / 100000};',
        '    if (inf_per100k_perweek >= ${tier3_threshold} && dyn.scratch["tier"] < 3) {',
        '        P.changes.ch[2].mode = Change::Bypass;',
        '        P.changes.ch[3].mode = Change::Assign;',
        '        dyn.scratch["tier"] = 3;',
        '        dyn.scratch["revt"] = t + 28;',
        '    }',
        '    else if (inf_per100k_perweek >= ${tier2_threshold} && (dyn.scratch["tier"] < 2',
        '             || (t >= dyn.scratch["revt"] && dyn.scratch["tier"] >= 2))) {',
        '        P.changes.ch[2].mode = Change::Assign;',
        '        P.changes.ch[3].mode = Change::Bypass;',
        '        dyn.scratch["tier"] = 2;',
        '        dyn.scratch["revt"] = t + 28;',
        '    }',
        '    else if (inf_per100k_perweek < ${tier2_threshold} && (dyn.scratch["tier"] < 1 || t >= dyn.scratch["revt"])) {',
        '        P.changes.ch[2].mode = Change::Bypass;',
        '        P.changes.ch[3].mode = Change::Bypass;',
        '        dyn.scratch["tier"] = 1;',
        '        dyn.scratch["revt"] = t + 28;',
        '     }',
        '}',
        'dyn.Obs(t, 0, 1, 0) = dyn.scratch["tier"];',
        .sep = "\n", .open = "${", .close = "}")
}

cpp_obsI_cb = function(cb_dates, cb_durations, behaviour)
{
    cb_durations = rep_len(cb_durations, length(cb_dates))
    ret = 'if (t == 0) { dyn.scratch["cb"] = 0; }'; # cb yes/no;
    
    for (k in seq_along(cb_dates)) {
        cb_date = cb_dates[k];
        cb_duration = cb_durations[k];
        
        cb_t0 = as.integer(ymd(cb_date) - ymd("2020-01-01"));
        cb_t1 = cb_t0 + cb_duration;
        
        ret = c(ret,
            glue::glue(
                'if (t == ${cb_t0}) {',
                switch(behaviour,
                    default      = 'dyn.scratch["cb"] = 1;',
                    stay_in_tier = 'dyn.scratch["revt"] = t + 42; dyn.scratch["cb"] = 1;',
                    go_to_tier_3 = 'dyn.scratch["revt"] = t + 42; dyn.scratch["cb"] = 1;'
                ),
                '}',
                'if (t == ${cb_t1}) {',
                '    dyn.scratch["cb"] = 0;',
                if (behaviour == "go_to_tier_3") 'P.changes.ch[2].mode = Change::Bypass; P.changes.ch[3].mode = Change::Assign; dyn.scratch["tier"] = 3;' else '',
                '}',
                .sep = "\n", .open = "${", .close = "}")
        )
    }
    
    ret = c(ret,
        'dyn.Obs(t, 0, 2, 0) = dyn.scratch["cb"];'
    )
    
    return (ret)
}

cpp_chgI_close_schools = function(cb_dates, cb_durations)
{
    cb_durations = rep_len(cb_durations, length(cb_dates))
    ret = c(
        'P.changes.ch[5].values = { vector<double>(8, x[18]) };',
        'P.changes.ch[5].times = { x[24] };'
    )
    
    for (k in seq_along(cb_dates)) {
        cb_date = cb_dates[k];
        cb_duration = cb_durations[k];
        
        cb_t0 = as.integer(ymd(cb_date) - ymd("2020-01-01"));
        cb_t1 = cb_t0 + cb_duration;
        
        ret = c(ret,
            '{',
            glue::glue('double school_close = {cb_t0}, school_open = {cb_t1};'),
            
            # Sep boost
            'P.changes.ch[5].values.push_back(vector<double>(8, 1.0));',
            'P.changes.ch[5].values.push_back(vector<double>(8, x[18]));',
            'P.changes.ch[5].times.push_back(school_close);',
            'P.changes.ch[5].times.push_back(school_open + 1);',
    
            # fitting of google mobility indices
            'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
            '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
            '        double t = P.changes.ch[k].times[i];',
            '        if (t >= school_close && t <= school_open) {',
            '            P.changes.ch[k].values[i][2] = 0;',
            '            P.changes.ch[k].values[i][6] = 0;',
            '        }',
            '    }',
            '}',
            '}'
        )
    }
    
    return (ret)
}

project = function(popset, tiers = FALSE, tier2 = NULL, tier3 = NULL, 
    cb_date = NA, cb_duration = 14, lockdown = NULL, close_schools = FALSE, cb_behaviour = "default", 
    waning_duration = -1, seasonality = 0, n_run = 100, expire = NULL)
{
    dynamicsO = list()
    for (p in popset) {
        # calculate population size
        pop_size = sum(parametersI[[p]]$pop[[1]]$size);
        
        # load user defined functions
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = c(
                        cpp_chgI(),
                        if (close_schools) cpp_chgI_close_schools(cb_date, cb_duration) else NULL
                    ),
                    cpp_loglikelihood = "",
                    cpp_observer = c(
                        cpp_obsI(P.death), 
                        if (!is.na(cb_date[1])) cpp_obsI_cb(cb_date, cb_duration, cb_behaviour) else NULL,
                        if (tiers) cpp_obsI_tiers("2020-10-14", 700, 2100, pop_size) else NULL,
                        if (waning_duration > 0) paste("if (t == 274) { P.pop[0].wn = vector<double>(16,", 1 / waning_duration, "); }") else NULL,
                        if (seasonality > 0) paste("if (t == 274) { P.pop[0].season_A[0] = ", seasonality, "; }") else NULL
                    )
                )
            )
        )
        
        # Set parameters for this run . . .
        paramsI = rlang::duplicate(parametersI[[p]])
        paramsI$time1 = "2021-03-31";
        
        # Expire mobility changes after date "expire"
        if (!is.null(expire)) {
            expire_t = as.numeric(ymd(expire) - ymd(paramsI$date0))
            index = which(paramsI$schedule[[1]]$times == expire_t)
            for (i in index:length(paramsI$schedule[[1]]$values)) {
                paramsI$schedule[[1]]$values[[i]] = paramsI$schedule[[1]]$values[[index]]
            }
            paramsI$schedule[[1]]$values
        }
        
        # adjust schedule1 for lockdown 
        if (!is.na(cb_date[1])) {
            for (k in seq_along(cb_date)) {
                for (i in seq_along(paramsI$schedule[[1]]$values)) {
                    today = ymd(paramsI$date0) + paramsI$schedule[[1]]$times[i];
                    if (today %between% c(ymd(cb_date[k]), ymd(cb_date[k]) + rep_len(cb_duration, length(cb_date))[k] - 1)) {
                        paramsI$schedule[[1]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + lockdown$res / 100); # res
                        paramsI$schedule[[1]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + lockdown$wor / 100); # work
                        paramsI$schedule[[1]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + lockdown$gro / 100); # groc
                        paramsI$schedule[[1]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + lockdown$ret / 100); # retail
                        paramsI$schedule[[1]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + lockdown$tra / 100); # transit
                    }
                }
            }
        }
        
        # contacts for tier 2
        for (i in seq_along(paramsI$schedule[[3]]$values)) {
            paramsI$schedule[[3]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + tier2$res / 100); # res
            paramsI$schedule[[3]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + tier2$wor / 100); # work
            paramsI$schedule[[3]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + tier2$gro / 100); # groc
            paramsI$schedule[[3]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + tier2$ret / 100); # retail
            paramsI$schedule[[3]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + tier2$tra / 100); # transit
        }
    
        # contacts for tier 3
        for (i in seq_along(paramsI$schedule[[4]]$values)) {
            paramsI$schedule[[4]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + tier3$res / 100); # res
            paramsI$schedule[[4]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + tier3$wor / 100); # work
            paramsI$schedule[[4]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + tier3$gro / 100); # groc
            paramsI$schedule[[4]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + tier3$ret / 100); # retail
            paramsI$schedule[[4]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + tier3$tra / 100); # transit
        }
        
        paramsI$processes[[length(paramsI$processes) + 1]] = list(
            source = "E",
            type = "multinomial",
            names = "infections",
            report = "ip",
            prob = matrix(1, nrow = 1, ncol = 16),
            delays = matrix(c(rep(0, 7 * 4 - 1), 1), nrow = 1, byrow = T)
        )
        
        # Sampling fits
        test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI), posteriorsI[[p]], n_run, seed = 0);

        test = rbindlist(test)
        test = test[, population := p]
        dynamicsO[[p]] = test
        print(p)
    }

    testX = rbindlist(dynamicsO)
    testX[, population := nhs_regions[population]]

    return (testX)
}

make_data = function(ld, sitreps)
{
    rbind(
        ld[, .(ValueType = "type28_death_inc_line", Geography = name,
            dmin = NA, d = as.Date(date), dmax = NA, ymin = NA, y = N, ymax = NA)],
        sitreps[, .(ValueType = "icu_prev", Geography = name,
            dmin = NA, d = as.Date(date), dmax = NA, ymin = NA, y = n_in_itu, ymax = NA)],
        sitreps[, .(ValueType = "hospital_prev", Geography = name,
            dmin = NA, d = as.Date(date), dmax = NA, ymin = NA, y = n_in_all_beds, ymax = NA)],
        sitreps[, .(ValueType = "hospital_inc", Geography = name,
            dmin = NA, d = as.Date(date), dmax = NA, ymin = NA, y = n_admitted_diagnosed, ymax = NA)]
    )
}

niceq = function(x)
{
    q = quantile(x, c(0.05, 0.5, 0.95));
    f = function(y) prettyNum(signif(y, 3), big.mark = ",")
    paste0(f(q[2]), " (", f(q[1]), " - ", f(q[3]), ")")
}

nicepc = function(x)
{
    q = quantile(x, c(0.05, 0.5, 0.95));
    f = function(y) round(y * 100, 0)
    paste0(f(q[2]), "% (", f(q[1]), " - ", f(q[3]), "%)")
}

england_pops = c(1, 3, 4, 5, 6, 9, 10)

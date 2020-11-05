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

popsize = NULL
for (i in seq_along(parametersI)) {
    if (!is.null(parametersI[[i]])) {
        popsize = rbind(popsize,
            data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = sum(parametersI[[i]]$pop[[1]]$size))
        )
    }
}
popsize = rbind(data.table(Geography = "England", population_size = sum(popsize$population_size)), popsize)
setkey(popsize, Geography)



# 2. Bringing in of interventions

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


summarize_projection = function(proj, from_date, popsize, scenario_output = NULL, to_date = NULL)
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
    
    post_peaks = proj[ymd("2020-01-01") + t >= from_date,
        .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
        by = .(population, run, t)][,
        .(all = max(all), icu = max(icu)),
        keyby = .(population, run)]
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
            totals[when == "after", .(indicator = "Deaths", value = niceq(deaths)), keyby = population],
            totals[when == "after", .(indicator = "Admissions", value = niceq(admissions)), keyby = population],
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

cpp_obsI_cb = function(cb_date, cb_duration, behaviour)
{
    cb_t0 = as.integer(ymd(cb_date) - ymd("2020-01-01"));
    cb_t1 = cb_t0 + cb_duration;
    
    glue::glue(
        'if (t == 0) { dyn.scratch["cb"] = 0; }', # cb yes/no
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
        'dyn.Obs(t, 0, 2, 0) = dyn.scratch["cb"];',
        .sep = "\n", .open = "${", .close = "}")
}

cpp_chgI_close_schools = function(cb_date, cb_duration)
{
    cb_t0 = as.integer(ymd(cb_date) - ymd("2020-01-01"));
    cb_t1 = cb_t0 + cb_duration;
    
    c(
        glue::glue('double school_close = {cb_t0}, school_open = {cb_t1};'),
        
        # Sep boost
        'P.changes.ch[5].values = { vector<double>(8, x[18]), vector<double>(8, 1.0), vector<double>(8, x[18])};',
        'P.changes.ch[5].times = { x[24], school_close, school_open + 1 };',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        double t = P.changes.ch[k].times[i];',
        '        if (t >= school_close && t <= school_open) {',
        '            P.changes.ch[k].values[i][2] = 0;',
        '            P.changes.ch[k].values[i][6] = 0;',
        '        }',
        '    }',
        '}'
    )
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
                        if (!is.na(cb_date)) cpp_obsI_cb(cb_date, cb_duration, cb_behaviour) else NULL,
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
        if (!is.na(cb_date)) {
            for (i in seq_along(paramsI$schedule[[1]]$values)) {
                today = ymd(paramsI$date0) + paramsI$schedule[[1]]$times[i];
                if (today %between% c(ymd(cb_date), ymd(cb_date) + cb_duration - 1)) {
                    paramsI$schedule[[1]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + lockdown$res / 100); # res
                    paramsI$schedule[[1]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + lockdown$wor / 100); # work
                    paramsI$schedule[[1]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + lockdown$gro / 100); # groc
                    paramsI$schedule[[1]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + lockdown$ret / 100); # retail
                    paramsI$schedule[[1]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + lockdown$tra / 100); # transit
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


tier2 = list(res = 0.55, wor = -0.36, gro = -1.41, ret = -1.99, tra = -1.32)
tier3 = list(res = 1.44, wor = -2.72, gro = -1.98, ret = -8.83, tra = -4.62)

lockdown0 = list(res = 0, wor = 0, gro = 0, ret = 0, tra = 0)

# Wales lockdown
#                 variable         V1
# 1:  grocery_and_pharmacy -19.626778
# 2:                 parks -44.768158
# 3:           residential   7.596776
# 4: retail_and_recreation -41.182127
# 5:      transit_stations -20.958923
# 6:            workplaces -21.720779
lockdownW = list(res = 7.59, wor = -21.72, gro = -19.62, ret = -41.18, tra = -20.96)

# NI lockdown
# IMPACT OF LOCKDOWN FOR NORTHERN IRELAND:
#                 variable          V1
# 1:  grocery_and_pharmacy  -0.7792208
# 2:                 parks   7.5283163
# 3:           residential   4.7287157
# 4: retail_and_recreation -13.8701299
# 5:      transit_stations  -9.7984127
# 6:            workplaces -14.4220779
lockdownN = list(res = 4.73, wor = -14.42, gro = -0.78, ret = -13.87, tra = -9.79)

# Scotland restrictions # OLD
# IMPACT OF LOCKDOWN FOR SCOTLAND:
#                 variable         V1
# 1:  grocery_and_pharmacy  0.5931034
# 2:           residential  2.0497126
# 3: retail_and_recreation -7.0221675
# 4:      transit_stations -6.6267492
# 5:            workplaces -8.2720801
lockdownS = list(res = 2.05, wor = -8.27, gro = 0.59, ret = -7.02, tra = -6.63)





# RUN PROJECTIONS

# BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
# Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
# Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd

# ldN4o = lockdown, like Northern Ireland, 4 weeks, schools open
proj_basic = project(england_pops, tiers = FALSE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_tiers = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldN4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownN, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldN4c = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownN, close_schools = TRUE,  cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4c = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = TRUE,  cb_behaviour = "default", expire = "2020-10-13")

plot_projection(list(proj_basic, proj_tiers), list("No tiers", "Tiers"), "2020-10-01", "Dark2")
ggsave("./figures/proj_tiers.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_tiers.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)
plot_projection(list(proj_tiers, proj_ldN4o, proj_ldN4c), list("Tiers only", "NI lockdown, schools open", "NI lockdown, schools closed"), "2020-10-01", "Set2")
ggsave("./figures/proj_ld_N.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_ld_N.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)
plot_projection(list(proj_tiers, proj_ldW4o, proj_ldW4c), list("Tiers only", "Wales lockdown, schools open", "Wales lockdown, schools closed"), "2020-10-01", "Set2")
ggsave("./figures/proj_ld_W.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_ld_W.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)

# Different duration of lockdown
proj_ldW1o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 6, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW2o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 13, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW3o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 20, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW5o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 34, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW6o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 41, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")

# Different timing of lockdown
proj_ldW4o_10_08 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-08", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_10_15 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-15", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_10_22 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-22", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_10_29 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-29", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_11_05 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_11_12 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-12", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_11_19 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-19", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")

# Sensitivity analyses with waning and seasonality
proj_tiers_w  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0)
proj_ld4Wo_w  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0)
proj_tiers_s  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = -1, seasonality = 0.1)
proj_ld4Wo_s  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = -1, seasonality = 0.1)
proj_tiers_ws = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0.1)
proj_ld4Wo_ws = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0.1)

plot_projection(list(proj_tiers_w, proj_ldW4o_w), list("Tiers + waning immunity", "Lockdown + waning immunity"), "2020-10-01", "Set1")
ggsave("./figures/proj_sensitivity_w.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_sensitivity_w.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)
plot_projection(list(proj_tiers_s, proj_ldW4o_s), list("Tiers + seasonality", "Lockdown + seasonality"), "2020-10-01", "Set1")
ggsave("./figures/proj_sensitivity_s.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_sensitivity_s.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)
plot_projection(list(proj_tiers_ws, proj_ldW4o_ws), list("Tiers + waning immunity + seasonality", "Lockdown + waning immunity + seasonality"), "2020-10-01", "Set1")
ggsave("./figures/proj_sensitivity_ws.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_sensitivity_ws.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)


# Summary tables
tbBa = summarize_projection(proj_basic, "2020-10-01", popsize)
fwrite(tbBa, "./figures/proj_basic.csv");
tbTi = summarize_projection(proj_tiers, "2020-10-01", popsize)
fwrite(tbTi, "./figures/proj_tiers.csv");
tbNo = summarize_projection(proj_ldN4o, "2020-10-01", popsize)
fwrite(tbNo, "./figures/proj_ldN4o.csv");
tbNc = summarize_projection(proj_ldN4c, "2020-10-01", popsize)
fwrite(tbNc, "./figures/proj_ldN4c.csv");
tbWo = summarize_projection(proj_ldW4o, "2020-10-01", popsize)
fwrite(tbWo, "./figures/proj_ldW4o.csv");
tbWc = summarize_projection(proj_ldW4c, "2020-10-01", popsize)
fwrite(tbWc, "./figures/proj_ldW4c.csv");

tbEngland = england_only(list(tbBa, tbTi, tbNo, tbNc, tbWo, tbWc), 
    c("Baseline", "Tiers only", "NI-style lockdown, schools open", "NI-style lockdown, schools closed", "Wales-style lockdown, schools open", "Wales-style lockdown, schools closed"))
fwrite(tbEngland, "./figures/table_england.csv");

tbs_tiers_w = summarize_projection(proj_tiers_w, "2020-10-01", popsize)
tbs_ldown_w = summarize_projection(proj_ldW4o_w, "2020-10-01", popsize)
tbs_tiers_s = summarize_projection(proj_tiers_s, "2020-10-01", popsize)
tbs_ldown_s = summarize_projection(proj_ldW4o_s, "2020-10-01", popsize)
tbs_tiers_ws = summarize_projection(proj_tiers_ws, "2020-10-01", popsize)
tbs_ldown_ws = summarize_projection(proj_ldW4o_ws, "2020-10-01", popsize)
tbSensitivity = england_only(list(tbTi, tbWo, tbs_tiers_s, tbs_ldown_s, tbs_tiers_w, tbs_ldown_w, tbs_tiers_ws, tbs_ldown_ws), 
    c("Tiers only", "Lockdown", "Tiers only + seasonality", "Lockdown + seasonality", "Tiers only + waning", "Lockdown + waning", "Tiers only + seasonality + waning", "Lockdown + seasonality + waning"))
fwrite(tbSensitivity, "./figures/table_sensitivity.csv");


arrs = function(proj) {
    w = proj[group == 1 & t %in% c(308, 310), obs0, by = .(run, population, t)]
    w = rbind(w, data.table(population = w[t == 310, population], obs0 = w[t == 308]$obs0 - w[t == 310]$obs0), fill = TRUE)
    
    w = w[, .(Rt = niceq(obs0)), by = .(population, t)]
    w[t == 308, when := "Pre-lockdown"]
    w[t == 310, when := "Post-lockdown"]
    w[is.na(t), when := "Reduction"]
    w[, when := factor(when, unique(when))]
    dcast(w, population ~ when, value.var = "Rt")
}

arrs(proj_ldN4o)
arrs(proj_ldN4c)
fwrite(arrs(proj_ldW4o), "./figures/R0_table_Wo.csv")
fwrite(arrs(proj_ldW4c), "./figures/R0_table_Wc.csv")

series_weeks = rbind(
    summarize_projection(proj_tiers, "2020-10-01", popsize, "No lockdown"),
    summarize_projection(proj_ldW1o, "2020-10-01", popsize, "1 week"),
    summarize_projection(proj_ldW2o, "2020-10-01", popsize, "2 weeks"),
    summarize_projection(proj_ldW3o, "2020-10-01", popsize, "3 weeks"),
    summarize_projection(proj_ldW4o, "2020-10-01", popsize, "4 weeks"),
    summarize_projection(proj_ldW5o, "2020-10-01", popsize, "5 weeks"),
    summarize_projection(proj_ldW6o, "2020-10-01", popsize, "6 weeks")
)
series_weeks[, scenario := factor(scenario, unique(scenario))]

series_timing = rbind(
    summarize_projection(proj_ldW4o_10_08, "2020-10-01", popsize, "8 Oct"),
    summarize_projection(proj_ldW4o_10_15, "2020-10-01", popsize, "15 Oct"),
    summarize_projection(proj_ldW4o_10_22, "2020-10-01", popsize, "22 Oct"),
    summarize_projection(proj_ldW4o_10_29, "2020-10-01", popsize, "29 Oct"),
    summarize_projection(proj_ldW4o_11_05, "2020-10-01", popsize, "5 Nov"),
    summarize_projection(proj_ldW4o_11_12, "2020-10-01", popsize, "12 Nov"),
    summarize_projection(proj_ldW4o_11_19, "2020-10-01", popsize, "19 Nov")
)
series_timing[, scenario := factor(scenario, unique(scenario))]

series_what = rbind(
    summarize_projection(proj_basic, "2020-10-01", popsize, "Nothing"),
    summarize_projection(proj_tiers, "2020-10-01", popsize, "Tiers"),
    summarize_projection(proj_ldN4o, "2020-10-01", popsize, "Ld N/o"),
    summarize_projection(proj_ldN4c, "2020-10-01", popsize, "Ld N/c"),
    summarize_projection(proj_ldW4o, "2020-10-01", popsize, "Ld W/o"),
    summarize_projection(proj_ldW4c, "2020-10-01", popsize, "Ld W/c")
)
series_what[, scenario := factor(scenario, unique(scenario))]

# series_weeks_to_dec31 = rbind(
#     summarize_projection(proj_tiers, "2020-10-01", popsize, "No lockdown", "2020-12-31"),
#     summarize_projection(proj_ldW1o, "2020-10-01", popsize, "1 week", "2020-12-31"),
#     summarize_projection(proj_ldW2o, "2020-10-01", popsize, "2 weeks", "2020-12-31"),
#     summarize_projection(proj_ldW3o, "2020-10-01", popsize, "3 weeks", "2020-12-31"),
#     summarize_projection(proj_ldW4o, "2020-10-01", popsize, "4 weeks", "2020-12-31"),
#     summarize_projection(proj_ldW5o, "2020-10-01", popsize, "5 weeks", "2020-12-31"),
#     summarize_projection(proj_ldW6o, "2020-10-01", popsize, "6 weeks", "2020-12-31")
# )
# series_weeks_to_dec31[, scenario := factor(scenario, unique(scenario))]

plot_indicator = function(series, indicator_name)
{
    ggplot(series[indicator == indicator_name]) +
        geom_pointrange(aes(x = population, ymin = lo, y = mid, ymax = hi, colour = scenario), position = position_dodge(width = 0.8), shape = 20, size = 0.0) +
        geom_linerange(aes(x = population, ymin = lo, y = mid, ymax = hi, colour = scenario), position = position_dodge(width = 0.8), size = 0.2) +
        labs(x = NULL, y = indicator_name) + ylim(0, NA)
}

plot_indicators_england = function(series, indicator_names, y_axis_title, unit = 1, legpos = "none", pal = "Dark2")
{
    ggplot(series[population == "England" & indicator %in% indicator_names]) +
        geom_ribbon(aes(x = scenario, ymin = lo / unit, ymax = hi / unit, fill = indicator, group = indicator), alpha = 0.5) +
        geom_line(aes(x = scenario, y = mid / unit, colour = indicator, group = indicator), size = 0.2) +
        labs(x = NULL, y = y_axis_title, fill = NULL, colour = NULL) + ylim(0, NA) +
        theme(panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2),
            legend.position = legpos) +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
}

plot_indicators_england_stack = function(series, indicator_names, y_axis_title, stack_order, legpos = "none", pal = "Set2")
{
    ser2 = rlang::duplicate(series[population == "England" & indicator %in% indicator_names])
    ser2[, indicator := factor(indicator, stack_order)]
    ser3 = rlang::duplicate(ser2)
    ser3 = ser3[order(scenario, rev(indicator))]
    ser3[, lo := lo - mid]
    ser3[, hi := hi - mid]
    ser3[, mid := cumsum(mid), by = scenario]
    ser3[, lo := mid + lo]
    ser3[, hi := mid + hi]

    ggplot() +
        geom_area(data = ser2, aes(x = scenario, y = mid, fill = indicator, group = indicator)) +
        geom_ribbon(data = ser3, aes(x = scenario, ymin = lo, ymax = hi, group = indicator), alpha = 0.25) +
        labs(x = NULL, y = y_axis_title, fill = NULL, colour = NULL) + ylim(0, NA) +
        theme(panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2),
            legend.position = legpos) +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
}

plot_indicator_regions = function(series, indicator_name, label = indicator_name, suffix = "\n(Thousands)", unit = 1000, legpos = "none", bigtitle = NULL)
{
    ser2 = rlang::duplicate(series[population != "England" & indicator == indicator_name])
    ser2[population == "East of England", population := "EE"]
    ser2[population == "London", population := "Ldn"]
    ser2[population == "Midlands", population := "Mlds"]
    ser2[population == "North East and Yorkshire", population := "NE&Y"]
    ser2[population == "North West", population := "NW"]
    ser2[population == "South East", population := "SE"]
    ser2[population == "South West", population := "SW"]
    ggplot(ser2) +
        geom_point(aes(x = population, y = mid / unit, colour = scenario), position = position_dodge(width = 0.8), shape = 20, size = 0.5) +
        geom_linerange(aes(x = population, ymin = lo / unit, ymax = hi / unit, colour = scenario), position = position_dodge(width = 0.8), size = 0.3) +
        labs(x = NULL, y = paste0(label, suffix), colour = NULL, title = bigtitle) + 
        ylim(0, NA) +
        theme(legend.position = legpos)
}

# Regional plots: Type of intervention
pr1  = plot_indicator_regions(series_what, "Deaths", legpos = "top", bigtitle = "Type of intervention")
pr2  = plot_indicator_regions(series_what, "Admissions")
pr3  = plot_indicator_regions(series_what, "Cases")
pr4  = plot_indicator_regions(series_what, "Infections")
pr5  = plot_indicator_regions(series_what, "Peak ICU beds")
pr6  = plot_indicator_regions(series_what, "Weeks of high ICU occupancy", "High ICU weeks", suffix = "", unit = 1)
pr7  = plot_indicator_regions(series_what, "Peak hospital beds", "Peak hosp beds")
pr8  = plot_indicator_regions(series_what, "Weeks of high hospital occupancy", "High hosp weeks", suffix = "", unit = 1)
pr9  = plot_indicator_regions(series_what, "Weeks in Tier 2", "Weeks in\nTier 2", suffix = "", unit = 1)
pr10 = plot_indicator_regions(series_what, "Weeks in Tier 3", "Weeks in\nTier 3", suffix = "", unit = 1)
pr11 = plot_indicator_regions(series_what, "Weeks in lockdown", "Weeks in\nlockdown", suffix = "", unit = 1)
ppr1 = plot_grid(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, pr9, pr10, pr11, ncol = 1, align = "v", rel_heights = c(1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
ggsave("./figures/regional_what.pdf", ppr1, width = 7, height = 21, units = "cm")

# Regional plots: Duration of lockdown
pr1  = plot_indicator_regions(series_weeks, "Deaths", legpos = "top", bigtitle = "Duration of lockdown")
pr2  = plot_indicator_regions(series_weeks, "Admissions")
pr3  = plot_indicator_regions(series_weeks, "Cases")
pr4  = plot_indicator_regions(series_weeks, "Infections")
pr5  = plot_indicator_regions(series_weeks, "Peak ICU beds")
pr6  = plot_indicator_regions(series_weeks, "Weeks of high ICU occupancy", "High ICU weeks", suffix = "", unit = 1)
pr7  = plot_indicator_regions(series_weeks, "Peak hospital beds", "Peak hosp beds")
pr8  = plot_indicator_regions(series_weeks, "Weeks of high hospital occupancy", "High hosp weeks", suffix = "", unit = 1)
pr9  = plot_indicator_regions(series_weeks, "Weeks in Tier 2", "Weeks in\nTier 2", suffix = "", unit = 1)
pr10 = plot_indicator_regions(series_weeks, "Weeks in Tier 3", "Weeks in\nTier 3", suffix = "", unit = 1)
pr11 = plot_indicator_regions(series_weeks, "Weeks in lockdown", "Weeks in\nlockdown", suffix = "", unit = 1)
ppr2 = plot_grid(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, pr9, pr10, pr11, ncol = 1, align = "v", rel_heights = c(1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
ggsave("./figures/regional_weeks.pdf", ppr2, width = 7, height = 21, units = "cm")

# Regional plots: Timing of lockdown
pr1  = plot_indicator_regions(series_timing, "Deaths", legpos = "top", bigtitle = "Timing of lockdown")
pr2  = plot_indicator_regions(series_timing, "Admissions")
pr3  = plot_indicator_regions(series_timing, "Cases")
pr4  = plot_indicator_regions(series_timing, "Infections")
pr5  = plot_indicator_regions(series_timing, "Peak ICU beds")
pr6  = plot_indicator_regions(series_timing, "Weeks of high ICU occupancy", "High ICU weeks", suffix = "", unit = 1)
pr7  = plot_indicator_regions(series_timing, "Peak hospital beds", "Peak hosp beds")
pr8  = plot_indicator_regions(series_timing, "Weeks of high hospital occupancy", "High hosp weeks", suffix = "", unit = 1)
pr9  = plot_indicator_regions(series_timing, "Weeks in Tier 2", "Weeks in\nTier 2", suffix = "", unit = 1)
pr10 = plot_indicator_regions(series_timing, "Weeks in Tier 3", "Weeks in\nTier 3", suffix = "", unit = 1)
pr11 = plot_indicator_regions(series_timing, "Weeks in lockdown", "Weeks in\nlockdown", suffix = "", unit = 1)
ppr3 = plot_grid(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, pr9, pr10, pr11, ncol = 1, align = "v", rel_heights = c(1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
ggsave("./figures/regional_timing.pdf", ppr3, width = 7, height = 21, units = "cm")

ppr = plot_grid(ppr1, ppr2, ppr3, nrow = 1, label_size = 8, labels = letters)
ggsave("./figures/regional.pdf", ppr, width = 21, height = 21, units = "cm", useDingbats = FALSE)
ggsave("./figures/regional.png", ppr, width = 21, height = 21, units = "cm")



# Comparative plots: Type of intervention
p0 = plot_cum_deaths(list(proj_basic, proj_tiers, proj_ldN4o, proj_ldN4c, proj_ldW4o, proj_ldW4c), 
    c("Nothing", "Tiers", "Ld N/o", "Ld N/c", "Ld W/o", "Ld W/c"), 
    "2020-10-01", ystart = 28, ydiff = -2.1, titl = "Type of intervention")
p7 = plot_indicators_england(series_what, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack(series_what, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp1 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/big_what.pdf", pp1, width = 7, height = 12, units = "cm")


# Comparative plots: Duration of lockdown
p0 = plot_cum_deaths(list(proj_tiers, proj_ldW1o, proj_ldW2o, proj_ldW3o, proj_ldW4o, proj_ldW5o, proj_ldW6o), 
    c("No lockdown", "1 week", "2 weeks", "3 weeks", "4 weeks", "5 weeks", "6 weeks"), 
    "2020-10-01", ystart = 25.9, ydiff = -1.7, titl = "Duration of lockdown")
p7 = plot_indicators_england(series_weeks, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack(series_weeks, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp2 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/big_duration.pdf", pp2, width = 7, height = 12, units = "cm")


# Comparative plots: Timing of lockdown
p0 = plot_cum_deaths(list(proj_ldW4o_10_08, proj_ldW4o_10_15, proj_ldW4o_10_22, proj_ldW4o_10_29, proj_ldW4o_11_05, proj_ldW4o_11_12, proj_ldW4o_11_19), 
    c(  "                  8 Oct", 
        "                      15 Oct", 
        "                          22 Oct", 
        "29 Oct", 
        "5 Nov", 
        "12 Nov", 
        "19 Nov"), 
    "2020-10-01", ystart = 21.6, ydiff = -1.35, titl = "Timing of lockdown")
p7 = plot_indicators_england(series_timing, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack(series_timing, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp3 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/big_timing.pdf", pp3, width = 7, height = 12, units = "cm")


pp = plot_grid(pp1, pp2, pp3, nrow = 1, label_size = 8, labels = letters)
ggsave("./figures/big.pdf", pp, width = 21, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./figures/big.png", pp, width = 21, height = 12, units = "cm")

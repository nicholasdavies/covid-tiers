# Run projections for tiers and lockdown scenarios
source("./R/projections_setup.R")
source("./R/define_tiers.R")


# RUN PROJECTIONS

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

## TEMP ##
proj_ld2Woa = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = c("2020-10-08", "2021-01-07"), cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ld2Wob = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = c("2020-10-22", "2021-01-07"), cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ld2Woc = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = c("2020-10-08", "2020-12-31"), cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ld2Wod = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = c("2020-10-22", "2020-12-31"), cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ld2Woe = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = c("2020-10-08", "2020-12-10"), cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ld2Wof = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = c("2020-10-01", "2020-12-03"), cb_duration = c(27, 13), lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")

plot_projection(list(proj_ldW4o, proj_ld2Woa, proj_ld2Wob, proj_ld2Woc, proj_ld2Wod, proj_ld2Woe, proj_ld2Wof), 
    list("One lockdown", "Two lockdowns a", "Two lockdowns b", "Two lockdowns c", "Two lockdowns d", "Two lockdowns e", "Two lockdowns f"), 
    "2020-10-01", "Set2")

plot_cum_deaths(list(proj_ldW4o, proj_ld2Woa, proj_ld2Wob, proj_ld2Woc, proj_ld2Wod, proj_ld2Woe, proj_ld2Wof), 
    c("1 Ld", "2 Ld a", "2 Ld b", "2 Ld c", "2 Ld d", "2 Ld e", "2 Ld f"), 
    "2020-10-01", ystart = 28, ydiff = -2.1, titl = "Type of intervention")
## / TEMP ##

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
proj_ldW4o_w  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0)
proj_tiers_s  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = -1, seasonality = 0.1)
proj_ldW4o_s  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = -1, seasonality = 0.1)
proj_tiers_ws = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0.1)
proj_ldW4o_ws = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
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
tbBa_before = summarize_projection(proj_basic, "2020-10-01", popsize, wh = "before")

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

R_ldN4o = arrs(proj_ldN4o)
R_ldN4c = arrs(proj_ldN4c)
R_ldW4o = arrs(proj_ldW4o)
R_ldW4c = arrs(proj_ldW4c)

combine_R = function(R_ldN4o, R_ldN4c) {
    R_ldN4 = cbind(R_ldN4o[, 1:3], R_ldN4c[, 3], R_ldN4o[, 4], R_ldN4c[, 4])
    names(R_ldN4)[c(3,5)] = paste0(names(R_ldN4)[c(3,5)], ", schools open")
    names(R_ldN4)[c(4,6)] = paste0(names(R_ldN4)[c(4,6)], ", schools closed")
    R_ldN4
}

fwrite(combine_R(R_ldN4o, R_ldN4c), "./figures/R0_table_N.csv")
fwrite(combine_R(R_ldW4o, R_ldW4c), "./figures/R0_table_W.csv")

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
    c("Baseline", "Tiers", "Ld N/o", "Ld N/c", "Ld W/o", "Ld W/c"), 
    "2020-10-01", ystart = 28, ydiff = -2.1, titl = "Type of intervention")
p1 = plot_icu(list(proj_basic, proj_tiers, proj_ldN4o, proj_ldN4c, proj_ldW4o, proj_ldW4c), 
    c("Baseline", "Tiers", "Ld N/o", "Ld N/c", "Ld W/o", "Ld W/c"), 
    "2020-10-01")
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
p1 = plot_icu(list(proj_tiers, proj_ldW1o, proj_ldW2o, proj_ldW3o, proj_ldW4o, proj_ldW5o, proj_ldW6o), 
    c("No lockdown", "1 week", "2 weeks", "3 weeks", "4 weeks", "5 weeks", "6 weeks"), 
    "2020-10-01")
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
p1 = plot_icu(list(proj_ldW4o_10_08, proj_ldW4o_10_15, proj_ldW4o_10_22, proj_ldW4o_10_29, proj_ldW4o_11_05, proj_ldW4o_11_12, proj_ldW4o_11_19), 
    c(  "                  8 Oct", 
        "                      15 Oct", 
        "                          22 Oct", 
        "29 Oct", 
        "5 Nov", 
        "12 Nov", 
        "19 Nov"), 
    "2020-10-01")
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

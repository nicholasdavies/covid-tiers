library(data.table)
library(ggplot2)
library(zoo)
library(mgcv)
library(lubridate)
library(stringr)
library(cowplot)
library(qs)

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

# Load Google Mobility data
gm = qread("./data/google_mobility_uk.qs")

# Other places' lockdowns / circuit breakers
gm_melt = function(gm) {
    g = melt(gm, id.vars = 1:8)
    g[, variable := str_remove_all(variable, "_percent_change_from_baseline")]
    return (g)
}

# Look at Northern Ireland and at Wales.

# Northern Ireland
ni = c(
"Antrim and Newtownabbey",
"Ards and North Down",
"Armagh City, Banbridge and Craigavon",
"Belfast",
"Causeway Coast and Glens",
"Derry and Strabane",
"Fermanagh and Omagh",
"Lisburn and Castlereagh",
"Mid and East Antrim",
"Mid Ulster",
"Newry, Mourne and Down")

gm_ni = gm_melt(gm[country_region_code == "GB" & sub_region_1 %in% ni])
gm_ni = gm_ni[, .(value = mean(value, na.rm = T)), keyby = .(variable, date)]
gm_ni[, rollval := rollmean(value, k = 7, fill = "extend"), by = variable]

# Plot NI google mobility data
pl_ni_mob = ggplot(gm_ni) +
    geom_line(aes(x = date, y = rollval, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-16")), linetype = "22", size = 0.3) +
    facet_wrap(~variable) +
    theme(legend.position = "none") +
    labs(x = "Date", y = "Rolling mean - Northern Ireland")

# Compare post-lockdown values to pre-lockdown values (2 week lookback)
compare = gm_ni[date %between% c(ymd("2020-10-17", "2020-10-30"))]
compare = merge(compare,
    gm_ni[date %between% c(ymd("2020-10-17", "2020-10-30") - 14), 
        .(date = date + 14, variable, baseline = value)],
    by = c("date", "variable"))
compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, mean(change), by = variable]
compare[, variable := str_replace_all(variable, "_", " ")]
compare[, variable := str_to_sentence(variable)]
compare[, variable := factor(variable)]

summary = compare[, mean(change), by = variable]

# IMPACT OF LOCKDOWN FOR NORTHERN IRELAND:
#                 variable          V1
# 1:  grocery_and_pharmacy  -0.7792208
# 2:                 parks   7.5283163
# 3:           residential   4.7287157
# 4: retail_and_recreation -13.8701299
# 5:      transit_stations  -9.7984127
# 6:            workplaces -14.4220779

# Visual comparison
pl_ni = ggplot(compare) +
    geom_hline(aes(yintercept = 0), linetype = "23") + 
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    geom_label(data = summary, aes(x = ymd("2020-10-23") + 0.5, y = 45, 
        label = paste0(ifelse(V1 > 0, "+", ""), round(V1, 2))), alpha = 0.75) +
    facet_wrap(~variable) +
    labs(x = "Date", y = "Change in Google Mobility index")



# Wales
iso3166 = fread("./data/county_iso3166_2_gb.csv")

wales_iso3166 = iso3166[gss_code %like% "^W", unique(iso_code)];
gm_w = gm[iso_3166_2_code %in% wales_iso3166, lapply(.SD, function(x) mean(x, na.rm = TRUE)), .SDcols = 9:14, by = .(date)]
gm_w = melt(gm_w, id.vars = "date")
gm_w[, variable := str_remove_all(variable, "_percent_change_from_baseline")]
gm_w[, rollval := rollmean(value, k = 7, fill = "extend"), by = variable]


# Plot Wales google mobility data
pl_w_mob = ggplot(gm_w) + 
    geom_line(aes(x = date, y = rollval, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-23")), linetype = "22", size = 0.3) +
    facet_wrap(~variable) +
    theme(legend.position = "none") +
    labs(x = "Date", y = "Rolling mean - Wales")

# Compare post-lockdown values to pre-lockdown values (2 week lookback)
compare = gm_w[date %between% c(ymd("2020-10-24", "2020-11-06"))]
compare = merge(compare,
    gm_w[date %between% c(ymd("2020-10-24", "2020-11-06") - 14), 
        .(date = date + 14, variable, baseline = value)],
    by = c("date", "variable"))
compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, mean(change), by = variable]
compare[, variable := str_replace_all(variable, "_", " ")]
compare[, variable := str_to_sentence(variable)]
compare[, variable := factor(variable)]

summary = compare[, mean(change), by = variable]


# IMPACT OF LOCKDOWN FOR WALES:
#                 variable         V1
# 1:  grocery_and_pharmacy -19.626778
# 2:                 parks -44.768158
# 3:           residential   7.596776
# 4: retail_and_recreation -41.182127
# 5:      transit_stations -20.958923
# 6:            workplaces -21.720779

# Visual comparison
pl_w = ggplot(compare) +
    geom_hline(aes(yintercept = 0), linetype = "23") + 
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    geom_label(data = summary, aes(x = ymd("2020-10-27"), y = -10, 
        label = paste0(ifelse(V1 > 0, "+", ""), round(V1, 2))), alpha = 0.75) +
    facet_wrap(~variable) +
    labs(x = "Date", y = "Change in Google Mobility index")

cowplot::plot_grid(pl_ni_mob, pl_ni, pl_w_mob, pl_w, nrow = 2, labels = letters, label_size = 10)
ggsave("./figures/lockdown_analysis_2.png", width = 24, height = 18, units = "cm")
ggsave("./figures/lockdown_analysis_2.pdf", width = 24, height = 18, units = "cm", useDingbats = F)

cowplot::plot_grid(pl_ni, pl_w, nrow = 2, labels = letters, label_size = 10)

ggsave("./figures/lockdown_analysis.pdf", width = 18, height = 18, units = "cm", useDingbats = F)
ggsave("./figures/lockdown_analysis.png", width = 18, height = 18, units = "cm")




# Scotland
iso3166 = fread("./data/county_iso3166_2_gb.csv")

scot_iso3166 = iso3166[gss_code %like% "^S" & iso_code != "", unique(iso_code)];
gm_s = gm[iso_3166_2_code %in% scot_iso3166, lapply(.SD, function(x) mean(x, na.rm = TRUE)), .SDcols = 9:14, by = .(date)]
gm_s = melt(gm_s, id.vars = "date")
gm_s[, variable := str_remove_all(variable, "_percent_change_from_baseline")]

# Plot Scotland google mobility data
ggplot(gm_s) + 
    geom_line(aes(x = date, y = rollmean(value, k = 7, fill = "extend"), colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-09"))) +
    geom_vline(aes(xintercept = ymd("2020-10-25"))) +
    facet_wrap(~variable) +
    theme(legend.position = "none")

# Compare post-lockdown values to pre-lockdown values (2 week lookback)
compare = gm_s[date %between% c(ymd("2020-10-10", "2020-10-24")) & variable != "parks"]
compare = merge(compare,
    gm_s[date %between% c(ymd("2020-10-10", "2020-10-24") - 14) & variable != "parks", 
        .(date = date + 14, variable, baseline = value)],
    by = c("date", "variable"))
compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, mean(change), by = variable]

# IMPACT OF LOCKDOWN FOR SCOTLAND:
#                 variable         V1
# 1:  grocery_and_pharmacy  0.5931034
# 2:           residential  2.0497126
# 3: retail_and_recreation -7.0221675
# 4:      transit_stations -6.6267492
# 5:            workplaces -8.2720801

# Visual comparison
ggplot(compare) +
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    facet_wrap(~variable)




# NEW ADDITION 16 Nov 2020: England
# England
iso3166 = fread("./data/county_iso3166_2_gb.csv")

en_iso3166 = iso3166[gss_code %like% "^E", unique(iso_code)];
gm_e = gm[iso_3166_2_code %in% en_iso3166, lapply(.SD, function(x) mean(x, na.rm = TRUE)), .SDcols = 9:14, by = .(date)]
gm_e = melt(gm_e, id.vars = "date")
gm_e[, variable := str_remove_all(variable, "_percent_change_from_baseline")]
gm_e[, rollval := rollmean(value, k = 7, fill = "extend"), by = variable]


# Plot England google mobility data
pl_e_mob = ggplot(gm_e) + 
    geom_line(aes(x = date, y = value, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-11-05")), linetype = "22", size = 0.3) +
    facet_wrap(~variable) +
    theme(legend.position = "none") +
    labs(x = "Date", y = "Rolling mean - England")

# Compare post-lockdown values to pre-lockdown values (2 week lookback)
compare = gm_e[date %between% c(ymd("2020-11-05", "2020-11-13"))]
compare = merge(compare,
    gm_e[date %between% c(ymd("2020-11-05", "2020-11-13") - 21), 
        .(date = date + 21, variable, baseline = value)],
    by = c("date", "variable"))
compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, mean(change), by = variable]
compare[, variable := str_replace_all(variable, "_", " ")]
compare[, variable := str_to_sentence(variable)]
compare[, variable := factor(variable)]

summary = compare[, mean(change), by = variable]


# IMPACT OF LOCKDOWN FOR ENGLAND:
#                 variable         V1

# Visual comparison
pl_e = ggplot(compare) +
    geom_hline(aes(yintercept = 0), linetype = "23") + 
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    geom_label(data = summary, aes(x = ymd("2020-11-07") + 0.5, y = -10, 
        label = paste0(ifelse(V1 > 0, "+", ""), round(V1, 2))), alpha = 0.75) +
    facet_wrap(~variable) +
    labs(x = "Date", y = "Change in Google Mobility index")

cowplot::plot_grid(pl_ni_mob, pl_ni, pl_w_mob, pl_w, nrow = 2, labels = letters, label_size = 10)
ggsave("./figures/lockdown_analysis_2.png", width = 24, height = 18, units = "cm")
ggsave("./figures/lockdown_analysis_2.pdf", width = 24, height = 18, units = "cm", useDingbats = F)


summ_all = cbind(summaryN, summaryW[, 2], summary[, 2])
names(summ_all) = c("Indicator", "NI change", "Wales change", "England change")
summ_all

library(data.table)
library(ggplot2)
library(zoo)
library(mgcv)
library(lubridate)
library(stringr)
library(cowplot)
library(qs)

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
gm_ni = gm_ni[, .(value = mean(value, na.rm = T)), keyby = .(date, variable)]

# Plot NI google mobility data
ggplot(gm_ni) + 
    geom_line(aes(x = date, y = rollmean(value, k = 7, fill = "extend"), colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-16"))) +
    facet_wrap(~variable) +
    theme(legend.position = "none")

# Compare post-lockdown values to pre-lockdown values (2 week lookback)
compare = gm_ni[date %between% c(ymd("2020-10-17", "2020-10-27")) & variable != "parks"]
compare = merge(compare,
    gm_ni[date %between% c(ymd("2020-10-17", "2020-10-27") - 14) & variable != "parks", 
        .(date = date + 14, variable, baseline = value)],
    by = c("date", "variable"))
compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, mean(change), by = variable]

# IMPACT OF LOCKDOWN FOR NORTHERN IRELAND:
#                 variable          V1
# 1:  grocery_and_pharmacy  -0.6363636
# 2:           residential   4.6299357
# 3: retail_and_recreation -14.1074380
# 4:      transit_stations  -8.7707071
# 5:            workplaces -13.3636364

# Visual comparison
ggplot(compare) +
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    facet_wrap(~variable)

# Wales
iso3166 = fread("./data/county_iso3166_2_gb.csv")

wales_iso3166 = iso3166[gss_code %like% "^W", unique(iso_code)];
gm_w = gm[iso_3166_2_code %in% wales_iso3166, lapply(.SD, function(x) mean(x, na.rm = TRUE)), .SDcols = 9:14, by = .(date)]
gm_w = melt(gm_w, id.vars = "date")
gm_w[, variable := str_remove_all(variable, "_percent_change_from_baseline")]

# Plot Wales google mobility data
ggplot(gm_w) + 
    geom_line(aes(x = date, y = rollmean(value, k = 7, fill = "extend"), colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-16"))) +
    facet_wrap(~variable) +
    theme(legend.position = "none")

# Compare post-lockdown values to pre-lockdown values (2 week lookback)
compare = gm_w[date %between% c(ymd("2020-10-23", "2020-10-27")) & variable != "parks"]
compare = merge(compare,
    gm_w[date %between% c(ymd("2020-10-23", "2020-10-27") - 7) & variable != "parks", 
        .(date = date + 7, variable, baseline = value)],
    by = c("date", "variable"))
compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, mean(change), by = variable]

# IMPACT OF LOCKDOWN FOR WALES:
#                 variable         V1
# 1:  grocery_and_pharmacy -16.086580
# 2:           residential   5.699122
# 3: retail_and_recreation -32.978788
# 4:      transit_stations -17.884917
# 5:            workplaces -18.981818


# Visual comparison
ggplot(compare) +
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    facet_wrap(~variable)

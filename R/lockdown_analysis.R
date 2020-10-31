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

iso3166 = fread("./data/county_iso3166_2_gb.csv")

ggplot(gm_melt(gm[country_region_code == "IE" & sub_region_1 == ""])) + 
    geom_line(aes(x = date, y = value, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-21")))

ggplot(gm_melt(gm[country_region_code == "GB" & sub_region_1 == "Edinburgh"])) + 
    geom_line(aes(x = date, y = value, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-09")))

wales_iso3166 = iso3166[gss_code %like% "^W", unique(iso_code)];
wales_gm = gm[iso_3166_2_code %in% wales_iso3166, lapply(.SD, function(x) mean(x, na.rm = TRUE)), .SDcols = 9:14, by = .(date)]
wales_gm = melt(wales_gm, id.vars = "date")
wales_gm[, variable := str_remove_all(variable, "_percent_change_from_baseline")]

ggplot(wales_gm) + 
    geom_line(aes(x = date, y = rollmean(value, 7, fill = "extend"), colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-23")))

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
gm_ni = gm_ni[, .(value = mean(value, na.rm = T)), by = .(date, variable)]

ggplot(gm_ni) + 
    geom_line(aes(x = date, y = value, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-16"))) +
    facet_wrap(~variable) +
    theme(legend.position = "none")

# Add GSS code
gm = gm[country_region_code == "GB"]
iso3166 = fread("./data/county_iso3166_2_gb.csv")
gm = merge(gm, iso3166[iso_code != "", .(iso_code, gss_code)], by.x = "iso_3166_2_code", by.y = "iso_code", all.x = TRUE)



# Google places
gp = rbind(
    gm[sub_region_1 != "", .(place = unique(sub_region_1), type = "sub_region_1")],
    gm[sub_region_2 != "", .(place = unique(sub_region_2), type = "sub_region_2")]
)

# Tiers data
tiers = fread("./data/tiers.txt")

# Interactive match finder
imatch = function(strings, table, table_attr)
{
    match_table = data.table(str = table, attr = table_attr);
    match_table[, match_str := str_replace_all(str, "&", "and")];
    
    result = NULL;
    
    for (string in strings)
    {
        match_table[, score := 0];
    
        s2 = str_replace_all(string, "&", "and");
        s2a = str_split(s2, " | ")[[1]]
        s2b = s2a[2];
        s2a = s2a[1];
        
        match_table[str %ilike% paste0("^", s2a, "$"), score := score + 100];
        match_table[str %ilike% paste0("^", s2b, "$"), score := score + 100];
        words = str_split(s2, " ")[[1]];
        for (w in words) {
            match_table[str %ilike% w, score := score + str_length(w)];
        }
        
        best = match_table[order(-score)][1:10, .(str, attr, score)];
        
        cat(paste0(string, ":\n"));
        cat(paste0(1:10, ": ", best$str, " (", best$attr, ")", collapse = "\n"));
        cat("\n0: No match\n");
        
        a = as.numeric(readline(prompt = "Which match is best? "));
        
        result = rbind(result, data.table(string = string, match = best$str[a], attr = best$attr[a]));
    }
    
    return (result)
}

# Create match table
# result = imatch(paste(tiers$reg1, "|", tiers$reg2), gp$place, gp$type)
# result[string == "North Yorkshire | York", match := "York"]
# result[string == "North Yorkshire | York", attr := "sub_region_1"]
# fwrite(result, "./data/tiers_match_table.txt")

# Update match table
result = fread("./data/tiers_match_table.txt")
tiers[, reg12 := paste(reg1, reg2, sep = " | ")]
missing = tiers[!reg12 %in% result[, string], reg12]
updates = imatch(missing, gp$place, gp$type);
result = rbind(result, updates);
fwrite(result, "./data/tiers_match_table.txt")

tiers = merge(tiers, result, by.x = "reg12", by.y = "string", all = T)

# Create tiers timeline
changes = tiers[, .(gp = match, type = attr, date = as.Date(date), tier)]
changes = rbind(gp[, .(gp = place, type, date = ymd("2020-10-14"), tier = 1)], changes)

# Remove first changes on same date
changes = changes[!duplicated(changes[, .(gp, type, date)], fromLast = TRUE)]

# Add "tier 0"
changes = rbind(gp[, .(gp = place, type, date = ymd("2020-09-01"), tier = 0)], changes)

# Put in order
changes = changes[order(date, gp, type)]

# Merge into original gm
gm1 = gm[sub_region_1 != "" & sub_region_2 == ""]
gm2 = gm[sub_region_1 != "" & sub_region_2 != ""]

gm1 = merge(gm1, changes[type == "sub_region_1", .(gp, date, tier)], by.x = c("sub_region_1", "date"), by.y = c("gp", "date"), all.x = TRUE);
gm2 = merge(gm2, changes[type == "sub_region_2", .(gp, date, tier)], by.x = c("sub_region_2", "date"), by.y = c("gp", "date"), all.x = TRUE);

gm = rbindlist(list(gm[sub_region_1 == ""], gm1, gm2), fill = TRUE);
gm = gm[order(country_region, sub_region_1, sub_region_2, date)];

# Remove "all UK"
gm = gm[sub_region_1 != ""]

# Fill out tiers
gm[, tier := na.locf(tier, na.rm = FALSE), by = .(sub_region_1, sub_region_2)]

# Remove pre Tier 0
gm = gm[!is.na(tier)]

# Restrict to England
regions = fread("./data/goog_regions.txt")
english_places = regions[gss_code %like% "^E", .(goog_name)]
gm = gm[sub_region_1 %in% english_places$goog_name]

# Analysis setup
gm[, full_region := paste(sub_region_1, "|", sub_region_2)]
gm[, weekday := factor(wday(date, label = TRUE), ordered = FALSE)]
gm[, weekday_t := wday(date)]
gm[, date_t := as.numeric(ymd(date) - ymd("2020-01-01"))]
gm[, week_t := week(date)]

indicators = names(gm)[names(gm) %like% "_percent_change"]


# ANALYSIS

regions = gm[, unique(full_region)]
gm2 = rlang::duplicate(gm)
for (col in 9:14) set(gm2, j = col, value = as.double(gm2[[col]]))

# Remove weekday effect from each trendline
results = NULL
for (reg in regions) {
    cat(".")
    for (ind in indicators) {
        model = try(gam(get(ind) ~ s(date_t) + s(weekday_t, bs = "cp", k = 7), 
            knots = list(weekday_t = c(0, 7)), data = gm2[full_region == reg]), silent = TRUE)
        if ("try-error" %in% class(model)) {
            next;
        }
        pp = predict(model, type = "terms")
        results = rbind(results,
            gm2[full_region == reg & !is.na(get(ind)), 
                .(sub_region_1, sub_region_2, date, date_t, weekday, tier, indicator = ind, 
                    original = get(ind), unweekday = get(ind) - pp[, "s(weekday_t)"])]
        )
    }
}

### Remove weekday effect (the old-fashioned way)
#results[, unweekday2 := rollmean(original, k = 7, fill = "extend"), by = .(sub_region_1, sub_region_2, indicator)]

# Remove national trend from trendlines
for (ind in indicators) {
    cat(".")
    model = gam(unweekday ~ s(date_t), data = results[indicator == ind])
    results[indicator == ind, detrend := unweekday - predict(model, newdata = .SD)];
}

# Factorise
results[, full_region := factor(paste(sub_region_1, sub_region_2))]
results[, tier := factor(paste("T", tier))]
results[, indicator := factor(indicator)]

# - detrended, compared to tier 0
baselineD = results[tier == "T 0", .(baseline = mean(tail(detrend, 7))), by = .(full_region, indicator)]
resultsD = merge(results, baselineD, by = c("full_region", "indicator"))
resD = resultsD[tier != "T 0", .(detrend = mean(detrend), baseline = mean(baseline), change = mean(detrend - baseline)), 
    by = .(full_region, sub_region_1, indicator, tier)]
ggplot(resD) +
    geom_jitter(aes(x = tier, y = change, colour = as.factor(tier)), width = 0.4, height = 0) +
    facet_wrap(~indicator, scales = "free_y") +
    theme(legend.position = "none")
summaryD = resD[, .(change = mean(change)), keyby = .(indicator, tier)]
summaryD[, T0 := first(change), by = indicator]
summaryD[, change - T0, keyby = .(indicator, tier)]

# Plots for the above
resultsD[, indic := str_remove_all(indicator, "_percent_change_from_baseline")]
resultsD[, indic := str_replace_all(indic, "_", " ")]
resultsD[, indic := str_to_sentence(indic)]
resultsD[, indic := factor(indic)]

theme_set(cowplot::theme_cowplot(font_size = 8) + theme(strip.background = element_blank()))

# Start with GM trends
ggplot(resultsD[indic != "Parks"]) +
    geom_line(aes(x = date, y = original, group = full_region), size = 0.1) +
    facet_wrap(~indic) +
    ylim(-80, 90) +
    labs(x = "Date", y = "Google Mobility indices")
ggsave("./figures/tiers_analysis_1.png", width = 25, height = 15, units = "cm")

# Subtract the weekday effect individually
ggplot(resultsD[indic != "Parks"]) +
    geom_line(aes(x = date, y = unweekday, group = full_region), size = 0.1) +
    facet_wrap(~indic) +
    ylim(-80, 90) +
    labs(x = "Date", y = "Google Mobility indices, minus weekday effect")
ggsave("./figures/tiers_analysis_2.png", width = 25, height = 15, units = "cm")

# Subtract the national trend
ggplot(resultsD[indic != "Parks"]) +
    geom_line(aes(x = date, y = detrend, group = full_region), size = 0.1) +
    facet_wrap(~indic) +
    ylim(-80, 90) +
    labs(x = "Date", y = "Google Mobility indices, minus weekday effect, detrended")
ggsave("./figures/tiers_analysis_3.png", width = 25, height = 15, units = "cm")

# Comparison
ggplot(resultsD[indic != "Parks"]) +
    geom_ribbon(data = data.table(x = ymd(c("2020-10-07", "2020-10-13"))), aes(x = x, ymin = -80, ymax = 90), fill = "#eeeeee") +
    geom_line(data = resultsD[indic != "Parks" & tier == "T 0"], 
        aes(x = date, y = detrend, group = paste(tier, full_region)), size = 0.1) +
    geom_line(data = resultsD[indic != "Parks" & tier != "T 0"], 
        aes(x = date, y = detrend, group = paste(tier, full_region), colour = tier), size = 0.5) +
    scale_colour_manual(values = c("#4499cc", "#ee99ee", "#ee4433")) +
    facet_wrap(~indic) +
    ylim(-80, 90) +
    theme(legend.position = c(0.8, 0.2)) +
    labs(x = "Date", y = "Google Mobility indices, minus weekday effect, detrended", colour = "Tier")
ggsave("./figures/tiers_analysis_4.png", width = 25, height = 15, units = "cm")

# Impact of tiers
individual = resultsD[tier != "T 0", .(detrend = mean(detrend), baseline = mean(original), change = mean(detrend - baseline)), 
    by = .(full_region, sub_region_1, indic, tier)]
summary = individual[, .(change = mean(change)), keyby = .(indic, tier)]
summary[, T0 := first(change), by = indic]
summary = summary[, .(change = change - T0), keyby = .(indic, tier)]
summary

ggplot(individual[indic != "Parks"]) +
    geom_violin(aes(x = tier, y = change, fill = as.factor(tier)), alpha = 0.25, colour = NA) +
    geom_jitter(aes(x = tier, y = change, colour = as.factor(tier)), width = 0.4, height = 0, size = 0.75) +
    geom_label(data = summary[indic != "Parks"], aes(x = tier, y = change, label = round(change, 2)), alpha = 0.75) +
    facet_wrap(~indic, scales = "free") +
    scale_colour_manual(values = c("#4499cc", "#ee99ee", "#ee4433"), aesthetics = c("colour", "fill")) +
    theme(legend.position = c(0.8, 0.2)) +
    labs(x = "Tier", y = "Change in Google Mobility index", fill = "Tier", colour = "Tier")
ggsave("./figures/tiers_analysis_5.png", width = 25, height = 15, units = "cm")

# Is there a change over time?
bytime = resultsD[tier != "T 0", .(detrend = detrend, baseline = mean(baseline), change = detrend - mean(baseline)), 
    keyby = .(full_region, sub_region_1, indic, tier, date, weekday)]
bytime[, days := 1:.N, by = .(full_region, sub_region_1, indic, tier)]

ggplot(bytime[indic != "Parks" & indic != "Residential"]) +
    geom_jitter(aes(x = days, y = change, colour = tier)) +
    #geom_smooth(aes(x = days, y = change, colour = tier, group = tier)) +
    facet_wrap(~indic)



# How does this compare to week on week changes?
resultsC = results[date >= "2020-09-02" & date <= "2020-10-20"]
resultsC[, week := factor(paste("W", as.numeric(as.Date(date) - ymd("2020-09-02")) %/% 7 + 1))]
baselineC = resultsC[week == "W 1", .(baseline = mean(unweekday)), by = .(full_region, indicator)]
resultsC = merge(resultsC, baselineC, by = c("full_region", "indicator"))
resultsC[, indic := str_remove_all(indicator, "_percent_change_from_baseline")]
resultsC[, indic := str_replace_all(indic, "_", " ")]
resultsC[, indic := str_to_sentence(indic)]
resultsC[, indic := factor(indic)]
individualC = resultsC[, .(unweekday = mean(unweekday), baseline = mean(baseline), change = mean(unweekday - baseline)), 
    by = .(full_region, sub_region_1, indic, week)]
summaryC = individualC[, .(change = mean(change)), keyby = .(indic, week)]
summaryC[, T0 := first(change), by = indic]
summaryC = summaryC[, .(change = change - T0), keyby = .(indic, week)]

ggplot(individualC[indic != "Parks"]) +
    geom_violin(aes(x = week, y = change, fill = week), alpha = 0.25, colour = NA, scale = "width") +
    geom_jitter(aes(x = week, y = change, colour = week), width = 0.4, height = 0, size = 0.75) +
    geom_label(data = summaryC[indic != "Parks"], aes(x = week, y = change, label = round(change, 2)), alpha = 0.75, size = 3) +
    facet_wrap(~indic, scales = "free") +
    theme(legend.position = c(0.8, 0.2)) +
    labs(x = "Week", y = "Change in Google Mobility index", fill = "Week", colour = "Week")
ggsave("./figures/tiers_analysis_6.png", width = 25, height = 15, units = "cm")





# OLD MATERIAL - For reference.




# Add lat & lon
latlon = fread("./data/iso-3166-2-latlon.txt")
gm = merge(gm, latlon[, .(code, lat, lon)], by.x = "iso_3166_2_code", by.y = "code")
gm = gm[!is.na(lat)]

# Add rural/urban
ruc = fread("./data/Rural_Urban_Classification_(2011)_for_Local_Authority_Districts_in_England/RUC_LAD_2011_EN_LU.csv")
ruc[, p_urban := `Total Urban population 2011` / `Total population 2011`]
wards = fread("~/Dropbox/uk_covid_data/fitting/data/Ward_to_Local_Authority_District_to_County_to_Region_to_Country_December_2019_Lookup_in_United_Kingdom.csv")
counties = wards[CTY19CD %like% "^E10", unique(CTY19CD)]
for (county in counties) {
    lads = wards[CTY19CD == county, unique(LAD19CD)]
    p_urban = ruc[LAD18CD %in% lads | LAD11CD %in% lads, sum(`Total Urban population 2011`) / sum(`Total population 2011`)]
    ruc = rbind(ruc,
        data.table(LAD11CD = county, p_urban = p_urban),
        fill = TRUE)
}
ruc = rbind(ruc,
    data.table(LAD11CD = "E10000009", 
        p_urban = ruc[LAD18NM %like% "Dorset", sum(`Total Urban population 2011`) / sum(`Total population 2011`)]),
    fill = TRUE)
gm = merge(gm, ruc[, .(LAD11CD, p_urban)], by.x = "gss_code", by.y = "LAD11CD", all.x = TRUE)





# Analysis . . .
models = list()
for (ind in indicators) {
    models[[ind]] = gam(get(ind) ~ as.factor(tier) + s(weekday_t, bs = "cp") + date_t:sub_region_1 + sub_region_1 + sub_region_2, data = gm)
}


for (ind in indicators) {
    cat(paste0(ind, ":\n"))
    s = summary(models[[ind]])
    cat(paste("Tier 1", "d = ", signif(s$p.coeff[2], 3), "P = ", signif(s$p.pv[2], 3), "\n", sep = "\t"))
    cat(paste("Tier 2", "d = ", signif(s$p.coeff[3], 3), "P = ", signif(s$p.pv[3], 3), "\n", sep = "\t"))
    cat(paste("Tier 3", "d = ", signif(s$p.coeff[4], 3), "P = ", signif(s$p.pv[4], 3), "\n", sep = "\t"))
    cat("\n")
}

# V2 of analysis...

models2 = list()
for (ind in indicators) {
    models2[[ind]] = gam(get(ind) ~ s(weekday_t, bs = "cp", k = 7) + s(date_t) + date_t:sub_region_1 + sub_region_1 + sub_region_2, data = gm[tier <= 1])
}

for (ind in indicators) {
    cat(paste(ind, "\n"))
    
    terms = names(models2[[ind]]$coefficients);
    sr1 = str_remove(terms[terms %like% "^sub_region_1"], "^sub_region_1")
    sr2 = str_remove(terms[terms %like% "^sub_region_2"], "^sub_region_2")
    sr1 = str_remove(sr1, ":.*$")
    sr2 = str_remove(sr2, ":.*$")
    sr1 = unique(sr1)
    sr2 = c(unique(sr2), "")
    
    gmo = gm[sub_region_1 %in% sr1 & sub_region_2 %in% sr2]
    
    gmo[, predicted := predict(models2[[ind]], newdata = gmo)]
    print(gmo[, .(diff = mean(get(ind) - predicted, na.rm = T)), keyby = tier])
    cat("\n")
}

# V3 of analysis...
regions = gm[, unique(full_region)]
results = NULL
for (reg in regions) {
    gmo = gm[full_region == reg & date >= "2020-09-15"]
    for (ind in indicators[-3]) {
        if (gmo[, sum(!is.na(get(ind)))] < 28) {
            next;
        }
        model = gam(get(ind) ~ s(weekday_t, bs = "cp", k = 7) + s(date_t, k = 4), knots = list(weekday_t = c(0, 7)), data = gmo[tier == 0])
        #model = gam(get(ind) ~ s(weekday_t, bs = "cp", k = 7) + date_t, knots = list(weekday_t = c(0, 7)), data = gmo[tier == 0])
        gmo[, predicted := predict(model, newdata = gmo)]
        
        results = rbind(results,
            gmo[, .(sub_region_1, sub_region_2, date, weekday, tier, indicator = ind, measured = get(ind), predicted)]
        )
    }
}

results[, median(measured - predicted, na.rm = T), keyby = .(indicator, tier)]

results[, mean(measured - predicted, na.rm = T), keyby = .(indicator, tier, weekday)][, mean(V1), keyby = .(indicator, tier)]

ggplot(results[, median(measured - predicted, na.rm = T), keyby = .(indicator, tier)]) + 
    geom_point(aes(tier, V1)) + 
    facet_wrap(~indicator)

ggplot(results) + 
    geom_violin(aes(as.factor(tier), measured - predicted)) + 
    facet_wrap(~indicator)

# by day of week . . .
ggplot(results[sub_region_1 != "Greater Manchester", median(measured - predicted, na.rm = T), keyby = .(indicator, tier, weekday)]) + 
    geom_line(aes(weekday, V1, group = as.factor(tier), colour = as.factor(tier))) + 
    geom_point(aes(weekday, V1, group = as.factor(tier), colour = as.factor(tier))) + 
    facet_wrap(~indicator)

reg
ggplot(results) +
    geom_line(aes(x = date, y = measured, colour = "measured")) +
    geom_line(aes(x = date, y = predicted, colour = "predicted")) +
    facet_wrap(~indicator)


# Direct weekday comparison...

gmm = melt(gm, id.vars = c("date", "full_region", "tier", "weekday"), measure.vars = 9:14)
baseline = gmm[tier <= 1, .(date = tail(date, 7), weekday = tail(weekday, 7), baseline = tail(value, 7)), by = .(full_region, variable)]
comparison = gmm[tier >= 2, .(comparison = mean(value, na.rm = T)), by = .(full_region, variable, weekday, tier)]
merged = merge(baseline, comparison, by = c("full_region", "variable", "weekday"))
merged = merged[!is.na(baseline)]

merged[, mean(comparison - baseline, na.rm = T), by = .(variable, weekday, tier)][, mean(V1), by = .(variable, tier)]

ggplot(merged) +
    geom_violin(aes(x = as.factor(tier), y = comparison - baseline)) +
    facet_wrap(~variable)

merged[, unique(tier)]



ind = indicators[5]
gmo = gm[sub_region_1 == "Nottingham" & date > "2020-09-23"]
models4 = gam(get(ind) ~ s(weekday_t, bs = "cp", k = 7) + s(date_t, k = 10), knots = list(weekday_t = c(0, 7)), data = gmo[tier == 0])
plot(models4)


gmo[, mean(get(ind) - predicted, na.rm = T), by = tier]


for (ind in indicators[4]) {
    cat(paste(ind, "\n"))
    
    terms = names(models3[[ind]]$coefficients);
    sr1 = str_remove(terms[terms %like% "^sub_region_1"], "^sub_region_1")
    sr2 = str_remove(terms[terms %like% "^sub_region_2"], "^sub_region_2")
    sr1 = str_remove(sr1, ":.*$")
    sr2 = str_remove(sr2, ":.*$")
    sr1 = unique(sr1)
    sr2 = c(unique(sr2), "")
    
    gmo = gm[sub_region_1 %in% sr1 & sub_region_2 %in% sr2]
    
    gmo[, predicted := predict(models3[[ind]], newdata = gmo)]
    print(gmo[, .(diff = mean(get(ind) - predicted, na.rm = T)), keyby = tier])
    cat("\n")
    
    print(
    ggplot(gmo) + 
        geom_line(aes(x = date, y = get(ind) - predicted, group = full_region, colour = tier)) + 
        geom_line(aes(x = date, y = predicted, group = full_region)) + 
        facet_wrap(~tier) + 
        geom_smooth(aes(x = date, y = get(ind) - predicted, group = tier), method = "lm")
    )
}


# Spatial analysis...

# TEST
# model = gam(get(ind) ~ s(weekday_t, bs = "cp", k = 7) + s(date_t) + s(lon, lat), 
#     knots = list(weekday_t = c(0, 7)), data = gm[tier < 1 & lat > 45 & full_region != "Rutland | "])
# model = gam(get(ind) ~ s(weekday_t, bs = "cp", k = 7) + s(date_t, lon, lat), 
#     knots = list(weekday_t = c(0, 7)), data = gm[tier < 1 & lat > 45 & full_region != "Rutland | "])
# model = gam(get(ind) ~ s(weekday_t, bs = "cp", k = 7) + s(date_t, lat), 
#     knots = list(weekday_t = c(0, 7)), data = gm[tier < 1 & lat > 45 & full_region != "Rutland | "])
# model = gam(get(ind) ~ s(date_t, lon, lat), 
#     data = gm[tier < 1 & lat > 45 & full_region != "Rutland | "])
# model = gam(get(ind) ~ s(weekday_t, bs = "cp", k = 7) + s(lon, lat, by = date_t), 
#     knots = list(weekday_t = c(0, 7)), data = gm[tier < 1 & lat > 45 & full_region != "Rutland | "])

modelsS = list()
for (ind in indicators) {
    cat(".")
    modelsS[[ind]] = gam(get(ind) ~ as.factor(tier) + s(weekday_t, bs = "cp", k = 7) + s(date_t, lon, lat) + full_region, 
        knots = list(weekday_t = c(0, 7)), data = gm[lat > 45 & full_region != "Rutland | "])
}
indicators
summary(modelsS[[6]])

results = NULL
for (ind in indicators[-3]) {
    cat(".")
    gm[lat > 45 & full_region != "Rutland | ", predicted := predict(modelsS[[ind]], newdata = gm[lat > 45 & full_region != "Rutland | "])]
    results = rbind(results,
        gm[, .(full_region, date, weekday, tier, indicator = ind, predicted, observed = get(ind))]
    )
}

ggplot(results) +
    geom_histogram(aes(x = observed - predicted)) +
    facet_wrap(~indicator)

ggplot(results[full_region %like% "London"]) +
    geom_line(aes(x = date, y = predicted, colour = "predicted")) +
    geom_line(aes(x = date, y = observed, colour = "observed")) +
    facet_grid(indicator ~ full_region)

results[, mean(observed - predicted, na.rm = T), keyby = .(indicator, tier)]

ggplot(gm) + geom_line(aes(x = date, y = predicted, group = full_region))

ggplot(gm) +
    geom_line(aes(x = date, y = residential_percent_change_from_baseline - predicted, colour = tier, group = full_region))

ggplot(gm[tier > 0]) +
    geom_point(aes(y = residential_percent_change_from_baseline, x = predicted, colour = as.factor(tier)))





resD = resultsD[tier != "T 0", .(detrend = mean(detrend), baseline = mean(baseline), change = mean(detrend - baseline)), 
    by = .(full_region, sub_region_1, indicator, tier)]
ggplot(resD) +
    geom_jitter(aes(x = tier, y = change, colour = as.factor(tier)), width = 0.4, height = 0) +
    facet_wrap(~indicator, scales = "free_y") +
    theme(legend.position = "none")
summaryD = resD[, .(change = mean(change)), keyby = .(indicator, tier)]
summaryD[, T0 := first(change), by = indicator]
summaryD[, change - T0, keyby = .(indicator, tier)]



ggplot(individual[indic != "Parks" & tier == "T 3"]) +
    geom_jitter(aes(x = tier, y = change, colour = sub_region_1), width = 0.2, height = 0) +
    facet_wrap(~indic) +
    theme(legend.position = c(0.7, 0.2)) +
    labs(x = "Tier", y = "Change in Google Mobility index")

summaryD = resD[, .(change = mean(change)), keyby = .(indicator, tier)]
summaryD[, T0 := first(change), by = indicator]
summaryD[, change - T0, keyby = .(indicator, tier)]

    scale_colour_manual(values = c("#4499cc", "#ee99ee", "#ee4433")) +


# - detrended, compared to tier 0/1
baselineD1 = results[tier %in% c("T 0", "T 1"), .(baseline = mean(head(tail(detrend, 10), 7))), by = .(full_region, indicator)]
resultsD1 = merge(results, baselineD1, by = c("full_region", "indicator"))
resD1 = resultsD1[tier != "T 0", mean(detrend - baseline), by = .(full_region, indicator, tier)]
ggplot(resD1) +
    geom_jitter(aes(x = tier, y = V1, colour = as.factor(tier)), width = 0.4, height = 0) +
    facet_wrap(~indicator, scales = "free_y") +
    theme(legend.position = "none")
summaryD1 = resD1[, mean(V1), keyby = .(indicator, tier)]
summaryD1


# - not detrended, compared to tier 0
baselineU = results[tier == "T 0", .(baseline = mean(tail(unweekday, 7))), by = .(full_region, indicator)]
resultsU = merge(results, baselineU, by = c("full_region", "indicator"))
resU = resultsU[tier != "T 0", mean(unweekday - baseline), by = .(full_region, sub_region_1, indicator, tier)]
ggplot(resU) +
    geom_jitter(aes(x = tier, y = V1, colour = as.factor(tier)), width = 0.4, height = 0) +
    facet_wrap(~indicator, scales = "free_y") +
    theme(legend.position = "none")
resU[, mean(V1), keyby = .(indicator, tier)]
summaryU = resU[, mean(V1), keyby = .(indicator, tier)]
summaryU[, T0 := first(V1), by = indicator]
summaryU[, V1 - T0, keyby = .(indicator, tier)]


# ugggghh 2
baseline = results[tier %in% c("T 0", "T 1"), .(baseline = mean(tail(unweekday, 7))), by = .(full_region, indicator)]
results = merge(results, baseline, by = c("full_region", "indicator"))
res3 = results[, mean(unweekday - baseline), by = .(full_region, indicator, tier)]
ggplot(res3) +
    geom_jitter(aes(x = tier, y = V1, colour = as.factor(tier)), width = 0.4, height = 0) +
    facet_wrap(~indicator, scales = "free_y") +
    theme(legend.position = "none")
res3[, mean(V1), keyby = .(indicator, tier)]




results = results[lon > -60]
ggplot(results) + geom_point(aes(lon, lat))

coeffs = NULL
for (ind in indicators) {
    cat(".")
    oldtier = rlang::duplicate(results$tier)

    #results[tier %in% c("T 0", "T 1"), tier := "T 0"]
    model = gam(unweekday ~ tier + s(lat, lon, date_t) + full_region, data = results[indicator == ind])
    results[indicator == ind, predicted := predict(model, newdata = .SD)]
    results[, tier := "T 0"]
    results[indicator == ind, predicted_tier0 := predict(model, newdata = .SD)]
    results[, tier := oldtier]
    
    coeffs = rbind(coeffs, 
        data.table(indicator = ind, tier = c(1, 2, 3), coeff = model$coefficients[2:4]))
}

plot(model)

coeffs
ggplot(results[indicator == "retail_and_recreation_percent_change_from_baseline"]) + 
    geom_line(aes(x = date_t, y = unweekday, colour = "unweekday")) + 
    geom_line(aes(x = date_t, y = predicted_tier0, colour = "predictedT0")) + 
    facet_wrap(~full_region) + cowplot::theme_cowplot(font_size = 6)

res2 = results[, mean(unweekday - predicted_tier0), by = .(full_region, tier, indicator)][tier != "T 0"]





results[, mean(unweekday), keyby = .(indicator, tier)][tier != "T 0"]

ggplot(res2) +
    geom_jitter(aes(x = tier, y = V1, colour = as.factor(tier)), width = 0.4, height = 0) +
    facet_wrap(~indicator, scales = "free_y") +
    theme(legend.position = "none")

results[, mean((unweekday - predicted)^2, na.rm = T), by = .(full_region, indicator)][order(indicator, -V1)]
View(results[, mean((unweekday - predicted)^2, na.rm = T), by = .(full_region, indicator)][order(indicator, -V1)])

ggplot(results) + geom_point(aes(x = original, y = unweekday, colour = weekday)) + facet_wrap(~indicator)
ggplot(results) + geom_point(aes(x = unweekday, y = predicted, colour = tier)) + facet_wrap(~indicator)
ggplot(results[tier != "T 0"]) + geom_point(aes(x = unweekday, y = predicted, colour = tier)) + facet_wrap(~indicator)

ggplot(results) + geom_point(aes(x = predicted, y = predicted_tier0))
ggplot(results) + geom_histogram(aes(unweekday - predicted)) + facet_wrap(~indicator)


#model = gam(unweekday ~ tier_f + s(date_t, by = full_region), data = results[indicator == indicators[1]])


plot(model)
indicators
summary(modelsT[[6]])


results = fread(
"context	tier	diff
retail_and_recreation	0	-0.05752478
retail_and_recreation	1	0.55583373
retail_and_recreation	2	-2.93082065
retail_and_recreation	3	-9.58820971
grocery_and_pharmacy	0	-0.06584278
grocery_and_pharmacy	1	0.64911179
grocery_and_pharmacy	2	-0.94165802
grocery_and_pharmacy	3	-1.82715398
parks	0	-0.4607227
parks	1	3.7967036
parks	2	-3.9297007
parks	3	-0.7824230
transit_stations	0	-0.07237979
transit_stations	1	0.70670155
transit_stations	2	-2.74575021
transit_stations	3	-6.97454232
workplaces	0	-0.03200341
workplaces	1	0.33949308
workplaces	2	0.10869918
workplaces	3	-7.70573322
residential	0	0.01192354
residential	1	-0.12863987
residential	2	0.12112563
residential	3	1.95181705")

results[, tier := as.character(paste0("tier", tier))]
results = dcast(results, context ~ tier)

results[, tier2 - tier1, by = context]
results[, tier3 - tier1, by = context]

# ^^ these put into tier 1 / tier 2 schedule

gm[tier >= 2, unique(full_region)]
plot(model)
summary(model)
?predict.gam


gm[, max_tier := max(tier), by = sub_region_2]
ggplot(gm) + geom_smooth(aes(x = date, y = retail_and_recreation_percent_change_from_baseline, colour = as.factor(max_tier), group = as.factor(max_tier)))







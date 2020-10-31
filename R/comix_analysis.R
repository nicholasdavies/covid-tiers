library(data.table)
library(stringr)
library(lubridate)
library(zoo)
library(ggplot2)
library(socialmixr)
library(mgcv)
library(qs)

# LOAD COMIX AND GOOGLE MOBILITY
participants = readRDS("../comix/data/clean_participants.rds")
contacts = readRDS("../comix/data/clean_contacts.rds")
gm = qread("./data/google_mobility_uk.qs")
gm = gm[country_region_code == "GB" & sub_region_1 == ""]

# Ensure google and comix data are aligned
participants <- participants[date <= max(gm$date)]
contacts <- contacts[date <= max(gm$date)]
gm <- gm[date <= max(participants$date)]


# CLEAN COMIX DATA

map_minutes_min = c(
    "Less than 5 minutes" = 0,
    "5 minutes or more, but less than 15 minutes" = 5,
    "15 minutes or more, but less than 1 hour" = 15,
    "1 hour or more, but less than 4 hours" = 60,
    "4 hours or more" = 240
)

map_minutes_max = c(
    "Less than 5 minutes" = 5,
    "5 minutes or more, but less than 15 minutes" = 15,
    "15 minutes or more, but less than 1 hour" = 60,
    "1 hour or more, but less than 4 hours" = 240,
    "4 hours or more" = 1440
)

map_type = c(
    "They are a family member who is not in my household" = "family_non_household",
    "They are someone I work with" = "work",
    "They are someone I go to school, college or university with" = "school",
    "They are a friend" = "friend",
    "They are family members not in our household" = "family_non_household",
    "They are a babysitter, childminder, or nanny" = NA,
    "They are someone they go to nursery, pre-school, school, college or university with" = "school",
    "They are friends" = "friend",
    "They are a client or customer" = "work",
    "They are my spouse/girlfriend/boyfriend/partner" = "partner",
    "They are someone they see at nursery, pre-school, school, college or university" = "school"
)

map_freq = c(
    "Never met them before" = "never met",
    "Less often than once per month" = "occasional",
    "About once per month" = "1 month",
    "Every 2-3 weeks" = "2-3 weeks",
    "About once or twice a week" = "3-7 days",
    "Every day or almost every day" = "1-2 days"
)

map_gender = c(
    "Female" = "female",
    "Male" = "male",
    "In another way" = "other"
)

# Map "Yes", "No", NA to logical
YesNoNA = function(x)
{
    ifelse(x == "Yes", TRUE,
        ifelse(x == "No", FALSE, NA))
}

# Map "Yes", "No", NA to 1/0/NA
YesNoNA_Ind = function(x)
{
    ifelse(x == "Yes", 1,
        ifelse(x == "No", 0, NA))
}

# Pull in participant information
parts = participants[, .(
    part_age, 
    part_age_group,
    part_gender
    ), by = .(part_id, panel, wave, date, survey_date)]

# Clean participant ages
parts[, part_age_min := part_age]
parts[, part_age_max := part_age]
parts[is.na(part_age) & !is.na(part_age_group), part_age_min := as.numeric(str_replace_all(part_age_group, "\\[|,[0-9]+\\)", ""))]
parts[is.na(part_age) & !is.na(part_age_group), part_age_max := as.numeric(str_replace_all(part_age_group, "\\[[0-9]+,|\\)", ""))]
parts[, part_age := NULL]
parts[, part_age_group := NULL]

# Clean contact ages
contacts[cnt_age == "Under 1", cont_age_min := 0]
contacts[cnt_age == "Under 1", cont_age_max := 1]
contacts[cnt_age %like% "^[0-9]+\\+$", cont_age_min := as.numeric(str_replace_all(cnt_age, "\\+", ""))]
contacts[cnt_age %like% "^[0-9]+\\+$", cont_age_max := 120]
contacts[cnt_age %like% "^[0-9]+-[0-9]+$", cont_age_min := as.numeric(str_replace_all(cnt_age, "-[0-9]+", ""))]
contacts[cnt_age %like% "^[0-9]+-[0-9]+$", cont_age_max := as.numeric(str_replace_all(cnt_age, "[0-9]+-", ""))]

# Merge contacts and participants
cpmerged = merge(parts, contacts, by = c("part_id", "panel", "wave", "date"), all = TRUE)


# Tidy data
contacts_clean = cpmerged[, .(
    part_id,
    part_age_min,
    part_age_max,
    part_gender = map_gender[part_gender],
    panel = str_remove(panel, "Panel "), 
    wave = as.numeric(str_remove(wave, "Wave ")), 
    date = date, 
    week = lubridate::week(date),
    weekday = factor(lubridate::wday(date, label = TRUE), ordered = FALSE),
    contacts = ifelse(is.na(cont_id) | suspected_non_contact == 1 | suspected_multiple_contact == 1, 0, 1),
    cont_id,
    cont_age_min,
    cont_age_max,
    cont_gender = map_gender[cnt_gender],
    c_home = YesNoNA_Ind(cnt_home),
    c_sport = YesNoNA_Ind(cnt_sport),
    c_outside_other = YesNoNA_Ind(cnt_outside_other),
    c_other = YesNoNA_Ind(cnt_otheryn),
    c_other_house = YesNoNA_Ind(cnt_other_house),
    c_work = YesNoNA_Ind(cnt_work),
    c_worship = YesNoNA_Ind(cnt_worship),
    c_public_transport = YesNoNA_Ind(cnt_public_transport),
    c_school = YesNoNA_Ind(cnt_school),
    c_supermarket = YesNoNA_Ind(cnt_supermarket),
    c_shop = YesNoNA_Ind(cnt_shop),
    c_leisure = YesNoNA_Ind(cnt_leisure),
    c_health_facility = YesNoNA_Ind(cnt_health_facility),
    c_salon = YesNoNA_Ind(cnt_salon),
    sum_c = 0,
    other_text = cnt_other_text,
    inside = YesNoNA(cnt_inside),
    outside = YesNoNA(cnt_outside),
    phys = as.logical(phys_contact),
    hhm = (hhm_contact_yn == "Yes"),
    relationship = ifelse(hhm_contact_yn == "Yes", "household", map_type[cnt_type]),
    freq = map_freq[cnt_frequency],
    minutes_min = ifelse(!is.na(cnt_total_mins), cnt_total_mins, map_minutes_min[cnt_total_time]),
    minutes_max = ifelse(!is.na(cnt_total_mins), cnt_total_mins, map_minutes_max[cnt_total_time]),
    prec_none = YesNoNA(cnt_precautions_none),
    prec_2m = YesNoNA(cnt_two_metres_plus),
    prec_1m = YesNoNA(cnt_one_metre_plus),
    prec_0m = YesNoNA(cnt_within_one_metre),
    prec_mask = YesNoNA(cnt_mask),
    prec_wash_before = YesNoNA(cnt_wash_before),
    prec_wash_after = YesNoNA(cnt_wash_after),
    prec_prefer_not = YesNoNA(cnt_precautions_prefer_not_to_say),
    flag = cnt_nickname_flag
)]

contacts_clean[, sum_c := nafill(rowSums(.SD, na.rm = TRUE), fill = 0), .SDcols = patterns("^c_")]
warning("Next line: check column indices!")
set(contacts_clean, i = which(contacts_clean[, contacts == 0]), j = 11:46, value = NA)
set(contacts_clean, i = which(contacts_clean[, contacts == 0]), j = 15:29, value = 0)
contacts_clean = contacts_clean[order(date, part_id, cont_id)]

contacts_clean[, total_c  := sum(contacts), by = .(date, part_id)]
contacts_clean[, d_home   := c_home == 1]
contacts_clean[, d_school := ifelse(part_age_max <= 18, c_school == 1 & c_home == 0, c_school == 1 & c_home == 0 & c_work == 0)]
contacts_clean[, d_work := ifelse(part_age_min >= 18, c_work == 1 & c_home == 0, c_work == 1 & c_home == 0 & c_school == 0)]
contacts_clean[, d_other  := contacts == 1 & c_home == 0 & c_school == 0 & c_work == 0]
contacts_clean[, e_home_non_hh := d_home == 1 & hhm == FALSE]
contacts_clean[, e_other_house := d_other == 1 & c_other_house == 1]
contacts_clean[contacts == 0, d_home := FALSE]
contacts_clean[contacts == 0, d_work := FALSE]
contacts_clean[contacts == 0, d_school := FALSE]
contacts_clean[contacts == 0, d_work18 := FALSE]
contacts_clean[contacts == 0, d_school18 := FALSE]
contacts_clean[contacts == 0, d_other := FALSE]
contacts_clean[contacts == 0, e_home_non_hh := FALSE]
contacts_clean[contacts == 0, e_other_house := FALSE]




# ANALYSIS BIT

num = contacts_clean[, .(home = sum(d_home), work = sum(d_work), school = sum(d_school), other = sum(d_other), 
    e_home_non_hh = sum(e_home_non_hh), e_other_house = sum(e_other_house)), 
    by = .(part_id, part_age = part_age_min, date)]
num = num[, lapply(.SD, function(x) pmin(50, x)), by = .(part_id, part_age, date), .SDcols = c("home", "work", "school", "other", "e_home_non_hh", "e_other_house")]
num = num[!is.na(home)]
num[, t := as.numeric(date - ymd("2020-01-01"))]

gm2 = rlang::duplicate(gm)
names(gm2) = str_replace(names(gm2), "_percent_change_from_baseline", "")
names(gm2) = str_replace(names(gm2), "_and", "")
num = merge(num, gm2[, 8:14], by = "date", all.x = T)

# Add POLYMOD
data(polymod)
poly_p = polymod$participants[country == "United Kingdom"]
poly_c = polymod$contacts
poly = merge(poly_p, poly_c, by = "part_id", all.x = TRUE);


poly[, d_home   := cnt_home == 1]
poly[, d_school := ifelse(part_age < 18, cnt_school == 1 & cnt_home == 0, cnt_school == 1 & cnt_home == 0 & cnt_work == 0)]
poly[, d_work   := ifelse(part_age >= 18, cnt_work == 1 & cnt_home == 0, cnt_work == 1 & cnt_home == 0 & cnt_school == 0)]
poly[, d_other  := (cnt_transport + cnt_leisure + cnt_otherplace) >= 1 & cnt_home == 0 & cnt_school == 0 & cnt_work == 0]
poly = poly[!is.na(d_home)] # only cuts out 6 contacts.

pnum = poly[, .(home = sum(d_home), work = sum(d_work), school = sum(d_school), other = sum(d_other)), by = .(part_id, part_age, date = ymd(sday_id))]
pnum = pnum[, lapply(.SD, function(x) pmin(50, x)), by = .(part_id, part_age, date), .SDcols = c("home", "work", "school", "other")]

pnum[, t := as.numeric(date - ymd("2006-01-01"))]
pnum[, retail_recreation := 0]
pnum[, grocery_pharmacy := 0]
pnum[, parks := 0]
pnum[, transit_stations := 0]
pnum[, workplaces := 0]
pnum[, residential := 0]

num[, study := "CoMix"]
pnum[, study := "POLYMOD"]

num = rbind(num, pnum, fill = TRUE)

## Add weekday weights
num[,weekday = factor(lubridate::wday(date, label = TRUE), ordered = FALSE)]
num[,weekday_weights := ifelse(weekday %in% c("Saturday", "Sunday"), 2,5)]

num[, retail_recreation := (100 + retail_recreation) * 0.01]
num[, grocery_pharmacy  := (100 + grocery_pharmacy ) * 0.01]
num[, parks             := (100 + parks            ) * 0.01]
num[, transit_stations  := (100 + transit_stations ) * 0.01]
num[, workplaces        := (100 + workplaces       ) * 0.01]
num[, residential       := (100 + residential      ) * 0.01]


# Un-oversample young people from POLYMOD
num = rbind(
    num[study == "CoMix"],
    num[study == "POLYMOD" & part_age <= 20][seq(0, .N, by = 2)],
    num[study == "POLYMOD" & part_age > 20]
)



# MODELS

asc = function(x, y0, y1, s0, s1)
{
    xx = s0 + x * (s1 - s0)
    h0 = exp(s0) / (1 + exp(s0))
    h1 = exp(s1) / (1 + exp(s1))
    h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0)
    y0 + (y1 - y0) * h
}

## Weighted means accounting for weekdays. 
wm <- function(x)   weighted.mean(x, weekday_weights)


# *** Home model
another = num[, .(home = wm(home), residential = mean(residential), transit = mean(transit_stations)), 
    by = .(week = ifelse(study == "CoMix", week(date) %/% 2 * 2, rep(0, length(date))), study)]

home_f = another[, mean(home), by = .(context = ifelse(study == "CoMix", "pandemic", "pre-pandemic"))]

days = data.table(x = 0:365)
plh = ggplot(another) + 
    geom_point(aes(x = week, y = home, colour = study)) + 
    geom_line(data = days, aes(x = x / 7, y = asc(x / 366, 3.875622, 1.545019, -79 * 0.6, 286 * 0.6))) +
    ylim(0, 4) + labs(x = "Week of 2020", y = "Home contacts", colour = "Study") +
    theme(legend.position = c(0.8, 0.8))

# Alternative home models
ggplot(another) + 
    geom_point(aes(x = week, y = home, colour = study)) + 
    geom_smooth(aes(x = week, y = home), method = "gam", formula = y ~ s(x))

ggplot(another) + 
    geom_point(aes(x = residential, y = home, colour = study)) + 
    geom_smooth(aes(x = residential, y = home), method = "gam", formula = y ~ s(x))

ggplot(another) + 
    geom_point(aes(x = transit, y = home, colour = study)) + 
    geom_smooth(aes(x = transit, y = home), method = "gam", formula = y ~ s(x))

# Alternative: Home, non household member
another = num[, .(home = mean(home), non_hh = mean(e_home_non_hh), residential = mean(residential), transit = mean(transit_stations)), 
    by = .(week = ifelse(study == "CoMix", week(date) %/% 2, rep(0, length(date))), study)]

#home_f = another[, mean(home), by = .(context = ifelse(study == "CoMix", "pre-pandemic", "pandemic"))]

ggplot(another) + 
    geom_point(aes(x = week, y = home, colour = study)) + 
    geom_smooth(aes(x = week, y = home), method = "gam", formula = y ~ s(x)) +
    geom_point(aes(x = week, y = non_hh, colour = study)) + 
    geom_smooth(aes(x = week, y = non_hh), method = "gam", formula = y ~ s(x))

ggplot(another) + 
    geom_point(aes(x = residential, y = non_hh, colour = study)) + 
    geom_smooth(aes(x = residential, y = non_hh), method = "gam", formula = y ~ s(x))

ggplot(another) + 
    geom_point(aes(x = transit, y = non_hh, colour = study)) + 
    geom_smooth(aes(x = transit, y = non_hh), method = "gam", formula = y ~ s(x))

ggplot(another) + 
    geom_point(aes(x = transit, y = home - non_hh, colour = study)) + 
    geom_smooth(aes(x = transit, y = home - non_hh), method = "gam", formula = y ~ s(x))

# Alternative: Other house
# TODO go back up and fix NAs in e_other_house . . .
another = num[, .(other = mean(other), other_house = mean(e_other_house, na.rm = T), 
    residential = mean(residential), retail = mean(retail_recreation), 
    grocery = mean(grocery_pharmacy), transit = mean(transit_stations)), 
    by = .(week = ifelse(study == "CoMix", week(date) %/% 2, rep(0, length(date))), study)]

another[, predictor := retail * 0.345 + transit * 0.445 + grocery * 0.210]

#home_f = another[, mean(home), by = .(context = ifelse(study == "CoMix", "pandemic", "pre-pandemic"))]

ggplot(another) + 
    geom_point(aes(x = week, y = other, colour = study)) + 
    geom_smooth(aes(x = week, y = other), method = "gam", formula = y ~ s(x)) +
    geom_point(aes(x = week, y = other_house, colour = study)) + 
    geom_smooth(aes(x = week, y = other_house), method = "gam", formula = y ~ s(x))

ggplot(another) + 
    geom_point(aes(x = residential, y = other_house, colour = study)) + 
    geom_smooth(aes(x = residential, y = other_house), method = "gam", formula = y ~ s(x))

ggplot(another) + 
    geom_point(aes(x = transit, y = other_house, colour = study)) + 
    geom_smooth(aes(x = transit, y = other_house), method = "gam", formula = y ~ s(x))

ggplot(another) + 
    geom_point(aes(x = predictor, y = other_house, colour = study)) + 
    geom_smooth(aes(x = predictor, y = other_house), method = "gam", formula = y ~ s(x))


# *** Work model
another = num[part_age %between% c(18, 65), .(work = mean(work), workplaces = mean(workplaces), transit = mean(transit_stations)), 
    by = .(week = ifelse(study == "CoMix", week(date) %/% 2, rep(0, length(date))), study)]

model = gam(work ~ s(workplaces), family = gaussian, data = another)
work_f = data.table(workplaces = seq(0, 1.25, by = 0.01));
work_f[, work := pmax(0.0, predict(model, work_f, type = "response"))]

plw = ggplot(another) + 
    geom_point(aes(x = workplaces, y = work, colour = study)) + 
    geom_line(data = work_f, aes(x = workplaces, y = work)) +
    ylim(0, 3.5) +
    labs(x = "Google Mobility\n'workplaces' visits", y = "Work contacts", colour = "Study") +
    theme(legend.position = "none")

# Alternative work models

model3 = gam(work ~ s(workplaces), family = gaussian, data = another)
model3 = gam(work ~ 0 + I(workplaces^2), family = gaussian, data = another)
model3 = gam(work ~ 0 + I(workplaces), family = gaussian, data = another)
summary(model3)
plot(model3)

# *** School model
another = num[part_age %between% c(0, 18), .(school = mean(school), residential = mean(residential)), 
    by = .(week = ifelse(study == "CoMix", week(date) %/% 2 * 2, rep(0, length(date))), study)]

days = data.table(x = 0:365)

pls = ggplot(another) + 
    geom_point(aes(x = week, y = school, colour = study)) + 
    geom_line(data = days, aes(x = x / 7, y = ifelse(x <= 81 | x >= 244, 5.67, 0))) +
    labs(x = "Week of 2020", y = "School contacts", colour = "Study") +
    theme(legend.position = "none")


# *** Other model
another = num[, .(other = mean(other), retail = mean(retail_recreation), 
    grocery = mean(grocery_pharmacy), transit = mean(transit_stations)), 
    by = .(week = ifelse(study == "CoMix", week(date) %/% 2, rep(0, length(date))), study)]
another[, predictor := retail * 0.345 + transit * 0.445 + grocery * 0.210] # See "optimisation" below for how this was arrived at

model = gam(other ~ s(predictor), family = gaussian, data = another)
other_f = data.table(predictor = seq(0, 1.25, by = 0.01));
other_f[, other := pmax(0.0, predict(model, other_f, type = "response"))]

plo = ggplot(another) + 
    geom_point(aes(x = predictor, y = other, colour = study)) + 
    geom_line(data = other_f, aes(x = predictor, y = other)) +
    xlim(0, 1.25) + ylim(0, 5) + labs(x = "Google Mobility weighted 'transit stations',\n'retail and recreation', and 'grocery and pharmacy' visits", y = "Other contacts", colour = "Study") +
    theme(legend.position = "none")

# Alternative "other" models
ggplot(another) + 
    geom_point(aes(x = retail * 0.333 + transit * 0.333 + grocery * 0.333, y = other, colour = study)) + 
    geom_smooth(aes(x = retail * 0.333 + transit * 0.333 + grocery * 0.333, y = other), method = "gam", formula = y ~ s(x), fullrange = T) +
    xlim(0, 1) + ylim(0, 4)

ggplot(another) + 
    geom_point(aes(x = retail, y = other, colour = study)) + 
    geom_smooth(aes(x = retail, y = other), method = "gam", formula = y ~ s(x), fullrange = T) +
    xlim(0, 1)

ggplot(another) + 
    geom_point(aes(x = grocery, y = other, colour = study)) + 
    geom_smooth(aes(x = grocery, y = other), method = "gam", formula = y ~ s(x), fullrange = T) +
    xlim(0, 1)

ggplot(another) + 
    geom_point(aes(x = transit, y = other, colour = study)) + 
    geom_smooth(aes(x = transit, y = other), method = "gam", formula = y ~ s(x), fullrange = T) +
    xlim(0, 1) + ylim(0, 4)

ggplot(another) + 
    geom_point(aes(x = week, y = other, colour = study)) + 
    geom_smooth(aes(x = week, y = other), method = "gam", formula = y ~ s(x), fullrange = T) +
    xlim(0, 24) + ylim(0, 4)

model3 = gam(other ~ 0 + transit + I(transit^2), family = gaussian, data = another)
model3 = gam(other ~ s(transit), family = gaussian, data = another)
summary(model3)
plot(model3)
pred = data.table(transit = (0:100) / 100)
pred = cbind(pred, predict(model3, pred, type = "response"))

ggplot(pred) + 
    geom_point(data = another, aes(x = transit, y = other, colour = study)) + 
    geom_line(aes(x = transit, y = exp(V2)))


# *** Plot
cowplot::plot_grid(plh, pls, plw, plo, nrow = 2, align = "hv")
ggsave("./figures/comix_analysis.png", width = 20, height = 15, units = "cm")


# OPTIMISING "other"
another = num[, .(other = mean(other), 
    retail = mean(retail_recreation), 
    grocery = mean(grocery_pharmacy), 
    transit = mean(transit_stations)),
    by = .(week = ifelse(study == "CoMix", week(date) %/% 2, rep(0, length(date))), study)]

objective = function(x)
{
    if (sum(x) > 1 || any(x < 0)) {
        return (Inf)
    }
    #model = gam(other ~ 0 + s(I(x[1] * retail + x[2] * grocery + x[3] * transit +
    #                            x[4] * parks + x[5] * workplaces + (1 - sum(x)) * residential)), family = gaussian, data = another)

    model = gam(other ~ 0 + s(I(x[1] * retail + x[2] * grocery + (1 - sum(x)) * transit)), family = gaussian, data = another)

    #model$aic
    -summary(model)$dev.expl
}

model = optim(rep(1/3, 2), objective, method = "Nelder-Mead", control = list(trace = 5))
model$par
1 - sum(model$par)


# OPTIMISING "work" (?)
another = num[, .(work = mean(work), 
    retail = mean(retail_recreation), 
    grocery = mean(grocery_pharmacy), 
    transit = mean(transit_stations),
    parks = mean(parks),
    workplaces = mean(workplaces),
    residential = mean(residential)), 
    by = .(week = ifelse(study == "CoMix", week(date) %/% 2, rep(0, length(date))), study)]

objective = function(x)
{
    if (sum(x) > 1 || any(x < 0)) {
        return (Inf)
    }
    model = gam(work  ~  0 + s(I(x[1] * retail + x[2] * grocery + x[3] * transit +
                               (1 - sum(x)) * workplaces)), family = gaussian, data = another)

    #model$aic
    -summary(model)$dev.expl
}

model2 = optim(c(0,0,0), objective, method = "Nelder-Mead", control = list(trace = 5))
model2$par
1 - sum(model2$par)

x = model2$par
model = gam(work ~ 0 + workplaces + I(workplaces * transit), family = gaussian, data = another)
model = gam(work ~ workplaces, family = gaussian, data = another)

ggplot(another[, .(x = 0.1784 * workplaces + 1.9290 * workplaces * transit, y = work)]) +
    geom_point(aes(x, y))

ggplot(another[, .(x = -0.6448 + 2.4115 * workplaces, y = work)]) +
    geom_point(aes(x, y))


# Interpolation curves
work_f[, paste(round(work / work[101], 3), collapse = ", ")]
other_f[, paste(round(other / other[101], 3), collapse = ", ")]



summary(model)

model = gam(school ~ s(t), family = ziP, method = "REML", data = num)
summary(model)
plot(model)

model = gam(other ~ s(retail_recreation) + s(grocery_pharmacy) + s(transit_stations), family = ziP, method = "REML", data = num)
model = gam(other ~ retail_recreation + grocery_pharmacy + transit_stations, family = ziP, method = "REML", data = num)
model = gam(other ~ retail_recreation + grocery_pharmacy + transit_stations, family = ziP, method = "REML", data = num)
summary(model)
plot(model)
gam.check(model)

pred = data.table(t = 83:271)
pred = cbind(pred, predict(model, pred, type = "response"))

ggplot(pred) + 
    geom_line(aes(x = ymd("2020-01-01") + t, y = V2 * 0.25 + 0.25), colour = "blue") +
    geom_line(aes(x = date, y = rollmean(1 + workplaces_percent_change_from_baseline / 100, 7, fill = NA)), data = gm) +
    ylim(0, NA)

# why so different?
ggplot() +
    geom_point(data = work_contacts, aes(x = t, y = n_cnt_work, colour = "work_contacts")) +
    geom_point(data = num, aes(x = t, y = work, colour = "num"))


num = melt(num, id.vars = c("part_id", "date", "t"))

ggplot(num, aes(x = date, y = value)) + 
    #geom_point() +
    geom_smooth(method = "gam") +
    facet_wrap(~variable)



# TODO: see what else from participants would be useful.

contacts_clean[, part_age := floor(runif(.N, part_age_min, part_age_max))]
contacts_clean[, cont_age := floor(runif(.N, cont_age_min, cont_age_max))]
ggplot(contacts_clean) + geom_bin2d(aes(x = part_age, y = cont_age), binwidth = 5) + viridis::scale_fill_viridis()






parts = participants[, .(part_age = first(part_age[!is.na(part_age)]), 
    part_age_group = first(part_age_group[!is.na(part_age_group)])), by = part_id]

# Impute missing participant ages
parts[, amin := part_age]
parts[, amax := part_age]
parts[is.na(part_age) & !is.na(part_age_group), amin := as.numeric(str_replace_all(part_age_group, "\\[|,[0-9]+\\)", ""))]
parts[is.na(part_age) & !is.na(part_age_group), amax := as.numeric(str_replace_all(part_age_group, "\\[[0-9]+,|\\)", ""))]
parts[, part_age := floor(runif(.N, amin, amax))]
parts = parts[, .(part_id, part_age)];

full = merge(contacts, parts[, .(part_id, part_age)], by = "part_id")
full[, table(cnt_age)]

# Impute contact ages
full[cnt_age == "Under 1", cnt_age_min := 0]
full[cnt_age == "Under 1", cnt_age_max := 1]
full[cnt_age %like% "^[0-9]+\\+$", cnt_age_min := as.numeric(str_replace_all(cnt_age, "\\+", ""))]
full[cnt_age %like% "^[0-9]+\\+$", cnt_age_max := 93]
full[cnt_age %like% "^[0-9]+-[0-9]+$", cnt_age_min := as.numeric(str_replace_all(cnt_age, "-[0-9]+", ""))]
full[cnt_age %like% "^[0-9]+-[0-9]+$", cnt_age_max := as.numeric(str_replace_all(cnt_age, "[0-9]+-", ""))]
full[, cnt_age := floor(runif(.N, cnt_age_min, pmin(93, cnt_age_max)))]

ggplot(full) + geom_bin2d(aes(x = part_age, y = cnt_age), binwidth = 5)

asc = function(x, y0, y1, s0, s1)
{
    xx = s0 + x * (s1 - s0);
    h0 = exp(s0) / (1 + exp(s0));
    h1 = exp(s1) / (1 + exp(s1));
    h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);
    y0 + (y1 - y0) * h;
}

# How to define home decrease?
asch = data.table(t = 0:365);
asch[, y := asc(t/365, 0, -75, -7.9*6, 28.6*6)]
asch[, date := ymd("2020-01-01") + t]

ggplot(gm) + 
    geom_line(aes(x = date, y = retail_and_recreation_percent_change_from_baseline)) +
    geom_line(data = asch, aes(x = date, y = y), colour = "red")



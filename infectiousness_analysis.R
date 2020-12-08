library(shmanipulate)
library(deSolve)
library(ggplot2)
library(data.table)
library(mcmc)
library(stringr)

# Infectiousness distribution from Ferretti et al
shift = -0.0747
scale = 1.8567
df = 3.3454
curve(dt((x - shift) / scale, df) / scale, -10, 10)

inf_curve = data.table(x = seq(-10, 10, by = 0.05))
inf_curve[, infectiousness := dt((x - shift) / scale, df) / scale]
inf_curve[, infectiousness := infectiousness / max(infectiousness)]

# Exploratory analysis
shmanipulate({
    inf_rand = data.table(person = 1:nsamp)
    inf_rand[, E := rgamma(.N, shape = Es, rate = Es / Em)]
    inf_rand[, P := rgamma(.N, shape = Ps, rate = Ps / Pm)]
    inf_rand[, S := rgamma(.N, shape = Ss, rate = Ss / Sm)]
    
    tost_summ = data.table(TOST = seq(-10, 10, by = 0.05))
    inf_summ = data.table(t = seq(0, 20, by = 0.05))
    
    for (i in 1:nrow(tost_summ)) {
        tost_summ[i, p_inf := ifelse(TOST < 0, inf_rand[, sum(P > -TOST)/.N], inf_rand[, sum(S > TOST)/.N])]
        inf_summ[i, p_inf := inf_rand[, sum(t >= E & t <= E + P + S)/.N]]
    }
    
    incubation = inf_rand[, mean(E + P)]

    plA = ggplot() +
        geom_line(data = inf_curve, aes(x, infectiousness), size = 2, colour = "red") +
        geom_line(data = tost_summ, aes(x = TOST + delta, y = p_inf)) +
        labs(x = "TOST", y = "Infectiousness")
    
    plB = ggplot() +
        geom_line(data = inf_summ, aes(x = t, y = p_inf)) +
        annotate("label", label = incubation, x = incubation, y = 0.5) +
        labs(x = "Time since exposure", y = "Infectiousness")
    
    cowplot::plot_grid(plA, plB)
}, Em = 1, Es = 1, Pm = 1, Ps = 1, Sm = 1, Ss = 1, delta = shift, nsamp = list(1000, 10000, 100000, 1000000))

Em = 1
Es = 1
Pm = 1
Ps = 1
Sm = 1
Ss = 1

inf_rand = data.table(person = 1:1000)
inf_rand[, E := rgamma(.N, shape = Es, rate = Es / Em)]
inf_rand[, P := rgamma(.N, shape = Ps, rate = Ps / Pm)]
inf_rand[, S := rgamma(.N, shape = Ss, rate = Ss / Sm)]

inf_summ = data.table(TOST = seq(-10, 10, by = 0.05))

for (i in 1:nrow(inf_summ)) {
    inf_summ[i, p_inf := ifelse(TOST < 0, inf_rand[, sum(P > -TOST)], inf_rand[, sum(S > TOST)])]
}


# Positivity function for integration.
# Assumes time to conversion is gamma distributed with mean mean_conv and shape shape_conv,
# and time from conversion to reversion is gamma distributed with mean mean_wane and shape shape_wane.
# Then the probability that an individual tests positive at time x days after infection is
#    integral from 0 to x ( dgamma(y | mean_conv, shape_conv) * (1 - pgamma(x - y | mean_wane, shape_wane)) dy
# Above, the first term is the probability that the invididual has converted at time y < x and the second
# term is the probability that they have not reverted by time x.
positive = function(y, x, mean_conv, shape_conv, mean_wane, shape_wane)
{
    rate_conv = shape_conv / mean_conv
    rate_wane = shape_wane / mean_wane
    dgamma(y, shape_conv, rate_conv) * (1 - pgamma(x - y, shape_wane, rate_wane))
}



# Approach 1: maximum likelihood, Nelder-Mead optimisation
objective = function(x) {
    mean_conv = x[1]
    shape_conv = x[2]
    mean_wane = x[3]
    shape_wane = x[4]
    
    if (any(x <= 0)) {
        return (10000)
    }
    
    d = data.table(t = setdiff(pos_curve$days_since_infection, 0), pos = 0)
    for (i in 1:nrow(d)) {
        d[i, pos := integrate(positive, 0, t, x = t, 
            mean_conv = mean_conv, shape_conv = shape_conv, mean_wane = mean_wane, shape_wane = shape_wane)[[1]]]
    }
    
    d = merge(pos_curve, d, by.x = "days_since_infection", by.y = "t")
    
    return (d[, mean((pos - median)^2)])
}

ans = optim(c(1, 1, 1, 1), objective, method = "Nelder-Mead", control = list(trace = 5, maxit = 1000))
ans
# Result: Conversion ~ gamma(mean = 2.64 [= ans$par[1]], shape = 6.07 [= ans$par[2]])
#         Waning ~ gamma(mean = 8.77 [= ans$par[3]], shape = 1.52 [= ans$par[4]])


# Approach 2: Posterior mean with MCMC
posterior = function(x) {
    mean_conv = x[1]
    shape_conv = x[2]
    mean_wane = x[3]
    shape_wane = x[4]
    
    if (any(x <= 0)) {
        return (-1e12)
    }
    
    d = data.table(t = setdiff(pos_curve$days_since_infection, 0), pos = 0)
    for (i in 1:nrow(d)) {
        d[i, pos := integrate(positive, 0, t, x = t, 
            mean_conv = mean_conv, shape_conv = shape_conv, mean_wane = mean_wane, shape_wane = shape_wane)[[1]]]
    }
    
    d = merge(pos_curve, d, by.x = "days_since_infection", by.y = "t")
    d[, ll := dsn(pos, p1, p2, p3, log = TRUE)]
    
    return (d[2:.N, sum(ll)])
}

# Initial exploration via MCMC
set.seed(42)
x_init = c(1, 1, 1, 1)
out = metrop(posterior, x_init, 1e3)
out$accept
# Acceptance probability too low. Rerunning with more iterations and lower scale.

# Rerunning -- scale found by experimenting...
out = metrop(out, scale = c(0.4, 1, 0.25, 0.25) * 0.7, nbatch = 1e4)
out$accept # 0.28 acceptance rate - good

# Results
d = as.data.table(out$batch)
d[, trial := 1:.N]
# Examine chains
ggplot(d) + 
    geom_line(aes(x = trial, y = V1, colour = "V1")) +
    geom_line(aes(x = trial, y = V2, colour = "V2")) +
    geom_line(aes(x = trial, y = V3, colour = "V3")) +
    geom_line(aes(x = trial, y = V4, colour = "V4"))
d = d[trial > 5000]
# Looks like decent mixing.

# Examine posterior distributions
ggplot(d) + geom_histogram(aes(x = V1))
ggplot(d) + geom_histogram(aes(x = V2))
ggplot(d) + geom_histogram(aes(x = V3))
ggplot(d) + geom_histogram(aes(x = V4))

# Get posterior means
d[, mean(V1)]
d[, mean(V2)]
d[, mean(V3)]
d[, mean(V4)]

# > d[, mean(V1)]
# [1] 2.76385
# > d[, mean(V2)]
# [1] 4.787206
# > d[, mean(V3)]
# [1] 8.470949
# > d[, mean(V4)]
# [1] 1.960072


# Figure for supplement.
library(cowplot)

theme_set(theme_cowplot(font_size = 8))

# use thinned posterior
d0 = d[seq(1, .N, by = 10)]

# return a single PCR positivity curve using parameters x
trace = function(x) {
    mean_conv = x[1]
    shape_conv = x[2]
    mean_wane = x[3]
    shape_wane = x[4]
    
    d = data.table(t = 0:30, pos = 0)
    for (i in 2:nrow(d)) {
        d[i, pos := integrate(positive, 0, t, x = t, 
            mean_conv = mean_conv, shape_conv = shape_conv, mean_wane = mean_wane, shape_wane = shape_wane)[[1]]]
    }

    return (d)
}

# Get posterior distribution of curves
for (i in 1:nrow(d0)) {
    tr = trace(d0[i, c(V1, V2, V3, V4)])
    d0[i, paste0("t", 0:30) := as.list(tr$pos)]
}

# Calculate quantiles
d1 = melt(d0[, .SD, .SDcols = -(1:4)], id.vars = "trial", variable.name = "time")
d1[, time := as.numeric(str_remove_all(time, "t"))]
d1 = d1[, as.list(quantile(value, c(0.025, 0.5, 0.975))), by = time]

pl_pos = ggplot() +
    geom_ribbon(data = d1, aes(x = time, ymin = `2.5%`, ymax = `97.5%`), fill = "#88eeee") +
    geom_line(data = d1, aes(x = time, y = `50%`), colour = "#449999") +
    geom_pointrange(data = pos_curve, aes(x = days_since_infection, y = median, ymin = lower, ymax = upper), fatten = 0.2) +
    labs(x = "Days since infection", y = "Probability of positive PCR test")

# Histograms for values
d2 = melt(d, id.vars = "trial")
d2[variable == "V1", variable := "(i) Days from infection\nto positivity"]
d2[variable == "V2", variable := "(ii) Shape parameter of\ngamma distribution for (i)"]
d2[variable == "V3", variable := "(iii) Days from positivity\nto loss of positivity"]
d2[variable == "V4", variable := "(iv) Shape parameter of\ngamma distribution for (iii)"]
d2[, variable := factor(variable, unique(variable))]

pl_hist = ggplot(d2) +
    geom_histogram(aes(x = value, y = after_stat(density), fill = variable), bins = 20) +
    facet_wrap(~variable, scales = "free") +
    theme(strip.background = element_blank(), legend.position = "none") +
    labs(x = "Value", y = "Posterior density")


pl = plot_grid(pl_pos, pl_hist, nrow = 1, labels = letters[1:2], label_size = 10)
ggsave("./figures/pcr_analysis.pdf", plot = pl, width = 18, height = 8, units = "cm", useDingbats = FALSE)
ggsave("./figures/pcr_analysis.png", plot = pl, width = 18, height = 8, units = "cm")


library(shmanipulate)
library(deSolve)
library(ggplot2)
library(data.table)
library(mcmc)
library(sn)

pos_curve = fread("~/Downloads/pos_curve__1_.csv")

solve_sn = function(lower, median, upper)
{
    objective = function(x) { 
        (psn(lower, x[1], x[2], x[3]) - 0.025) ^ 2 + 
        (psn(median, x[1], x[2], x[3]) - 0.5) ^ 2 + 
        (psn(upper, x[1], x[2], x[3]) - 0.975) ^ 2 
    }
    optim(c(1, 1, 1), objective, method = "Nelder-Mead")
}

for (i in 1:nrow(pos_curve)) {
    ans = solve_sn(pos_curve[i, lower], pos_curve[i, median], pos_curve[i, upper])
    pos_curve[i, p1 := ans$par[1]]
    pos_curve[i, p2 := ans$par[2]]
    pos_curve[i, p3 := ans$par[3]]
}


positive = function(y, x, mean_conv, shape_conv, mean_wane, shape_wane)
{
    rate_conv = shape_conv / mean_conv
    rate_wane = shape_wane / mean_wane
    dgamma(y, shape_conv, rate_conv) * (1 - pgamma(x - y, shape_wane, rate_wane))
}

# Exploratory analysis
shmanipulate({
    dt = data.table(t = 0:30)
    for (i in 1:nrow(dt)) {
        dt[i, pos := integrate(positive, 0, t, x = t, 
            mean_conv = conv, shape_conv = shape_conv, mean_wane = wane, shape_wane = shape_wane)[[1]]]
    }

    ggplot() +
        geom_ribbon(data = pos_curve, aes(x = days_since_infection, ymin = lower, ymax = upper), fill = "#ddddff") +
        geom_line(data = pos_curve, aes(x = days_since_infection, y = median), colour = "#9999bb") +
        geom_line(data = dt, aes(x = t, y = pos), colour = "#ff0000", size = 2)
        
}, conv = c(1), shape_conv = c(1), wane = c(1), shape_wane = c(1))


# Approach 1: maximum likelihood, Nelder-Mead optimisation
objective = function(x) {
    mean_conv = x[1]
    shape_conv = x[2]
    mean_wane = x[3]
    shape_wane = x[4]
    
    if (any(x <= 0)) {
        return (10000)
    }
    
    d = data.table(t = 0:30, pos = 0)
    for (i in 2:nrow(d)) {
        d[i, pos := integrate(positive, 0, t, x = t, 
            mean_conv = mean_conv, shape_conv = shape_conv, mean_wane = mean_wane, shape_wane = shape_wane)[[1]]]
    }
    
    d = merge(pos_curve, d, by.x = "days_since_infection", by.y = "t")
    
    return (d[, mean((pos - median)^2)])
}

ans = optim(c(1, 1, 1, 1), objective, method = "Nelder-Mead", control = list(trace = 5, maxit = 1000))
ans
# Result: Conversion ~ gamma(mean = 2.64, shape = 6.07)
#         Waning ~ gamma(mean = 8.77, shape = 1.52)



# Approach 2: Posterior mean with MCMC
posterior = function(x) {
    mean_conv = x[1]
    shape_conv = x[2]
    mean_wane = x[3]
    shape_wane = x[4]
    
    if (any(x <= 0)) {
        return (-1e12)
    }
    
    d = data.table(t = 0:30, pos = 0)
    for (i in 2:nrow(d)) {
        d[i, pos := integrate(positive, 0, t, x = t, 
            mean_conv = mean_conv, shape_conv = shape_conv, mean_wane = mean_wane, shape_wane = shape_wane)[[1]]]
    }
    
    d = merge(pos_curve, d, by.x = "days_since_infection", by.y = "t")
    d[, ll := dsn(pos, p1, p2, p3, log = TRUE)]
    
    return (d[2:.N, sum(ll)])
}

# Initial exploration
set.seed(42)
x_init = c(1, 1, 1, 1)
out = metrop(posterior, x_init, 1e3)
out$accept
# Acceptance probability too low. Rerunning with more iterations and lower scale.

# Rerunning -- scale found by experimenting...
out = metrop(out, scale = c(0.4, 1, 0.25, 0.25) * 0.7, nbatch = 1e4)
out$accept # 0.28 - good

# Results
d = as.data.table(out$batch)
d[, trial := 1:.N]
# Examine chains
ggplot(d) + 
    geom_line(aes(x = trial, y = V1, colour = "V1")) +
    geom_line(aes(x = trial, y = V2, colour = "V2")) +
    geom_line(aes(x = trial, y = V3, colour = "V3")) +
    geom_line(aes(x = trial, y = V4, colour = "V4"))
d = d[trial > 3000]


ggplot(d) + geom_histogram(aes(x = V1))
ggplot(d) + geom_histogram(aes(x = V2))
ggplot(d) + geom_histogram(aes(x = V3))
ggplot(d) + geom_histogram(aes(x = V4))

d[, mean(V1)]
d[, mean(V2)]
d[, mean(V3)]
d[, mean(V4)]







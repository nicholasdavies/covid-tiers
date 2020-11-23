library(ggplot2)
library(data.table)
library(shmanipulate)

# approximating the sum of two beta variables


p = 1
q = 10

a = 3
A = 10

b = 9
B = 12

pp = p / (p + q)
qq = q / (p + q)



h = pp^2 * (a*(A-a)/(A^2*(A+1))) + qq^2 * (b*(B-b)/(B^2*(B+1)))

K = ((2 * pp * a / A + 2 * qq * b / B - 1)^2 + 4 * h - 1) / (-4*h)
k1 = K * (1 - sqrt(1 - 4 * h * (K + 1))) / 2
k2 = K * (1 + sqrt(1 - 4 * h * (K + 1))) / 2


test_data = data.table(true = p * rbeta(10000, a, A - a) + q * rbeta(10000, b, B - b),
    approx1 = (p + q) * rbeta(10000, k1, K - k1),
    approx2 = (p + q) * rbeta(10000, k2, K - k2)
    )

ggplot(test_data) +
    geom_freqpoly(aes(true, colour = "true")) +
    geom_freqpoly(aes(approx1, colour = "approx1")) +
    geom_freqpoly(aes(approx2, colour = "approx2"))





# approximating the sum of two beta variables
shmanipulate({
    p = p1
    q = p2
    
    a = alpha1
    A = alpha1 + beta1
    
    b = alpha2
    B = alpha2 + beta2
    
    pp = p / (p + q)
    qq = q / (p + q)
    
    h = pp^2 * (a*(A-a)/(A^2*(A+1))) + qq^2 * (b*(B-b)/(B^2*(B+1)))
    
    K = ((2 * pp * a / A + 2 * qq * b / B - 1)^2 + 4 * h - 1) / (-4*h)
    k1 = K * (1 - sqrt(1 - 4 * h * (K + 1))) / 2
    k2 = K * (1 + sqrt(1 - 4 * h * (K + 1))) / 2
    
    # which k to use? Compare true mean to approximated mean... that is the one we want
    true_mean = p * a / A + q * b / B;
    k1_mean = (p + q) * k1 / K;
    k2_mean = (p + q) * k2 / K;
    if (abs(k1_mean - true_mean) < abs(k2_mean - true_mean)) {
        k = k1;
    } else {
        k = k2;
    }
    
    illust = data.table(x = seq(0, p + q, by = (p + 1) / 1000))
    # illust[, approx1 := dbeta(x / (p + q), k1, K - k1) / (p + q)]
    # illust[, approx2 := dbeta(x / (p + q), k2, K - k2) / (p + q)]
    illust[, approx := dbeta(x / (p + q), k, K - k) / (p + q)]

    rs = data.table(x = p * rbeta(100000, a, A - a) + q * rbeta(100000, b, B - b))
    
    ggplot(illust) +
        geom_histogram(data = rs, aes(x = x, y = after_stat(density))) +
        geom_line(aes(x, approx, colour = "approx"))

}, p1 = c(0, 10), alpha1 = c(0, 10), beta1 = c(0, 10), p2 = c(0, 10), alpha2 = c(0, 10), beta2 = c(0, 10))


p = 1
q = 10

a = 3
A = 10

b = 9
B = 12

pp = p / (p + q)
qq = q / (p + q)



h = pp^2 * (a*(A-a)/(A^2*(A+1))) + qq^2 * (b*(B-b)/(B^2*(B+1)))

K = ((2 * pp * a / A + 2 * qq * b / B - 1)^2 + 4 * h - 1) / (-4*h)
k1 = K * (1 - sqrt(1 - 4 * h * (K + 1))) / 2
k2 = K * (1 + sqrt(1 - 4 * h * (K + 1))) / 2


test_data = data.table(true = p * rbeta(10000, a, A - a) + q * rbeta(10000, b, B - b),
    approx1 = (p + q) * rbeta(10000, k1, K - k1),
    approx2 = (p + q) * rbeta(10000, k2, K - k2)
    )

ggplot(test_data) +
    geom_freqpoly(aes(true, colour = "true")) +
    geom_freqpoly(aes(approx1, colour = "approx1")) +
    geom_freqpoly(aes(approx2, colour = "approx2"))



# approximating the sum of two gamma variables
shmanipulate({
    A = (alpha1/beta1 + alpha2/beta2)^2 / (alpha1/beta1^2 + alpha2/beta2^2)
    B = (alpha1/beta1 + alpha2/beta2) / (alpha1/beta1^2 + alpha2/beta2^2)

    illust = data.table(x = seq(0, 10 * A/B, by = 0.01 * A/B))
    illust[, approx := dgamma(x, shape = A, rate = B)]

    rs = data.table(x = rgamma(100000, shape = alpha1, rate = beta1) + rgamma(100000, shape = alpha2, rate = beta2))
    
    ggplot(illust) +
        geom_histogram(data = rs, aes(x = x, y = after_stat(density))) +
        geom_line(aes(x, approx, colour = "approx"))

}, alpha1 = c(0, 10), beta1 = c(0, 10), alpha2 = c(0, 10), beta2 = c(0, 10))

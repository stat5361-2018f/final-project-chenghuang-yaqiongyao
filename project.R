#####################################################################
## project computing
#####################################################################


## generate bivariate copula data
## https://www.r-bloggers.com/copulas-made-easy/
## Gaussiam Copula
## one rv is from gamma(1,2)
## one rv is from chisq(1)

require(mvtnorm)
n <- 1000
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
set.seed(123)
dat <- rmvnorm(n, sigma = sigma)
dat <- pnorm(dat)
x <- qgamma(dat[,1], shape = 1, scale = 2)
y <- qchisq(dat[,2], df = 10)
plot(x, y)


## generate random number with empirical cdf

sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
set.seed(456)
dat <- rmvnorm(n, sigma = sigma)
dat <- pnorm(dat)

x_ecdf <- quantile(x, dat[,1])
y_ecdf <- quantile(y, dat[,2])

par(mfrow = c(1,2))
plot(x, y, xlim = c(0, 15), ylim = c(0, 35))
plot(x_ecdf, y_ecdf, xlim = c(0, 15), ylim = c(0, 35))


## generate random number from KDE
## use uniform kernel since it gives us closed form solution
## check the influence of bandwidth to the random number generating
sample.kernel <- function(n, x){
  bw <- density(x)$bw
  sample(x, n, replace = TRUE) + runif(n, -bw, bw)
}


x_kde <- sample.kernel(n, x)
y_kde <- sample.kernel(n, y)

par(mfrow = c(2,2))
plot(x, y, xlim = c(0, 15), ylim = c(0, 35))
plot(x_ecdf, y_ecdf, xlim = c(0, 15), ylim = c(0, 35))
plot(x_kde, y_kde, xlim = c(0, 15), ylim = c(0, 35))

## use close form solution of uniform kernel

sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
set.seed(456)
dat <- rmvnorm(n, sigma = sigma)
dat <- pnorm(dat)

bw_x <- density(x)$bw
bw_y <- density(y)$bw
x_kde_c <- (2*n*bw_x*dat[,1] - n*bw_x)/n + sample(x, n, replace = TRUE)
y_kde_c <- (2*n*bw_y*dat[,2] - n*bw_y)/n + sample(y, n, replace = TRUE)

par(mfrow = c(2,2))
plot(x, y, xlim = c(0, 15), ylim = c(0, 35))
plot(x_ecdf, y_ecdf, xlim = c(0, 15), ylim = c(0, 35))
plot(x_kde, y_kde, xlim = c(0, 15), ylim = c(0, 35))
plot(x_kde_c, y_kde_c, xlim = c(0, 15), ylim = c(0, 35))







#####################################################################
## project computing
#####################################################################


## univariate
n <- 100
set.seed(123)
u <- runif(100, 0, 1)
x <- as.numeric(u < 0.3) * rgamma(100, shape = 1, scale = 5) + 
  as.numeric(u > 0.3) * rgamma(100, shape = 5, scale = 1)
hist(x, breaks = 50)
plot(density(x))

true.density <- function(x){
  0.3 * dgamma(x, shape = 1, scale = 5) + 0.7 * dgamma(x, shape = 5, scale = 1)
}

hist(x, breaks = 50)
curve(true.density(x), add = TRUE)




## generate bivariate copula data
## https://www.r-bloggers.com/copulas-made-easy/
## Gaussiam Copula
## one rv is from gamma(1,2)
## one rv is from chisq(1)

require(mvtnorm)
n <- 100
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
set.seed(098)
dat <- rmvnorm(n, sigma = sigma)
dat <- pnorm(dat)
x <- qgamma(dat[,1], shape = 1, scale = 2)
y <- qchisq(dat[,2], df = 10)
plot(x, y)


## generate random number with empirical cdf
n <- 1000
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
## check the influence of bandwidth to the random number generating
sample.kernel <- function(n, x, y, adj){
  bw_x <- adj * density(x)$bw
  bw_y <- adj * density(y)$bw
  ind <- sample(1:100, n, replace = TRUE)
  x_kde <- x[ind] + runif(n, -bw_x, bw_x)
  y_kde <- y[ind] + runif(n, -bw_y, bw_y)
  cbind(x_kde, y_kde)
}
adj = 1
xy_kde <- sample.kernel(n, x, y, adj)

par(mfrow = c(2,2))
plot(x, y, xlim = c(0, 15), ylim = c(0, 35))
plot(x_ecdf, y_ecdf, xlim = c(0, 15), ylim = c(0, 35))
plot(xy_kde[,1], xy_kde[,2], xlim = c(0, 15), ylim = c(0, 35))








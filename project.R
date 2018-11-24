#####################################################################
## project computing
#####################################################################



## univariate example

n <- 100
set.seed(17)
x <- c(rnorm(n), rnorm(n, 4, 1/4), rnorm(n, 8, 1/4))

hist(x)


rdens <- function(n, density=z, data=x, kernel="gaussian") {
  width <- z$bw                              # Kernel width
  rkernel <- function(n) runif(n, -width, width)  # Kernel sampler
  sample(x, n, replace=TRUE) + rkernel(n)    # Here's the entire algorithm
}

dx <- function(x) (dnorm(x) + dnorm(x, 4, 1/4) + dnorm(x, 8, 1/4))/3
#
# Compute a kernel density estimate.
# It returns a kernel width in $bw as well as $x and $y vectors for plotting.
#
z <- density(x, bw=0.15, kernel="gaussian")
#
# Sample from the KDE.
#
system.time(y <- rdens(3*n, z, x)) # Millions per second
#
# Plot the sample.
#
h.density <- hist(y, breaks=60, plot=FALSE)
#
# Plot the KDE for comparison.
#
h.sample <- hist(x, breaks=h.density$breaks, plot=FALSE)
#
# Display the plots side by side.
#
histograms <- list(Sample=h.sample, Density=h.density)
y.max <- max(h.density$density) * 1.25
par(mfrow=c(1,2))
for (s in names(histograms)) {
  h <- histograms[[s]]
  plot(h, freq=FALSE, ylim=c(0, y.max), col="#f0f0f0", border="Gray",
       main=paste("Histogram of", s))
  curve(dx(x), add=TRUE, col="Black", lwd=2, n=501) # Underlying distribution
  lines(z$x, z$y, col="Red", lwd=2)                 # KDE of data
  
}
par(mfrow=c(1,1))








# use close form
# cannot use close form



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

## use close form solution of uniform kernel

# sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
# set.seed(456)
# dat <- rmvnorm(n, sigma = sigma)
# dat <- pnorm(dat)
# 
# adj <- 1
# ind <- sample(1:100, n, replace = TRUE)
# bw_x <- adj * density(x)$bw
# bw_y <- adj * density(y)$bw
# x_kde_c <- (2*n*bw_x*dat[,1] - n*bw_x)/n + x[ind]
# y_kde_c <- (2*n*bw_y*dat[,2] - n*bw_y)/n + y[ind]




# par(mfrow = c(2,2))
# plot(x, y, xlim = c(0, 15), ylim = c(0, 35))
# plot(x_ecdf, y_ecdf, xlim = c(0, 15), ylim = c(0, 35))
# plot(xy_kde[,1], xy_kde[,2], xlim = c(0, 15), ylim = c(0, 35))
# plot(x_kde_c, y_kde_c, xlim = c(0, 15), ylim = c(0, 35))







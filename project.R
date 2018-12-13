#####################################################################
## project computing
#####################################################################


## univariate
## generate data
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

hist(x, breaks = 50, probability = TRUE)
curve(true.density(x), add = TRUE)


## empirical CDF with kernel function
N <- 100   # number of samples

find.int.point <- function(x){
  y <- sort(x)
  bw <- density(x)$bw
  interval <- cbind(y-bw/2, y+bw/2)
  int.point <- NULL
  for (i in 1:dim(interval)[1]) {
    num.s <- sum((interval[1:i, 1] <= interval[i, 1]) & 
                   (interval[1:i, 2] >= interval[i, 1]))
    num.e <- sum((interval[(i):dim(interval)[1], 1] <= interval[i, 2]) & 
                   (interval[(i):dim(interval)[1], 2] >= interval[i, 2])) - 1
    num.m <- matrix(c(interval[i, 1], interval[i, 2], num.s, num.e), nrow = 2)
    int.point <- rbind(int.point, num.m)
  }
  int.point <- int.point[order(int.point[,1]),]
  sum((int.point[2:dim(int.point)[1], 1] - int.point[1:(dim(int.point)[1]-1),1])*
        int.point[1:(dim(int.point)[1]-1), 2]*(1/bw)) == 100
  int.point <- cbind(int.point, c((int.point[2:dim(int.point)[1], 1] - int.point[1:(dim(int.point)[1]-1),1])*
                                    int.point[1:(dim(int.point)[1]-1), 2]*(1/bw)*(1/length(x)), 0))
  int.point <- cbind(int.point, c(0,cumsum(int.point[,3])[-dim(int.point)[1]]))
  return(int.point)
}

sample.int <- function(int.point, u){
  num.int <- sapply(1:length(u), function(i) sum(int.point[,4] < u[i]))
  sample <- int.point[num.int, 1] + (u - int.point[num.int, 4])/(int.point[num.int, 2]*(1/bw))
  return(sample)
}

sample.kernel <- function(x, u){
  sample.int(find.int.point(x), u)
}

N <- 10000
set.seed(123)
u <- runif(N, 0, 1)
hist(u)
sample <- sample.kernel(x, u)
hist(sample, breaks = 20, probability = TRUE)
lines(density(sample), col = "blue")
lines(density(x), col = "red")



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








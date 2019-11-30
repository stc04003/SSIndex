#######################################################################
## Load package
#######################################################################

library(survival)
library(SSIndex)
library(methods)

#######################################################################
## Load function
#######################################################################

## Test only
do <- function(n, model, frailty = FALSE, type1 = FALSE, B = 200) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty, type1, offset = 5)
    p <- 2
    bi <- cbind(cos(seq(0, pi, length.out = 150)), sin(seq(0, pi, length.out = 150)))
    k0.tmp <- getk0s(dat, bi)
    k0 <- c(k0.tmp[1,], -k0.tmp[1,])
    k02 <- c(k0.tmp[2,], -k0.tmp[2,])
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    n <- length(unique(dat$id))
    getBootk <- function(dat) {
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
	dat0$Time[duplicated(dat0$Time) & dat0$Time < max(dat0$Time)] <-
            dat0$Time[duplicated(dat0$Time) & dat0$Time < max(dat0$Time)] +
            abs(rnorm(sum(duplicated(dat0$Time) & dat0$Time < max(dat0$Time)), sd = .0001))
        dat0$id <- rep(1:n, mm[ind] + 1)
        k0B.tmp <- getk0s(dat0, bi)
        k0B <- c(k0B.tmp[1,], -k0B.tmp[1,]) - k0
        k02B <- c(k0B.tmp[2,], -k0B.tmp[2,]) - k02
	c(max(k0B), max(k02B))
    }
    tmp <- replicate(B, getBootk(dat))
    c(mean(max(k0) > tmp[1,]), 
      mean(max(k02) > tmp[2,]))
}

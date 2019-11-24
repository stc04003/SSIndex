#######################################################################
## Load package
#######################################################################

library(survival)
library(SSIndex)
library(methods)

#######################################################################
## Load function
#######################################################################

sdOut <- function(x) {
    rm <- which(x %in% boxplot(x, plot = FALSE)$out)
    if (length(rm) > 0) x <- x[-rm]
    sd(x)
}

do <- function(n, model, frailty = FALSE, B = 200) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    p <- 2
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .05)), p)))
    r <- rowSums(tmp * tmp)
    keep <- which(r < 1 & r > 0)
    bi <- tmp[keep,] / sqrt(r[keep])
    k0.tmp <- sapply(1:NROW(bi), function(x) getk04(dat, bi[x,]))
    k0 <- k0.tmp[2,]
    k02 <- k0.tmp[1,]
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    n <- length(unique(dat$id))
    getBootk <- function(dat) {
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
	dat0$Time[duplicated(dat0$Time) & dat0$Time < max(dat0$Time)] <-
            dat0$Time[duplicated(dat0$Time) & dat0$Time < max(dat0$Time)] +
            abs(rnorm(sum(duplicated(dat0$Time) & dat0$Time < max(dat0$Time)), sd = .0001))
        dat0$id <- rep(1:n, mm[ind] + 1)
        fit0 <- gsm(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
                    data = dat0, shp.ind = FALSE, bIni = fit$b00, rIni = fit$r00)
        k0B.tmp <- sapply(1:NROW(bi), function(x) getk04(dat0, bi[x,]))
        k0B <- k0B.tmp[2,] - k0
        k02B <- k0B.tmp[1,] - k02
        ## c(max(k0B), max(k02B), fit0$b0, fit0$b00, fit0$r0, fit0$r00)
	c(fit0$b0, fit0$b00, fit0$r0, fit0$r00, max(k0B), max(k02B))
    }
    tmp <- replicate(B, getBootk(dat))
    c(fit$b0, fit$b00, fit$r0, fit$r00,
    apply(tmp[1:8,], 1, sd), 
    apply(tmp[1:8,], 1, sdOut), 
    mean(max(k0) > tmp[9,]), 
    mean(max(k02) > tmp[10,]))
}


library(parallel)

cl <- makePSOCKcluster(8)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do"))
invisible(clusterExport(NULL, "sdOut"))
invisible(clusterEvalQ(NULL, library(SSIndex)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(methods)))

s1 <- parSapply(NULL, 1:200, function(z) tryCatch(do(200, "M2"), error = function(e) do(200, "M2")))
s2 <- parSapply(NULL, 1:200, function(z) tryCatch(do(200, "M3"), error = function(e) do(200, "M3")))
s3 <- parSapply(NULL, 1:200, function(z) tryCatch(do(200, "M4"), error = function(e) do(200, "M4")))
s4 <- parSapply(NULL, 1:200, function(z) tryCatch(do(200, "M5"), error = function(e) do(200, "M5")))

stopCluster(cl)

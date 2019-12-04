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

getk05 <- function(dat, b) {
    time0 <- proc.time()
    if (any(dat$m == 0)) {
        tmp <- subset(dat, m == 0)
        tmp$Time <- 0
        tmp$event <- 1
        tmp$status <- 0
        dat <- rbind(dat, tmp)
    }
    dat <- dat[order(dat$id, dat$Time),]
    tid <- subset(dat, event == 1)$id
    tij <- subset(dat, event == 1)$Time
    yi <- subset(dat, event == 0)$Time
    m0 <- pmin(subset(dat, event == 1)$m, 1)
    Xij0 <- as.matrix(subset(dat, event == 0)[, grep("x", names(dat))])
    ## print(proc.time() - time0) ##### 
    tijyi0 <- outer(tij, yi, "<=")
    bxSgn0 <- outer(c(Xij0 %*% b), c(Xij0 %*% b), function(x, y) sign(x - y))
    ## print(proc.time() - time0) ##### 
    Xij <- Xij0[tid,]
    tijyi <- tijyi0[,tid]
    bxSgn <- bxSgn0[tid, tid]
    tikSgn <- outer(tij, tij, function(x, y) sign(x - y))
    ## print(proc.time() - time0) ##### 
    Nij <- do.call(rbind, lapply(split(tijyi0, tid), function(x)
        colSums(matrix(x, ncol = ncol(tijyi0)))))
    ## print(proc.time() - time0) ##### 
    k01 <- sum((m0 %*% t(m0)) * bxSgn * tikSgn * tijyi * t(tijyi)) / n / (n - 1)
    ## print(proc.time() - time0) ##### 
    k02 <- sum(bxSgn0 * (Nij - t(Nij))) / n / (n - 1)
    ## print(proc.time() - time0) ##### 
    c(k01, k02)
}


getk06 <- function(dat, b) {
    time0 <- proc.time()
    if (any(dat$m == 0)) {
        tmp <- subset(dat, m == 0)
        tmp$Time <- 0
        tmp$event <- 1
        tmp$status <- 0
        dat <- rbind(dat, tmp)
    }
    n <- length(unique(dat$id))
    dat <- dat[order(dat$id, dat$Time),]
    tid <- subset(dat, event == 1)$id
    tij <- subset(dat, event == 1)$Time
    yi <- subset(dat, event == 0)$Time
    m0 <- pmin(subset(dat, event == 1)$m, 1)
    Xij0 <- as.matrix(subset(dat, event == 0)[, grep("x|y", names(dat))])
    ## tijyi0 <- outer(tij, yi, "<=")
    ## bxSgn0 <- outer(c(Xij0 %*% fit$b0), c(Xij0 %*% fit$b0), function(x, y) sign(x - y))
    tijyi0 <- matrix(.C("outerC", as.double(tij), as.double(yi), as.integer(length(tij)),
                        as.integer(length(yi)), result = double(length(tij) * length(yi)),
                        PACKAGE = "SSIndex")$result, length(tij))
    bxSgn0 <- matrix(.C("outerCsign", as.double(Xij0 %*% b), as.double(Xij0 %*% b),
                        as.integer(n), as.integer(n), result = double(n * n),
                        PACKAGE = "SSIndex")$result, n)
    Xij <- Xij0[tid,]
    tijyi <- tijyi0[,tid]
    bxSgn <- bxSgn0[tid, tid]
    ## tikSgn <- outer(tij, tij, function(x, y) sign(x - y))
    tikSgn <- matrix(.C("outerCsign", as.double(tij), as.double(tij),
                    as.integer(length(tij)), as.integer(length(tij)),
                    result = double(length(tij) * length(tij)),
                    PACKAGE = "SSIndex")$result, length(tij))   
    Nij <- do.call(rbind, lapply(split(tijyi0, tid), function(x)
        colSums(matrix(x, ncol = ncol(tijyi0)))))
    k01 <- sum((m0 %*% t(m0)) * bxSgn * tikSgn * tijyi * t(tijyi)) / n / (n - 1)
    k02 <- sum(bxSgn0 * (Nij - t(Nij))) / n / (n - 1)
    c(k01, k02)
}


## Test only
do <- function(n, model, frailty = FALSE, type1 = FALSE, B = 200) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty, type1)
    p <- 2
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .1)), p)))
    r <- rowSums(tmp * tmp)
    keep <- which(r < 1 & r > 0)
    bi <- tmp[keep,] / sqrt(r[keep])
    k0.tmp <- getk0s(dat, bi) ## sapply(1:NROW(bi), function(x) getk04(dat, bi[x,]))
    k0 <- k0.tmp[1,]
    k02 <- k0.tmp[2,]
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
        k0B <- k0B.tmp[1,] - k0
        k02B <- k0B.tmp[2,] - k02
	c(max(k0B), max(k02B))
    }
    tmp <- replicate(B, getBootk(dat))
    c(mean(max(k0) > tmp[1,]), 
      mean(max(k02) > tmp[2,]))
}

set.seed(1)
system.time(print(do(200, "M2", B = 50)))

set.seed(1)
system.time(print(do(400, "M2", B = 50)))

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

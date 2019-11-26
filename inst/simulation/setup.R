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

## Test only
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
        k0B.tmp <- sapply(1:NROW(bi), function(x) getk04(dat0, bi[x,]))
        k0B <- k0B.tmp[2,] - k0
        k02B <- k0B.tmp[1,] - k02
	c(max(k0B), max(k02B))
    }
    tmp <- replicate(B, getBootk(dat))
    c(fit$b0, fit$b00, fit$r0, fit$r00,
    mean(max(k0) > tmp[9,]), 
    mean(max(k02) > tmp[10,]))
}

system.time(print(do(200, "M2")))

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


## ------------------------------------------------------------------------------------
## Debug
## ------------------------------------------------------------------------------------

n <- 5
## n <- 200
model <- "M2"
frailty  <- FALSE

## set.seed(961)
set.seed(9)

dat <- simDat(n, model, frailty)
## as.data.frame(dat)
fit <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
           data = dat, shp.ind = FALSE)
p <- 2
tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .05)), p)))
r <- rowSums(tmp * tmp)
keep <- which(r < 1 & r > 0)
bi <- tmp[keep,] / sqrt(r[keep])
## system.time(k0.tmp <- sapply(1:NROW(bi), function(x) getk04(dat, bi[x,])))

rev(getk04(dat, fit$b0))
getk05(dat, fit$b0)

rev(getk04(dat, fit$b0)) / getk05(dat, fit$b0)
system.time(print(getk04(dat, fit$b0))) ## 
system.time(print(getk05(dat, fit$b0))) ## 
system.time(print(getk06(dat, fit$b0))) ## 

microbenchmark::microbenchmark(getk04(dat, fit$b0),
                               getk05(dat, fit$b0),
                               getk06(dat, fit$b0))

e


## debug
## adding Tij = 0 for those with m = 0

if (any(dat$m == 0)) {
    tmp <- subset(dat, m == 0)
    tmp$Time <- 0
    tmp$event <- 1
    tmp$status <- 0
    dat <- rbind(dat, tmp)
}
dat <- dat[order(dat$id, dat$Time),]
## as.data.frame(dat)

tid <- subset(dat, event == 1)$id
tij <- subset(dat, event == 1)$Time
yi <- subset(dat, event == 0)$Time
m0 <- pmin(subset(dat, event == 1)$m, 1)
Xij0 <- as.matrix(subset(dat, event == 0)[, grep("x|y", names(dat))])

## tijyi0 <- outer(tij, yi, "<=")
tijyi0 <- matrix(.C("outerC", as.double(tij), as.double(yi), as.integer(length(tij)),
                    as.integer(length(yi)), result = double(length(tij) * length(yi)),
                    PACKAGE = "SSIndex")$result, length(tij))

## bxSgn0 <- outer(c(Xij0 %*% fit$b0), c(Xij0 %*% fit$b0), function(x, y) sign(x - y))
bxSgn0 <- matrix(.C("outerCsign", as.double(Xij0 %*% fit$b0), as.double(Xij0 %*% fit$b0),
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
Nij <- do.call(rbind, lapply(split(tijyi0, tid), function(x) colSums(matrix(x, ncol = ncol(tijyi0)))))
    ## apply(matrix(x, ncol = ncol(tijyi0)), 2, cumsum)))

## k0
sum((m0 %*% t(m0)) * bxSgn * tikSgn * tijyi * t(tijyi)) / n / (n - 1)

## k02
sum(bxSgn0 * (Nij - t(Nij))) / n / (n - 1)


mm <- aggregate(event ~ id, dat, sum)[, 2]
midx <- c(0, cumsum(mm)[-length(mm)])

sum(mm)

## k0 matrix
tmp0 <- (tikSgn * tijyi * t(tijyi))[outer(tid, tid, "<")]
tmp2 <- tikSgn * tijyi * t(tijyi) * outer(tid, tid, "<")
    
mat <- matrix(0, sum(mm), sum(mm))
mat[t(outer(tid, tid, "<"))] <- 
    .C("k0Mat", as.integer(n), as.integer(mm), as.integer(midx), 
       as.double(tij), as.double(yi), result = double(sum(outer(tid, tid, "<"))),
       PACKAGE = "SSIndex")$result
mat <- t(mat)
identical(mat, tmp2)

## k02 matrix
mat2 <- matrix(0, n, n)
mat2[t(outer(1:n, 1:n, "<"))] <- 
    .C("k02Mat", as.integer(n), as.integer(mm), as.integer(midx), 
       as.double(tij), as.double(yi), result = double(sum(1:(n - 1))),
       PACKAGE = "SSIndex")$result
mat2 <- t(mat2)
mat2 <- mat2 - t(mat2)


Xij0 %*% fit$b0

getk07 <- function(dat, bi) {
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
    Xij <- Xij0[tid,]
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    midx <- c(0, cumsum(mm)[-length(mm)])
    mat <- matrix(0, sum(mm), sum(mm))
    mat[t(outer(tid, tid, "<"))] <- 
        .C("k0Mat", as.integer(n), as.integer(mm), as.integer(midx), 
           as.double(tij), as.double(yi), result = double(sum(outer(tid, tid, "<"))),
           PACKAGE = "SSIndex")$result
    mat <- t(mat)
    mat2 <- matrix(0, n, n)
    mat2[t(outer(1:n, 1:n, "<"))] <- 
        .C("k02Mat", as.integer(n), as.integer(mm), as.integer(midx), 
           as.double(tij), as.double(yi), result = double(sum(1:(n - 1))),
           PACKAGE = "SSIndex")$result
    mat2 <- t(mat2)
    mat2 <- mat2 - t(mat2)
    ## need betas here
    bxSgn0 <- outer(c(Xij0 %*% fit$b0), c(Xij0 %*% fit$b0), function(x, y) sign(x - y))
    k0 <- sapply(1:NROW(bi), function(x) {
        xb <- c(Xij0 %*% bi[x,])
        bxSgn0 <- outer(xb, xb, function(x, y) sign(x - y))
        c(sum((m0 %*% t(m0)) * bxSgn0[tid, tid] * mat),
          sum(bxSgn0 * mat2))
    })
    return(k0 / n / (n - 1))
}

getk07(dat, bi)

debug(do)
do(200, "M3")



xb <- c(Xij0 %*% bi[1,])
bxSgn0 <- outer(xb, xb, function(x, y) sign(x - y))
c(sum((m0 %*% t(m0)) * bxSgn0[tid, tid] * mat) / n / (n - 1),
  sum(bxSgn0 * mat2) / n / (n - 1))

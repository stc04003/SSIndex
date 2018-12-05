library(GSM)

dim(simDat(100, "M1"))
dim(simDat(100, "M2"))
dim(simDat(100, "M3"))
dim(simDat(100, "M4"))

c(7, 24) / 25

do <- function(n, model, shp.ind = FALSE) {
    dat <- simDat(n, model)
    unlist(gsm(dat, shp.ind))
}

do(200, "M1")
do(200, "M2")
do(200, "M3")

sim1 <- t(replicate(100, do(200, "M1")))
sim2 <- t(replicate(100, do(200, "M2")))
sim3 <- t(replicate(100, do(200, "M3")))
sim4 <- t(replicate(100, do(200, "M4")))

summary(sim1)
summary(sim2)
summary(sim3)
summary(sim4)

Douglas06(simDat(100, "M2")) ## confirms data generation

set.seed(1);round(do(200, "M1"), 3) # 0.986 -0.168  0.180  0.984  0.181 -0.031  0.022  0.120 
set.seed(1);round(do(200, "M2"), 3) # -0.598 -0.801  0.641  0.767 -1.308 -1.752  0.128  0.154 
set.seed(1);round(do(200, "M3"), 3) # -0.650 -0.760 -0.554 -0.832 -2.989 -3.496 -1.780 -2.674 
set.seed(1);round(do(200, "M4"), 3) # -0.609 -0.793  0.328  0.945 -0.807 -1.050  0.050  0.143 

## Parallel computing on 1000 replications
library(parallel)
library(xtable)

sim1 <- sim1.2 <- sim2 <- sim3 <- sim4 <- NULL
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do')))
invisible(clusterEvalQ(NULL, library(GSM)))
set.seed(1)
sim1 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do(200, "M1"))), 8))
set.seed(1)
sim1.2 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do(200, "M1", TRUE))), 8))
sim2 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do(200, "M2"))), 8))
sim3 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do(200, "M3"))), 8))
sim4 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do(200, "M4"))), 8))
stopCluster(cl)

makeTab <- function(dat) {
    cbind(colMeans(dat)[1:4], apply(dat, 2, sd)[1:4])
}

tab <- cbind(makeTab(sim1), makeTab(sim1.2), makeTab(sim2), makeTab(sim3), makeTab(sim4))

print(xtable(tab, digits = 3), math.style.negative = TRUE)

## -------------------------------------------------------------------------------------
## Testing \beta_0 = 0
## -------------------------------------------------------------------------------------
library(GSM)

set.seed(1)
dat <- simDat(n = 100, model = "M1")

system.time(f1 <- gsm(dat, shp.ind = TRUE))
system.time(f2 <- gsm(dat, shp.ind = FALSE))
system.time(f3 <- gsm(dat, shp.ind = "test"))

do <- function(n, model, shp.ind = "test") {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model)
    set.seed(seed)
    tmp <- tryCatch(gsm(dat, shp.ind), error = function(e) NA)
    if (is.character(shp.ind)) {
        if (any(is.na(tmp)))
            return(c(rep(NA, 5), seed))
        else
            return(c(tmp$b0, tmp$r0, tmp$d / sd(tmp$dstar), seed))
    }
    return(c(tmp$b0, tmp$r0))
}

set.seed(1);round(do(100, "M1"), 3)
set.seed(1);round(do(100, "M2"), 3)
set.seed(1);round(do(100, "M3"), 3)
set.seed(1);round(do(100, "M4"), 3)

system.time(print(do(300, "M1")))
system.time(print(do(300, "M1", FALSE)))
system.time(print(do(300, "M2")))
system.time(print(do(300, "M3")))
system.time(print(do(300, "M4")))

library(parallel)
library(xtable)

sim1 <- sim2 <- sim3 <- NULL
## cl <- makePSOCKcluster(8)
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do')))
invisible(clusterEvalQ(NULL, library(GSM)))
sim1 <- t(parSapply(NULL, 1:100, function(z) do(300, "M1")))
sim1 <- t(parSapply(NULL, 1:1000, function(z)
    tryCatch(do(300, "M1", FALSE), error = function(e) rep(NA, 4))))
sim2 <- t(parSapply(NULL, 1:100, function(z) do(300, "M2")))
sim3 <- t(parSapply(NULL, 1:100, function(z) do(300, "M3")))
stopCluster(cl)

sum(abs(sim1[,5]) > qnorm(.975)) / nrow(sim1)
sum(abs(sim2[,5]) > qnorm(.975)) / nrow(sim2)
sum(abs(sim3[,5]) > qnorm(.975)) / nrow(sim3)

## ------------------------------------------------------------------------------------
## Testing \beta_0 = 0 with grid approach
## ------------------------------------------------------------------------------------
#'
#' Testing for 2 dimensional
#' 
#' @param B is the bootstrap size
#' @param len is the length of segments used for uniform sample
library(GSM)

do <- function(n, model, B = 200) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model)
    bi <- seq(0, pi, length = 100)
    b0 <- getb0(dat)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    getBootk <- function(dat) {
        n <- length(unique(dat$id))
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        ## max(sapply(1:100, function(x) getk0(dat0, c(cos(bi[x]), sin(bi[x])))))
        max(sapply(1:100, function(x) getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x]))
    }
    tmp <- replicate(B, getBootk(dat))
    list(k = getk0(dat, b0$bhat),
         kStar = k0, kBoot = tmp)
    ## c(b0$bhat, max(k0), tmp)
}

do2 <- function(n, model, B = 200) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model)
    bi <- seq(0, 2 * pi, length = 100)
    b0 <- getb0(dat)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    getBootk <- function(dat) {
        n <- length(unique(dat$id))
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        max(sapply(1:length(bi), function(x) getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x]))
    }
    tmp <- replicate(B, getBootk(dat))
    1 * (max(k0) > quantile(tmp, .95))
}

do3 <- function(n, model, B = 200) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model)
    bi <- seq(0, pi, length = 100)
    b0 <- getb0(dat)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    getBootk <- function(dat) {
        n <- length(unique(dat$id))
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        tmp <- sapply(1:length(bi), function(x) getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x])
        max(tmp, -tmp)
    }
    tmp <- replicate(B, getBootk(dat))
    1 * (max(k0, -k)0 > quantile(tmp, .95))
}

system.time(print(do2(100, "M1")))

library(parallel)
library(xtable)

sim150 <- sim250 <- sim350 <- sim450 <- NULL
sim1100 <- sim2100 <- sim3100 <- sim4100 <- NULL
sim1200 <- sim2200 <- sim3200 <- sim4200 <- NULL
## cl <- makePSOCKcluster(8)
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do', 'do2')))
invisible(clusterEvalQ(NULL, library(GSM)))

sim150 <- t(parSapply(NULL, 1:1000, function(z) do2(100, "M1")))
sim250 <- t(parSapply(NULL, 1:1000, function(z) do2(100, "M2")))
sim350 <- t(parSapply(NULL, 1:1000, function(z) do2(100, "M3")))
sim450 <- t(parSapply(NULL, 1:1000, function(z) do2(100, "M4")))

save(sim150, file = "sim150.RData")
save(sim250, file = "sim250.RData")
save(sim350, file = "sim350.RData")
save(sim450, file = "sim450.RData")

sim1100 <- t(parSapply(NULL, 1:1000, function(z) do2(100, "M1")))
sim2100 <- t(parSapply(NULL, 1:1000, function(z) do2(100, "M2")))
sim3100 <- t(parSapply(NULL, 1:1000, function(z) do2(100, "M3")))
sim4100 <- t(parSapply(NULL, 1:1000, function(z) do2(100, "M4")))

save(sim1100, file = "sim1100.RData")
save(sim2100, file = "sim2100.RData")
save(sim3100, file = "sim3100.RData")
save(sim4100, file = "sim4100.RData")

sim1200 <- t(parSapply(NULL, 1:1000, function(z) do2(200, "M1")))
sim2200 <- t(parSapply(NULL, 1:1000, function(z) do2(200, "M2")))
sim3200 <- t(parSapply(NULL, 1:1000, function(z) do2(200, "M3")))
sim4200 <- t(parSapply(NULL, 1:1000, function(z) do2(200, "M4")))
stopCluster(cl)

save(sim1200, file = "sim1200.RData")
save(sim2200, file = "sim2200.RData")
save(sim3200, file = "sim3200.RData")
save(sim4200, file = "sim4200.RData")

#######

sumSim <- function(n, model) {
    fname <- paste("output-", n, "-", model, "-test", sep = "")
    if (file.exists(fname)) dat <- read.table(fname)
    else stop("file name does not exist")
    ## dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, nrow(dat) / 50))))
    dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, 6))))
    mean(dat[,5] > qnorm(.95))
}

sumSim(100, "M1")
sumSim(500, "M1")

sumSim(100, "M2")
sumSim(500, "M2")

sumSim(100, "M3")
sumSim(500, "M3")
sumSim(1000, "M3")
sumSim(2000, "M3")


################################################################################################################################################
################################################################################################################################################

set.seed(2)
dat <- simDat(400, "M4")
n <- length(unique(dat$id))
mm <- aggregate(event ~ id, dat, sum)[,2]
tij <- subset(dat, event == 1)$t
yi <- subset(dat, event == 0)$t
midx <- c(0, cumsum(mm)[-length(mm)])
X <- as.matrix(subset(dat, event == 0, select = c(x1, x2)))
p <- ncol(X)
d <- dstar <- NULL

lik <- function(theta) {
    ## reduce the number of parameter by polar cordinate system
    bhat <- c(sin(theta),cos(theta))
    ## The estimating equation Sn needs Yi even for the m = 0's
    xb <- X %*% bhat
    ## h and h2 depend on the data, here I use some reasonable value
    h <- 0.25
    h2 <- .2
    
    ## This calculate R(C_i,x,\beta)
    Rhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
           as.double(xb), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "GSM")$result,
        X %*% bhat, yi))
    Rhat <- ifelse(is.na(Rhat), 0, Rhat) #### assign 0/0, Inf/Inf to 0
    Rhat1 <- rep(Rhat, times = table(dat$id)-1)
    
    Xbij <- as.matrix(dat[,4:5]) %*% (bhat)
    Xbij <- Xbij[dat$event==1]
    
    ## This calculate r(T_{ik},x,\beta)
    rhat <- unlist(mapply(FUN = function(x,t)
        .C("shapeFun3", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
           as.double(xb), as.double(x), as.double(t), 
           as.double(h), as.double(h2),
           result = double(1), PACKAGE = "GSM")$result,
        Xbij,tij))
    
    ## This calculate R(T_{ik},x,\beta)
    Rhat2 <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
           as.double(xb), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "GSM")$result,
        Xbij,tij))
    
    -sum(log(rhat)) + sum(Rhat2 - Rhat1)
}


## theta takes value on the interval [0,pi] for identifiability, 
res <- optimize(lik, interval = c(0,pi))
## to avoid local minimum, we probably need to 
## search the minimum on subintervals, and see which is the global optimum
res1 <- optimize(lik, interval = c(0,pi/2))
res2 <- optimize(lik, interval = c(pi/2,pi))

## Estimated \beta
c(sin(res$minimum), cos(res$minimum))
## [1] 0.6317432 0.7751777
  
  

## old 
## lik2 <- function(theta)
## {
##   bhat <- c(sin(theta),cos(theta))
##   #### The estimating equation Sn needs Yi even for the m = 0's
##   xb <- X %*% bhat
##   #### h <- 1.06 * sd(xb) * n^-.2
##   h <- 0.25
##   h2 <- .2
##   Fhat <- unlist(mapply(FUN = function(x, y)
##     .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
##        as.double(xb), as.double(x), as.double(y), as.double(h), 
##        result = double(1), PACKAGE = "GSM")$result,
##     X %*% bhat, yi))
##   Fhat <- ifelse(is.na(Fhat), 0, Fhat) #### assign 0/0, Inf/Inf to 0
##   Fhat <- exp(-Fhat)
##   Fhat1 <- rep(Fhat, times = table(dat$id)-1)
##   Xbij <- as.matrix(dat[,4:5]) %*% (bhat)
##   Xbij <- Xbij[dat$event==1]
##   Fhat2 <- unlist(mapply(FUN = function(x,t)
##     .C("shapeFun3", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
##        as.double(xb), as.double(x), as.double(t), 
##        as.double(h), as.double(h2),
##        result = double(1), PACKAGE = "GSM")$result,
##     Xbij,tij))
##   Fhat3 <- unlist(mapply(FUN = function(x, y)
##     .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
##        as.double(xb), as.double(x), as.double(y), as.double(h), 
##        result = double(1), PACKAGE = "GSM")$result,
##     Xbij,tij))
##   Fhat22 <- exp(-Fhat3)*Fhat2
##   -sum(log(Fhat22))+sum(log(Fhat1))
## }
## 

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

do <- function(n, model) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model)
    ## bi <- seq(0, pi, length = 100)
    bi <- seq(0, pi / 2, length = 100)
    b0 <- getb0(dat)
    k0 <- sapply(bi, function(x) -getk0(dat, c(cos(x), sin(x))))
    getBootk <- function(dat) {
        n <- length(unique(dat$id))
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        mm <- aggregate(event ~ id, dat, sum)[, 2]
        dat0$id <- rep(1:n, mm[ind] + 1)
        max(sapply(1:100, function(x) -getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x]))
    }
    tmp <- replicate(200, getBootk(dat))
    c(b0$bhat, max(k0), tmp)
}


do2 <- function(n, model) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model)
    b0 <- getb0(dat)
    bi <- seq(0, pi, length = 100)
    tmp <- sapply(bi, function(x) -getk0(dat, c(sin(x), cos(x))))
    c(b0$bhat, -getk0(dat, b0$bhat), tmp)
}

system.time(foo1 <- do(1000, "M1"))
foo1[3]
summary(foo1[4:203])

system.time(foo4 <- do(100, "M4"))

system.time(print(foo1 <- do(100, "M1")))
system.time(print(foo4 <- do(100, "M4")))


library(parallel)
library(xtable)

sim1 <- sim2 <- sim3 <- sim4 <- NULL
cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do', 'do2')))
invisible(clusterEvalQ(NULL, library(GSM)))
sim1 <- t(parSapply(NULL, 1:200, function(z) do(100, "M1")))
sim2 <- t(parSapply(NULL, 1:200, function(z) do(100, "M2")))
sim3 <- t(parSapply(NULL, 1:200, function(z) do(100, "M3")))
sim4 <- t(parSapply(NULL, 1:200, function(z) do(100, "M4")))
stopCluster(cl)

apply(sim1, 1, function(x) mean(x[3] < x[4:503]))
apply(sim2, 1, function(x) mean(x[3] < x[4:503]))
apply(sim3, 1, function(x) mean(x[3] < x[4:503]))
apply(sim4, 1, function(x) mean(x[3] < x[4:503]))

apply(sim1, 1, function(x) mean(x[3] > x[4:503]))

apply(sim1, 1, function(x) mean(x[3] > quantile(x[4:503], .95)))
apply(sim4, 1, function(x) mean(x[3] > quantile(x[4:503], .95)))


apply(sim1, 1, function(x) mean(x[3] < x[4:203]))
apply(sim2, 1, function(x) mean(x[3] < x[4:203]))
apply(sim3, 1, function(x) mean(x[3] < x[4:203]))
apply(sim4, 1, function(x) mean(x[3] < x[4:203]))

apply(sim1, 1, function(x) mean(x[3] > quantile(x[4:203], .95)))
apply(sim2, 1, function(x) mean(x[3] > quantile(x[4:203], .95)))
apply(sim3, 1, function(x) mean(x[3] > quantile(x[4:203], .95)))
apply(sim4, 1, function(x) mean(x[3] > quantile(x[4:203], .95)))

apply(sim1, 1, function(x) mean(x[3] > quantile(x[4:203] - x[3], .95)))
apply(sim2, 1, function(x) mean(x[3] > quantile(x[4:203] - x[3], .95)))
apply(sim3, 1, function(x) mean(x[3] > quantile(x[4:203] - x[3], .95)))
apply(sim4, 1, function(x) mean(x[3] > quantile(x[4:203] - x[3], .95)))

apply(sim1, 1, function(x) mean(x[3] > quantile(x[4:203], .95)))
apply(sim2, 1, function(x) mean(x[3] > quantile(x[4:203], .95)))
apply(sim3, 1, function(x) mean(x[3] > quantile(x[4:203], .95)))
apply(sim4, 1, function(x) mean(x[3] > quantile(x[4:203], .95)))

apply(sim1, 1, function(x) mean(x[3] < quantile(x[4:203] - x[3], .05)))
apply(sim2, 1, function(x) mean(x[3] < quantile(x[4:203] - x[3], .05)))
apply(sim3, 1, function(x) mean(x[3] < quantile(x[4:203] - x[3], .05)))
apply(sim4, 1, function(x) mean(x[3] < quantile(x[4:203] - x[3], .05)))

e
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

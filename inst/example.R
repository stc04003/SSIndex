library(GSM)
library(reReg)
## example

set.seed(1)
dat <- simDat(100, "M1")
gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = FALSE)$b0
gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = TRUE)
gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = "test")

set.seed(1)
dat <- simDat(100, "M1", TRUE)
gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = FALSE)$b0
gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = TRUE)
gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = "test")

head(dat)

## simulation and checks

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

####################################################################################3
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
    bi <- seq(0, 2 * pi, length = 100)
    fit <- gsm(dat)
    fit.indep <- gsm(dat, TRUE)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    getBootk <- function(dat) {
        n <- length(unique(dat$id))
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        fit0 <- gsm(dat0)
        fit1 <- gsm(dat0, TRUE)
        c(max(sapply(1:length(bi), function(x) getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x])),
          fit0$b0, fit0$r0, fit1$b0, fit1$r0)
    }
    tmp <- replicate(B, getBootk(dat))
    c(fit$b0, 1 * (max(k0) > quantile(tmp[1,], .95)), fit$r0, fit.indep$r0,
      sqrt(diag(var(t(tmp[2:3,])))),
      sqrt(diag(var(t(tmp[4:5,])))),
      sqrt(diag(var(t(tmp[8:9,])))))
}

do.r0 <- function(n, model) {
    dat <- simDat(n, model)
    gsm(dat, TRUE)$r0
}

r050M1 <- replicate(1000, do.r0(50, "M1"))
r050M2 <- replicate(1000, do.r0(50, "M2"))
r050M3 <- replicate(1000, do.r0(50, "M3"))
r050M4 <- replicate(1000, do.r0(50, "M4"))

r0100M1 <- replicate(500, do.r0(100, "M1"))
r0100M2 <- replicate(500, do.r0(100, "M2"))
r0100M3 <- replicate(500, do.r0(100, "M3"))
r0100M4 <- replicate(500, do.r0(100, "M4"))

sim150 <- t(parSapply(NULL, 1:1000, function(z) do2(50, "M1")))
sim250 <- t(parSapply(NULL, 1:1000, function(z) do2(50, "M2")))
sim350 <- t(parSapply(NULL, 1:1000, function(z) do2(50, "M3")))
sim450 <- t(parSapply(NULL, 1:1000, function(z) do2(50, "M4")))

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

sumSim <- function(n, model, indB0 = "test") {
    fname <- paste(c("results", n, model, B), collapse = "-")
    if (file.exists(fname)) dat <- read.table(fname)
    else stop("file name does not exist")
    ## dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, nrow(dat) / 50))))
    dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, 13))))
    if (is.logical(indB0)) {
        pwr <- mean(dat[,3])
        if (indB0) {
            PE.b0 <- c(0, 0)
            ESE.b0 <- c(0, 0)
            ASE.b0 <- c(0, 0)
            PE.r0 <- eval(parse(text = paste("rowMeans(r0", n, model, ")", sep = "")))
            ESE.r0 <- eval(parse(text = paste("apply(r0", n, model, ", 1, sd)", sep = "")))
            ASE.r0 <- colMeans(dat[,12:13])
        } else {
            PE.b0 <- colMeans(dat[,1:2])
            ESE.b0 <- apply(dat[,1:2], 2, sd)
            ASE.b0 <- colMeans(dat[,6:7])
            PE.r0 <- colMeans(dat[,4:5])
            ESE.r0 <- apply(dat[,4:5], 2, sd)
            ASE.r0 <- colMeans(dat[,8:9])
        }
    } else {
        PE.b0 <- colMeans(dat[,1:2])
        ESE.b0 <- apply(dat[,1:2], 2, sd)
        ASE.b0 <- colMeans(dat[,6:7])
        PE.r0 <- colMeans(dat[,4:5])
        ESE.r0 <- apply(dat[,4:5], 2, sd) * (1 - mean(dat[,3])) +
            apply(eval(parse(text = paste("t(r0", n, model, ")", sep = ""))), 2, sd) * mean(dat[,3]) 
        ASE.r0 <- colMeans(dat[,3] * dat[,8:9] + (1 - dat[,3]) * dat[,12:13])
    }
    cbind(PE.b0, ESE.b0, ASE.b0, PE.r0, ESE.r0, ASE.r0)
}

sumPwr <- function(n, model) {
    fname <- paste(c("results", n, model, B), collapse = "-")
    if (file.exists(fname)) dat <- read.table(fname)
    else stop("file name does not exist")
    ## dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, nrow(dat) / 50))))
    dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, 13))))
    mean(dat[,3])
}

tabPwr <- rbind(c(sumPwr(50, "M1"), sumPwr(100, "M1"), sumPwr(200, "M1"), sumPwr(500, "M1")),
                c(sumPwr(50, "M2"), sumPwr(100, "M2"), sumPwr(200, "M2"), sumPwr(500, "M2")),
                c(sumPwr(50, "M3"), sumPwr(100, "M3"), sumPwr(200, "M3"), sumPwr(500, "M3")),
                c(sumPwr(50, "M4"), sumPwr(100, "M4"), sumPwr(200, "M4"), sumPwr(500, "M4")))

makeTab <- function(n) {
    cbind(rbind(sumSim(n, "M1", FALSE)[,1:3], sumSim(n, "M1", FALSE)[,4:6]),
          rbind(sumSim(n, "M1", TRUE)[,1:3], sumSim(n, "M1", TRUE)[,4:6]),
          rbind(sumSim(n, "M2", FALSE)[,1:3], sumSim(n, "M2", FALSE)[,4:6]),
          rbind(sumSim(n, "M3", FALSE)[,1:3], sumSim(n, "M3", FALSE)[,4:6]),
          rbind(sumSim(n, "M4", FALSE)[,1:3], sumSim(n, "M4", FALSE)[,4:6]))
}

tab <- rbind(makeTab(50), makeTab(100), makeTab(200), makeTab(500))

makeTab2 <- function(n) {
    rbind(rbind(sumSim(n, "M1", FALSE)[,1:3], sumSim(n, "M1", FALSE)[,4:6]),
          rbind(sumSim(n, "M1", TRUE)[,1:3], sumSim(n, "M1", TRUE)[,4:6]),
          rbind(sumSim(n, "M2", FALSE)[,1:3], sumSim(n, "M2", FALSE)[,4:6]),
          rbind(sumSim(n, "M3", FALSE)[,1:3], sumSim(n, "M3", FALSE)[,4:6]),
          rbind(sumSim(n, "M4", FALSE)[,1:3], sumSim(n, "M4", FALSE)[,4:6]))
}

tab <- cbind(makeTab2(50), makeTab2(100), makeTab2(200), makeTab2(500))

library(xtable)
print(xtable(tab, digits = 3), math.style.negative = TRUE, include.rownames=FALSE)

sumSim(100, "M3")
sumSim(500, "M3")
sumSim(1000, "M3")
sumSim(2000, "M3")


sumSim2 <- function(n, model, indB0 = "test") {
    fname <- paste(c("results", n, model, B), collapse = "-")
    if (file.exists(fname)) dat <- read.table(fname)
    else stop("file name does not exist")
    ## dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, nrow(dat) / 50))))
    dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, 13))))
    if (is.logical(indB0)) {
        pwr <- mean(dat[,3])
        if (indB0) {
            PE.b0 <- colMeans(dat[,1:2])
            ESE.b0 <- apply(dat[,1:2], 2, sd)
            ASE.b0 <- colMeans(dat[,6:7])
            PE.r0 <- colMeans(dat[,4:5])
            ESE.r0 <- apply(dat[,4:5], 2, sd)
            ASE.r0 <- colMeans(dat[,12:13])
        } else {
            PE.b0 <- colMeans(dat[,1:2])
            ESE.b0 <- apply(dat[,1:2], 2, sd)
            ASE.b0 <- colMeans(dat[,6:7])
            PE.r0 <- colMeans(dat[,4:5])
            ESE.r0 <- apply(dat[,4:5], 2, sd)
            ASE.r0 <- colMeans(dat[,8:9])
        }
    } else {
        PE.b0 <- colMeans(dat[,1:2])
        ESE.b0 <- apply(dat[,1:2], 2, sd)
        ASE.b0 <- colMeans(dat[,6:7])
        PE.r0 <- colMeans(dat[,4:5])
        ESE.r0 <- apply(dat[,4:5], 2, sd)
        ASE.r0 <- colMeans(dat[,3] * dat[,8:9] + (1 - dat[,3]) * dat[,12:13])
    }
    cbind(PE.b0, ESE.b0, ASE.b0, PE.r0, ESE.r0, ASE.r0)
}




####################################################################################3
## Parallel computing on 1000 replications
library(parallel)
library(xtable)

do <- function(n, model, B = 200, frailty = FALSE) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    fit.indep <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
                     data = dat, shp.ind = TRUE)
    bi <- seq(0, 2 * pi, length = 100)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    getBootk <- function(dat) {
        n <- length(unique(dat$id))
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        fit0 <- gsm(reSurv(time1 = t, id = id, event = event, status = status) ~ x1 + x2,
                    data = dat0, shp.ind = FALSE)
        fit1 <- gsm(reSurv(time1 = t, id = id, event = event, status = status) ~ x1 + x2,
                    data = dat0, shp.ind = TRUE)
        c(max(sapply(1:length(bi), function(x) getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x])),
          fit0$b0, fit0$r0, fit1$b0, fit1$r0)
    }
    tmp <- replicate(B, getBootk(dat))
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta, (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0,
      1 * (max(k0) > quantile(tmp[1,], .95)), fit$r0, fit.indep$r0,
      sqrt(diag(var(t(tmp[2:3,])))),
      sqrt(diag(var(t(tmp[4:5,])))),
      ## sqrt(diag(var(t(tmp[6:7,])))),
      sqrt(diag(var(t(tmp[8:9,])))))
}

do2 <- function(n, model, B = 200, frailty = FALSE) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    fit.indep <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
                     data = dat, shp.ind = TRUE)
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta, (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0, fit$r0, fit.indep$r0)
}

do(100, "M1", B = 50)
do(50, "M4", B = 50)
do2(100, "M1", B = 50)
do2(50, "M4", B = 50)

sim150 <- sim250 <- sim350 <- sim450 <- NULL
sim1100 <- sim2100 <- sim3100 <- sim4100 <- NULL
sim1200 <- sim2200 <- sim3200 <- sim4200 <- NULL

sim1502 <- sim2502 <- sim3502 <- sim4502 <- NULL
sim11002 <- sim21002 <- sim31002 <- sim41002 <- NULL
sim12002 <- sim22002 <- sim32002 <- sim42002 <- NULL

cl <- makePSOCKcluster(8)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do', 'do2')))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

sim150 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(50, "M1"))), 6))
sim250 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(50, "M2"))), 6))
sim350 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(50, "M3"))), 6))
sim450 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(50, "M4"))), 6))

sim1100 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(100, "M1"))), 6))
sim2100 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(100, "M2"))), 6))
sim3100 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(100, "M3"))), 6))
sim4100 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(100, "M4"))), 6))

sim1200 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(200, "M1"))), 6))
sim2200 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(200, "M2"))), 6))
sim3200 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(200, "M3"))), 6))
sim4200 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(200, "M4"))), 6))

sim1502 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(50, "M1", TRUE))), 6))
sim2502 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(50, "M2", TRUE))), 6))
sim3502 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(50, "M3", TRUE))), 6))
sim4502 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(50, "M4", TRUE))), 6))

sim11002 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(100, "M1", TRUE))), 6))
sim21002 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(100, "M2", TRUE))), 6))
sim31002 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(100, "M3", TRUE))), 6))
sim41002 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(100, "M4", TRUE))), 6))

sim12002 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(200, "M1", TRUE))), 6))
sim22002 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(200, "M2", TRUE))), 6))
sim32002 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(200, "M3", TRUE))), 6))
sim42002 <- t(matrix(unlist(parLapply(NULL, 1:1e3, function(z) do2(200, "M4", TRUE))), 6))

stopCluster(cl)

makeTab2 <- function(dat) {
    cbind(colMeans(dat), apply(dat, 2, sd))
}

tab <- rbind(cbind(makeTab2(sim150), makeTab2(sim1100), makeTab2(sim1200)),
             cbind(makeTab2(sim250), makeTab2(sim2100), makeTab2(sim2200)),
             cbind(makeTab2(sim350), makeTab2(sim3100), makeTab2(sim3200)),
             cbind(makeTab2(sim450), makeTab2(sim4100), makeTab2(sim4200)),
             cbind(makeTab2(sim1502), makeTab2(sim11002), makeTab2(sim12002)),
             cbind(makeTab2(sim2502), makeTab2(sim21002), makeTab2(sim22002)),
             cbind(makeTab2(sim3502), makeTab2(sim31002), makeTab2(sim32002)),
             cbind(makeTab2(sim4502), makeTab2(sim41002), makeTab2(sim42002)))
             
print(xtable(tab, digits = 3), math.style.negative = TRUE, include.rownames = FALSE)

## X ~ truncated normal
## Without frailty
## True (0, 0), (0.28, 0.96)
0.005 & 0.692 & $-$0.038 & 0.681 & $-$0.023 & 0.694 \\ 
0.023 & 0.722 & $-$0.005 & 0.732 & $-$0.013 & 0.720 \\ 
0.253 & 0.130 & 0.262 & 0.087 & 0.265 & 0.060 \\ 
0.958 & 0.039 & 0.960 & 0.040 & 0.962 & 0.015 \\ 
0.251 & 0.132 & 0.261 & 0.078 & 0.264 & 0.053 \\ 
0.958 & 0.033 & 0.962 & 0.021 & 0.963 & 0.014 \\
## True (0.6, 0.8), (0.6, 0.8)
$-$0.542 & 0.346 & $-$0.580 & 0.214 & $-$0.596 & 0.122 \\ 
$-$0.707 & 0.296 & $-$0.765 & 0.179 & $-$0.788 & 0.094 \\ 
0.585 & 0.161 & 0.597 & 0.106 & 0.594 & 0.082 \\ 
0.786 & 0.119 & 0.788 & 0.104 & 0.796 & 0.083 \\ 
0.608 & 0.178 & 0.610 & 0.098 & 0.610 & 0.063 \\ 
0.759 & 0.151 & 0.782 & 0.082 & 0.788 & 0.055 \\
## True (0.6, 0.8), (0.6, 0.8)
$-$0.552 & 0.294 & $-$0.600 & 0.157 & $-$0.594 & 0.098 \\ 
$-$0.744 & 0.236 & $-$0.775 & 0.121 & $-$0.795 & 0.074 \\ 
$-$0.577 & 0.205 & $-$0.595 & 0.116 & $-$0.597 & 0.071 \\ 
$-$0.775 & 0.156 & $-$0.790 & 0.096 & $-$0.797 & 0.054 \\ 
$-$0.582 & 0.196 & $-$0.598 & 0.122 & $-$0.601 & 0.071 \\ 
$-$0.774 & 0.154 & $-$0.787 & 0.089 & $-$0.794 & 0.054 \\
## True (0.6, 0.8), (0.28, 0.96)
$-$0.535 & 0.344 & $-$0.570 & 0.223 & $-$0.592 & 0.113 \\ 
$-$0.701 & 0.322 & $-$0.771 & 0.174 & $-$0.793 & 0.090 \\ 
0.224 & 0.235 & 0.241 & 0.159 & 0.264 & 0.102 \\ 
0.942 & 0.088 & 0.956 & 0.041 & 0.958 & 0.043 \\ 
0.224 & 0.224 & 0.247 & 0.144 & 0.266 & 0.090 \\ 
0.946 & 0.071 & 0.957 & 0.046 & 0.959 & 0.025 \\

## With frailty
## True (0, 0), (0.28, 0.96)
$-$0.029 & 0.703 & $-$0.027 & 0.700 & $-$0.003 & 0.691 \\ 
$-$0.005 & 0.712 & $-$0.020 & 0.714 & $-$0.011 & 0.724 \\ 
0.253 & 0.136 & 0.261 & 0.078 & 0.265 & 0.064 \\ 
0.957 & 0.046 & 0.962 & 0.022 & 0.962 & 0.015 \\ 
0.250 & 0.124 & 0.259 & 0.081 & 0.263 & 0.054 \\ 
0.960 & 0.031 & 0.962 & 0.022 & 0.963 & 0.014 \\ 
## True (0.6, 0.8), (0.6, 0.8)
$-$0.527 & 0.348 & $-$0.587 & 0.208 & $-$0.596 & 0.125 \\ 
$-$0.716 & 0.299 & $-$0.765 & 0.166 & $-$0.787 & 0.099 \\ 
0.588 & 0.147 & 0.590 & 0.126 & 0.590 & 0.083 \\ 
0.783 & 0.139 & 0.789 & 0.115 & 0.797 & 0.095 \\ 
0.596 & 0.189 & 0.613 & 0.098 & 0.609 & 0.067 \\ 
0.768 & 0.142 & 0.779 & 0.092 & 0.787 & 0.069 \\
## True (0.6, 0.8), (0.6, 0.8)
$-$0.563 & 0.313 & $-$0.602 & 0.162 & $-$0.595 & 0.098 \\ 
$-$0.707 & 0.292 & $-$0.769 & 0.139 & $-$0.795 & 0.073 \\ 
$-$0.579 & 0.202 & $-$0.596 & 0.116 & $-$0.597 & 0.076 \\ 
$-$0.771 & 0.171 & $-$0.789 & 0.099 & $-$0.797 & 0.057 \\ 
$-$0.580 & 0.200 & $-$0.599 & 0.113 & $-$0.601 & 0.074 \\ 
$-$0.775 & 0.149 & $-$0.787 & 0.093 & $-$0.794 & 0.056 \\
## True (0.6, 0.8), (0.28, 0.96)
$-$0.540 & 0.363 & $-$0.573 & 0.219 & $-$0.591 & 0.111 \\ 
$-$0.690 & 0.318 & $-$0.771 & 0.175 & $-$0.795 & 0.084 \\ 
0.225 & 0.246 & 0.243 & 0.153 & 0.262 & 0.112 \\ 
0.936 & 0.115 & 0.957 & 0.040 & 0.957 & 0.058 \\ 
0.227 & 0.233 & 0.247 & 0.141 & 0.267 & 0.087 \\ 
0.941 & 0.088 & 0.958 & 0.035 & 0.959 & 0.023 \\ 


makeTab2(sim150)
makeTab2(sim250)
makeTab2(sim350)
makeTab2(sim450)

makeTab2(sim1100)
makeTab2(sim2100)
makeTab2(sim3100)
makeTab2(sim4100)

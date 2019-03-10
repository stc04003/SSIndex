library(GSM)
library(reReg)
## example


dat <- simDat(100, "M1")
gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = FALSE)
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

<<<<<<< HEAD
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

do <- function(n, model, B = 200) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model)
    fit <- gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = FALSE)
    fit.indep <- gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat, shp.ind = TRUE)
    bi <- seq(0, 2 * pi, length = 100)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    getBootk <- function(dat) {
        n <- length(unique(dat$id))
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        fit0 <- gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat0, shp.ind = FALSE)
        fit1 <- gsm(reSurv(t, id, event, status) ~ x1 + x2, data = dat0, shp.ind = TRUE)
        c(max(sapply(1:length(bi), function(x) getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x])),
          fit0$b0, fit0$r0, fit1$b0, fit1$r0)
    }
    tmp <- replicate(B, getBootk(dat))
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta, (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0, 1 * (max(k0) > quantile(tmp[1,], .95)), fit$r0, fit.indep$r0,
      sqrt(diag(var(t(tmp[2:3,])))),
      sqrt(diag(var(t(tmp[4:5,])))),
      ## sqrt(diag(var(t(tmp[6:7,])))),
      sqrt(diag(var(t(tmp[8:9,])))))
}

do(100, "M1", B = 50)

sim1 <- sim2 <- sim3 <- sim4 <- NULL
cl <- makePSOCKcluster(8)
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

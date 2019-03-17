library(GSM)
library(reReg)
## example

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
## March 10
####################################################################################3
library(parallel)
library(xtable)

## March 17, 2019
do <- function(n, model, frailty = FALSE) {
    B <- 200
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    bi <- seq(0, 2 * pi, length = 100)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    k02 <- sapply(bi, function(x) getk02(dat, c(cos(x), sin(x)), fit$Fhat))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    n <- length(unique(dat$id))
    getBootk <- function(dat) {
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        fit0 <- gsm(reSurv(time1 = t, id = id, event = event, status = status) ~ x1 + x2,
                    data = dat0, shp.ind = FALSE)
        c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
          max(sapply(1:length(bi), function(x)
              getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x])),
          max(sapply(1:length(bi), function(x)
              getk02(dat0, c(cos(bi[x]), sin(bi[x])), fit$Fhat) - k02[x])))
    }
    tmp <- replicate(B, getBootk(dat))
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta, (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0, fit$b00, fit$r0, fit$r00, 
      sqrt(diag(var(t(tmp[1:2,])))),
      sqrt(diag(var(t(tmp[3:4,])))),
      sqrt(diag(var(t(tmp[5:6,])))),
      sqrt(diag(var(t(tmp[7:8,])))), 
      1 * (max(k0) > quantile(tmp[9,], .95)),
      1 * (max(k0) > quantile(tmp[10,], .95)))
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

system.time(print(do(100, "M1", TRUE)))
do(50, "M4", FALSE)
do2(100, "M1", FALSE)
do2(50, "M4", FALSE)

cl <- makePSOCKcluster(8)
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do', 'do2')))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

runPara <- function(model, n, frailty) {
    obj <- paste(c(model, n, frailty), collapse = "")
    eval(parse(text = paste(obj, " <- NULL")))
    toRun <- paste(obj, " <- t(matrix(unlist(parLapply(NULL, 1:500, function(z) do(",
                   n, ",'", model, "',", frailty, "))), 13))", sep = "")
    ## toRun <- paste(obj, " <- parSapply(NULL, 1:500, function(z) do(",
    ##                n, ",'", model, "',", frailty, "))", sep = "")
    ptm <- proc.time()
    eval(parse(text = toRun))
    print(proc.time() - ptm)
    fname <- paste(obj, ".RData", sep = "")
    eval(parse(text = paste("save(", obj, ",file = '", fname, "')", sep = "")))
}

runPara("M1", 50, TRUE)
runPara("M1", 100, TRUE)
runPara("M1", 200, TRUE)

runPara("M2", 50, TRUE)
runPara("M2", 100, TRUE)
runPara("M2", 200, TRUE)

runPara("M3", 50, TRUE)
runPara("M3", 100, TRUE)
runPara("M3", 200, TRUE)

runPara("M4", 50, TRUE)
runPara("M4", 100, TRUE)
runPara("M4", 200, TRUE)

runPara("M1", 50, FALSE)
runPara("M1", 100, FALSE)
runPara("M1", 200, FALSE)

runPara("M2", 50, FALSE)
runPara("M2", 100, FALSE)
runPara("M2", 200, FALSE)

runPara("M3", 50, FALSE)
runPara("M3", 100, FALSE)
runPara("M3", 200, FALSE)

runPara("M4", 50, FALSE)
runPara("M4", 100, FALSE)
runPara("M4", 200, FALSE)

stopCluster(cl)

makeTab2 <- function(dat) {
    cbind(colMeans(dat), apply(dat, 2, sd))
}




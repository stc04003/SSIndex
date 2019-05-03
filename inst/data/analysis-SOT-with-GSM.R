############################################################################################
## Load package and data
## SOT data with kidney transplant
############################################################################################

library(parallel)
library(GSM)
library(reReg)
library(survival)
load("dat.SOT.RData")

## exract baseline matrix
base.SOT <- dat.SOT[cumsum(lapply(with(dat.SOT, split(id, id)), length)),]
## scale age
base.SOT$scaleAge <- scale(base.SOT$age)
base.SOT$age01 <- 1 * (base.SOT$age > 65)
base.SOT$m <- aggregate(event ~ id, dat.SOT, sum)[,2]

## put the scaled age back
dat.SOT$scaleAge <- base.SOT$scaleAge[dat.SOT$id]
dat.SOT$age01 <- base.SOT$age01[dat.SOT$id]
dat.SOT$m <- base.SOT$m[dat.SOT$id]

## event plots
with(dat.SOT, plot(reSurv(Time, id, event, status))) ## Only has 6 death
with(dat.SOT, plotEvents(reSurv(Time, id, event, status) ~ race))

############################################################################################
## Baseline information
############################################################################################

dim(base.SOT)
head(dat.SOT)

## median follow-up times in days
survfit(Surv(Time, rep(1, nrow(base.SOT))) ~ 1, data = base.SOT)

## total number of infection episodes
sum(dat.SOT$event)

## number of deaths
sum(base.SOT$status)

## age
summary(base.SOT$age)

## race
table(base.SOT$race)
summary(base.SOT$race)

## HLA
table(base.SOT$HLA.incomp)
summary(base.SOT$HLA.incomp)

dat0 <- dat.SOT[,c(1:4, 10, 5:9)]
head(dat0)
###############################################################################################
## Analysis with GSM (Sinica paper)
## With scaled age
###############################################################################################

pValShape <- function(fname, B = 100, dat0 = dat0) {
    xNames <- attr(terms(fname), "term.labels")
    p <- length(attr(terms(fname), "term.labels"))
    dat1 <- dat0
    xCol <- as.numeric(sapply(xNames, function(x) which(names(dat0) == x)))
    colnames(dat1)[xCol] <- paste("x", 1:p, sep = "")
    dat1 <- dat1[,c(1:5, xCol)]
    ## head(dat1)
    ## bi <- as.matrix(expand.grid(rep(list(seq(0, 2 * pi, length = 100)), p - 1)))
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .1)), p)))
    r <- apply(tmp, 1, function(z) sqrt(sum(z^2)))
    bi <- (tmp / r)[r < 1 & r > 0,]
    k0 <- sapply(1:NROW(bi), function(x) getk0(dat1, bi[x,]))
    getBootK <- function(dat) {
        n <- length(unique(dat$id))
        mm <- aggregate(event ~ id, dat, length)[, 2]
        ind <- sample(1:n, n, TRUE)
        datB <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        datB$id <- rep(1:n, mm[ind])
        rownames(datB) <- NULL
        datB1 <- datB
        xCol <- as.numeric(sapply(xNames, function(x) which(names(datB) == x)))
        colnames(datB1)[xCol] <- paste("x", 1:p, sep = "")
        datB1 <- datB1[,c(1:5, xCol)]
        kb <- max(sapply(1:NROW(bi), function(x) getk0(datB1, bi[x,]) - k0[x]))
            ## getk0(datB1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1)) - k0[x]))
        max(kb)
    }
    cl <- makePSOCKcluster(8)
    ## cl <- makePSOCKcluster(16)
    setDefaultCluster(cl)
    invisible(clusterExport(cl, c("bi", "k0", "fname", "dat0", "xNames", "p", "getBootK"),
                            environment()))
    invisible(clusterEvalQ(NULL, library(GSM)))
    invisible(clusterEvalQ(NULL, library(reReg)))
    system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootK(dat0))) 
    stopCluster(cl)
    mean(max(k0) > tmp)
}

system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + race, 50, dat0))) # 0.64
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + HLA.incomp, 50, dat0))) # 0.28
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + HLA.incomp + race, 50, dat0))) # 0.64

################################################################################################################
## With CMV
load("SOT_kidney.RData")
kid <- kid[complete.cases(kid),]
mm <- aggregate(event ~ id, kid, sum)[, 2]
kid$id <- rep(1:length(unique(kid$id)), mm + 1)
kid$m <- unlist(aggregate(event ~ id, kid, function(x) rep(sum(x), length(x)))[,2])
kid <- kid[,c(1:4, 11, 5:10)]
head(kid)

system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + CMV, 50, kid))) # 0.34
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + HLA, 100, kid))) # 0.65
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + CMV + HLA, 100, kid))) # 0.7
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + CMV + race, 100, kid))) # 0.3
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + CMV + diabetes, 100, kid))) # 0.5
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + CMV + hypertension, 50, kid))) # 0.38
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + CMV + diabetes + hypertension, 50, kid))) # 0.62
system.time(print(pValShape(reSurv(Time, id, event, status) ~ CMV + diabetes + hypertension, 50, kid))) # 0.94

system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + HLA + race, 100, kid))) # 0.5
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + HLA + diabetes, 100, kid))) # 0.67
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + HLA + hypertension, 50, kid))) # 0.54
system.time(print(pValShape(reSurv(Time, id, event, status) ~ scaleAge + HLA + diabetes + hypertension, 50, kid))) # 0.

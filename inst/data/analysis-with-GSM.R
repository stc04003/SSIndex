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


###############################################################################################
## Analysis with reReg (Sinica paper)
###############################################################################################

fit <- gsm(reSurv(Time, id, event, status) ~ scaleAge + race + HLA.incomp, data = dat.SOT)

## Custom function for this data set, need to generalize this later...
dat.SOT2 <- dat.SOT
names(dat.SOT2)[c(8, 6, 7)] <- c("x1", "x2", "x3")
dat.SOT2 <- dat.SOT2[,c(1:5, 8, 6:7, 9:10)]
head(dat.SOT2)

bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 50), seq(0, 2 * pi, length = 50)))
k0 <- sapply(1:NROW(bi), function(x) getk0(dat.SOT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))))
k02 <- sapply(1:NROW(bi), function(x) getk02(dat.SOT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0))

B <- 1000
mm <- aggregate(event ~ id, dat.SOT, sum)[, 2]
n <- length(unique(dat.SOT$id))
getBootk <- function(dat.SOT) {
    ind <- sample(1:n, replace = TRUE)
    dat.SOT0 <- dat.SOT[unlist(sapply(ind, function(x) which(dat.SOT$id %in% x))),]
    dat.SOT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.SOT0) <- NULL
    fit0 <- gsm(reSurv(time1 = Time, id = id, event = event, status = status) ~ scaleAge + race + HLA.incomp,
                data = dat.SOT0, shp.ind = FALSE)
    names(dat.SOT0)[c(8, 6, 7)] <- c("x1", "x2", "x3")
    dat.SOT0 <- dat.SOT0[,c(1:5, 8, 6:7, 9:10)]    
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.SOT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.SOT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit0$Fhat0) - k02[x])))
}

set.seed(1)    
system.time(print(getBootk(dat.SOT)))

cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "getBootk"))
invisible(clusterExport(NULL, c("n", "mm", "dat.SOT", "bi", "k0", "k02")))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))
set.seed(1)
system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootk(dat.SOT))) ## 450 seconds for B = 400; 2202 seconds for B = 1000
stopCluster(cl)

1 * (max(k0) > quantile(tmp[13,], .95)) ## 0
1 * (max(k02) > quantile(tmp[14,], .95)) ## 0

sqrt(diag(var(t(tmp[1:3,]))))
sqrt(diag(var(t(tmp[1:3 + 3,]))))
sqrt(diag(var(t(tmp[1:3 + 3 * 2,]))))
sqrt(diag(var(t(tmp[1:3 + 3 * 3,]))))
mean(max(k0) > tmp[13,]) ## .760
mean(max(k02) > tmp[14,]) ## .955

## > sqrt(diag(var(t(tmp[1:3,]))))
## [1] 0.4034014 0.4172543 0.4465839
## > sqrt(diag(var(t(tmp[1:3 + 3,]))))
## [1] 0.4035078 0.4146305 0.4365111
## > sqrt(diag(var(t(tmp[1:3 + 3 * 2,]))))
## [1] 0.12219159 0.21317145 0.03398074
## > sqrt(diag(var(t(tmp[1:3 + 3 * 3,]))))
## [1] 0.12968396 0.21964168 0.04322522

## tmp <- replicate(B, getBootk(dat))

system.time(print(getBootk(dat.SOT)))

e

###############################################################################################
## Analysis with reReg (Sinica paper)
## Tab 3, upper panel, SOT:
###############################################################################################
set.seed(1)
## 0's as initial value
fit.SOT0 <- reReg(reSurv(Time, id, event, status) ~ scaleAge + race + HLA.incomp,
             data = dat.SOT, method = "sc.XCYH", se = "resampling", control = list(B = 500))
summary(fit.SOT0)


fit.SOT.HW <- reReg(reSurv(Time, id, event, status) ~ scaleAge + race + HLA.incomp,
                data = dat.SOT, method = "cox.HW", se = "bootstrap", control = list(B = 500))
summary(fit.SOT.HW)

## HW as initial value
fit.SOT2 <- reReg(reSurv(Time, id, event, status) ~ scaleAge + race + HLA.incomp,
                  data = dat.SOT, method = "sc.XCYH", se = "resampling",
                  control = list(B = 500, b0 = coef(fit.SOT.HW)[1:3]))
summary(fit.SOT2)

## Gehan weight for a and HW for b
fit.SOT3 <- reReg(reSurv(Time, id, event, status) ~ scaleAge + race + HLA.incomp,
                  data = dat.SOT, method = "sc.XCYH", se = "resampling",
                  control = list(B = 500, b0 = coef(fit.SOT.HW)[1:3], eqType = "Gehan"))
summary(fit.SOT3)

fit.SOT4 <- reReg(reSurv(Time, id, event, status) ~ scaleAge + race + HLA.incomp,
                  data = dat.SOT, method = "sc.XCYH", se = "resampling",
                  control = list(B = 500, a0 = coef(fit.SOT3)[1:3], b0 = coef(fit.SOT.HW)[1:3]))
summary(fit.SOT4)

## scale-age^2

summary(reReg(reSurv(Time, id, event, status) ~ scaleAge + I(scaleAge^2) + race + HLA.incomp,
              data = dat.SOT, method = "sc.XCYH", se = "resampling", control = list(B = 500)))

summary(reReg(reSurv(Time, id, event, status) ~ scaleAge + I(scaleAge^2) + race + HLA.incomp + scaleAge:race + scaleAge:HLA.incomp,
              data = dat.SOT, method = "sc.XCYH", se = "resampling", control = list(B = 500)))


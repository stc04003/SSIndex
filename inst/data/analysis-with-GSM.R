############################################################################################
## Load package and data
## SOT data with kidney transplant
############################################################################################

library(GSM)
library(reReg)
library(survival)
load("dat.SOT.RData")

## exract baseline matrix
base.SOT <- dat.SOT[cumsum(lapply(with(dat.SOT, split(id, id)), length)),]
## scale age
base.SOT$scaleAge <- scale(base.SOT$age)
base.SOT$age01 <- 1 * (base.SOT$age > 65)

## put the scaled age back
dat.SOT$scaleAge <- base.SOT$scaleAge[dat.SOT$id]
dat.SOT$age01 <- base.SOT$age01[dat.SOT$id]

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

gsm(reSurv(Time, id, event, status) ~ scaleAge + race + HLA.incomp, data = dat.SOT)
    
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


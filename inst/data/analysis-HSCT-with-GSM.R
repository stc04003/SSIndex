############################################################################################
## Load package and data
## SOT data with kidney transplant
############################################################################################

library(parallel)
library(GSM)
library(reReg)
library(survival)

load("dat.HSCT.RData")
head(dat.HSCT)
summary(dat.HSCT) ## Missing values at CMVstatus

## exract baseline matrix
base.HSCT <- dat.HSCT[cumsum(lapply(with(dat.HSCT, split(id, id)), length)),]
## scale age
base.HSCT$scaleAge <- scale(base.HSCT$age)
base.HSCT$age01 <- 1 * (base.HSCT$age > 65)
base.HSCT$m <- aggregate(event ~ id, dat.HSCT, sum)[,2]
## put the scaled age back
dat.HSCT$scaleAge <- base.HSCT$scaleAge[dat.HSCT$id]
dat.HSCT$age01 <- base.HSCT$age01[dat.HSCT$id]
dat.HSCT$m <- base.HSCT$m[dat.HSCT$id]
## make gender 0-1
dat.HSCT$gender <- dat.HSCT$gender - 1

###############################################################################################
## Analysis with scaled age
###############################################################################################

fname <- reSurv(Time, id, event, status) ~ scaleAge + allo + race
fit <- gsm(fname, data = dat.HSCT)

## Custom function for this data set, need to generalize this later...

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("scaleAge", "allo", "race"), function(x) which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(10, 7, 6)] <- c("x1", "x2", "x3")
dat.HSCT2 <- dat.HSCT2[,c(1:5, 10, 7, 6, 8:9, 11:12)]
head(dat.HSCT2)

bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 50), seq(0, 2 * pi, length = 50)))
k0 <- sapply(1:NROW(bi), function(x) getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))))
k02 <- sapply(1:NROW(bi), function(x) getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0))


B <- 1000
mm <- aggregate(event ~ id, dat.HSCT, sum)[, 2]
n <- length(unique(dat.HSCT$id))
getBootk <- function(dat.HSCT) {
    ind <- sample(1:n, replace = TRUE)
    dat.HSCT0 <- dat.HSCT[unlist(sapply(ind, function(x) which(dat.HSCT$id %in% x))),]
    dat.HSCT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.HSCT0) <- NULL
    fit0 <- gsm(fname, data = dat.HSCT0, shp.ind = FALSE)
    names(dat.HSCT0)[c(10, 7, 6)] <- c("x1", "x2", "x3")
    dat.HSCT0 <- dat.HSCT0[,c(1:5, 10, 7, 6, 8:9, 11:12)]
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit0$Fhat0) - k02[x])))
}

set.seed(1)    
system.time(print(getBootk(dat.HSCT)))

## cl <- makePSOCKcluster(16)
cl <- makePSOCKcluster(8)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "getBootk"))
invisible(clusterExport(NULL, c("n", "mm", "dat.HSCT", "bi", "k0", "k02", "fname")))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

set.seed(1)
system.time(tmp <- parSapply(NULL, 1:200, function(z) getBootk(dat.HSCT))) ## B = 200 takes 1097 seconds.
stopCluster(cl)


1 * (max(k0) > quantile(tmp[13,], .95)) ## 0; don't reject H0
1 * (max(k02) > quantile(tmp[14,], .95)) ## 0


###############################################################################################
## age01, allo, gender
###############################################################################################

fname <- reSurv(Time, id, event, status) ~ age01 + allo + gender
fit <- gsm(fname, data = dat.HSCT)

## Custom function for this data set, need to generalize this later...

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("age01", "allo", "gender"), function(x) which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(11, 7, 8)] <- c("x1", "x2", "x3")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 12, 11, 7, 8)]
head(dat.HSCT2)

bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 50), seq(0, 2 * pi, length = 50)))
k0 <- sapply(1:NROW(bi), function(x) getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))))
k02 <- sapply(1:NROW(bi), function(x) getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0))
max(k0)
max(k02)

B <- 1000
mm <- aggregate(event ~ id, dat.HSCT, sum)[, 2]
n <- length(unique(dat.HSCT$id))
getBootk <- function(dat.HSCT) {
    ind <- sample(1:n, replace = TRUE)
    dat.HSCT0 <- dat.HSCT[unlist(sapply(ind, function(x) which(dat.HSCT$id %in% x))),]
    dat.HSCT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.HSCT0) <- NULL
    fit0 <- gsm(fname, data = dat.HSCT0, shp.ind = FALSE)
    names(dat.HSCT0)[c(11, 7, 8)] <- c("x1", "x2", "x3")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 12, 11, 7, 8)]
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit0$Fhat0) - k02[x])))
}

set.seed(1)    
system.time(print(getBootk(dat.HSCT)))

cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "getBootk"))
invisible(clusterExport(NULL, c("n", "mm", "dat.HSCT", "bi", "k0", "k02", "fname")))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

set.seed(1)
system.time(tmp <- parSapply(NULL, 1:300, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)


1 * (max(k0) > quantile(tmp[13,], .95)) ## 1; reject H0; shape dependence
1 * (max(k02) > quantile(tmp[14,], .95)) ## 0

mean(max(k0) > tmp[13,]) 
mean(max(k02) > tmp[14,])

summary(k0)
summary(k02)


max(k0)
max(k02)
summary(tmp[13,])
summary(tmp[14,])

###############################################################################################
## Analysis with CMV status; there are missing values
###############################################################################################

dat.HSCT$id[which(is.na(dat.HSCT$CMVstatus))]
table(dat.HSCT$id[which(is.na(dat.HSCT$CMVstatus))])
sapply(unique(dat.HSCT$id[which(is.na(dat.HSCT$CMVstatus))]), function(x) sum(dat.HSCT$id == x))


dat.HSCT <- dat.HSCT[complete.cases(dat.HSCT),]
dat.HSCT$id <- rep(1:length(unique(dat.HSCT$id)), aggregate(event ~ id, dat.HSCT, sum)[, 2] + 1)
head(dat.HSCT)



fname <- reSurv(Time, id, event, status) ~ allo + gender + CMVstatus + scaleAge
fit <- gsm(fname, data = dat.HSCT)
str(fit)

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("allo", "gender", "CMVstatus", "scaleAge"), function(x)
    which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[7:10] <- c("x1", "x2", "x3", "x4")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 7:10, 12)]
head(dat.HSCT2)

b0 <- seq(0, 2 * pi, length = 50)
bi <- as.matrix(expand.grid(b0, b0, b0))
k0 <- sapply(1:NROW(bi), function(x) getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))))
k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0))
max(k0)
max(k02)


B <- 1000
mm <- aggregate(event ~ id, dat.HSCT, sum)[, 2]
n <- length(unique(dat.HSCT$id))
getBootk <- function(dat.HSCT) {
    ind <- sample(1:n, replace = TRUE)
    dat.HSCT0 <- dat.HSCT[unlist(sapply(ind, function(x) which(dat.HSCT$id %in% x))),]
    dat.HSCT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.HSCT0) <- NULL
    fit0 <- gsm(fname, data = dat.HSCT0, shp.ind = FALSE)
    names(dat.HSCT0)[7:10] <- c("x1", "x2", "x3", "x4")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 7:10, 12)]
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit0$Fhat0) - k02[x])))
}

set.seed(1)    
system.time(print(getBootk(dat.HSCT)))

cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "getBootk"))
invisible(clusterExport(NULL, c("n", "mm", "dat.HSCT", "bi", "k0", "k02", "fname")))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

set.seed(1)
system.time(tmp <- parSapply(NULL, 1:200, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)


1 * (max(k0) > quantile(tmp[13,], .95)) ## 0; fail to reject H0; shape independence
1 * (max(k02) > quantile(tmp[14,], .95)) ## 0

mean(max(k0) > tmp[13,])
## [1] 0.025
mean(max(k02) > tmp[14,])
## [1] 1


###############################################################################################
## Analysis with CMV status; there are missing values
## binary age
###############################################################################################
head(dat.HSCT)

fname <- reSurv(Time, id, event, status) ~ allo + gender + CMVstatus + age01
fit <- gsm(fname, data = dat.HSCT)
str(fit)


dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("allo", "gender", "CMVstatus", "age01"), function(x)
    which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(7:9, 11)] <- c("x1", "x2", "x3", "x4")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 7:9, 11, 12)]
head(dat.HSCT2)

b0 <- seq(0, 2 * pi, length = 50)
bi <- as.matrix(expand.grid(b0, b0, b0))
k0 <- sapply(1:NROW(bi), function(x) getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))))
k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0))
max(k0)
max(k02)


B <- 1000
mm <- aggregate(event ~ id, dat.HSCT, sum)[, 2]
n <- length(unique(dat.HSCT$id))
getBootk <- function(dat.HSCT) {
    ind <- sample(1:n, replace = TRUE)
    dat.HSCT0 <- dat.HSCT[unlist(sapply(ind, function(x) which(dat.HSCT$id %in% x))),]
    dat.HSCT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.HSCT0) <- NULL
    fit0 <- gsm(fname, data = dat.HSCT0, shp.ind = FALSE)
    names(dat.HSCT0)[c(7:9, 11)] <- c("x1", "x2", "x3", "x4")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 7:9, 11, 12)]
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit0$Fhat0) - k02[x])))
}

set.seed(1) 
system.time(print(getBootk(dat.HSCT)))

cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "getBootk"))
invisible(clusterExport(NULL, c("n", "mm", "dat.HSCT", "bi", "k0", "k02", "fname")))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

set.seed(1)
system.time(tmp <- parSapply(NULL, 1:100, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)



1 * (max(k0) > quantile(tmp[13,], .95)) ## 0
1 * (max(k02) > quantile(tmp[14,], .95)) ## 1

mean(max(k0) > tmp[13,])
## [1] 0.025
mean(max(k02) > tmp[14,])
## [1] 1



###############################################################################################
## Analysis with CMV status; there are missing values
###############################################################################################

dat.HSCT$id[which(is.na(dat.HSCT$CMVstatus))]
table(dat.HSCT$id[which(is.na(dat.HSCT$CMVstatus))])
sapply(unique(dat.HSCT$id[which(is.na(dat.HSCT$CMVstatus))]), function(x) sum(dat.HSCT$id == x))


dat.HSCT <- dat.HSCT[complete.cases(dat.HSCT),]
dat.HSCT$id <- rep(1:length(unique(dat.HSCT$id)), aggregate(event ~ id, dat.HSCT, sum)[, 2] + 1)
head(dat.HSCT)



fname <- reSurv(Time, id, event, status) ~ allo + gender + CMVstatus + scaleAge
fit <- gsm(fname, data = dat.HSCT)
str(fit)

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("allo", "gender", "CMVstatus", "scaleAge"), function(x)
    which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[7:10] <- c("x1", "x2", "x3", "x4")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 7:10, 12)]
head(dat.HSCT2)

b0 <- seq(0, 2 * pi, length = 50)
bi <- as.matrix(expand.grid(b0, b0, b0))
k0 <- sapply(1:NROW(bi), function(x) getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))))
k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0))
max(k0)
max(k02)


B <- 1000
mm <- aggregate(event ~ id, dat.HSCT, sum)[, 2]
n <- length(unique(dat.HSCT$id))
getBootk <- function(dat.HSCT) {
    ind <- sample(1:n, replace = TRUE)
    dat.HSCT0 <- dat.HSCT[unlist(sapply(ind, function(x) which(dat.HSCT$id %in% x))),]
    dat.HSCT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.HSCT0) <- NULL
    fit0 <- gsm(fname, data = dat.HSCT0, shp.ind = FALSE)
    names(dat.HSCT0)[7:10] <- c("x1", "x2", "x3", "x4")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 7:10, 12)]
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit0$Fhat0) - k02[x])))
}

set.seed(1)    
system.time(print(getBootk(dat.HSCT)))

cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "getBootk"))
invisible(clusterExport(NULL, c("n", "mm", "dat.HSCT", "bi", "k0", "k02", "fname")))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

set.seed(1)
system.time(tmp <- parSapply(NULL, 1:200, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)


1 * (max(k0) > quantile(tmp[13,], .95)) ## 0; fail to reject H0; shape independence
1 * (max(k02) > quantile(tmp[14,], .95)) ## 0

mean(max(k0) > tmp[13,])
## [1] 0.025
mean(max(k02) > tmp[14,])
## [1] 1


###############################################################################################
## Analysis with CMV status; there are missing values
## binary age
###############################################################################################
head(dat.HSCT)

fname <- reSurv(Time, id, event, status) ~ allo + gender + CMVstatus 
fit <- gsm(fname, data = dat.HSCT)
str(fit)


dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("allo", "gender", "CMVstatus"), function(x)
    which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(7:9)] <- c("x1", "x2", "x3")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 7:9, 12)]
head(dat.HSCT2)

b0 <- seq(0, 2 * pi, length = 50)
bi <- as.matrix(expand.grid(b0, b0))
k0 <- sapply(1:NROW(bi), function(x) getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))))
k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0))
max(k0)
max(k02)


B <- 1000
mm <- aggregate(event ~ id, dat.HSCT, sum)[, 2]
n <- length(unique(dat.HSCT$id))
getBootk <- function(dat.HSCT) {
    ind <- sample(1:n, replace = TRUE)
    dat.HSCT0 <- dat.HSCT[unlist(sapply(ind, function(x) which(dat.HSCT$id %in% x))),]
    dat.HSCT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.HSCT0) <- NULL
    fit0 <- gsm(fname, data = dat.HSCT0, shp.ind = FALSE)
    names(dat.HSCT0)[c(7:9)] <- c("x1", "x2", "x3")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 7:9, 12)]
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit0$Fhat0) - k02[x])))
}

set.seed(1) 
system.time(print(getBootk(dat.HSCT)))

cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "getBootk"))
invisible(clusterExport(NULL, c("n", "mm", "dat.HSCT", "bi", "k0", "k02", "fname")))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

set.seed(1)
system.time(tmp <- parSapply(NULL, 1:100, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)



1 * (max(k0) > quantile(tmp[13,], .95)) ## 1
1 * (max(k02) > quantile(tmp[14,], .95)) ## 1

mean(max(k0) > tmp[13,])
## [1] 0.95
mean(max(k02) > tmp[14,])
## [1] 1


###############################################################################################
## age01, allo, gender, CMV
###############################################################################################

fname <- reSurv(Time, id, event, status) ~ allo + gender + CMVstatus + age01
dat.HSCT <- dat.HSCT[complete.cases(dat.HSCT),]
dat.HSCT$id <- rep(1:length(unique(dat.HSCT$id)), aggregate(event ~ id, dat.HSCT, sum)[, 2] + 1)

fit <- gsm(fname, data = dat.HSCT)
str(fit)

## Custom function for this data set, need to generalize this later...

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("age01", "allo", "gender", "CMVstatus"), function(x) which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(7:9, 11)] <- c("x1", "x2", "x3", "x4")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 7:9, 11, 12)]
head(dat.HSCT2)

bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 50), seq(0, 2 * pi, length = 50), seq(0, 2 * pi, length = 50)))
k0 <- sapply(1:NROW(bi), function(x) getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))))
k02 <- sapply(1:NROW(bi), function(x) getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0))
max(k0)
max(k02)

B <- 1000
mm <- aggregate(event ~ id, dat.HSCT, sum)[, 2]
n <- length(unique(dat.HSCT$id))
getBootk <- function(dat.HSCT) {
    ind <- sample(1:n, replace = TRUE)
    dat.HSCT0 <- dat.HSCT[unlist(sapply(ind, function(x) which(dat.HSCT$id %in% x))),]
    dat.HSCT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.HSCT0) <- NULL
    fit0 <- gsm(fname, data = dat.HSCT0, shp.ind = FALSE)
    names(dat.HSCT0)[c(7:9, 11)] <- c("x1", "x2", "x3", "x4")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 7:9, 11, 12)]
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.HSCT0, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit0$Fhat0) - k02[x])))
}

set.seed(1)    
system.time(print(getBootk(dat.HSCT)))

cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "getBootk"))
invisible(clusterExport(NULL, c("n", "mm", "dat.HSCT", "bi", "k0", "k02", "fname")))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

set.seed(1)
system.time(tmp <- parSapply(NULL, 1:300, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)


1 * (max(k0) > quantile(tmp[13,], .95)) 
1 * (max(k02) > quantile(tmp[14,], .95)) 

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
head(dat.HSCT)

###############################################################################################
## Analysis with scaled age
###############################################################################################

fname <- reSurv(Time, id, event, status) ~ race + allo + scaleAge
fit <- gsm(fname, data = dat.HSCT)
str(fit)

## Custom function for this data set, need to generalize this later...

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("race", "allo", "scaleAge"), function(x) which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(6:7, 10)] <- c("x1", "x2", "x3")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 6:7, 10, 12)]
head(dat.HSCT2)

bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 100), seq(0, 2 * pi, length = 100)))
system.time(k0 <- sapply(1:NROW(bi), function(x)
    getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)))))
system.time(k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0)))


B <- 1000
mm <- aggregate(event ~ id, dat.HSCT, sum)[, 2]
n <- length(unique(dat.HSCT$id))
getBootk <- function(dat.HSCT) {
    ind <- sample(1:n, replace = TRUE)
    dat.HSCT0 <- dat.HSCT[unlist(sapply(ind, function(x) which(dat.HSCT$id %in% x))),]
    dat.HSCT0$id <- rep(1:n, mm[ind] + 1)
    rownames(dat.HSCT0) <- NULL
    fit0 <- gsm(fname, data = dat.HSCT0, shp.ind = FALSE)
    names(dat.HSCT0)[c(6:7, 10)] <- c("x1", "x2", "x3")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 6:7, 10, 12)]
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
system.time(tmp <- parSapply(NULL, 1:100, function(z) getBootk(dat.HSCT))) 
stopCluster(cl)


1 * (max(k0) > quantile(tmp[13,], .95)) ## 0; don't reject H0
1 * (max(k02) > quantile(tmp[14,], .95)) ## 0
mean(max(k0) > tmp[13,]) ## .74
mean(max(k02) > tmp[14,]) ## .93

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

bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 50),
                            seq(0, 2 * pi, length = 50),
                            seq(0, 2 * pi, length = 50)))
system.time(k0 <- sapply(1:NROW(bi), function(x)
    getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)))))
system.time(k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0)))
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
system.time(tmp <- parSapply(NULL, 1:200, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)


1 * (max(k0) > quantile(tmp[nrow(tmp) - 1,], .95)) ## 0
1 * (max(k02) > quantile(tmp[nrow(tmp),], .95)) ## 1 


###############################################################################################

fname <- reSurv(Time, id, event, status) ~ allo + gender + scaleAge
fit <- gsm(fname, data = dat.HSCT)
str(fit)

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("scaleAge", "allo", "gender"), function(x) which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(7, 8, 10)] <- c("x1", "x2", "x3")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 7, 8, 10, 12)]
head(dat.HSCT2)

bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 100), seq(0, 2 * pi, length = 100)))
system.time(k0 <- sapply(1:NROW(bi), function(x)
    getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)))))
system.time(k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0)))
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
    names(dat.HSCT0)[c(7, 8, 10)] <- c("x1", "x2", "x3")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 7, 8, 10, 12)]
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


1 * (max(k0) > quantile(tmp[13,], .95)) ## 0; reject H0; shape dependence
1 * (max(k02) > quantile(tmp[14,], .95)) ## 0

mean(max(k0) > tmp[13,])  ## ~ .9
mean(max(k02) > tmp[14,])

summary(k0)
summary(k02)

###############################################################################################

fname <- reSurv(Time, id, event, status) ~ allo + scaleAge
fit <- gsm(fname, data = dat.HSCT)
str(fit)

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("scaleAge", "allo"), function(x) which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(7, 10)] <- c("x1", "x2")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 7, 10, 12)]
head(dat.HSCT2)

## bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 100), seq(0, 2 * pi, length = 100)))
bi <- seq(0, 2 * pi, length = 100)

system.time(k0 <- sapply(1:NROW(bi), function(x)
    getk0(dat.HSCT2, cumprod(c(1, sin(bi[x])) * c(cos(bi[x]), 1)))))
system.time(k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x])) * c(cos(bi[x]), 1)), fit$Fhat0)))
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
    names(dat.HSCT0)[c(7, 10)] <- c("x1", "x2")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 7, 10, 12)]
    c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
      max(sapply(1:NROW(bi), function(x)
          getk0(dat.HSCT0, cumprod(c(1, sin(bi[x])) * c(cos(bi[x]), 1))) - k0[x])),
      max(sapply(1:NROW(bi), function(x)
          getk02(dat.HSCT0, cumprod(c(1, sin(bi[x])) * c(cos(bi[x]), 1)), fit0$Fhat0) - k02[x])))
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
system.time(tmp <- parSapply(NULL, 1:500, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)


1 * (max(k0) > quantile(tmp[9,], .95)) ## 0; reject H0; shape dependence
1 * (max(k02) > quantile(tmp[10,], .95)) ## 0

mean(max(k0) > tmp[9,]) 
mean(max(k02) > tmp[10,])

summary(k0)
summary(k02)


######################################################################################################
head(dat.HSCT)

fname <- reSurv(Time, id, event, status) ~ race + allo + gender + scaleAge

fit <- gsm(fname, data = dat.HSCT)
str(fit)

dat.HSCT2 <- dat.HSCT
as.numeric(sapply(c("scaleAge", "allo", "gender", "race"), function(x) which(names(dat.HSCT2) == x)))
names(dat.HSCT2)[c(6:8, 10)] <- c("x1", "x2", "x3", "x4")
dat.HSCT2 <- dat.HSCT2[,c(1:4, 6:8, 10, 12)]
head(dat.HSCT2)

bi <- as.matrix(expand.grid(seq(0, 2 * pi, length = 50),
                            seq(0, 2 * pi, length = 50),
                            seq(0, 2 * pi, length = 50)))

system.time(k0 <- sapply(1:NROW(bi), function(x)
    getk0(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)))))
system.time(k02 <- sapply(1:NROW(bi), function(x)
    getk02(dat.HSCT2, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0)))
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
    names(dat.HSCT0)[c(6:8, 10)] <- c("x1", "x2", "x3", "x4")
    dat.HSCT0 <- dat.HSCT0[,c(1:4, 6:8, 10, 12)]
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


1 * (max(k0) > quantile(tmp[nrow(tmp) - 1,], .95)) ## 0
1 * (max(k02) > quantile(tmp[nrow(tmp),], .95)) ## 1 


############################################################################################
## Load package and data
## HSCT cohort from csv
############################################################################################

library(tidyverse)
library(parallel)
library(GSM)
library(reReg)
library(survival)

dat.hsct <- read.csv("HSCT.csv", sep=",", header=TRUE)
dat.hsct <- subset(dat.hsct, !is.na(hsct_type) & Time > 0)
dat.hsct$id <- rep(1:length(unique(dat.hsct$id)),
                   with(dat.hsct, unlist(lapply(split(id, id), length))))

## hsct, prepare variables:
dat.hsct$allo <- dat.hsct$hsct_type
dat.hsct$lym <- dat.hsct$disease_lym___0 + dat.hsct$disease_lym___1
dat.hsct$race0 <- ifelse(dat.hsct$race == 1, 1, 0)
dat.hsct$race1 <- ifelse(dat.hsct$race == 2, 1, 0)
base.hsct <- subset(dat.hsct, select = -c(Time, ser.inf.type, event, status))
base.hsct <- base.hsct[cumsum(aggregate(gender ~ id, data = base.hsct, length)[,2]),]
base.hsct$scaleAge <- scale(base.hsct$age)
## base.hsct$scaleAge <- base.hsct$age / max(base.hsct$age)
dat.hsct$scaleAge <- base.hsct$scaleAge[dat.hsct$id]
dat.hsct$scaleAge2 <- pmax(0, base.hsct$scaleAge[dat.hsct$id])


dat.hsct$heme1 <- ifelse(dat.hsct$heme_state == 1, 1, 0)
dat.hsct$heme2 <- ifelse(dat.hsct$heme_state == 2, 1, 0)
dat.hsct$cmv1 <- ifelse(dat.hsct$sero_cmv_hsct == 0, 1, 0)
dat.hsct$cmv2 <- ifelse(dat.hsct$sero_cmv_hsct == 0 & dat.hsct$dsero_cmv_hsct == 0, 1, 0)

head(dat.hsct)
dim(dat.hsct)
summary(dat.hsct)

## dat.hsct is the original data
## dat0 is what we will use
dat.hsct$m <- rep(aggregate(event ~ id, dat.hsct, sum)[, 2], aggregate(Time ~ id, dat.hsct, length)[, 2])
dat0 <- dat.hsct %>% select(id, Time, event, status, m, 
                            age, scaleAge, scaleAge2, race0, allo, gender, lym, agvhd,
                            heme1, heme2, cmv1, cmv2)
summary(dat0)
dim(dat0)
head(dat0)

dat00 <- subset(dat0, event == 0)
summary(dat00)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

fname <- reSurv(Time, id, event, status) ~ scaleAge + race0 + allo
head(dat0)

pVal <- function(fname, B = 100, dat0 = dat0) {
    xNames <- attr(terms(fname), "term.labels")
    p <- length(attr(terms(fname), "term.labels"))
    fit <- gsm(fname, data = dat0)
    ## str(fit)
    dat1 <- dat0
    xCol <- as.numeric(sapply(xNames, function(x) which(names(dat0) == x)))
    colnames(dat1)[xCol] <- paste("x", 1:p, sep = "")
    dat1 <- dat1[,c(1:5, xCol)]
    ## head(dat1)
    ## bi <- as.matrix(expand.grid(rep(list(seq(0, 2 * pi, length = 200)), p - 1)))
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .1)), p)))
    r <- apply(tmp, 1, function(z) sqrt(sum(z^2)))
    bi <- (tmp / r)[r < 1 & r > 0,]
    k0 <- sapply(1:NROW(bi), function(x) getk0(dat1, bi[x,]))
    k02 <- sapply(1:NROW(bi), function(x) getk02(dat1, bi[x,], fit$Fhat0))
    getBootK <- function(dat) {
        n <- length(unique(dat$id))    
        mm <- aggregate(event ~ id, dat, length)[, 2]
        ind <- sample(1:n, n, TRUE)
        datB <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        datB$id <- rep(1:n, mm[ind])
        rownames(datB) <- NULL
        fitB <- gsm(fname, dat = datB)
        datB1 <- datB
        xCol <- as.numeric(sapply(xNames, function(x) which(names(datB) == x)))
        colnames(datB1)[xCol] <- paste("x", 1:p, sep = "")
        datB1 <- datB1[,c(1:5, xCol)]
        kb <- max(sapply(1:NROW(bi), function(x) getk0(datB1, bi[x,]) - k0[x]))
        kb2 <- max(sapply(1:NROW(bi), function(x) getk02(datB1, bi[x,], fitB$Fhat0) - k02[x]))
        c(max(kb), max(kb2), fitB$b0, fitB$b00, fitB$r0, fitB$r00)
    }
    cl <- makePSOCKcluster(8)
    ## cl <- makePSOCKcluster(16)
    setDefaultCluster(cl)
    invisible(clusterExport(cl, c("bi", "k0", "k02", "fname", "dat0", "xNames", "p", "getBootK"),
                            environment()))
    invisible(clusterEvalQ(NULL, library(GSM)))
    invisible(clusterEvalQ(NULL, library(reReg)))
    system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootK(dat0))) 
    stopCluster(cl)
    list(h1 = mean(max(k0) > tmp[1,]), h2 = mean(max(k02) > tmp[2,]),
         coef11 = fit$b0, coef12 = fit$b00,
         coef21 = fit$r0, coef22 = fit$r00,
         se11 = apply(tmp[3:5,], 1, sd),
         se12 = apply(tmp[6:8,], 1, sd),
         se21 = apply(tmp[9:11,], 1, sd),
         se22 = apply(tmp[12:14,], 1, sd))      
}


pValShape <- function(fname, B = 100, dat0 = dat0) {
    xNames <- attr(terms(fname), "term.labels")
    p <- length(attr(terms(fname), "term.labels"))
    dat1 <- dat0
    xCol <- as.numeric(sapply(xNames, function(x) which(names(dat0) == x)))
    colnames(dat1)[xCol] <- paste("x", 1:p, sep = "")
    dat1 <- dat1[,c(1:5, xCol)]
    dat1 <- dat1[complete.cases(dat1),]
    ## head(dat1)
    ## bi <- as.matrix(expand.grid(rep(list(seq(0, 2 * pi, length = 100)), p - 1)))
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .05)), p)))
    r <- apply(tmp, 1, function(z) sqrt(sum(z^2)))
    bi <- (tmp / r)[r < 1 & r > 0,]
    k0 <- sapply(1:NROW(bi), function(x) getk0(dat1, bi[x,]))
    getBootK <- function(dat) {
        n <- length(unique(dat$id))
        mm <- aggregate(event ~ id, dat, length)[, 2]
        ind <- sample(1:n, n, TRUE)
        datB <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        datB$id <- rep(1:n, mm[ind])
        datB <- datB[complete.cases(datB),]
        rownames(datB) <- NULL
        datB1 <- datB
        xCol <- as.numeric(sapply(xNames, function(x) which(names(datB) == x)))
        colnames(datB1)[xCol] <- paste("x", 1:p, sep = "")
        datB1 <- datB1[,c(1:5, xCol)]
        kb <- max(sapply(1:NROW(bi), function(x) getk0(datB1, bi[x,]) - k0[x]))
            ## getk0(datB1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1)) - k0[x]))
        max(kb)
    }
    ## cl <- makePSOCKcluster(8)
    cl <- makePSOCKcluster(16)
    setDefaultCluster(cl)
    invisible(clusterExport(cl, c("bi", "k0", "fname", "dat0", "xNames", "p", "getBootK"),
                            environment()))
    invisible(clusterEvalQ(NULL, library(GSM)))
    invisible(clusterEvalQ(NULL, library(reReg)))
    system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootK(dat0))) 
    stopCluster(cl)
    mean(max(k0) > tmp)
}

system.time(f1 <- pValShape(reSurv(Time, id, event, status) ~ allo + heme1 + scaleAge + scaleAge2, 100, dat0)) ## .03
system.time(f1 <- pValShape(reSurv(Time, id, event, status) ~ allo + heme2 + scaleAge + scaleAge2, 100, dat0)) ## .13
system.time(f1 <- pValShape(reSurv(Time, id, event, status) ~ allo + cmv1 + scaleAge + scaleAge2, 100, dat0)) ## .19
system.time(f1 <- pValShape(reSurv(Time, id, event, status) ~ allo + cmv2 + scaleAge + scaleAge2, 100, dat0)) ## .19

system.time(f2 <- pValShape(reSurv(Time, id, event, status) ~ allo + heme1 + cmv2 + scaleAge + scaleAge2, 100, dat0)) ## .11

system.time(f2 <- pValShape(reSurv(Time, id, event, status) ~ allo + heme2 + cmv2 + scaleAge + scaleAge2, 100, dat0)) ## .11



plotEvents(reSurv(Time, id, event, status) ~ allo, data = dat0)
plotEvents(reSurv(Time, id, event, status) ~ cmv1, data = dat0)

## Black and white event plots
n <- length(unique(dat0$id))

base <- subset(dat0, event == 0)
rownames(base) <- NULL
datPlot <- NULL
for (i in 1:n) {                                    
    tmp <- split(dat0$Time, dat0$id)[[i]]
    mi <- length(tmp) - 1
    datPlot <- rbind(datPlot, cbind(base$id[i], tmp, c(rep(0, mi), base$status[i]),
                                    c(rep(1, mi), 0), base$cmv1[i], base$allo[i]))
    ## datPlot <- rbind(datPlot, cbind(base$id[i], c(tmp, base$Time[i]), c(rep(0, mi),
    ## base$status[i]),
    ##                                 c(rep(1, mi), 0), base$cmv1[i], base$allo[i]))
}
colnames(datPlot) <- c("id", "T", "status", "event", "cmv", "allo")
datPlot <- data.frame(datPlot)
datPlot <- datPlot[datPlot$T > 0,]
datPlot <- datPlot[datPlot$event + datPlot$status > 0,]
## datPlot$type <- ifelse(datPlot$allo == 1, "allogeneic", "autologous")
datPlot$type <- ifelse(datPlot$cmv == 1, "positive", "negative")
datPlot$cmv <- ifelse(datPlot$cmv == 1, "positive", "negative")
datPlot$type <- as.factor(datPlot$type)
datPlot$cmv <- as.factor(datPlot$cmv)
datPlot$event <- as.factor(datPlot$event)
datPlot$status <- as.factor(datPlot$status)
datPlot$id0 <- 0

base <- rbind(base, base)
base$Time[1:n] <- 0
base <- base[order(base$id),]
rownames(base) <- NULL
## base$type <- ifelse(base$allo == 1, "allogeneic", "autologous")
base$type <- ifelse(base$cmv1 == 1, "positive", "negative")
base$cmv <- ifelse(base$cmv1 == 1, "positive", "negative")
base$type <- as.factor(base$type)
base$cmv <- as.factor(base$cmv)
base <- base[order(base$type, base$id), ]

datPlot$type <- factor(datPlot$type, levels = c("positive", "negative")
base$type <- factor(base$type, levels = c("positive", "negative")
## datPlot$type <- factor(datPlot$type, levels = c("allogeneic", "autologous"))
## base$type <- factor(base$type, levels = c("allogeneic", "autologous"))

## re-sort for allo vs auto
base <- base[base$Time > 0,]
base <- base[order(base$type, base$Time),] ##, decreasing = TRUE),]
base$id0 <- 1:n
base <- rbind(base, base)
base$Time[1:n] <- 0
base <- base[order(base$id0),]
rownames(base) <- NULL
## head(base, 30)
base$id0[which(base$allo == 1)] <- base$id0[which(base$allo == 1)] - 36
datPlot$id0 <- 0
for (i in 1:nrow(datPlot)) datPlot$id0[i] <- base$id0[min(which(base$id == datPlot$id[i]))]

ggplot(base, aes(x = Time, y = id0, group = id0)) +
geom_line(color = "gray55", size = 1.5) +
facet_grid(type ~ ., scales = "free_y", space = "free_y") +
geom_point(data = datPlot,
           aes(x = T, y = id0, shape = event:status), size = 1.5, alpha = .7) +
scale_shape_manual(values = c(16, 4), name = "Event types:", labels = c("Death", "Infection")) +
theme_bw() + labs(x = "Time in days", y = "Subjects") +
theme(legend.position="none", axis.text.y = element_blank())

ggsave("HSCT-allo.pdf")

head(base)
## re-sort for cmv
base <- base[base$Time > 0,]
base <- base[order(base$type, base$Time),] ##, decreasing = TRUE),]
base$id0 <- 1:n
base <- rbind(base, base)
base$Time[1:n] <- 0
base <- base[order(base$id0),]
rownames(base) <- NULL
head(base, 30)
base$id0[which(base$cmv == 1)] <- base$id0[which(base$cmv == 1)] - 36
datPlot$id0 <- 0
for (i in 1:nrow(datPlot)) datPlot$id0[i] <- base$id0[min(which(base$id == datPlot$id[i]))]


ggplot(base, aes(x = Time, y = id0, group = id0)) +
geom_line(color = "gray55", size = 1.5) +
facet_grid(type ~ ., scales = "free_y", space = "free_y") +
geom_point(data = datPlot,
           aes(x = T, y = id0, shape = event:status), size = 1.5, alpha = .7) +
scale_shape_manual(values = c(16, 4), name = "Event types:", labels = c("Death", "Infection")) +
theme_bw() + labs(x = "Time in days", y = "Subjects") +
theme(legend.position="none", axis.text.y = element_blank())

ggsave("HSCT-cmv.pdf")

#########################################################################################
## event plot by age



n <- length(unique(dat0$id))
base <- subset(dat0, event == 0)
base <- base[order(base$age),]
base$id <- 1:n
rownames(base) <- datPlot <- NULL
datPlot <- subset(dat0,  select = c("id", "Time", "event", "status", "m", "age", "scaleAge"))
datPlot <- datPlot[order(datPlot$age),]
datPlot$id <- rep(1:n, c(1, diff(which(datPlot$event == 0))))
base <- rbind(base, base)
base$Time[1:n] <- 0
base <- base[order(base$id),]
datPlot <- subset(datPlot, event + status > 0)
datPlot$event <- as.factor(datPlot$event)
datPlot$status <- as.factor(datPlot$status)

                  
ggplot(base, aes(x = Time, y = id, group = id)) + geom_line(color = "gray55", size = 1.5) +
geom_point(data = datPlot,
           aes(x = Time, y = id, shape = event:status), size = 1.5, alpha = .7) +
scale_shape_manual(values = c(16, 4), name = "Event types:", labels = c("Death", "Infection")) +
theme_bw() + labs(x = "Time in days", y = "Subjects") +
theme(legend.position="none", axis.text.y = element_blank()) +
geom_hline(yintercept = 52, color = "red")

ggsave("HSCT-age.pdf")

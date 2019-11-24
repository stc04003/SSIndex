############################################################################################
## Load package and data
## SOT data with kidney transplant
############################################################################################

library(parallel)
library(SSIndex)
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
## List of 8
##  $ b0   : num [1:3] 0.475 0.8781 0.0581
##  $ r0   : num [1:3] 0.47 0.874 0.124
##  $ b00  : num [1:3] 0.475 0.8781 0.0581
##  $ r00  : num [1:3] 0.3721 0.9277 0.0302
##  $ d    : NULL
##  $ dstar: NULL
##  $ Fhat : num [1:164] 0.826 1 0.585 0.945 0.227 ...
##  $ Fhat0: num [1:164] 0.829 0.609 0.589 0.943 0.231 ...

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

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
invisible(clusterEvalQ(NULL, library(SSIndex)))

set.seed(1)
system.time(tmp <- parSapply(NULL, 1:200, function(z) getBootk(dat.HSCT))) ## 
stopCluster(cl)


1 * (max(k0) > quantile(tmp[nrow(tmp) - 1,], .95)) ## 0
1 * (max(k02) > quantile(tmp[nrow(tmp),], .95)) ## 1 


############################################################################################
## Load package and data
## HSCT cohort from csv
## Used in the paper
############################################################################################

library(tidyverse)
library(parallel)
library(SSIndex)
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

fname <- reSurv(Time, id, event, status) ~ allo + scaleAge + scaleAge2 + gender + race0 + cmv1
## fname <- reSurv(Time, id, event, status) ~ scaleAge + scaleAge2 + allo + cmv1 + gender + race0
fit <- gsm(fname, data = dat0)
str(fit)
tab <- data.frame(b0 = fit$b0, r0 = fit$r0, b00 = fit$b00, r00 = fit$r0)
names(tab) <- c("b", "r", "b.smooth", "r.smooth")
rownames(tab) <- c("allogeneic", "age", "age2", "gender", "white", "cmv")
tab
##                     b           r   b.smooth    r.smooth
## allogeneic  0.7563023  0.82574611  0.7707468  0.82574611
## age        -0.2160166 -0.02609171 -0.2094977 -0.02609171
## age2        0.5446405  0.12314392  0.5329246  0.12314392
## gender      0.1124970  0.12469424  0.1008363  0.12469424
## whtie       0.2402366 -0.43123901  0.2457787 -0.43123901
## cmv1        0.1197548  0.31746247  0.0864659  0.31746247

fit$Fhat0
datFhat <- data.frame(Time = dat00$Time, Fhat = fit$Fhat0)
datFhat <- datFhat[order(datFhat$Time),]

library(ggplot2)

## Ploting F(t, a) with a = \bar{X} %*% beta
ggplot(datFhat, aes(x = Time, y = Fhat)) + geom_step() +
    labs(x = "Time", y = expression(F(t, paste(hat(beta)^T, " ", bar(X)))), title = "")

## ggsave("FtaAT0.pdf")
## ggsave("Fta.pdf")

## Ploting F(t, a) with fixed t at each Yi and varying a

i <- 111
ggplot(data.frame(Time = fit$xb[order(fit$xb)], Fhat = fit$Fhat0[[i]][order(fit$xb)]),
       aes(x = Time, y = Fhat)) + geom_step()
summary(fit$xb)

fhat <- do.call(rbind, fit$Fhat0)

for (i in 1:164) {
    plot(fit$xb, fhat[i,], cex = .5, pch = 19)
    ## plot(dat00$Time, fhat[,i], cex = .5, pch = 19)
    Sys.sleep(.5)
}

e

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
    invisible(clusterEvalQ(NULL, library(SSIndex)))
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
    invisible(clusterEvalQ(NULL, library(SSIndex)))
    system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootK(dat0))) 
    stopCluster(cl)
    mean(max(k0) > tmp)
}

system.time(f1 <- pVal(fname, 100, dat0))
f1

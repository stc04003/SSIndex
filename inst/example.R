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

makeTab <- function(dat) {
    cbind(colMeans(dat)[1:4], apply(dat, 2, sd)[1:4])
}

tab <- cbind(makeTab(sim1), makeTab(sim1.2), makeTab(sim2), makeTab(sim3), makeTab(sim4))

print(xtable(tab, digits = 3), math.style.negative = TRUE)


## Testing \beta_0 = 0
library(tidyverse)

dat <- simDat(100, "M1")
dat <- simDat(100, "M2")

gsm(dat, FALSE)
gsm(dat, TRUE)
gsm(dat, "test")

system.time(f1 <- gsm(simDat(100, "M1"), "test"))
system.time(f2 <- gsm(simDat(100, "M1"), "test"))
system.time(f6 <- gsm(simDat(100, "M1"), "test"))
system.time(f5 <- gsm(simDat(100, "M1"), "test"))

system.time(f3 <- gsm(simDat(100, "M2"), "test"))
system.time(f4 <- gsm(simDat(100, "M2"), "test"))

e


sqrt(50) * with(f1, d / sd(dstar))
sqrt(50) * with(f2, d / sd(dstar))
sqrt(50) * with(f3, d / sd(dstar))
sqrt(50) * with(f4, d / sd(dstar))
sqrt(50) * with(f5, d / sd(dstar))
sqrt(50) * with(f6, d / sd(dstar))

library(GSM)

set.seed(2)
system.time(f1 <- gsm(simDat(100, "M1"), "test"))

f1$d
sd(f1$dstar)
sqrt(50) * with(f1, d / sd(dstar))

set.seed(2)
system.time(f1 <- gsm(simDat(100, "M1"), "test"))

debug(gsm)
debug(getd)


Ft0 <- exp(-unlist(mapply(FUN = function(x, y)
    .C("shapeFun", 
       as.integer(n2), as.integer(mm2), as.integer(midx2), as.double(tij2), 
       as.double(yi2), as.double(X2 %*% tilde.b), as.double(x), 
       as.double(y), result = double(1), PACKAGE = "GSM")$result, 
    rep(.6, length(u)), u)))

Ft <- exp(-.C("shapeFun2", as.integer(n2), as.integer(mm2), as.integer(midx2), 
              as.double(tij2), as.double(yi2), as.double(X2 %*% tilde.b), 
              as.double(.6), as.double(u), as.integer(length(u)),
              result = double(length(u)),
              PACKAGE = "GSM")$result)


##

do <- function(n, model) {
    dat <- simDat(n, model)
    tmp <- gsm(dat, "test")
    c(tmp$b0, tmp$r0, sqrt(round(n / 2)) * tmp$d / sd(tmp$dstar))
}

do(100, "M2")

library(parallel)
library(xtable)

sim1 <- sim1.2 <- sim2 <- sim3 <- sim4 <- NULL
cl <- makePSOCKcluster(8)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do')))
invisible(clusterEvalQ(NULL, library(GSM)))
sim2 <- t(matrix(unlist(parLapply(NULL, 1:20, function(z) do(100, "M2"))), 5))
sim1 <- t(matrix(unlist(parLapply(NULL, 1:200, function(z) do(100, "M1"))), 5))

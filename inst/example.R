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

set.seed(1);round(do(100, "M1"), 3) # 0.000       0.000       0.341       0.940       0.946 2655087.000
set.seed(1);round(do(100, "M2"), 3) # -0.670      -0.742       0.614       0.789      -2.980 2655087.000
set.seed(1);round(do(100, "M3"), 3) # 0.000       0.000      -0.593      -0.805      -0.740 2655087.000
set.seed(1);round(do(100, "M4"), 3) # 0.165      -0.986       0.335       0.942      -2.485 2655087.000

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
## Debug area
## ------------------------------------------------------------------------------------

dat <- simDat(100, "M4")
gsm(dat)

gsm(dat, "test")
debug(getd)

system.time(print(do(1000, "M1")))
system.time(print(do(1000, "M4")))

plot(xb, tmp1, "p", cex = .5, pch = 16, ylim = range(c(tmp1, tmp2, tmp3)))
points(xb, tmp2, col = 2, cex = .5, pch = 16)
points(xb, tmp3, col = 3, cex = .5, pch = 16)

summary(tmp1)
summary(tmp2)

summary(xb)

h <- 0.1682727

h <- .16
z <- 1 ## which xb to evaluate

Ri <- .C("shapeFun2", as.integer(n2), as.integer(mm2), as.integer(midx2), 
         as.double(tij2), as.double(yi2), as.double(xb), 
         as.double(z), as.double(h), result = double(sum(mm2)), PACKAGE = "GSM")$result
Ft <- exp(-colSums((Ri * outer(tij2, u, ">="))))

plot(u, Ft, 's')
lines(u, pbeta(u, 2, 1 + exp(z)), col = 2, lwd = 2)
lines(u, bcdf(u, exp(z)), col = 3, lwd = 1)

bcdf <- function(t, x) {
    ((x + 1) * (1 - t)^(x + 2) - (x + 2) * (1 - t)^(x + 1) + 1)
}


x <- exp(rnorm(1))
t0 <- sort(runif(1e4))
summary(bcdf(t0, x) - pbeta(t0, 2, 1 + x))

2 / (3 + x)
sum((1 - bcdf(t0, x)) * diff(c(0, t0)))
sum(t0 * diff(c(0, bcdf(t0, x))))

## new test
library(GSM)
set.seed(123)

do <- function(n, model) {
    dat <- simDat(n, model) 
    f1 <- gsm(dat)
    k0 <- getk(dat, f1$b0)
    kb <- replicate(1000, boot.k(dat, f1$b0))
    k0 / sd(kb)
}

library(parallel)

f1 <- f2 <- f3 <- f4 <- NULL

cl <- makePSOCKcluster(8)
## cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do')))
invisible(clusterEvalQ(NULL, library(GSM)))

f1 <- c(parSapply(NULL, 1:1000, function(z) do(200, "M1")))
f2 <- c(parSapply(NULL, 1:1000, function(z) do(100, "M1")))

f5 <- c(parSapply(NULL, 1:1000, function(z) do(200, "M4")))
f3 <- c(parSapply(NULL, 1:1000, function(z) do(200, "M2")))
f4 <- c(parSapply(NULL, 1:1000, function(z) do(200, "M3")))
f2 <- c(parSapply(NULL, 1:1000, function(z) do(400, "M1")))

stopCluster(cl)

mean(f1 > -qnorm(.95))
mean(f2 > -qnorm(.95))
mean(f3 > -qnorm(.95))
mean(f4 > -qnorm(.95))
mean(f5 > -qnorm(.95))

mean(f1 > qnorm(.95))
mean(f2 > qnorm(.95))
mean(f3 > qnorm(.95))
mean(f4 > qnorm(.95))
mean(f5 > qnorm(.95))

mean(f1 * sqrt(200) > qnorm(.95))
mean(f2 * sqrt(400) > qnorm(.95))
mean(f3 * sqrt(200) > qnorm(.95))
mean(f4 * sqrt(200) > qnorm(.95))
mean(f5 * sqrt(200) > qnorm(.95))


dat <- simDat(1000, "M1")
b0 <- gsm(dat)$b0
getk(dat, b0)
bb <- replicate(1000, boot.k(dat, b0))

summary(bb)
mean(bb < 0)

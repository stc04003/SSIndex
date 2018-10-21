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

sim1 <- sim1.2 <- sim2 <- sim3 <- NULL
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

tilde.mu <- function(x) {
    Ft <- exp(-unlist(mapply(FUN = function(x, y)
        .C("shapeFun", 
           as.integer(n2), as.integer(mm2), as.integer(midx2), as.double(tij2), 
           as.double(yi2), as.double(X2 %*% tilde.b), as.double(x), as.double(y), 
           result = double(1), PACKAGE = "GSM")$result, rep(x, length(u)), u)))
    sum(diff(c(0, Ft)) * u)
}

b.test <- function(dat) {
    n <- length(unique(dat$id))
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    tij <- subset(dat, event == 1)$t
    yi <- subset(dat, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = c(x1, x2)))
    p <- ncol(X)
    ind <- sample(1:n)[1:round(n/2)]
    dat1 <- subset(dat, id %in% ind)
    dat2 <- subset(dat, !(id %in% ind))
    tilde.b <- gsm(dat1)$b0
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    tij2 <- subset(dat2, event == 1)$t
    yi2 <- subset(dat2, event == 0)$t
    midx2 <- c(0, cumsum(mm2)[-length(mm2)])
    X2 <- as.matrix(subset(dat2, event == 0, select = c(x1, x2)))
    u <- unique(sort(c(tij, yi)))
    xb <- X %*% tilde.b
    d1 <- sapply(xb[xb <= median(xb)], function(z) tilde.mu(z))
    d2 <- sapply(xb[xb > median(xb)], function(z) tilde.mu(z))
    return((sum(d1) - sum(d2)) / n)
}

system.time(d <- replicate(200, b.test(dat)))

library(SSIndex)

####################################################################################3
## Parallel computing on 1000 replications
## March 10
####################################################################################3
library(parallel)
library(xtable)

## March 17, 2019
do <- function(n, model, frailty = FALSE) {
    B <- 5
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    bi <- seq(0, 2 * pi, length = 100)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    k02 <- sapply(bi, function(x) getk02(dat, c(cos(x), sin(x)), fit$Fhat0))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    n <- length(unique(dat$id))
    getBootk <- function(dat) {
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        dat0$Time[duplicated(dat0$Time) & dat0$Time < max(dat0$Time)] <-
            dat0$Time[duplicated(dat0$Time) & dat0$Time < max(dat0$Time)] +
            rnorm(sum(duplicated(dat0$Time) & dat0$Time < max(dat0$Time)), sd = .001)
        fit0 <- gsm(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
                    data = dat0, shp.ind = FALSE)
        c(fit0$b0, fit0$b00, fit0$r0, fit0$r00,
          max(sapply(1:length(bi), function(x)
              getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x])),
          max(sapply(1:length(bi), function(x)
              getk02(dat0, c(cos(bi[x]), sin(bi[x])), fit0$Fhat0) - k02[x])))
    }
    tmp <- replicate(B, getBootk(dat))
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta,
    ## (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0, fit$b00, fit$r0, fit$r00, 
      sqrt(diag(var(t(tmp[1:2,])))),
      sqrt(diag(var(t(tmp[3:4,])))),
      sqrt(diag(var(t(tmp[5:6,])))),
      sqrt(diag(var(t(tmp[7:8,])))), 
      1 * (max(k0) > quantile(tmp[9,], .95)),
      1 * (max(k02) > quantile(tmp[10,], .95)))
}

## No bootstrap
do2 <- function(n, model, B = 200, frailty = FALSE) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    fit.indep <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
                     data = dat, shp.ind = TRUE)
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta, (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0, fit$r0, fit.indep$r0)
}

set.seed(1)
system.time(print(do(50, "M1", TRUE)))
## 0.2703377 -0.9627656  0.2703377 -0.9627656 -0.9068044  0.4215517 -0.9068044 
## 0.4215517  0.6043511  0.7671864  0.6481810  0.7593999  0.4532586  0.3688905 
##                             95%        95% 
## 0.4492857  0.3486952  0.0000000  0.0000000 

set.seed(1); do(50, "M4", FALSE)
## -0.58236436 -0.81292789 -0.58898479 -0.80814412  0.20725537  0.97828687 
##  0.20725537  0.97828687  0.05281645  0.03351502  0.06749574  0.04468926 
##                                                         95%         95% 
##  0.04674488  0.01047292  0.06570734  0.01518660  1.00000000  1.00000000 

set.seed(1); do2(100, "M1", FALSE)
## [1] 0.9207841 0.3900727 0.3896230 0.9209744 0.3308185 0.9436944
set.seed(1); do2(50, "M4", FALSE)
## [1] -0.5823644 -0.8129279  0.2072554  0.9782869  0.2094693  0.9778152

set.seed(1)
system.time(print(do2(50, "M1", TRUE)))
## 0.7747495 0.6322683 0.3468791 0.9379099 0.2405921 0.9706263

cl <- makePSOCKcluster(8)
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do', 'do2')))
invisible(clusterEvalQ(NULL, library(SSIndex)))

runPara <- function(model, n, frailty) {
    obj <- paste(c(model, n, frailty), collapse = "")
    eval(parse(text = paste(obj, " <- NULL")))
    toRun <- paste(obj, " <- t(matrix(unlist(parLapply(NULL, 1:200, function(z) do(",
                   n, ",'", model, "',", frailty, "))), 18))", sep = "")
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



##################
## Test hypothesis

n <- 100
model <- "M2"
frailty  <- FALSE

## March 17, 2019
do <- function(n, model, frailty = FALSE) {
    B <- 200
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    bi <- seq(0, 2 * pi, length = 100)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))    
    k02 <- sapply(bi, function(x) getk02(dat, c(cos(x), sin(x)), fit$Fhat0))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    n <- length(unique(dat$id))
    getBootk <- function(dat) {
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        fit0 <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
                    data = dat0, shp.ind = FALSE)        
        c(max(sapply(1:length(bi), function(x)
            getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x])),
          max(sapply(1:length(bi), function(x)
            getk02(dat0, c(cos(bi[x]), sin(bi[x])), fit0$Fhat0) - k02[x])))
    }
    tmp <- replicate(B, getBootk(dat))
    c(1 * (max(k0) > quantile(tmp[1,], .95)),
      1 * (max(k02) > quantile(tmp[2,], .95)))
}

system.time(print(do(100, "M2")))


cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, 'do'))
invisible(clusterEvalQ(NULL, library(SSIndex)))
set.seed(1)
f <- parSapply(NULL, 1:50, function(z) do(100, "M2"))
stopCluster(cl)
## rowMeans(f)
##    95%  95% 
##    0.26 0.96

dat <- simDat(200, "M2", TRUE)
head(dat)

pValShape(reSurv(Time, id, event, status) ~ x1 + x2, dat0 = dat)

do <- function(fname, B = 100, dat0 = dat0) {
    p <- length(attr(terms(fname), "term.labels"))
    dat1 <- dat0
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .1)), p)))
    r <- apply(tmp, 1, function(z) sqrt(sum(z^2)))
    bi <- (tmp / r)[r < 1 & r > 0,]
    k0 <- sapply(1:NROW(bi), function(x) getk0(dat1, bi[x,]))
    ## k02 <- sapply(1:NROW(bi), function(x) getk02(dat1, bi[x,], fit$Fhat0))
    getBootK <- function(dat) {
        n <- length(unique(dat$id))
        mm <- aggregate(event ~ id, dat, length)[, 2]
        ind <- sample(1:n, n, TRUE)
        datB <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        datB$id <- rep(1:n, mm[ind])
        datB <- datB[complete.cases(datB),]
        rownames(datB) <- NULL
        datB1 <- datB
        kb <- max(sapply(1:NROW(bi), function(x) getk0(datB1, bi[x,]) - k0[x]))
        ## getk0(datB1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1)) - k0[x]))
        max(kb)
    }
    cl <- makePSOCKcluster(8)
    ## cl <- makePSOCKcluster(16)
    setDefaultCluster(cl)
    invisible(clusterExport(cl, c("bi", "k0", "fname", "dat0", "p", "getBootK"),
                            environment()))
    invisible(clusterEvalQ(NULL, library(GSM)))
    invisible(clusterEvalQ(NULL, library(reReg)))
    system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootK(dat0))) 
    stopCluster(cl)
    mean(max(k0) > tmp)
}


do <- function(n, model, frailty = FALSE) {
    B <- 200
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    p <- 2
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .05)), p)))
    r <- apply(tmp, 1, function(z) sqrt(sum(z^2)))
    bi <- (tmp / r)[r < 1 & r > 0,]
    k0 <- sapply(1:NROW(bi), function(x) getk0(dat, bi[x,]))
    k02 <- sapply(1:NROW(bi), function(x) getk02(dat, bi[x,], fit$Fhat0))    
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    n <- length(unique(dat$id))
    getBootk <- function(dat) {
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        fit0 <- gsm(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
                    data = dat0, shp.ind = FALSE)
        k0B <- sapply(1:NROW(bi), function(x) getk0(dat0, bi[x,]) - k0[x])
        k02B <- sapply(1:NROW(bi), function(x) getk02(dat0, bi[x,], fit0$Fhat0) - k02[x])
        c(max(k0B), max(k02B), fit0$b0, fit0$b00, fit0$r0, fit0$r00)
    }
    tmp <- replicate(B, getBootk(dat))
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta, (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0, fit$b00, fit$r0, fit$r00,
      mean(k0 > quantile(tmp[1,], .95)),
      mean(k02 > quantile(tmp[2,], .95)),
      diag(var(t(tmp[3:10,]))))
}

system.time(f <- do(200, "M3", frailty = FALSE))

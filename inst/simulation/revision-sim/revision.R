#######################################################################
## Load package
#######################################################################

library(survival)
library(SSIndex)
library(methods)
library(tidyverse)
library(Rcpp)


source('pietg/fpigsim.R')
source('pietg/fisher.R')
library(locfit)
library(MAVE)
library(EDR)
library(MASS)
sourceCpp("pietg/SSE.cpp")
sourceCpp("pietg/ESE.cpp")
sourceCpp("pietg/LSE.cpp")
sourceCpp("pietg/spline.cpp")


#######################################################################
## Load function
#######################################################################


## Assume two dimentional
## Assume everything is length 476
sourceCpp(code = '
  #include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]
  using namespace arma;
  // [[Rcpp::export]]
double getb(double theta, double h, arma::vec tij, arma::vec yi, arma::mat X) {
    arma::vec b = zeros<vec>(2);
    arma::vec J = zeros<vec>(2);
    int n = yi.n_elem;
    double out = 0;
    b[0] = sin(theta);
    b[1] = cos(theta);
    J[0] = cos(theta);
    J[1] = -sin(theta);
    arma::vec xb = X * b;
    arma::vec xJ = X * J;
    // Rcpp::Rcout << xJ;
    double nu = 0;
    double de = 0;
    double Kh = 0;
    for (int i = 0; i < n; i++) {
      nu = 0;
      de = 0;
      for (int j = 0; j < n; j++) {
        if(yi[j] >= tij[i] & tij[j] <= tij[i]) {
          double dx = (xb[j] - xb[i]) / h;
          Kh = 0;
          if (dx <= 1 & dx >= -1) Kh = 0.75 * (1 - dx * dx) / h;
          de += Kh;
          nu += Kh * xJ[j];
        }
      }
      out += xJ[i] - nu / de;
    }
    return(out);
}')

sdOut <- function(x) {
    rm <- which(x %in% boxplot(x, plot = FALSE)$out)
    if (length(rm) > 0) x <- x[-rm]
    sd(x)
}

## Test only
do <- function(n, model, frailty = FALSE, type1 = FALSE, B = 200) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty, type1)
    p <- 2
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .1)), p)))
    r <- rowSums(tmp * tmp)
    keep <- which(r < 1 & r > 0)
    bi <- tmp[keep,] / sqrt(r[keep])
    k0.tmp <- getk0s(dat, bi) ## sapply(1:NROW(bi), function(x) getk04(dat, bi[x,]))
    k0 <- k0.tmp[1,]
    k02 <- k0.tmp[2,]
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    n <- length(unique(dat$id))
    getBootk <- function(dat) {
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
	dat0$Time[duplicated(dat0$Time) & dat0$Time < max(dat0$Time)] <-
            dat0$Time[duplicated(dat0$Time) & dat0$Time < max(dat0$Time)] +
            abs(rnorm(sum(duplicated(dat0$Time) & dat0$Time < max(dat0$Time)), sd = .0001))
        dat0$id <- rep(1:n, mm[ind] + 1)
        k0B.tmp <- getk0s(dat0, bi)
        k0B <- k0B.tmp[1,] - k0
        k02B <- k0B.tmp[2,] - k02
	c(max(k0B), max(k02B))
    }
    tmp <- replicate(B, getBootk(dat))
    c(mean(max(k0) > tmp[1,]), 
      mean(max(k02) > tmp[2,]))
}

pietg <- function(formula, data, B = 100, bIni = NULL, rIni = NULL) {
    ## Extract vectors
    Call <- match.call()
    if (missing(data)) obj <- eval(formula[[2]], parent.frame()) 
    if (!missing(data)) obj <- eval(formula[[2]], data) 
    ## if (!is.reSurv(obj)) stop("Response must be a reSurv object")
    formula[[2]] <- NULL
    if (formula == ~ 1) {
        DF <- cbind(obj$reDF[, -5], zero=0)
    } else {
        if (!missing(data)) DF <- cbind(obj$reDF[,-5], model.matrix(formula, data))
        if (missing(data)) DF <- cbind(obj$reDF[,-5], model.matrix(formula, parent.frame()))
        DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    }
    DF <- DF[order(DF$id, DF$Time), ]
    DF$id <- rep(1:length(unique(DF$id)), table(DF$id))
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    dat <- DF %>% add_column(m = rep(m, m + 1), .after = 4)
    ## computation
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id)) ## n = 80
    mm <- aggregate(event ~ id, dat0, sum)[, 2] ## 80 by 1
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time ## 476 by 1
    yi <- subset(dat0, event == 0)$Time ## 80 by 1
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status))) ## 80 by 2
    h <- .5
    ## Different ways to find root; choose the best root
    b1 <- BB::spg(par = 0, fn = function(x)
        getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),])^2,
        lower = 0, upper = pi, quiet = TRUE, alertConvergence = FALSE)$par
    b2 <- optimize(f = function(x)
        getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),])^2,
        interval = c(0, pi))$minimum
    gg <- seq(0, pi, length.out = 21)
    g0 <- sapply(gg, function(x)
        getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),]))
    b3 <- sapply(which(diff(sign(g0)) != 0), function(x)
        uniroot(f = function(x) getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),]),
                interval = c(gg[x], gg[x + 1]))$root)
    b <- c(b1, b2, b3)
    b <- b[which.min(sapply(b, function(x)
        abs(getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),]))))]
    ## c(sin(b), cos(b))
    bhat <- c(sin(b), cos(b))
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
           as.double(tij), as.double(yi),
           as.double(X %*% bhat), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "SSIndex")$result,
        X %*% bhat, yi))
    Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
    Fhat <- exp(-Fhat)
    LSE <- ComputeLSE(X, mm / Fhat)$alpha
    SSE <- ComputeSSE(X, mm / Fhat)$alpha
    ESE <- ComputeESE(X, mm / Fhat)$alpha
    spline <- Compute_spline(X, mm / Fhat)$alpha
    MAVE <- as.numeric(mave.compute(X, mm / Fhat, method = "meanmave")$dir[[1]])
    list(b0 = bhat, LSE = LSE, SSE = SSE, ESE = ESE, spline = spline, MAVE = MAVE)
}


seed <- sample(1:1e7, 1)
set.seed(seed)
dat <- simDat(200, "M2", TRUE, FALSE)
head(dat, 20)

## Implementation given dat; generated by simDat()

fit <- gsm(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
           data = dat, shp.ind = FALSE)

fit <- pietg(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
             data = dat)

do <- function(n, model, frailty = FALSE, type1 = FALSE) {
    dat <- simDat(n, model, frailty, type1)
    fit <- pietg(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
                 data = dat)
    unlist(fit)
}

M2n200 <- replicate(500, do(200, "M2"))
M3n200 <- replicate(500, do(200, "M3"))
M4n200 <- replicate(500, do(200, "M4"))
M5n200 <- replicate(500, do(200, "M5"))
save(M2n200, file = "M2n200.RData")
save(M3n200, file = "M3n200.RData")
save(M4n200, file = "M4n200.RData")
save(M5n200, file = "M5n200.RData")

M2n400 <- replicate(500, do(400, "M2"))
M3n400 <- replicate(500, do(400, "M3"))
M4n400 <- replicate(500, tryCatch(do(400, "M4"), error = function(e) do(400, "M4")))
M5n400 <- replicate(500, do(400, "M5"))
save(M2n400, file = "M2n400.RData")
save(M3n400, file = "M3n400.RData")
save(M4n400, file = "M4n400.RData")
save(M5n400, file = "M5n400.RData")



M2n200c <- replicate(500, do(200, "M2", TRUE))
M3n200c <- replicate(500, do(200, "M3", TRUE))
M4n200c <- replicate(500, do(200, "M4", TRUE))
M5n200c <- replicate(500, do(200, "M5", TRUE))
save(M2n200c, file = "M2n200c.RData")
save(M3n200c, file = "M3n200c.RData")
save(M4n200c, file = "M4n200c.RData")
save(M5n200c, file = "M5n200c.RData")

M2n400c <- replicate(500, do(400, "M2", TRUE))
M3n400c <- replicate(500, do(400, "M3", TRUE))
M4n400c <- replicate(500, do(400, "M4", TRUE))
M5n400c <- replicate(500, do(400, "M5", TRUE))
save(M2n400c, file = "M2n400c.RData")
save(M3n400c, file = "M3n400c.RData")
save(M4n400c, file = "M4n400c.RData")
save(M5n400c, file = "M5n400c.RData")


load("M2n200c.RData")
load("M3n200c.RData")
load("M4n200c.RData")
load("M5n200c.RData")
load("M2n400c.RData")
load("M3n400c.RData")
load("M4n400c.RData")
load("M5n400c.RData")
load("M2n200.RData")
load("M3n200.RData")
load("M4n200.RData")
load("M5n200.RData")
load("M2n400.RData")
load("M3n400.RData")
load("M4n400.RData")
load("M5n400.RData")

sumTab <- function(dat) {
    bhat <- c(mean(sin(dat[1,])), mean(cos(dat[1,])))
    lse <- rowMeans(dat[2:3,])
    sse <- rowMeans(dat[4:5,])
    ese <- rowMeans(dat[6:7,])
    mave <- rowMeans(dat[10:11,])
    cbind(bhat, lse, sse, ese, mave)
}

sumTab(M2n200)
sumTab(M3n200)
sumTab(M4n200)
sumTab(M5n200)

sumTab(M2n400)
sumTab(M3n400)
sumTab(M4n400)
sumTab(M5n400)

## true value
## b = c(.6, .8)
## gamma = c(.28, .96)


plot(sin(M2n200[1,]))
plot(sin(M3n200[1,]))
plot(sin(M4n200[1,]))
plot(sin(M5n200[1,]))

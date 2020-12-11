#######################################################################
## Load package
#######################################################################

library(survival)
library(SSIndex)
library(methods)
library(tidyverse)
library(Rcpp)
library(xtable)

source('pietg/fpigsim.R')
source('pietg/fisher.R')
library(locfit)
library(MAVE)
library(EDR)
library(MASS)
sourceCpp("pietg/SSE.cpp")
sourceCpp("pietg/SSE2.cpp")
sourceCpp("pietg/ESE.cpp")
sourceCpp("pietg/ESE2.cpp")
sourceCpp("pietg/LSE.cpp")
sourceCpp("pietg/LSE2.cpp")
sourceCpp("pietg/spline.cpp")

#######################################################################
## Load function
#######################################################################

sourceCpp(code = '
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] 
double M(arma::vec& beta,
         arma::mat& Z,
         arma::vec& T2,
         arma::vec& C2,
         double h,
         double ht)
{
  int n2 = T2.n_elem;
  arma::mat kmat(n2, n2);
  arma::mat bzmat(n2, n2);
  arma::mat imat(n2, n2);
  kmat.zeros();
  bzmat.zeros();
  imat.zeros();
  arma::vec bz = Z * beta;
  for(int i = 0; i < n2; i++)
  {
    for(int j = 0; j < i-1; j++)
    {
      imat(i,j) = T2(j) <= T2(i) && T2(i) <= C2(j);
      double x = (T2(i)-T2(j))/ht;
      if(x < 1 && x>-1)
      {
        kmat(i,j) = 0.75*(1-x*x)/ht;
      }
      x = (bz(i)-bz(j))/h;
      if(x < 1 && x>-1)
      {
        bzmat(i,j) = 0.75*(1-x*x)/h;
      }
    }
  }
  kmat = kmat + kmat.t();
  bzmat = bzmat + bzmat.t();
  imat = imat + imat.t();
  arma::vec ov(n2);
  ov.ones();
  kmat.diag() = ov*0.75/ht;
  bzmat.diag() = ov*0.75/h;
  imat.diag() = ov;
  double m1 = arma::sum(log(arma::sum(kmat % bzmat))) - arma::sum( log(arma::sum(kmat % imat)) );
  double m2 = 0;
  for(int i = 0; i < n2; i++)
  {
    for(int j = 0; i < n2; i++)
    {
      double dij = 0;
      double nij = kmat(i,j)*imat(i,j);
      for(int l = 0; l < n2; l++)
      {
        dij = dij + kmat(l,j)*imat(i,l);
      }
      if(dij > 0)
      {
        m2 += nij/dij;
      }
    }
  }
  return m1 - m2;
}')

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


## Both EE approaches
## pietg <- function(formula, data, B = 100, bIni = NULL, rIni = NULL) {
##     ## Extract vectors
##     Call <- match.call()
##     if (missing(data)) obj <- eval(formula[[2]], parent.frame()) 
##     if (!missing(data)) obj <- eval(formula[[2]], data) 
##     ## if (!is.reSurv(obj)) stop("Response must be a reSurv object")
##     formula[[2]] <- NULL
##     if (formula == ~ 1) {
##         DF <- cbind(obj$reDF[, -5], zero=0)
##     } else {
##         if (!missing(data)) DF <- cbind(obj$reDF[,-5], model.matrix(formula, data))
##         if (missing(data)) DF <- cbind(obj$reDF[,-5], model.matrix(formula, parent.frame()))
##         DF <- DF[,-which(colnames(DF) == "(Intercept)")]
##     }
##     DF <- DF[order(DF$id, DF$Time), ]
##     DF$id <- rep(1:length(unique(DF$id)), table(DF$id))
##     m <- aggregate(event ~ id, data = DF, sum)[,2]
##     dat <- DF %>% add_column(m = rep(m, m + 1), .after = 4)
##     ## computation
##     dat0 <- subset(dat, m > 0)
##     n <- length(unique(dat0$id)) ## n = 80
##     mm <- aggregate(event ~ id, dat0, sum)[, 2] ## 80 by 1
##     dat0$id <- rep(1:n, mm + 1)
##     tij <- subset(dat0, event == 1)$Time ## 476 by 1
##     yi <- subset(dat0, event == 0)$Time ## 80 by 1
##     midx <- c(0, cumsum(mm)[-length(mm)])
##     X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status))) ## 80 by 2
##     h <- 1.06 * sd(X %*% c(.6, .8)) * nrow(X)^-.2
##     ## h <- .5
##     ## Different ways to find root; choose the best root
##     b1 <- BB::spg(par = 0, fn = function(x)
##         getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),])^2,
##         lower = 0, upper = pi, quiet = TRUE, alertConvergence = FALSE)$par
##     b2 <- optimize(f = function(x)
##         getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),])^2,
##         interval = c(0, pi))$minimum
##     gg <- seq(0, pi, length.out = 21)
##     g0 <- sapply(gg, function(x)
##         getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),]))
##     b3 <- sapply(which(diff(sign(g0)) != 0), function(x)
##         uniroot(f = function(x) getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),]),
##                 interval = c(gg[x], gg[x + 1]))$root)
##     b <- c(b1, b2, b3)
##     b <- b[which.min(sapply(b, function(x)
##         abs(getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),]))))]
##     ## c(sin(b), cos(b))
##     X0 <- as.matrix(subset(dat, status == 1, select = -c(Time, id, m, event, status)))
##     bhat <- c(sin(b), cos(b))
##     Fhat <- unlist(mapply(FUN = function(x, y)
##         .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
##            as.double(tij), as.double(yi),
##            as.double(X %*% bhat), as.double(x), as.double(y), as.double(h), 
##            result = double(1), PACKAGE = "SSIndex")$result,
##         X0 %*% bhat, subset(dat, status == 1)$Time))
##     Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
##     Fhat <- exp(-Fhat)    
##     LSE <- ComputeLSE(X0, m / Fhat)$alpha
##     SSE <- ComputeSSE(X0, m / Fhat)$alpha
##     ESE <- ComputeESE(X0, m / Fhat)$alpha
##     spline <- Compute_spline(X0, m / Fhat)$alpha
##     MAVE <- as.numeric(mave.compute(X0, m / Fhat, method = "meanmave")$dir[[1]])
##     MAVE <- abs(MAVE)
##     suppressWarnings(EFM <- fisher(X0, m / Fhat, 1, mymodel = "none")$root)
##     est1 <- list(b0 = bhat, LSE = LSE, SSE = SSE, ESE = ESE, spline = spline, MAVE = MAVE, EFM = EFM)
##     bhat <- LSE <- SSE <- ESE <- spline <- MAVE <- EFM <- NULL
##     ## Use beta from M-approach
##     dat$Y <- rep(dat$Time[dat$status == 1], m + 1)
##     T2 <- dat$Time[dat$status == 0]
##     C2 <- dat$Y[dat$status == 0]
##     Z <- as.matrix(subset(dat, status == 0)[,c("x1", "x2")])
##     b2 <- optimize(f = function(x) -M(c(sin(x), cos(x)), Z, T2, C2, 0.2, 0.5),
##                    interval = c(0, pi))$minimum
##     bhat <- c(sin(b2), cos(b2))
##     Fhat <- unlist(mapply(FUN = function(x, y)
##         .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
##            as.double(tij), as.double(yi),
##            as.double(X %*% bhat), as.double(x), as.double(y), as.double(h), 
##            result = double(1), PACKAGE = "SSIndex")$result,
##         X0 %*% bhat, subset(dat, status == 1)$Time))
##     Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
##     Fhat <- exp(-Fhat)
##     LSE <- ComputeLSE(X0, m / Fhat)$alpha
##     SSE <- ComputeSSE(X0, m / Fhat)$alpha
##     ESE <- ComputeESE(X0, m / Fhat)$alpha
##     spline <- Compute_spline(X0, m / Fhat)$alpha
##     MAVE <- as.numeric(mave.compute(X0, m / Fhat, method = "meanmave")$dir[[1]])
##     MAVE <- abs(MAVE)
##     suppressWarnings(EFM <- fisher(X0, m / Fhat, 1, mymodel = "none")$root)
##     est2 <- list(b0 = bhat, LSE = LSE, SSE = SSE, ESE = ESE, spline = spline, MAVE = MAVE, EFM = EFM)
##     list(est1 = est1, est2 = est2)    
## }

## Only M approach
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
    ## h <- .5
    ## Different ways to find root; choose the best root    
    ## Use beta from M-approach
    dat$Y <- rep(dat$Time[dat$status == 1], m + 1)
    T2 <- dat$Time[dat$status == 0]
    C2 <- dat$Y[dat$status == 0]
    Z <- as.matrix(subset(dat, status == 0)[,c("x1", "x2")])
    h <- 1.06 * sd(Z %*% c(.6, .8)) * nrow(Z)^-.2
    ht <- 1.06 * sd(T2) * length(T2)^-.2
    ## print(c(h, ht))    
    b2 <- optimize(f = function(x) -M(c(sin(x), cos(x)), Z, T2, C2, h, ht),
                   interval = c(0, pi))$minimum
    bhat <- c(sin(b2), cos(b2))
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id)) 
    mm <- aggregate(event ~ id, dat0, sum)[, 2] 
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status, Y)))     
    X0 <- as.matrix(subset(dat, status == 1, select = -c(Time, id, m, event, status, Y)))    
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
           as.double(tij), as.double(yi),
           as.double(X %*% bhat), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "SSIndex")$result,
        X0 %*% bhat, subset(dat, status == 1)$Time))
    Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
    Fhat <- exp(-Fhat)
    LSE <- ComputeLSE(X0, m / Fhat)$alpha
    SSE <- ComputeSSE(X0, m / Fhat)$alpha
    ESE <- ComputeESE(X0, m / Fhat)$alpha
    spline <- Compute_spline(X0, m / Fhat)$alpha
    MAVE <- as.numeric(mave.compute(X0, m / Fhat, method = "meanmave")$dir[[1]])
    MAVE <- abs(MAVE)
    suppressWarnings(EFM <- fisher(X0, m / Fhat, 1, mymodel = "none")$root)
    list(b0 = bhat, LSE = LSE, SSE = SSE, ESE = ESE, spline = spline, MAVE = MAVE, EFM = EFM)
}

do <- function(n, model, frailty = FALSE, type1 = FALSE, noCen = FALSE) {
    dat <- simDat(n, model, frailty, type1, noCen = noCen)
    fit <- pietg(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
                 data = dat)
    unlist(fit)
}

set.seed(1); do(200, "M2", TRUE)
set.seed(10); do(200, "M2", TRUE) ## -1
set.seed(10); do(200, "M2", TRUE, noCen = TRUE)
set.seed(5); matrix(do(200, "M2", TRUE), ncol = 2)

matrix(do(200, "M2", TRUE, noCen = TRUE), ncol = 2)
f <- replicate(10, do(200, "M2", TRUE, noCen = TRUE))

f <- NULL
for (i in 1:10) {
    f <- cbind(f, do(200, "M2", TRUE, noCen = TRUE))
}

## Yifei code
set.seed(55)
dat <- simDat(200, "M2", TRUE)
m <- dat$m[dat$status == 1]
Y <- dat$Time[dat$status == 1]
dat$Y <- rep(Y, m+1)
T2 <- dat$Time[dat$status == 0]
C2 <- dat$Y[dat$status == 0]
Z <- as.matrix(dat[,c("x1","x2")])
thetaall <- seq(from = 0, to = pi, length = 100)
mall <- rep(0, 100)
for(i in 1:100) {
    theta <- thetaall[i]
    beta <- c(sin(theta),cos(theta))
    mall[i] <- M(beta, Z, T2, C2, 0.2, 0.5)
}
plot(thetaall, mall, 'l')


M2n200 <- replicate(500, do(200, "M2"))
M3n200 <- replicate(500, do(200, "M3"))
M4n200 <- replicate(500, do(200, "M4"))
M5n200 <- replicate(500, do(200, "M5"))
save(M2n200, file = "M2n200.RData")
save(M3n200, file = "M3n200.RData")
save(M4n200, file = "M4n200.RData")
save(M5n200, file = "M5n200.RData")

M2n200noCen <- replicate(500, do(200, "M2", noCen = TRUE))
M3n200noCen <- replicate(500, do(200, "M3", noCen = TRUE))
M4n200noCen <- replicate(500, do(200, "M4", noCen = TRUE))
M5n200noCen <- replicate(500, do(200, "M5", noCen = TRUE))
save(M2n200noCen, file = "M2n200noCen.RData")
save(M3n200noCen, file = "M3n200noCen.RData")
save(M4n200noCen, file = "M4n200noCen.RData")
save(M5n200noCen, file = "M5n200noCen.RData")

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

M2n200cnoCen <- replicate(500, do(200, "M2", TRUE, noCen = TRUE))
M3n200cnoCen <- replicate(500, do(200, "M3", TRUE, noCen = TRUE))
M4n200cnoCen <- replicate(500, do(200, "M4", TRUE, noCen = TRUE))
M5n200cnoCen <- replicate(500, do(200, "M5", TRUE, noCen = TRUE))
save(M2n200cnoCen, file = "M2n200cnoCen.RData")
save(M3n200cnoCen, file = "M3n200cnoCen.RData")
save(M4n200cnoCen, file = "M4n200cnoCen.RData")
save(M5n200cnoCen, file = "M5n200cnoCen.RData")

M2n400c <- replicate(500, tryCatch(do(400, "M2", TRUE), error = function(e) do(400, "M2", TRUE)))
M3n400c <- replicate(500, do(400, "M3", TRUE))
M4n400c <- replicate(500, do(400, "M4", TRUE))
M5n400c <- replicate(500, do(400, "M5", TRUE))
save(M2n400c, file = "M2n400c.RData")
save(M3n400c, file = "M3n400c.RData")
save(M4n400c, file = "M4n400c.RData")
save(M5n400c, file = "M5n400c.RData")


M6n200 <- replicate(500, tryCatch(do(200, "M6", FALSE), error = function(e) do(200, "M6", FALSE)))
save(M6n200, file = "M6n200.RData")
M6n200c <- replicate(500, tryCatch(do(200, "M6", TRUE), error = function(e) do(200, "M6", TRUE)))
save(M6n200c, file = "M6n200c.RData")
M6n400 <- replicate(500, tryCatch(do(400, "M6", FALSE), error = function(e) do(400, "M6", FALSE)))
save(M6n400, file = "M6n400.RData")
M6n400c <- replicate(500, tryCatch(do(400, "M6", TRUE), error = function(e) do(400, "M6", TRUE)))
save(M6n400c, file = "M6n400c.RData")

load("revision-sim/M2n200c.RData")
load("revision-sim/M3n200c.RData")
load("revision-sim/M4n200c.RData")
load("revision-sim/M5n200c.RData")
load("revision-sim/M2n400c.RData")
load("revision-sim/M3n400c.RData")
load("revision-sim/M4n400c.RData")
load("revision-sim/M5n400c.RData")
load("revision-sim/M2n200.RData")
load("revision-sim/M3n200.RData")
load("revision-sim/M4n200.RData")
load("revision-sim/M5n200.RData")
load("revision-sim/M2n400.RData")
load("revision-sim/M3n400.RData")
load("revision-sim/M4n400.RData")
load("revision-sim/M5n400.RData")
load("revision-sim/M6n200.RData")
load("revision-sim/M6n400.RData")
load("revision-sim/M6n200c.RData")
load("revision-sim/M6n400c.RData")

load("M2n200noCen.RData")
load("M3n200noCen.RData")
load("M4n200noCen.RData")
load("M5n200noCen.RData")

matrix(rowMeans(M2n200), ncol = 2) - c(.6, .8)
matrix(rowMeans(M3n200), ncol = 2) - c(.6, .8)
matrix(rowMeans(M4n200), ncol = 2) - c(.6, .8, rep(c(.28, .96), 6))
matrix(rowMeans(M5n200), ncol = 2) - c(.6, .8)

matrix(rowMeans(M2n200), ncol = 1) - c(.6, .8)
matrix(rowMeans(M3n200), ncol = 1) - c(.6, .8)
matrix(rowMeans(M4n200), ncol = 1) - c(.6, .8, rep(c(.28, .96), 6))
matrix(rowMeans(M5n200), ncol = 1) - c(.6, .8)

matrix(rowMeans(M2n200noCen), ncol = 2) - c(.6, .8)
matrix(rowMeans(M3n200noCen), ncol = 2) - c(.6, .8)
matrix(rowMeans(M4n200noCen), ncol = 2) - c(.6, .8, rep(c(.28, .96), 6))
matrix(rowMeans(M5n200noCen), ncol = 2) - c(.6, .8)

cbind((matrix(rowMeans(M2n200noCen), ncol = 2) - c(.6, .8))[1:12,],
      matrix(rowMeans(M2n200), ncol = 1) - c(.6, .8))      
cbind((matrix(rowMeans(M3n200noCen), ncol = 2) - c(.6, .8))[1:12,],
      matrix(rowMeans(M3n200), ncol = 1) - c(.6, .8))      
cbind((matrix(rowMeans(M4n200noCen), ncol = 2) - c(.6, .8))[1:12,],
      matrix(rowMeans(M4n200), ncol = 1) - c(.6, .8))      
cbind((matrix(rowMeans(M5n200noCen), ncol = 2) - c(.6, .8))[1:12,],
      matrix(rowMeans(M5n200), ncol = 1) - c(.6, .8))      

M2n200[10:11,] <- abs(M2n200[10:11,])
M2n200c[10:11,] <- abs(M2n200c[10:11,])
M2n400[10:11,] <- abs(M2n400[10:11,])
M2n400c[10:11,] <- abs(M2n400c[10:11,])
M3n200[10:11,] <- abs(M3n200[10:11,])
M3n200c[10:11,] <- abs(M3n200c[10:11,])
M3n400[10:11,] <- abs(M3n400[10:11,])
M3n400c[10:11,] <- abs(M3n400c[10:11,])
M4n200[10:11,] <- abs(M4n200[10:11,])
M4n200c[10:11,] <- abs(M4n200c[10:11,])
M4n400[10:11,] <- abs(M4n400[10:11,])
M4n400c[10:11,] <- abs(M4n400c[10:11,])
M5n200[10:11,] <- abs(M5n200[10:11,])
M5n200c[10:11,] <- abs(M5n200c[10:11,])
M5n400[10:11,] <- abs(M5n400[10:11,])
M5n400c[10:11,] <- abs(M5n400c[10:11,])


sumTab <- function(dat) {
    bhat <- c(mean(sin(dat[1,])), mean(cos(dat[1,])))
    lse <- rowMeans(dat[2:3,])
    sse <- rowMeans(dat[4:5,])
    ese <- rowMeans(dat[6:7,])
    mave <- rowMeans(dat[10:11,])
    cbind(bhat, lse, sse, ese, mave)
}

sumTab <- function(dat) {
    rowMeans(dat)
}

sumTab(M2n200)
sumTab(M3n200)
sumTab(M4n200)
sumTab(M5n200)

## true value
## b = c(.6, .8)
## gamma = c(.28, .96)


plot(sin(M2n200[1,]))
plot(sin(M3n200[1,]))
plot(sin(M4n200[1,]))
plot(sin(M5n200[1,]))


print(xtable(cbind(abs(rowMeans(M2n200) - c(.6, .8)),
             apply(M2n200, 1, sd),
             abs(rowMeans(M2n200c) - c(.6, .8)),
             apply(M2n200c, 1, sd)), digits = 3), include.rownames = FALSE)

print(xtable(cbind(abs(rowMeans(M3n200) - c(.6, .8)),
             apply(M3n200, 1, sd),
             abs(rowMeans(M3n200c) - c(.6, .8)),
             apply(M3n200c, 1, sd)), digits = 3), include.rownames = FALSE)

print(xtable(cbind(abs(rowMeans(M5n200) - c(.6, .8)),
             apply(M5n200, 1, sd),
             abs(rowMeans(M5n200c) - c(.6, .8)),
             apply(M5n200c, 1, sd)), digits = 3), include.rownames = FALSE)

print(xtable(cbind(abs(rowMeans(M4n200) - c(.6, .8, rep(c(.28, .96), 5))), 
             apply(M4n200, 1, sd),
             abs(rowMeans(M4n200c) - c(.6, .8, rep(c(.28, .96), 5))), 
             apply(M4n200c, 1, sd)), digits = 3), include.rownames = FALSE)


## n = 400


print(xtable(cbind(abs(rowMeans(M2n400) - c(.6, .8)),
             apply(M2n400, 1, sd),
             abs(rowMeans(M2n400c) - c(.6, .8)),
             apply(M2n400c, 1, sd)), digits = 3), include.rownames = FALSE)

print(xtable(cbind(abs(rowMeans(M3n400) - c(.6, .8)),
             apply(M3n400, 1, sd),
             abs(rowMeans(M3n400c) - c(.6, .8)),
             apply(M3n400c, 1, sd)), digits = 3), include.rownames = FALSE)

print(xtable(cbind(abs(rowMeans(M5n400) - c(.6, .8)),
             apply(M5n400, 1, sd),
             abs(rowMeans(M5n400c) - c(.6, .8)),
             apply(M5n400c, 1, sd)), digits = 3), include.rownames = FALSE)

print(xtable(cbind(abs(rowMeans(M4n400) - c(.6, .8, rep(c(.28, .96), 5))), 
             apply(M4n400, 1, sd),
             abs(rowMeans(M4n400c) - c(.6, .8, rep(c(.28, .96), 5))), 
             apply(M4n400c, 1, sd)), digits = 3), include.rownames = FALSE)

## ###############################################################################
## Print for paper
## ###############################################################################

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
    ## "b01"     "b02"     "LSE1"    "LSE2"    "SSE1"    "SSE2"    "ESE1"   
    ## "ESE2"    "spline1" "spline2" "MAVE1"   "MAVE2"   "EFM1"    "EFM2"   
    rowMeans(dat)
}

rowMeans(M2n200)

## ###############################################################################
## Flip sign for LSE, SSE, and ESE
## ###############################################################################

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
    ## h <- .5
    ## Different ways to find root; choose the best root    
    ## Use beta from M-approach
    dat$Y <- rep(dat$Time[dat$status == 1], m + 1)
    T2 <- dat$Time[dat$status == 0]
    C2 <- dat$Y[dat$status == 0]
    Z <- as.matrix(subset(dat, status == 0)[,c("x1", "x2")])
    h <- 1.06 * sd(Z %*% c(.6, .8)) * nrow(Z)^-.2
    ht <- 1.06 * sd(T2) * length(T2)^-.2
    ## print(c(h, ht))    
    b2 <- optimize(f = function(x) -M(c(sin(x), cos(x)), Z, T2, C2, h, ht),
                   interval = c(0, pi))$minimum
    bhat <- c(sin(b2), cos(b2))
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id)) 
    mm <- aggregate(event ~ id, dat0, sum)[, 2] 
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status, Y)))     
    X0 <- as.matrix(subset(dat, status == 1, select = -c(Time, id, m, event, status, Y)))    
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
           as.double(tij), as.double(yi),
           as.double(X %*% bhat), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "SSIndex")$result,
        X0 %*% bhat, subset(dat, status == 1)$Time))
    Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
    Fhat <- exp(-Fhat)
    LSE <- ComputeLSE(X0, m / Fhat)$alpha
    SSE <- ComputeSSE(X0, m / Fhat)$alpha
    ESE <- ComputeESE(X0, m / Fhat)$alpha
    spline <- Compute_spline(X0, m / Fhat)$alpha
    MAVE <- as.numeric(mave.compute(X0, m / Fhat, method = "meanmave")$dir[[1]])
    MAVE <- abs(MAVE)
    suppressWarnings(EFM <- fisher(X0, m / Fhat, 1, mymodel = "none")$root)
    list(b0 = bhat, LSE = LSE, SSE = SSE, ESE = ESE, spline = spline, MAVE = MAVE, EFM = EFM)
}

library(SSIndex)

n <- 10000
sum(simDat(n, "M2", TRUE)$event) / n ## 6.9931
sum(simDat(n, "M2", FALSE)$event) / n ## 8.4387

sum(simDat(n, "M3", TRUE)$event) / n ## 1.8758
sum(simDat(n, "M3", FALSE)$event) / n ## 2.1208

sum(simDat(n, "M5", TRUE)$event) / n ## 2.27
sum(simDat(n, "M5", FALSE)$event) / n ## 2.1

sum(simDat(n, "M4", TRUE)$event) / n ## 6.2866
sum(simDat(n, "M4", FALSE)$event) / n ## 6.3037


d <- simDat(1e4, "M5", TRUE)
summary(unlist(lapply(split(d$event, d$id), sum)))
sum(d$event) / 1e4
summary(d$Time)

x0 <- rexp(1)
b0 <- seq(-10, 10, .01)


rh <- lam(x0, b0) * exp(-Lam(x0, b0)) / (1 - exp(-Lam(x0, b0)))
plot(b0, rh, 'l')

rh2 <- lam(x0, b0) / Lam(x0, b0)
plot(b0, rh2, 'l')
lines(b0, rh, col = 2)


library(SSIndex)

## M1 
lam <- function(x, b) 1 / (1 + x) + .5 * exp(b)
Lam <- function(x, b) (2 * log(1 + x) + x * exp(b))/2

## M2
lam <- function(x, b) exp(-.5 * x * exp(b))
Lam <- function(x, b) ((1 - exp(-x * exp(b)/2)) * 2 * exp(-b))

## M3
lam <- function(x, b) .5 * log(exp(b) + 1) * (1 + .5 * x)^(log(exp(b) + 1) - 1)
Lam <- function(x, b) .5 * (1 + .5 * x)^log(exp(b) + 1)
## lam <- function(x, b)
##     5 * exp(b / 2) * (1 + x) ^ (exp(b / 2) / 2)
## Lam <- function(x, b)
##     10 * ((1 + x)^ (exp(b / 2) / 2) - 1)

## M4
lam <- function(x, b) dbeta(x, 2, 1 + exp(b))
Lam <- function(x, b, r) 4 * (pbeta(x, 2, 1 + exp(b)) * exp(r))



library(SSIndex)

xb <- matrix(rnorm(1e6), ncol = 2) %*% c(.6, .8)
summary(xb)
summary(exp(xb))
summary(exp(-xb))
summary(1.2^(xb))
summary(log(abs(xb)))
summary(exp(xb - 2))
summary(exp(xb - 3))

summary(replicate(100, max(10 * (1 + rexp(1)) ^ (-exp(xb) / 2) - 10)))
summary(replicate(100, max(10 * (1 + rexp(1)) ^ (log(abs(xb))) - 10)))

summary(replicate(100, min(5 * (1 + rexp(1)) ^ (exp(xb - 2)) - 5)))

t0 <- seq(0, 10, .001)
plot(t0, .1 * ((1 + t0)^exp(4) - 1), 'l')

n <- 100
sum(simDat(n, "M5", TRUE)$event) / n 
sum(simDat(n, "M5", FALSE)$event) / n 

d <- simDat(1e4, "M5", TRUE)
summary(unlist(lapply(split(d$event, d$id), sum)))
sum(d$event) / 1e4
summary(d$Time)

do <- function(n) {
    d <- simDat(n, "M5", TRUE)
    f <- gsm(reSurv(Time, id, event, status) ~ x1 + x2, data = d)
    c(f$b0, f$b00)
}

do(200)

foo <- replicate(100, do(200))
summary(t(foo))
apply(foo, 1, mean)
apply(foo, 1, sd)

#######################################################################
## Check direction
#######################################################################
library(SSIndex)

fo <- function(n, m) {
    d <- simDat(n, m, FALSE)
    d$id <- d$id / 100
    f <- gsm(reSurv(Time, id, event, status) ~ x1 + x2, data = d)
    c(f$b0, f$b00, f$r0, f$r00)
}

set.seed(1); fo(200, "M2") ## ++
set.seed(1); fo(200, "M3") ## --
set.seed(1); fo(200, "M5") ## ++
set.seed(1); fo(200, "M4") ## -+


d <- simDat(200, "M4", FALSE)
head(d)
summary(d)
f <- gsm(reSurv(Time, id, event, status) ~ x1 + x2, data = d)

debug(gsm)



d2 <- simDat(200, "M3", FALSE)
f <- gsm(reSurv(Time, id, event, status) ~ x1 + x2, data = d2)

str(f)

#######################################################################
## Reading from ESE.cpp and modify lines
#######################################################################
root <- "https://raw.githubusercontent.com/pietg/single_index/master/Comparisons/"
txtSSE <- readLines(paste0(root, "SSE.cpp"))
txtSSE[46] <- "List ComputeSSE2(NumericMatrix X, NumericVector y)"
txtSSE[91] <- "    beta1 = golden(0,2*pi,criterion);"
sourceCpp(code = paste0(txtSSE, collapse="\n"))

txtESE <- readLines("pietg/ESE.cpp")
txtESE[47] <- "List ComputeESE2(NumericMatrix X, NumericVector y)"
txtESE[93] <- "    beta1 = golden(0,2*pi,criterion);"
sourceCpp(code = paste0(txtESE, collapse="\n"))

txtLSE <- readLines(paste0(root, "LSE.cpp"))
txtLSE[44] <- "List ComputeLSE2(NumericMatrix X, NumericVector y)"
txtLSE[90] <- "    beta1 = golden(0,2*pi,criterion);"
sourceCpp(code = paste0(txtLSE, collapse="\n"))

do1 <- function(n, m) {
    dat <- simDat(n, model = m)
    m <- unlist(lapply(split(dat$Time, dat$id), function(x) sum(!is.na(x))))
    X0 <- as.matrix(subset(dat, select = c(x1, x2))[!duplicated(dat$id),])
    list(LSE = ComputeLSE2(X0, m)$alpha,
         SSE = ComputeSSE2(X0, m)$alpha,
         ESE = ComputeESE2(X0, m)$alpha)
}

set.seed(1); do1(200)

pietg3 <- function(formula, data, B = 100, bIni = NULL, rIni = NULL) {
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
    ## h <- .5
    ## Different ways to find root; choose the best root    
    ## Use beta from M-approach
    dat$Y <- rep(dat$Time[dat$status == 1], m + 1)
    T2 <- dat$Time[dat$status == 0]
    C2 <- dat$Y[dat$status == 0]
    Z <- as.matrix(subset(dat, status == 0)[,c("x1", "x2")])
    h <- 1.06 * sd(Z %*% c(.6, .8)) * nrow(Z)^-.2
    ht <- 1.06 * sd(T2) * length(T2)^-.2
    ## print(c(h, ht))
    bhat <- c(.6, .8)
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id)) 
    mm <- aggregate(event ~ id, dat0, sum)[, 2] 
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status, Y)))     
    X0 <- as.matrix(subset(dat, status == 1, select = -c(Time, id, m, event, status, Y)))    
    LSE <- ComputeLSE2(X0, m)$alpha
    SSE <- ComputeSSE2(X0, m)$alpha
    ESE <- ComputeESE2(X0, m)$alpha
    list(LSE = LSE, SSE = SSE, ESE = ESE)
}

do3 <- function(n, model, frailty = FALSE, type1 = FALSE, noCen = FALSE) {
    dat <- simDat(n, model, frailty, type1, noCen = noCen)
    fit <- pietg3(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
                 data = dat)
    unlist(fit)
}

set.seed(1); f1 <- replicate(50, do3(200, "M2"))
set.seed(1); f2 <- replicate(50, do3(200, "M3"))
set.seed(1); f3 <- replicate(50, do3(200, "M5"))
set.seed(1); f4 <- replicate(50, do3(200, "M4"))


library(ggplot2)
library(reshape2)
tmp <- do.call(rbind, foo[1:50 * 3])
ggplot(melt(tmp), aes(x = value, group = Var2)) +
    geom_boxplot() + coord_flip()


boxplot(tmp)
qplot(tmp[,1], geom = "boxplot")
qplot(tmp[,2], geom = "boxplot")


## ####################################################################################

pietg2 <- function(formula, data, B = 100, bIni = NULL, rIni = NULL) {
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
    ## h <- .5
    ## Different ways to find root; choose the best root    
    ## Use beta from M-approach
    dat$Y <- rep(dat$Time[dat$status == 1], m + 1)
    T2 <- dat$Time[dat$status == 0]
    C2 <- dat$Y[dat$status == 0]
    Z <- as.matrix(subset(dat, status == 0)[,c("x1", "x2")])
    h <- 1.06 * sd(Z %*% c(.6, .8)) * nrow(Z)^-.2
    ht <- 1.06 * sd(T2) * length(T2)^-.2
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id)) 
    mm <- aggregate(event ~ id, dat0, sum)[, 2] 
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status, Y)))     
    X0 <- as.matrix(subset(dat, status == 1, select = -c(Time, id, m, event, status, Y)))    
    ## SSE1 <- ComputeSSE(X0, m)$alpha
    ## SSE2 <- ComputeSSE2(X0, m)$alpha
    ## c(SSE1, SSE2)
    ESE1 <- ComputeESE(X0, m)$alpha
    ESE2 <- ComputeESE2(X0, m)$alpha
    c(ESE1, ESE2)
}

do2 <- function(n, model, frailty = FALSE, type1 = FALSE, noCen = FALSE) {
    dat <- simDat(n, model, frailty, type1, noCen = noCen)
    fit <- pietg2(reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2,
                  data = dat)
    unlist(fit)
}

set.seed(1); do2(200, "M2", noCen = TRUE)
set.seed(1); do2(200, "M3", noCen = TRUE)

do2(200, "M1", noCen = TRUE)
do2(200, "M2", noCen = TRUE)
do2(200, "M3", noCen = TRUE)
f <- replicate(100, do2(200, "M3", noCen = TRUE))
summary(t(f))

M2n200noCen <- replicate(500, do2(200, "M2", noCen = TRUE))
M3n200noCen <- replicate(500, do2(200, "M3", noCen = TRUE))
M4n200noCen <- replicate(500, do2(200, "M4", noCen = TRUE))
M5n200noCen <- replicate(500, do2(200, "M5", noCen = TRUE))
M2n200cnoCen <- replicate(500, do2(200, "M2", TRUE, noCen = TRUE))
M3n200cnoCen <- replicate(500, do2(200, "M3", TRUE, noCen = TRUE))
M4n200cnoCen <- replicate(500, do2(200, "M4", TRUE, noCen = TRUE))
M5n200cnoCen <- replicate(500, do2(200, "M5", TRUE, noCen = TRUE))

optESE200 <- list(M2n200noCen = M2n200noCen, M3n200noCen = M3n200noCen,
                  M4n200noCen = M4n200noCen, M5n200noCen = M5n200noCen,
                  M2n200cnoCen = M2n200cnoCen, M3n200cnoCen = M3n200cnoCen,
                  M4n200cnoCen = M4n200cnoCen, M5n200cnoCen = M5n200cnoCen)
save(optESE200, file = "optESE200.RData")


M2n400noCen <- replicate(500, do2(400, "M2", noCen = TRUE))
M3n400noCen <- replicate(500, do2(400, "M3", noCen = TRUE))
M4n400noCen <- replicate(500, do2(400, "M4", noCen = TRUE))
M5n400noCen <- replicate(500, do2(400, "M5", noCen = TRUE))
M2n400cnoCen <- replicate(500, do2(400, "M2", TRUE, noCen = TRUE))
M3n400cnoCen <- replicate(500, do2(400, "M3", TRUE, noCen = TRUE))
M4n400cnoCen <- replicate(500, do2(400, "M4", TRUE, noCen = TRUE))
M5n400cnoCen <- replicate(500, do2(400, "M5", TRUE, noCen = TRUE))

optESE400 <- list(M2n400noCen = M2n400noCen, M3n400noCen = M3n400noCen,
                  M4n400noCen = M4n400noCen, M5n400noCen = M5n400noCen,
                  M2n400cnoCen = M2n400cnoCen, M3n400cnoCen = M3n400cnoCen,
                  M4n400cnoCen = M4n400cnoCen, M5n400cnoCen = M5n400cnoCen)
save(optESE400, file = "optESE400.RData")


do2(200, "M3")
foo <- replicate(500, do2(200, "M3"))
rowMeans(foo)


#######################################################################
## Load package
#######################################################################

library(survival)
library(SSIndex)
library(tidyverse)
library(Rcpp)

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
sourceCpp(code = '
  #include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]
  using namespace arma;
  // [[Rcpp::export]]
double getb2(double theta,
             double h,
             arma::vec tij,
             arma::vec yi,
             arma::mat X,
             arma::vec J,
             arma::vec W) {
    arma::vec b = zeros<vec>(2);
    int n = yi.n_elem;
    double out = 0;
    b[0] = sin(theta);
    b[1] = cos(theta);
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
        if(yi[j] >= tij[i] && tij[j] <= tij[i]) {
          double dx = (xb[j] - xb[i]) / h;
          Kh = 0;
          if (dx < 1 && dx > -1) Kh = 0.75 * (1 - dx * dx) / h;
          de += Kh;
          nu += Kh * xJ[j];
        }
      }
      out += W[i] * (xJ[i] - nu / de);
    // Rcpp::Rcout << nu / de << std::endl;
    }
    return(out);
}')

sourceCpp(code = '
  #include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]
  using namespace arma;
  // [[Rcpp::export]]
double getb(double theta,
            double h,
            arma::vec tij,
            arma::vec yi,
            arma::mat X,
            arma::vec W) {
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
        if(yi[j] >= tij[i] && tij[j] <= tij[i]) {
          double dx = (xb[j] - xb[i]) / h;
          Kh = 0;
          if (dx < 1 && dx > -1) Kh = 0.75 * (1 - dx * dx) / h;
          de += Kh; 
          nu += Kh * xJ[j];
        }
      }
      out += W[i] * (xJ[i] - nu / de);
      // Rcpp::Rcout << xJ[i] - nu / de << std::endl;
    }
    return(out);
}')

betaEst <- function(formula, data) {
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
    y <- subset(dat, event == 0)$Time
    W <- sapply(tij, function(x) sum(x <= y)) / length(y)
    ## s <- survfit(Surv(y, y < max(y)) ~ 1)
    ## W <- approx(s$time, s$surv, tij, yleft = 1, yright = min(s$surv))$y
    ## h <- 1.06 * sd(X %*% c(.6, .8)) * nrow(X)^-.2
    h <- .5
    ## Different ways to find root; choose the best root
    b1 <- BB::spg(par = 0, fn = function(x)
        getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], W = W)^2,
        lower = 0, upper = pi, quiet = TRUE, alertConvergence = FALSE)$par
    b2 <- optimize(f = function(x)
        getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], W = W)^2,
        interval = c(0, pi))$minimum
    gg <- seq(0, pi, length.out = 21)
    g0 <- sapply(gg, function(x)
        getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], W = W))
    b3 <- sapply(which(diff(sign(g0)) != 0), function(x)
        uniroot(f = function(x) getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)],
                                     X = X[rep(1:n, mm),], W = W),
                interval = c(gg[x], gg[x + 1]))$root)
    b <- c(b1, b2, b3)
    b <- b[which.min(sapply(b, function(x)
        abs(getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)],
                 X = X[rep(1:n, mm),], W = W))))]
    b
}

betaEst2 <- function(formula, data) {
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
    h <- 1.06 * sd(X %*% c(.6, .8)) * nrow(X)^-.2
    y <- subset(dat, event == 0)$Time
    W <- sapply(tij, function(x) sum(x <= y)) / length(y)
    ## h <- .5
    ## Different ways to find root; choose the best root
    b1 <- BB::spg(par = 0, fn = function(x)
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(1, 0), W = W)^2 +
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(0, 1), W = W)^2, 
        lower = 0, upper = pi, quiet = TRUE, alertConvergence = FALSE)$par
    b2 <- optimize(f = function(x)
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(1, 0), W = W)^2 +
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(0, 1), W = W)^2,
        interval = c(0, pi))$minimum
    b <- c(b1, b2)
    keep <- which.min(sapply(b, function(x)
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(1, 0), W = W)^2 +
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(0, 1), W = W)^2))
    b[keep]
}

eeCurve <- function(formula, data) {
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
    y <- subset(dat, event == 0)$Time
    W <- sapply(tij, function(x) sum(x <= y)) / length(y)    
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status))) ## 80 by 2
    h <- 1.06 * sd(X %*% c(.6, .8)) * nrow(X)^-.2
    ## h <- .5
    t0 <- seq(0, 2 * pi, length.out = 1000)
    e0 <- sapply(t0, function(x)
        getb(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], W = W))
    e0
}

eeCurve2 <- function(formula, data) {
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
    y <- subset(dat, event == 0)$Time
    W <- sapply(tij, function(x) sum(x <= y)) / length(y)
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status))) ## 80 by 2
    h <- 1.06 * sd(X %*% c(.6, .8)) * nrow(X)^-.2
    ## h <- .5
    t0 <- seq(0, 2 * pi, length.out = 1000)
    e0 <- sapply(t0, function(x)
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(1, 0), W = W)^2 +
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(0, 1), W = W)^2)
    e1 <- sapply(t0, function(x)
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(1, 0), W = W))
    e2 <- sapply(t0, function(x)
        getb2(x, h = h, tij = tij, yi = yi[rep(1:n, mm)], X = X[rep(1:n, mm),], J = c(0, 1), W = W))
    list(e0 = e0, e1 = e1, e2 = e2)
}

do <- function(n, model, frailty = FALSE, type1 = FALSE) {
    dat <- simDat(n, model, frailty, type1)
    fm <- reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2
    b1 <- betaEst(fm, data = dat)
    m <- dat$m[dat$status == 1]
    Y <- dat$Time[dat$status == 1]
    dat$Y <- rep(Y, m+1)
    T2 <- dat$Time[dat$status == 0]
    C2 <- dat$Y[dat$status == 0]
    Z <- as.matrix(subset(dat, status == 0)[,c("x1", "x2")])
    b2 <- optimize(f = function(x) -M(c(sin(x), cos(x)), Z, T2, C2, 0.2, 0.5),
                   interval = c(0, pi))$minimum
    c(b1, b2)
}

doEE <- function(n, model, frailty = FALSE, type1 = FALSE) {
    dat <- simDat(n, model, frailty, type1)
    fm <- reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2
    b1 <- betaEst(fm, data = dat)
    b2 <- betaEst2(fm, data = dat)
    c(sin(b1), cos(b1), sin(b2), cos(b2))
}

set.seed(1); doEE(200, "M2", TRUE) # good
set.seed(10); doEE(200, "M2", TRUE) # bad
set.seed(55); doEE(200, "M2", TRUE) # bad

set.seed(10)
set.seed(1)
dat <- simDat(200, "M2", TRUE, FALSE)
fm <- reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2
betaEst(fm, data = dat)
betaEst2(fm, data = dat)
t0 <- seq(0, 2 * pi, length.out = 1000)
e1 <- eeCurve(fm, data = dat)
e2 <- eeCurve2(fm, data = dat)

## pdf("normal.pdf", width = 20, height = 10)
## pdf("bad.pdf", width = 20, height = 10)
par(mfrow = c(2, 2))
plot(t0, e1 / 200, 'l', main = "Phi is Jacobian; vertical line is true value")
abline(h = 0)
abline(v = asin(.6), col = 2)
plot(t0, e2$e0 / 200, 'l', main = "Phi is (x1, x2)")
abline(v = asin(.6), col = 2)
plot(t0, e2$e1 / 200, 'l', ylim = c(min(e2$e1, e2$e2) / 200, max(e2$e1, e2$e2) / 200))
lines(t0, e2$e2 / 200, col = 2)
legend(5, -.1, c("EE1", "EE2"), lty = 1, col = 1:2)
hist(with(dat, unlist(lapply(split(m, id), unique))), breaks = 20)
## dev.off()

f <- matrix(NA, 100, 4)
for (i in 1:100) {
    set.seed(i) 
    f[i,] <- doEE(200, "M2", TRUE)
}
colMeans(f)

doM <- function(n, model, frailty = FALSE, type1 = FALSE) {
    dat <- simDat(n, model, frailty, type1)
    m <- dat$m[dat$status == 1]
    Y <- dat$Time[dat$status == 1]
    dat$Y <- rep(Y, m+1)
    T2 <- dat$Time[dat$status == 0]
    C2 <- dat$Y[dat$status == 0]
    Z <- as.matrix(subset(dat, status == 0)[,c("x1", "x2")])
    b2 <- optimize(f = function(x) -M(c(sin(x), cos(x)), Z, T2, C2, 0.2, 0.5),
                   interval = c(0, pi))$minimum
    c(sin(b2), cos(b2))
}

f <- replicate(200, tryCatch(do(200, "M2", TRUE), error = function(e)
    tryCatch(do(200, "M2", TRUE), error = function(ee)
        tryCatch(do(200, "M2", TRUE), error = function(eee) do(200, "M2", TRUE))))) 

c(mean(sin(f[1,])), mean(cos(f[1,])))
c(mean(sin(f[2,])), mean(cos(f[2,])))
c(sd(sin(f[1,])), sd(cos(f[1,])))
c(sd(sin(f[2,])), sd(cos(f[2,])))

rowMeans(f)
summary(t(f))
asin(.6)

## #################################################################################
## M approach only
## #################################################################################

M2n200 <- replicate(500, doM(200, "M2", FALSE))
M2n200c <- replicate(500, doM(200, "M2", TRUE))
M3n200 <- replicate(500, doM(200, "M3", FALSE))
M3n200c <- replicate(500, doM(200, "M3", TRUE))
M4n200 <- replicate(500, doM(200, "M4", FALSE))
M4n200c <- replicate(500, doM(200, "M4", TRUE))
M5n200 <- replicate(500, doM(200, "M5", FALSE))
M5n200c <- replicate(500, doM(200, "M5", TRUE))
n200M <- list(M2n200 = M2n200, M2n200c = M2n200c,
              M3n200 = M3n200, M3n200c = M3n200c,
              M4n200 = M4n200, M4n200c = M4n200c,
              M5n200 = M5n200, M5n200c = M5n200c)
## save(n200M, file = "n200M.RData")

M2n400 <- replicate(500, doM(400, "M2", FALSE))
M2n400c <- replicate(500, doM(400, "M2", TRUE))
M3n400 <- replicate(500, doM(400, "M3", FALSE))
M3n400c <- replicate(500, doM(400, "M3", TRUE))
M4n400 <- replicate(500, doM(400, "M4", FALSE))
M4n400c <- replicate(500, doM(400, "M4", TRUE))
M5n400 <- replicate(500, doM(400, "M5", FALSE))
M5n400c <- replicate(500, doM(400, "M5", TRUE))
n400M <- list(M2n400 = M2n400, M2n400c = M2n400c,
              M3n400 = M3n400, M3n400c = M3n400c,
              M4n400 = M4n400, M4n400c = M4n400c,
              M5n400 = M5n400, M5n400c = M5n400c)
## save(n400M, file = "n400M.RData")

load("n200M.RData")
load("n400M.RData")

n200 <- lapply(n200M, function(x) cbind(abs(rowMeans(x) - 3:4 / 5), apply(x, 1, sd)))
tab200 <- rbind(cbind(n200[[1]], n200[[2]]),
                cbind(n200[[3]], n200[[4]]),
                cbind(n200[[7]], n200[[8]]),
                cbind(n200[[5]], n200[[6]]))

## n400 <- lapply(n400M, function(x) cbind(rowMeans(x), apply(x, 1, sd)))
n400 <- lapply(n400M, function(x) cbind(abs(rowMeans(x) - 3:4 / 5), apply(x, 1, sd)))
tab400 <- rbind(cbind(n400[[1]], n400[[2]]),
                cbind(n400[[3]], n400[[4]]),
                cbind(n400[[7]], n400[[8]]),
                cbind(n400[[5]], n400[[6]]))



## #################################################################################
## M approach only
## #################################################################################
set.seed(1); dat1 <- simDat(200, "M2", TRUE, FALSE)
set.seed(10); dat2 <- simDat(200, "M2", TRUE, FALSE)

with(subset(dat1, event == 0), table(x1 > 0, x2 > 0))
with(subset(dat2, event == 0), table(x1 > 0, x2 > 0))
summary(subset(dat1, event == 0)$x1)
summary(subset(dat1, event == 0)$x2)
summary(subset(dat2, event == 0)$x1)
summary(subset(dat2, event == 0)$x2)

betaEst(fm, data = dat1)
betaEst(fm, data = dat2)

dat3 <- dat2
dat3$x1 <- dat3$x1 / max(abs(dat3$x1))
dat3$x2 <- dat3$x2 / max(abs(dat3$x2))
betaEst(fm, data = dat3)

summary(subset(dat3, event == 0)$x1)
summary(subset(dat3, event == 0)$x2)

dat4 <- dat2
dat4$x1 <- with(dat4, (x1 - min(x1)) / (max(x1) - min(x1)))
dat4$x2 <- with(dat4, (x2 - min(x2)) / (max(x2) - min(x2)))
summary(subset(dat4, event == 0)$x1)
summary(subset(dat4, event == 0)$x2)
betaEst(fm, data = dat4)
                
with(dat1, Recur(Time, id, event))
with(dat1, Recur(Time, id, event, status))

library(reReg)
library(gridExtra)

grid.arrange(plotEvents(Recur(Time, id, event) ~ 1, data = dat1,
                        main = "Good example", legend.positive = "none"),
             plotEvents(Recur(Time, id, event) ~ 1, data = dat2,
                        main = "Bad example", legend.positive = "none"),
             ncol = 2)
## ggsave("eventPlots.pdf")

set.seed(55); dat <- simDat(200, "M2", TRUE, FALSE)
plotEvents(Recur(Time, id, event) ~ 1, data = dat, main = "", legend.positive = "none")


set.seed(10)
dat <- simDat(200, "M2", TRUE, FALSE)
fm <- reSurv(time1 = Time, id = id, event = event, status = status) ~ 
    x1 + x2
betaEst(fm, data = dat)
betaEst(fm, subset(dat, id %in% which(subset(dat, status == 1)$Time == 10)))
betaEst(fm, subset(dat, id %in% which(subset(dat, status == 1)$Time < 10)))

ss <- sapply(1:200, function(x) betaEst(fm, subset(dat, id <= x)))
asin(.6)
betaEst(fm, subset(dat, id <= 49))
betaEst(fm, subset(dat, id <= 50))
subset(dat, id == 50)

betaEst(fm, subset(dat, !(id %in% c(50:75))))

ss <- sapply(1:200, function(x) betaEst(fm, subset(dat, !(id %in% x))))

plot(1:200, ss, 'l')

dim(dat)
dim(subset(dat, id %in% which(subset(dat, status == 1)$Time == 10)))
dim(subset(dat, id %in% which(subset(dat, status == 1)$Time < 10)))



bb <- matrix(NA, 500, 4)
for (i in 1:500) {
    set.seed(i)
    bb[i,] <- doEE(100, "M2", TRUE)
}

set.seed(1); doEE(100, "M2", TRUE)
set.seed(1); dim(simDat(100, "M2", TRUE, FALSE))


doEE(50, "M2", TRUE)

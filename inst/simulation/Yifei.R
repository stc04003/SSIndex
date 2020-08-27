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
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[, 2] 
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time ## like T2
    yi <- subset(dat0, event == 0)$Time ## like C2
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status))) 
    y <- subset(dat, event == 0)$Time
    W <- sapply(tij, function(x) sum(x <= y)) / length(y)
    h <- 1.06 * sd(X %*% c(.6, .8)) * nrow(X)^-.2
    ## h <- .5
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

doEE <- function(n, model, frailty = FALSE, type1 = FALSE) {
    dat <- simDat(n, model, frailty, type1)
    fm <- reSurv(time1 = Time, id = id, event = event, status = status) ~ x1 + x2
    b1 <- betaEst(fm, data = dat)
    c(sin(b1), cos(b1))
}

## The "M2" in code is the "M1" in paper
set.seed(1); doEE(200, "M2", TRUE) # good
set.seed(10); doEE(200, "M2", TRUE) # bad

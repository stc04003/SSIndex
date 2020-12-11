#' Main estimation function
#'
#' Fits general rate model. The formula takes `reSurv` from `reReg`.
#'
#' @param formula a formula object, with the response on the left of a '~' operator.
#' @param data an optional argument for data.
#' @param shp.ind a character string indicating either to proceed with shape independence assumption. 
## #' @param B a numerical value specifying number of bootstrap in shape independence testing.
#' @param bIni initial value for estimating beta
#' @param rIni initial value for estimating gamma
#' @param opt how many optimization methods to try?
#' 1: optim
#' 2: 1 + spg
#' Default: 1 when p = 1, 2 when p > 1
#' @param smoothFirst estimate the smooth estimator first and use it to guide?
#' Default: FALSE when p = 1,  TRUE when p > 1
#' 
#' @importFrom BB spg
#' @importFrom tibble add_column
#' @importFrom dplyr %>%
#' 
#' @useDynLib SSIndex, .registration = TRUE
#' @export
#' 
ssfit <- function(formula, data, shp.ind = FALSE, bIni = NULL, rIni = NULL,
                optMethod = NULL, smoothFirst = NULL) {
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
    ## The first 5 columns of `DF` are fixed at id, time, event, status, m, followed by covariates.
    ## id Time event status m ...covariates...
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    dat <- DF %>% add_column(m = rep(m, m + 1), .after = 4)
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[,2]
    tij <- subset(dat, event == 1)$Time
    yi <- subset(dat, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = -c(Time, id, m, event, status)))
    p <- ncol(X)
    if (is.null(smoothFirst)) smoothFirst <- ifelse(p > 2, TRUE, FALSE)
    ## assuming data is generated from simDat
    ## estimate \beta first
    if (is.logical(shp.ind) && !shp.ind) {
        tmp <- getb0(dat, bIni, optMethod, smoothFirst)
        bhat1 <- tmp$bhat1
        bhat2 <- tmp$bhat2
    }
    if (is.logical(shp.ind) && shp.ind){
        bhat1 <- bhat2 <- double(2)
    }
    ## which bhat to use
    bhat <- bhat2
    ## The estimating equation Sn needs Yi even for the m = 0's
    xb <- X %*% bhat
    ## h <- 1.06 * sd(xb) * n^-.2
    h <- 2.78 * sd(xb) * n^(-1/3)
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
           as.double(tij), as.double(yi),
           as.double(xb), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "SSIndex")$result,
        X %*% bhat, yi))
    Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
    Fhat <- exp(-Fhat)
    ## #####################################################################################
    ## ## F(t, a) with a = \bar{X} %*% \beta
    ## Fhat0 <- unlist(mapply(FUN = function(x, y)
    ##     .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
    ##        as.double(xb), as.double(x), as.double(y), as.double(h), 
    ##        result = double(1), PACKAGE = "SSIndex")$result,
    ##     t(replicate(nrow(X), colMeans(X))) %*% bhat, yi))
    ## #####################################################################################
    ## F(t, a), with a = X %*% bhat,
    ## ## outputs a list with yi being each yi
    ## Fhat0 <- NULL
    ## for (i in 1:length(yi)) {
    ##     tmp <- unlist(mapply(FUN = function(x, y)
    ##         .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
    ##            as.double(tij), as.double(yi),
    ##            as.double(xb), as.double(x), as.double(y), as.double(h), 
    ##            result = double(1), PACKAGE = "SSIndex")$result,
    ##         xb, rep(yi[i], nrow(xb))))
    ##     tmp <- ifelse(is.na(tmp), 0, tmp) ## assign 0/0, Inf/Inf to 0
    ##     Fhat0[[i]] <- exp(-tmp)
    ## }
    Sn <- function(r) {
        r <- cumprod(c(1, sin(r))) * c(cos(r), 1)
        xr <- X %*% r
        ## -.C("shapeEq", as.integer(n), as.double(xr), as.double(mm / Fhat),
        ##     result = double(1), PACKAGE = "SSIndex")$result
        -.C("sizeEq2", as.integer(n), as.double(xr), as.double(mm / Fhat), as.double(yi),
            result = double(1), PACKAGE = "SSIndex")$result
    }
    Sn2 <- function(r) {
        r <- cumprod(c(1, sin(r))) * c(cos(r), 1)
        xr <- X %*% r
        -.C("sizeEqSmooth", as.integer(n), as.integer(p), as.double(xr), as.double(X),
            as.double(diag(p)), as.double(mm / Fhat),
            result = double(1), PACKAGE = "SSIndex")$result        
    }
    if (is.null(rIni)) rIni <- rep(1 / sqrt(p), p)
    if (smoothFirst) {
        rhat1 <- optSolver(fn = Sn2, par = rIni, optMethod = optMethod)
        rhat2 <- optSolver(fn = Sn, par = rhat1, optMethod = optMethod)
    } else {
        rhat1 <- NULL
        rhat2 <- optSolver(fn = Sn, par = rIni, optMethod = optMethod)
    }
    info <- list(mm = mm, midx = midx, tij = tij, yi = yi,
                 xb = as.numeric(xb), xg = as.numeric(X %*% rhat2), h = h)
    ## b0 and r0 are smooth estimator
    ## b00 and r00 are nonsmooth
    list(b0 = -1 * bhat1, r0 = rhat1, b00 = -1 * bhat2, r00 = rhat2, Fhat = Fhat) ## , info = info)
}


#' Function to get beta_0 estiamte
#' @noRd
getb0 <- function(dat, bIni, optMethod, smoothFirst) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status)))
    Cn <- function(b) {
        b <- cumprod(c(1, sin(b))) * c(cos(b), 1)
        -.C("rank", as.integer(n), as.integer(mm), as.integer(midx),
            as.double(tij), as.double(yi), as.double(X %*% b),
            result = double(1), PACKAGE = "SSIndex")$result
    }
    p <- ncol(X)
    Cn2 <- function(b) {
        b <- cumprod(c(1, sin(b))) * c(cos(b), 1)
        -.C("rankSmooth", as.integer(n), as.integer(p), as.integer(mm), as.integer(midx),
            as.double(solve(t(X) %*% X)), ## as.double(diag(p)), 
            as.double(tij), as.double(yi), as.double(X %*% b), as.double(X), 
            result = double(1), PACKAGE = "SSIndex")$result
    }
    dCn2 <- function(b) {
        b <- cumprod(c(1, sin(b))) * c(cos(b), 1)
        .C("drankSmooth", as.integer(n), as.integer(p), as.integer(mm), as.integer(midx),
            as.double(diag(p)), 
            as.double(tij), as.double(yi), as.double(X %*% b), as.double(X), 
            result = double(1), PACKAGE = "SSIndex")$result
    }
    if (is.null(bIni)) bIni <- rep(1 / sqrt(p), p)
    if (smoothFirst) {
        bhat1 <- optSolver(fn = Cn2, par = bIni, optMethod = optMethod)
        bhat2 <- optSolver(fn = Cn, par = bhat1, optMethod = optMethod)
    } else {
        bhat1 <- NULL
        bhat2 <- optSolver(fn = Cn, par = bIni, optMethod = optMethod)
    }
    list(bhat1 = bhat1, bhat2 = bhat2)
}

#' Function to get kappa for testing H0: beta0 = 0
#' @importFrom dplyr filter
#' @noRd
getk0 <- function(dat, b) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(dat0 %>% filter(event == 0) %>% select(starts_with("x")))
    ## X <- as.matrix(subset(dat0, event == 0, select = c(x1, x2)))
    .C("kappa", as.integer(n), as.integer(mm), as.integer(midx),
       as.double(tij), as.double(yi), as.double(X %*% b),
       result = double(1), PACKAGE = "SSIndex")$result
}

#' Function to get kappa for testing H0: beta0 = 0 & gamma0 = 0
#' @noRd
#'
#' @importFrom dplyr select
#' @importFrom tidyselect starts_with
getk02 <- function(dat, b, Fhat) {
    mm <- aggregate(event ~ id, dat, sum)[,2]
    dat0 <- subset(dat, event < 1)
    n <- nrow(dat0)
    rownames(dat0) <- NULL
    dat0$id <- 1:n
    X <- as.matrix(dat0 %>% select(starts_with("x")))
    .C("kappa2", as.integer(n), as.double(X %*% b), as.double(mm / Fhat),
       result = double(1), PACKAGE = "SSIndex")$result
}


#' Revised function to get kappa for testing H0: beta0 = 0 & gamma0 = 0
#' @noRd
getk03 <- function(dat, b) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status)))
    .C("kappa3", as.integer(n), as.integer(mm), as.integer(midx), 
       as.double(X %*% b), as.double(tij), as.double(yi),
       result = double(1), PACKAGE = "SSIndex")$result
}



#' Revised function to get kappa for testing H0: beta0 = 0 & gamma0 = 0;
#' give both shape and rate test statistics
#' @noRd
getk04 <- function(dat, b) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = -c(Time, id, m, event, status)))
    .C("kappa4", as.integer(n), as.integer(mm), as.integer(midx), 
       as.double(X %*% b), as.double(tij), as.double(yi),
       result = double(2), PACKAGE = "SSIndex")$result
}

#' Function to compute Fhat given data and beta
getFhat <- function(formula, data, b) {
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
    ## The first 5 columns of `DF` are fixed at id, time, event, status, m, followed by covariates.
    ## id Time event status m ...covariates...
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    dat <- DF %>% add_column(m = rep(m, m + 1), .after = 4)
    ## assuming data is generated from simDat
    ## estimate \beta first
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[,2]
    tij <- subset(dat, event == 1)$Time
    yi <- subset(dat, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = -c(Time, id, m, event, status)))
    p <- ncol(X)
    d <- dstar <- NULL
    if (ncol(X) != length(b)) stop("The length of beta does not equal to the number of covariate")
    xb <- X %*% b
    ## h <- 1.06 * sd(xb) * n^-.2
    h <- 2.78 * sd(xb) * n^(-1/3)
    t0 <- seq(0, max(yi), length.out = 500)
    ## f <- sapply(1:length(xb), function(i)
    ##     unlist(mapply(FUN = function(x, y)
    ##         .C("shapeFun", as.integer(n), 
    ##            as.integer(mm), as.integer(midx), as.double(tij), as.double(yi), 
    ##            as.double(xb), as.double(x), as.double(y), as.double(h), 
    ##            result = double(1), PACKAGE = "SSIndex")$result, xb[i], t0)))    
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
           as.double(tij), as.double(yi),
           as.double(xb), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "SSIndex")$result,
        xb, t0))
    Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
    Fhat <- exp(-Fhat)
    Shat <- 1 - Fhat
    mu <- diff(c(0, t0)) * Shat
    return(list(Time = t0, Fhat = Fhat, mu = mu))
}

## To minimize
optSolver <- function(fn, par, 
                      interval = c(0, 2 * pi),
                      optMethod = NULL, ...) {
    p <- length(par)
    if (is.null(optMethod)) optMethod <- ifelse(p <= 2, 1, 3)
    r0 <- trans(par)
    tmp1 <- tmp2 <- NULL
    if (length(r0) <= 1) {
        tmp1 <- optimize(f = function(z) fn(z, ...), interval = interval)
        res <- tmp1$minimum %% (2 * pi)
    } else {
        tmp1 <- optim(par = r0 %% (2 * pi), fn = function(z) fn(z, ...),
                      control = list(warn.1d.NelderMead = FALSE))
        res <- tmp1$par %% (2 * pi)
    }
    if (optMethod >= 2) {
        tmp2 <- spg(par = r0 %% (2 * pi), fn = function(z) fn(z, ...),
                    quiet = TRUE, control = list(trace = FALSE))
        if (c(tmp1$objective, tmp1$value) > tmp2$value) res <- tmp2$par %% (2 * pi)
    }
    cumprod(c(1, sin(res))) * c(cos(res), 1)
}

trans <- function(beta) {
  j <- stheta <- 1
  lb <- length(beta)
  theta <- rep(0, lb-1)
  while(j < lb) {
    theta[j] <- acos(beta[j] / stheta[j])
    stheta[j+1] <- stheta[j] * sin(theta[j])
    j <- j + 1
  }
  if(sin(theta[j-1]) * (beta[j]/stheta[j])<0){ 
    theta[j-1] <- -abs(theta[j-1])
  }
  theta
}

#' For old codes
#' @export
#' @noRd
gsm <- function(formula, data, shp.ind = FALSE, B = 100, bIni = NULL, rIni = NULL,
                optMethod = NULL, smoothFirst = NULL) {
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
    ## The first 5 columns of `DF` are fixed at id, time, event, status, m, followed by covariates.
    ## id Time event status m ...covariates...
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    dat <- DF %>% add_column(m = rep(m, m + 1), .after = 4)
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[,2]
    tij <- subset(dat, event == 1)$Time
    yi <- subset(dat, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = -c(Time, id, m, event, status)))
    p <- ncol(X)
    if (is.null(smoothFirst)) smoothFirst <- ifelse(p > 2, TRUE, FALSE)
    ## assuming data is generated from simDat
    ## estimate \beta first
    if (is.logical(shp.ind) && !shp.ind) {
        tmp <- getb0(dat, bIni, optMethod, smoothFirst)
        bhat1 <- tmp$bhat1
        bhat2 <- tmp$bhat2
    }
    if (is.logical(shp.ind) && shp.ind){
        bhat1 <- bhat2 <- double(2)
    }
    d <- dstar <- NULL
    if (!is.logical(shp.ind) && shp.ind == "test") {
        ind <- sample(1:n)[1:round(n/2)]
        dat1 <- subset(dat, id %in% ind)
        dat2 <- subset(dat, !(id %in% ind))
        tilde.b <- getb0(dat1, bIni, optMethod, smoothFirst)$bhat1
        n2 <- length(unique(dat2$id))
        mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
        d <- getd(dat2, tilde.b)
        dstar <- replicate(B, boot.d(dat2, tilde.b))
        if (abs(d / sd(dstar)) > qnorm(.975)) {
            tmp <- getb0(dat, bIni, optMethod, smoothFirst)
            bhat1 <- tmp$bhat1
            bhat2 <- tmp$bhat2
        }
        else bhat1 <- bhat2 <- double(p)
    }
    ## which bhat to use
    bhat <- bhat2
    ## The estimating equation Sn needs Yi even for the m = 0's
    xb <- X %*% bhat
    ## h <- 1.06 * sd(xb) * n^-.2
    h <- 2.78 * sd(xb) * n^(-1/3)
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
           as.double(tij), as.double(yi),
           as.double(xb), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "SSIndex")$result,
        X %*% bhat, yi))
    Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
    Fhat <- exp(-Fhat)
    ## F(t, a) with a = \bar{X} %*% \beta
    Fhat0 <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", 
           as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), 
           as.double(yi), as.double(X %*% double(p)), as.double(x), 
           as.double(y), as.double(h), result = double(1), PACKAGE = "SSIndex")$result, 
        X %*% double(p), yi))
    Fhat0 <- ifelse(is.na(Fhat0), 0, Fhat0)
    Fhat0 <- exp(-Fhat0)
    ## #####################################################################################
    ## ## F(t, a) with a = \bar{X} %*% \beta
    ## Fhat0 <- unlist(mapply(FUN = function(x, y)
    ##     .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
    ##        as.double(xb), as.double(x), as.double(y), as.double(h), 
    ##        result = double(1), PACKAGE = "SSIndex")$result,
    ##     t(replicate(nrow(X), colMeans(X))) %*% bhat, yi))
    ## #####################################################################################
    ## F(t, a), with a = X %*% bhat,
    ## ## outputs a list with yi being each yi
    ## Fhat0 <- NULL
    ## for (i in 1:length(yi)) {
    ##     tmp <- unlist(mapply(FUN = function(x, y)
    ##         .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx),
    ##            as.double(tij), as.double(yi),
    ##            as.double(xb), as.double(x), as.double(y), as.double(h), 
    ##            result = double(1), PACKAGE = "SSIndex")$result,
    ##         xb, rep(yi[i], nrow(xb))))
    ##     tmp <- ifelse(is.na(tmp), 0, tmp) ## assign 0/0, Inf/Inf to 0
    ##     Fhat0[[i]] <- exp(-tmp)
    ## }
    Sn <- function(r) {
        r <- cumprod(c(1, sin(r))) * c(cos(r), 1)
        xr <- X %*% r
        ## -.C("shapeEq", as.integer(n), as.double(xr), as.double(mm / Fhat),
        ##     result = double(1), PACKAGE = "SSIndex")$result
        -.C("sizeEq2", as.integer(n), as.double(xr), as.double(mm / Fhat), as.double(yi),
            result = double(1), PACKAGE = "SSIndex")$result
    }
    Sn2 <- function(r) {
        r <- cumprod(c(1, sin(r))) * c(cos(r), 1)
        xr <- X %*% r
        -.C("sizeEqSmooth", as.integer(n), as.integer(p), as.double(xr), as.double(X),
            as.double(diag(p)), as.double(mm / Fhat),
            result = double(1), PACKAGE = "SSIndex")$result        
    }
    if (is.null(rIni)) rIni <- rep(1 / sqrt(p), p)
    if (smoothFirst) {
        rhat1 <- optSolver(fn = Sn2, par = rIni, optMethod = optMethod)
        rhat2 <- optSolver(fn = Sn, par = rhat1, optMethod = optMethod)
    } else {
        rhat1 <- NULL
        rhat2 <- optSolver(fn = Sn, par = rIni, optMethod = optMethod)
    }
    info <- list(mm = mm, midx = midx, tij = tij, yi = yi,
                 xb = as.numeric(xb), xg = as.numeric(X %*% rhat2), h = h)
    ## b0 and r0 are smooth estimator
    ## b00 and r00 are nonsmooth
    list(b0 = -1 * bhat1, r0 = rhat1, b00 = -1 * bhat2, r00 = rhat2,
         d = d, dstar = dstar,
         Fhat = Fhat, Fhat0 = Fhat0) ## , info = info)
}

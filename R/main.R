#' Main estimation function
#'
#' Fits general rate model. The formula takes `reSurv` from `reReg`.
#'
#' @param formula a formula object, with the response on the left of a '~' operator.
#' @param data an optional argument for data.
#' @param shp.ind a character string indicating either to proceed with shape independence assumption. 
#' @param B a numerical value specifying number of bootstrap in shape independence testing.
#' 
#' @importFrom BB spg
#' @importFrom reReg reSurv is.reSurv
#' @importFrom tibble add_column
#' @importFrom dplyr %>%
#' 
#' @useDynLib GSM, .registration = TRUE
#' @export
#' 

gsm <- function(formula, data, shp.ind = FALSE, B = 100) {
    ## Extract vectors
    Call <- match.call()
    if (missing(data)) obj <- eval(formula[[2]], parent.frame()) 
    if (!missing(data)) obj <- eval(formula[[2]], data) 
    if (!is.reSurv(obj)) stop("Response must be a reSurv object")
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
    if (is.logical(shp.ind) && !shp.ind) {
        tmp <- getb0(dat)
        bhat1 <- tmp$bhat1
        bhat2 <- tmp$bhat2
    }
    if (is.logical(shp.ind) && shp.ind){
        bhat1 <- bhat2 <- double(2)
    }
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[,2]
    tij <- subset(dat, event == 1)$Time
    yi <- subset(dat, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = -c(Time, id, m, event, status)))
    p <- ncol(X)
    d <- dstar <- NULL
    if (!is.logical(shp.ind) && shp.ind == "test") {
        ind <- sample(1:n)[1:round(n/2)]
        dat1 <- subset(dat, id %in% ind)
        dat2 <- subset(dat, !(id %in% ind))
        tilde.b <- getb0(dat1)$bhat1
        n2 <- length(unique(dat2$id))
        mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
        d <- getd(dat2, tilde.b)
        dstar <- replicate(B, boot.d(dat2, tilde.b))
        if (abs(d / sd(dstar)) > qnorm(.975)) {
            tmp <- getb0(dat)
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
    h <- 2.78 * sd(xb) * n^-.2
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
           as.double(xb), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "GSM")$result,
        X %*% bhat, yi))
    Fhat <- ifelse(is.na(Fhat), 0, Fhat) ## assign 0/0, Inf/Inf to 0
    Fhat <- exp(-Fhat)
    Fhat0 <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
           as.double(X %*% double(p)), as.double(x), as.double(y), as.double(h), 
           result = double(1), PACKAGE = "GSM")$result,
        X %*% double(p), yi))    
    Fhat0 <- ifelse(is.na(Fhat0), 0, Fhat0) ## assign 0/0, Inf/Inf to 0
    Fhat0 <- exp(-Fhat0)
    Sn <- function(r) {
        r <- cumprod(c(1, sin(r))) * c(cos(r), 1)
        xr <- X %*% r
        -.C("shapeEq", as.integer(n), as.double(xr), as.double(mm / Fhat),
            result = double(1), PACKAGE = "GSM")$result
    }
    Sn2 <- function(r) {
        r <- cumprod(c(1, sin(r))) * c(cos(r), 1)
        xr <- X %*% r
        -.C("shapeEqSmooth", as.integer(n), as.integer(p), as.double(xr), as.double(X),
            as.double(diag(p)), as.double(mm / Fhat),
            result = double(1), PACKAGE = "GSM")$result        
    }
    if (p <= 2) {
        tmp1 <- spg(par = acos(1 / sqrt(p)), fn = Sn2, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optimize(f = Sn2, interval = c(-10, 10))
        if (tmp1$value < tmp2$objective) rhat1 <- tmp1$par %% (2 * pi)
        else rhat1 <- tmp2$minimum %% (2 * pi)
        tmp1 <- spg(par = rhat1, fn = Sn, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optimize(f = Sn, interval = c(rhat1 - pi / 2, rhat1 + pi / 2))
        if (tmp1$value < tmp2$objective) rhat2 <- tmp1$par %% (2 * pi)
        else rhat2 <- tmp2$minimum %% (2 * pi)
    }
    if (p > 2) {
        tmp1 <- spg(par = rep(acos(1 / sqrt(p)), p - 1),
                    fn = Sn2, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optim(par = rep(acos(1 / sqrt(p)), p - 1), fn = Sn2)
        if (tmp1$value < tmp2$value) rhat1 <- tmp1$par %% (2 * pi)
        else rhat1 <- tmp2$par %% (2 * pi)
        tmp1 <- spg(par = rhat1, fn = Sn, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optim(par = rhat1, fn = Sn)
        if (tmp1$value < tmp2$value) rhat2 <- tmp1$par %% (2 * pi)
        else rhat2 <- tmp2$par %% (2 * pi)        
    }
    rhat1 <- cumprod(c(1, sin(rhat1))) * c(cos(rhat1), 1)
    rhat2 <- cumprod(c(1, sin(rhat2))) * c(cos(rhat2), 1)
    list(b0 = bhat2, r0 = rhat1, b00 = bhat1, r00 = rhat2, d = d, dstar = dstar, Fhat = Fhat, Fhat0 = Fhat0)
}

#' Function to get beta_0 estiamte
#' @noRd
#' @export
getb0 <- function(dat) {
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
            result = double(1), PACKAGE = "GSM")$result
    }
    p <- ncol(X)
    Cn2 <- function(b) {
        b <- cumprod(c(1, sin(b))) * c(cos(b), 1)
        -.C("rankSmooth", as.integer(n), as.integer(p), as.integer(mm), as.integer(midx),
            as.double(diag(p)), 
            as.double(tij), as.double(yi), as.double(X %*% b), as.double(X), 
            result = double(1), PACKAGE = "GSM")$result
    }
    ## which one gives the absolution min?
    if (p <= 2) {
        ## Solve the induced smoothing version first, then un-smoothed
        tmp1 <- spg(par = acos(1 / sqrt(p)), fn = Cn2, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optimize(f = Cn2, interval = c(-10, 10))
        if (tmp1$value < tmp2$objective) bhat1 <- tmp1$par %% (2 * pi)
        else bhat1 <- tmp2$minimum %% (2 * pi)
        ## bhat1 is smooth version, bhat2 is unsmooth version
        tmp1 <- spg(par = bhat1, fn = Cn, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optimize(f = Cn, interval = c(bhat1 - pi / 2, bhat1 + pi / 2))
        if (tmp1$value < tmp2$objective) bhat2 <- tmp1$par %% (2 * pi)
        else bhat2 <- tmp2$minimum %% (2 * pi)
    }
    if (p > 2) {
        tmp1 <- spg(par = rep(acos(1 / sqrt(p)), p - 1),
                    fn = Cn2, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optim(par = rep(acos(1 / sqrt(p)), p - 1), fn = Cn2)
        if (tmp1$value < tmp2$value) bhat1 <- tmp1$par %% (2 * pi)
        else bhat1 <- tmp2$par %% (2 * pi)
        tmp1 <- spg(par = bhat1, fn = Cn, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optim(par = bhat1, fn = Cn)
        if (tmp1$value < tmp2$value) bhat2 <- tmp1$par %% (2 * pi)
        else bhat2 <- tmp2$par %% (2 * pi)
    }
    bhat1 <- cumprod(c(1, sin(bhat1))) * c(cos(bhat1), 1)
    bhat2 <- cumprod(c(1, sin(bhat2))) * c(cos(bhat2), 1)
    list(bhat1 = bhat1, bhat2 = bhat2)
}

#' Function to get kappa for testing H0: beta0 = 0
#' @importFrom dplyr filter
#' @noRd
#' @export
getk0 <- function(dat, b) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$t
    yi <- subset(dat0, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(dat0 %>% filter(event == 0) %>% select(starts_with("x")))
    ## X <- as.matrix(subset(dat0, event == 0, select = c(x1, x2)))
    .C("kappa", as.integer(n), as.integer(mm), as.integer(midx),
       as.double(tij), as.double(yi), as.double(X %*% b),
       result = double(1), PACKAGE = "GSM")$result
}

#' Function to get kappa for testing H0: beta0 = 0 & gamma0 = 0
#' @noRd
#'
#' @importFrom dplyr select
#' @importFrom tidyselect starts_with
#' @export
getk02 <- function(dat, b, Fhat) {
    mm <- aggregate(event ~ id, dat, sum)[,2]
    dat0 <- subset(dat, event < 1)
    n <- nrow(dat0)
    rownames(dat0) <- NULL
    dat0$id <- 1:n
    X <- as.matrix(dat0 %>% select(starts_with("x")))
    ## subset(dat0, select = c(x1, x2)))
    .C("kappa2", as.integer(n), as.double(X %*% b), as.double(mm / Fhat),
       result = double(1), PACKAGE = "GSM")$result
}

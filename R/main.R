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
        bhat <- tmp$bhat
        bhat0 <- tmp$bhat0
    }
    if (is.logical(shp.ind) && shp.ind){
        bhat <- bhat0 <- double(2)
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
        tilde.b <- getb0(dat1)$bhat
        n2 <- length(unique(dat2$id))
        mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
        d <- getd(dat2, tilde.b)
        dstar <- replicate(B, boot.d(dat2, tilde.b))
        if (abs(d / sd(dstar)) > qnorm(.975)) {
            tmp <- getb0(dat)
            bhat <- tmp$bhat
            bhat0 <- tmp$bhat0
        }
        else bhat <- bhat0 <- double(2)
    }
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
    Sn <- function(r) {
        r <- cumprod(c(1, sin(r))) * c(cos(r), 1)
        xr <- X %*% r
        h <- 1.06 * sd(xr) * n^-.2 
        -.C("shapeEq", as.integer(n), as.double(xr), as.double(mm / Fhat),
            result = double(1), PACKAGE = "GSM")$result
    }
    if (p <= 2) {
        tmp1 <- spg(par = double(p - 1), fn = Sn, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optimize(f = Sn, interval = c(-10, 10))
        if (tmp1$value < tmp2$objective) rhat <- rhat0 <- tmp1$par
        else rhat <- rhat0 <- tmp2$minimum
    }
    if (p > 2) {
        tmp1 <- spg(par = double(2), fn = Sn, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optim(par = double(2), fn = Sn)
        if (tmp1$value < tmp2$value) rhat <- rhat0 <- tmp1$par
        else rhat <- rhat0 <- tmp2$par
    }
    rhat <- cumprod(c(1, sin(rhat))) * c(cos(rhat), 1)
    ## rhat <- rhat / sqrt(sum(rhat^2))
    ## rhat <- optimize(f = Sn, interval = c(-10, 10))$minimum
    list(b0 = bhat, r0 = rhat, b00 = bhat0, r00 = rhat0, d = d, dstar = dstar)
}

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
    ## which one gives the absolution min?
    p <- ncol(X)
    if (p <= 2) {
        tmp1 <- spg(par = double(p - 1), fn = Cn, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optimize(f = Cn, interval = c(-10, 10))
        if (tmp1$value < tmp2$objective) bhat <- bhat0 <- tmp1$par
        else bhat <- bhat0 <- tmp2$minimum
    }
    if (p > 2) {
        tmp1 <- spg(par = double(p - 1), fn = Cn, quiet = TRUE, control = list(trace = FALSE))
        tmp2 <- optim(par = double(p - 1), fn = Cn)
        if (tmp1$value < tmp2$value) bhat <- bhat0 <- tmp1$par
        else bhat <- bhat0 <- tmp2$par
    }
    bhat <- cumprod(c(1, sin(bhat))) * c(cos(bhat), 1)
    ## bhat <- bhat / sqrt(sum(bhat^2))
    list(bhat = bhat, bhat0 = bhat0)
}

#' @export
getk0 <- function(dat, b) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$t
    yi <- subset(dat0, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = c(x1, x2)))
    Cn <- function(b) {
        .C("kappa", as.integer(n), as.integer(mm), as.integer(midx),
            as.double(tij), as.double(yi), as.double(X %*% b),
            result = double(1), PACKAGE = "GSM")$result
    }
    Cn(b)
}

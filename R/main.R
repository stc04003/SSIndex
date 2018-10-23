#' Main estimation function
#'
#' This function assumes the input data is generated from \code{simDat}.
#' The suage can be generated to take formula form.
#'
#' @param dat data frame generated from \code{simDat}.
#' @param shp.ind a character string indicating either to proceed with shape independence assumption. 
#' @param B a numerical value specifying number of bootstrap in shape independence testing.
#' 
#' @importFrom BB spg
#'
#' @useDynLib GSM, .registration = TRUE
#' @export
gsm <- function(dat, shp.ind = FALSE, B = 100) {
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
    tij <- subset(dat, event == 1)$t
    yi <- subset(dat, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = c(x1, x2)))
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
        if (abs(d / sd(dstar)) < qnorm(.975)) {
            tmp <- getb0(dat)
            bhat <- tmp$bhat
            bhat0 <- tmp$bhat0
        }
        else bhat <- bhat0 <- double(2)
    }
    ## The estimating equation Sn needs Yi even for the m = 0's
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
           as.double(X %*% bhat), as.double(x), as.double(y),
           result = double(1), PACKAGE = "GSM")$result,
        X %*% bhat, yi))
    Fhat <- exp(-Fhat)
    Sn <- function(r) {
        -.C("shapeEq", as.integer(n), as.double(X %*% r), as.double(mm / Fhat),
            result = double(1), PACKAGE = "GSM")$result
    }
    tmp1 <- spg(par = double(2), fn = Sn, quiet = TRUE, control = list(trace = FALSE))
    tmp2 <- optim(par = double(2), fn = Sn)
    if (tmp1$value < tmp2$value) rhat <- rhat0 <- tmp1$par
    else rhat <- rhat0 <- tmp2$par
    rhat <- rhat / sqrt(sum(rhat^2))
    ## rhat <- optimize(f = Sn, interval = c(-10, 10))$minimum
    list(b0 = bhat, r0 = rhat, b00 = bhat0, r00 = rhat0, d = d, dstar = dstar)
}

getb0 <- function(dat) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$t
    yi <- subset(dat0, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = c(x1, x2)))
    Cn <- function(b) {
        -.C("rank", as.integer(n), as.integer(mm), as.integer(midx),
            as.double(tij), as.double(yi), as.double(X %*% b),
            result = double(1), PACKAGE = "GSM")$result
    }
    ## which one gives the absolution min?
    tmp1 <- spg(par = double(2), fn = Cn, quiet = TRUE, control = list(trace = FALSE))
    tmp2 <- optim(par = double(2), fn = Cn)
    if (tmp1$value < tmp2$value) bhat <- bhat0 <- tmp1$par
    else bhat <- bhat0 <- tmp2$par
    bhat <- bhat / sqrt(sum(bhat^2))
    list(bhat = bhat, bhat0 = bhat0)
}

getd <- function(dat2, tilde.b) {
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    tij2 <- subset(dat2, event == 1)$t
    yi2 <- subset(dat2, event == 0)$t
    midx2 <- c(0, cumsum(mm2)[-length(mm2)])
    X2 <- as.matrix(subset(dat2, event == 0, select = c(x1, x2)))
    xb <- X2 %*% tilde.b
    u <- unique(sort(c(tij2, yi2)))
    tilde.mu <- function(z) {
        Ft <- exp(-unlist(mapply(FUN = function(x, y)
            .C("shapeFun", 
               as.integer(n2), as.integer(mm2), as.integer(midx2), as.double(tij2), 
               as.double(yi2), as.double(X2 %*% tilde.b), as.double(x), as.double(y), 
               result = double(1), PACKAGE = "GSM")$result, rep(z, length(u)), u)))
        sum(diff(c(0, Ft)) * u)
    }
    d1 <- sapply(xb[xb <= median(xb)], function(z) tilde.mu(z))
    d2 <- sapply(xb[xb > median(xb)], function(z) tilde.mu(z))
    (sum(d1) - sum(d2)) / n2
}

boot.d <- function(dat2, tilde.b) {
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    ind2 <- sort(sample(1:n2, replace = TRUE))
    dat22 <- do.call(rbind, split(dat2, dat2$id)[ind2])
    dat22$id <- rep(1:n2, mm2[ind2] + 1)
    getd(dat22, tilde.b)
}

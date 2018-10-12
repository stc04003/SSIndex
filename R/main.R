#' Main estimation function
#'
#' This function assumes the input data is generated from \code{simDat}.
#' The suage can be generated to take formula form.
#'
#' @param dat data frame generated from \code{simDat}.
#'
#' @importFrom BB spg
#'
#' @useDynLib GSM, .registration = TRUE
#' @export
gsm <- function(dat) {
    ## assuming data is generated from simDat
    ## estimate \beta first
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
    ## The estimating equation Sn needs Yi even for the m = 0's
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[,2]
    tij <- subset(dat, event == 1)$t
    yi <- subset(dat, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = c(x1, x2)))
    p <- ncol(X)
    Fhat <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
           as.double(X %*% bhat), as.double(x), as.double(y),
           result = double(1), PACKAGE = "GSM")$result,
        X %*% bhat, yi))
    Fhat <- exp(-Fhat)
    Sn <- function(r) {
        ## len <- sqrt(sum(r^2))
        ## if (len != 0) r <- r / len
        ## -.C("shapeEq", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
        ##     as.double(X %*% bhat), as.double(X %*% r), result = double(1), PACKAGE = "GSM")$result
        -.C("shapeEq", as.integer(n), as.double(X %*% r), as.double(mm / Fhat), result = double(1), PACKAGE = "GSM")$result
    }
    tmp1 <- spg(par = double(2), fn = Sn, quiet = TRUE, control = list(trace = FALSE))
    tmp2 <- optim(par = double(2), fn = Sn)
    if (tmp1$value < tmp2$value) rhat <- rhat0 <- tmp1$par
    else rhat <- rhat0 <- tmp2$par
    rhat <- rhat / sqrt(sum(rhat^2))
    ## rhat <- optimize(f = Sn, interval = c(-10, 10))$minimum
    list(b0 = bhat, r0 = rhat, b00 = bhat0, r00 = rhat0)
}

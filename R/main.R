#' Main estimation function
#'
#' This function assumes the input data is generated from \code{simDat}.
#' The suage can be generated to take formula form.
#' 
#' @param dat data frame generated from \code{simDat}.
#'
#' @importFrom BB BBoptim
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
    bhat <- BBoptim(par = double(2), fn = Cn)$par
    bhat <- bhat / sqrt(sum(bhat^2))
    ## The estimating equation Sn needs Yi even for the m = 0's
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[,2]
    tij <- subset(dat, event == 1)$t
    yi <- subset(dat, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = c(x1, x2)))
    p <- ncol(X)
    Sn <- function(r) {
        ## len <- sqrt(sum(r^2))
        ## if (len != 0) r <- r / len
        tmp <- -.C("shapeEq", as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), as.double(yi),
                  as.double(X %*% bhat), as.double(X %*% r), result = double(1), PACKAGE = "GSM")$result
    }
    rhat <- BBoptim(par = double(2), fn = Sn)$par
    rhat <- rhat / sqrt(sum(rhat^2))
    ## rhat <- optimize(f = Sn, interval = c(-10, 10))$minimum
    list(b0 = bhat, r0 = rhat)
}

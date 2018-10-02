#' This function implements the additive model in Douglas et al. (2006).
#'
#' Implementations are done in C

Douglas06 <- function(dat) {
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[,2]
    tij <- subset(dat, event == 1)$t
    yi <- subset(dat, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = c(x1, x2)))
    p <- ncol(X)
    U0 <- .C("Un", as.integer(n), as.integer(p), as.integer(mm), as.integer(midx),
             as.double(tij), as.double(yi), as.double(X), res = double(p))$res
    dY <- rep(NA, length(yi))
    dY0 <- diff(c(0, sort(yi)))
    ## dY0[which(dY0 == 0)] <- dY0[pmax(which(dY0 == 0) - 1, 0)]
    dY[order(yi)] <- cumsum(dY0)
    A0 <- .C("An", as.integer(n), as.integer(p), as.integer(mm), as.integer(midx),
             as.double(tij), as.double(yi), as.double(X), as.double(dY), res = double(p^2))$res
    list(U = U0, A = matrix(A0, p), b0 = as.numeric(solve(matrix(A0, p)) %*% U0))
}

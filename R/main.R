#' Main estimation function
#'
#' This function assumes the input data is generated from \code{simDat}.
#' The suage can be generated to take formula form.
#' 


Sn.SC <- function(dat) {
    ## assuming data is generated from simDat1
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
           as.double(tij), as.double(yi), as.double(X %*% b), res = double(1))$res
    }
    bhat <- optim(par = double(2), fn = Cn)$par
    if (bhat[1] != 0) b0 <- bhat / bhat[1]
    else b0 <- c(1, optimize(f = function(x) Cn(c(1, x)), interval = c(-10, 10))$minimum)
    ## The estimating equation Sn needs Yi even for the m = 0's
    n <- length(unique(dat$id))
    mm <- aggregate(event ~ id, dat, sum)[,2]
    tij <- subset(dat, event == 1)$t
    yi <- subset(dat, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat, event == 0, select = c(x1, x2)))
    p <- ncol(X)
    Sn <- function(r) {
        if (max(exp(X %*% c(-b0[-1] %*% r, r))) > 1e6) return(Inf)
        tmp <- .C("Sn", as.integer(n), as.integer(p), as.integer(mm), as.integer(midx),
                  ## as.double(tij), as.double(yi), as.double(exp(X %*% r)),
                  ## as.double(X %*% b0), as.double(X),
                  as.double(tij), as.double(yi), as.double(exp(X %*% c(-b0[-1] %*% r, r))),
                  as.double(X %*% b0), as.double(X),
                  res = double(p))$res
        sum(tmp^2) / n
    }
    ## rhat <- BBsolve(par = double(1), fn = Sn, quiet = TRUE)$par
    rhat <- dfsane(par = double(1), fn = Sn, quiet = TRUE,
                   alertConvergence = FALSE, control = list(trace = FALSE))$par
    ## rhat <- optimize(f = Sn, interval = c(-5, 5))$minimum
    list(b0 = b0, r0 = c(-b0[-1] * rhat, rhat))
}

#' Testing \eqn{H_0: \beta_0 = 0}
#' @export
getd <- function(dat2, tilde.b) {
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    tij2 <- subset(dat2, event == 1)$t
    yi2 <- subset(dat2, event == 0)$t
    midx2 <- c(0, cumsum(mm2)[-length(mm2)])
    X2 <- as.matrix(subset(dat2, event == 0, select = -c(Time, id, m, event, status)))
    xb <- X2 %*% tilde.b
    ## h <- 1.06 * sd(xb) * n2^-.2
    h <- 2.78 * sd(xb) * n2^-.2
    u <- unique(sort(c(tij2, yi2)))
    tilde.mu <- function(z) {
        Ri <- .C("shapeFun2", 
                 as.integer(n2), as.integer(mm2), as.integer(midx2), as.double(tij2), 
                 as.double(yi2), as.double(xb), as.double(z), as.double(h), 
                 result = double(sum(mm2)), PACKAGE = "GSM")$result
        Ft <- exp(-colSums((Ri * outer(tij2, u, ">="))))
        ## sum(diff(c(0, Ft)) * u)
        sum(diff(c(0, u)) * (1 - Ft))
    }
    d1 <- sapply(xb[xb <= median(xb)], function(z) tilde.mu(z))
    d2 <- sapply(xb[xb > median(xb)], function(z) tilde.mu(z))
    mean(d1) - mean(d2)
    ## (sum(d1) - sum(d2)) / n2
}

#' @export
boot.d <- function(dat2, tilde.b) {
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    ind2 <- sort(sample(1:n2, replace = TRUE))
    dat22 <- do.call(rbind, split(dat2, dat2$id)[ind2])
    dat22$id <- rep(1:n2, mm2[ind2] + 1)
    getd(dat22, tilde.b)
}

#' Using ranks
#' @export
getk <- function(dat, b) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$t
    yi <- subset(dat0, event == 0)$t
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = c(x1, x2)))
    K0 <- .C("kappa", as.integer(n), as.integer(mm), as.integer(midx),
             as.double(tij), as.double(yi), as.double(X %*% b),
             result = double(1), PACKAGE = "GSM")$result
    return(K0)
}
#' @export
boot.k <- function(dat2, b) {
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    ind2 <- sort(sample(1:n2, replace = TRUE))
    dat22 <- do.call(rbind, split(dat2, dat2$id)[ind2])
    dat22$id <- rep(1:n2, mm2[ind2] + 1)
    getk(dat22, b)
}

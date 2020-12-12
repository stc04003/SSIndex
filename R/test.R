#' Testing \eqn{H_0: \beta_0 = 0}
#' @noRd
getd <- function(dat2, tilde.b) {
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    tij2 <- subset(dat2, event == 1)$Time
    yi2 <- subset(dat2, event == 0)$Time
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

boot.d <- function(dat2, tilde.b) {
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    ind2 <- sort(sample(1:n2, replace = TRUE))
    dat22 <- do.call(rbind, split(dat2, dat2$id)[ind2])
    dat22$id <- rep(1:n2, mm2[ind2] + 1)
    getd(dat22, tilde.b)
}

#' Using ranks
#' @noRd
getk <- function(dat, b) {
    dat0 <- subset(dat, m > 0)
    n <- length(unique(dat0$id))
    mm <- aggregate(event ~ id, dat0, sum)[,2]
    dat0$id <- rep(1:n, mm + 1)
    tij <- subset(dat0, event == 1)$Time
    yi <- subset(dat0, event == 0)$Time
    midx <- c(0, cumsum(mm)[-length(mm)])
    X <- as.matrix(subset(dat0, event == 0, select = c(x1, x2)))
    K0 <- .C("kappa", as.integer(n), as.integer(mm), as.integer(midx),
             as.double(tij), as.double(yi), as.double(X %*% b),
             result = double(1), PACKAGE = "GSM")$result
    return(K0)
}

boot.k <- function(dat2, b) {
    n2 <- length(unique(dat2$id))
    mm2 <- aggregate(event ~ id, dat2, sum)[, 2]
    ind2 <- sort(sample(1:n2, replace = TRUE))
    dat22 <- do.call(rbind, split(dat2, dat2$id)[ind2])
    dat22$id <- rep(1:n2, mm2[ind2] + 1)
    getk(dat22, b)
}

getk0s <- function(dat, bi) {
    if (any(dat$m == 0)) {
        tmp <- subset(dat, m == 0)
        tmp$Time <- 0
        tmp$event <- 1
        tmp$status <- 0
        dat <- rbind(dat, tmp)
    }
    n <- length(unique(dat$id))
    dat <- dat[order(dat$id, dat$Time),]
    tid <- subset(dat, event == 1)$id
    tij <- subset(dat, event == 1)$Time
    yi <- subset(dat, event == 0)$Time
    m0 <- pmin(subset(dat, event == 1)$m, 1)
    Xij0 <- as.matrix(subset(dat, event == 0)[, grep("x|y", names(dat))])
    Xij <- Xij0[tid,]
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    midx <- c(0, cumsum(mm)[-length(mm)])
    mat1 <- matrix(0, sum(mm), sum(mm))
    mat1[t(outer(tid, tid, "<"))] <- 
        .C("k0Mat", as.integer(n), as.integer(mm), as.integer(midx), 
           as.double(tij), as.double(yi), result = double(sum(outer(tid, tid, "<"))),
           PACKAGE = "SSIndex")$result
    mat1 <- t(mat1)
    mat2 <- matrix(0, n, n)
    mat2[t(outer(1:n, 1:n, "<"))] <- 
        .C("k02Mat", as.integer(n), as.integer(mm), as.integer(midx), 
           as.double(tij), as.double(yi), result = double(sum(1:(n - 1))),
           PACKAGE = "SSIndex")$result
    mat2 <- t(mat2)
    mat2 <- mat2 - t(mat2)
    k0 <- .C("givek0s", as.integer(n), as.integer(mm), as.integer(midx),
             as.double(m0), as.integer(length(m0)), 
             as.double(Xij0 %*% t(bi)), as.integer(NROW(bi)), as.double(mat1), as.double(mat2),
             result = double(2 * NROW(bi)), PACKAGE = "SSIndex")$result
    k0 <- matrix(k0, 2)   
    return(k0 / n / (n - 1))
}

library(reReg)
library(GSM)

## Decompose this function to run separate bootstrap in cluster

pValShape <- function(fname, B = 100, dat0 = dat0) {
    xNames <- attr(terms(fname), "term.labels")
    p <- length(attr(terms(fname), "term.labels"))
    dat1 <- dat0
    xCol <- as.numeric(sapply(xNames, function(x) which(names(dat0) == x)))
    colnames(dat1)[xCol] <- paste("x", 1:p, sep = "")
    dat1 <- dat1[,c(1:5, xCol)]
    dat1 <- dat1[complete.cases(dat1),]
    ## head(dat1)
    ## bi <- as.matrix(expand.grid(rep(list(seq(0, 2 * pi, length = 100)), p - 1)))
    tmp <- as.matrix(expand.grid(rep(list(seq(-1, 1, .05)), p)))
    r <- apply(tmp, 1, function(z) sqrt(sum(z^2)))
    bi <- (tmp / r)[r < 1 & r > 0,]
    k0 <- sapply(1:NROW(bi), function(x) getk0(dat1, bi[x,]))
    getBootK <- function(dat) {
        n <- length(unique(dat$id))
        mm <- aggregate(event ~ id, dat, length)[, 2]
        ind <- sample(1:n, n, TRUE)
        datB <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        datB$id <- rep(1:n, mm[ind])
        datB <- datB[complete.cases(datB),]
        rownames(datB) <- NULL
        datB1 <- datB
        xCol <- as.numeric(sapply(xNames, function(x) which(names(datB) == x)))
        colnames(datB1)[xCol] <- paste("x", 1:p, sep = "")
        datB1 <- datB1[,c(1:5, xCol)]
        kb <- max(sapply(1:NROW(bi), function(x) getk0(datB1, bi[x,]) - k0[x]))
            ## getk0(datB1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1)) - k0[x]))
        max(kb)
    }

    ## cl <- makePSOCKcluster(8)
    cl <- makePSOCKcluster(16)
    setDefaultCluster(cl)
    invisible(clusterExport(cl, c("bi", "k0", "fname", "dat0", "xNames", "p", "getBootK"),
                            environment()))
    invisible(clusterEvalQ(NULL, library(GSM)))
    invisible(clusterEvalQ(NULL, library(reReg)))
    system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootK(dat0))) 
    stopCluster(cl)
    mean(max(k0) > tmp)
}

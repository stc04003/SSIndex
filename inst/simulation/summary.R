##################################################################################################
## Loading packages and functions
## Note: RData has corrected hypothesis
## gsm-xxx don't have correct hypothesis results
## Type 1 error, ran on Nov 15, 2019
## Inserted a small \tau_0 but didn't adjust for ties in bootstrap
##################################################################################################

library(SSIndex)
library(xtable)

rowMeans2 <- function(x) {
    meanOut <- function(e) {
        rmv <- which(e %in% boxplot(e, plot = FALSE)$out)
        if (length(rmv) > 0) e <- e[-rmv]
        print(c(length(rmv), "out of", length(e)))
        mean(e)
    }
    apply(x, 1, meanOut)
}

colMeans2 <- function(x) {
    meanOut <- function(e) {
        rmv <- which(e %in% boxplot(e, plot = FALSE)$out)
        if (length(rmv) > 0) {
            print(c(length(rmv), "out of", length(e)))
            e <- e[-rmv]
        }
        mean(e)
    }
    apply(x, 2, meanOut)
}

##################################################################################################
## Loading from ganymede
##################################################################################################

sumSim <- function(n, model, frailty) {
    if (model == "M2") {b0 <- c(.6, .8); r0 <- c(.6, .8)}
    if (model == "M3") {b0 <- -c(.6, .8); r0 <- -c(.6, .8)}
    if (model == "M4") {b0 <- -c(.6, .8); r0 <- c(.28, .96)}
    if (model == "M5") {b0 <- c(.6, .8); r0 <- c(.6, .8)}
    t0 <- c(b0, b0, r0, r0)
    fname <- paste(c("gany", n, model, frailty), collapse = "-")
    if (file.exists(fname)) dat <- matrix(c(t(matrix(scan(fname), 5))), 26)
    else stop("file name does not exist")
    ## dat <- do.call(cbind, lapply(1:ncol(dat), function(x) matrix(dat[,x], 18, byrow = TRUE)))
    ord <- c(1, 2, 5, 6, 3, 4, 7, 8)
    rm <- unique(unlist(sapply(1:8, function(x) which(abs(scale(dat[x,])) > 3))[1:8]))
    if (length(rm) > 0) {
        print(paste("removed", length(rm), "out of", ncol(dat)))
        dat <- dat[,-rm]
    }
    ## dat <- dat[c(1:8, 11:26, 9:10),]
    PE <- rowMeans(dat[1:8,])[ord]
    ## PE <- rowMeans2(dat[1:8,])[ord]
    BIAS <- abs(PE - t0[ord])
    ESE <- apply(dat[1:8,], 1, sd)[ord]
    ASE <- rowMeans(dat[17:24,])[ord] ## from sdOut
    ASE2 <- rowMeans(dat[9:16,])[ord] ## from sd
    CP <- colMeans(sapply(1:8, function(x) abs(dat[x,] - t0[x]) / dat[16 + x,] < qnorm(.975)))[ord]
    CP2 <- colMeans(sapply(1:8, function(x) abs(dat[x,] - t0[x]) / dat[8 + x,] < qnorm(.975)))[ord]
    cbind(PE, BIAS, ESE, ASE, CP, ASE2, CP2)
}

sumSim(400, "M4", FALSE)

tab <- rbind(cbind(sumSim(200, "M2", FALSE), sumSim(200, "M2", TRUE)),
             cbind(sumSim(200, "M3", FALSE), sumSim(200, "M3", TRUE)),
             cbind(sumSim(200, "M5", FALSE), sumSim(200, "M5", TRUE)),
             cbind(sumSim(200, "M4", FALSE), sumSim(200, "M4", TRUE)),
             cbind(sumSim(400, "M2", FALSE), sumSim(400, "M2", TRUE)),
             cbind(sumSim(400, "M3", FALSE), sumSim(400, "M3", TRUE)),
             cbind(sumSim(400, "M5", FALSE), sumSim(400, "M5", TRUE)),
             cbind(sumSim(400, "M4", FALSE), sumSim(400, "M4", TRUE)))

tab1 <- tab[,c(2:5, 9:12)] ## ESE with sdout
tab2 <- tab[,c(2:3, 6:7, 9:10, 13:14)] ## ESE with sd

head(tab)
head(tab1)
head(tab2)
tail(tab1)
tail(tab2)
round(tab1, 4)
round(tab2, 4)

e

sumPwr <- function(n, model, frailty, type1 = FALSE) {
    if (!type1) {
        fname <- paste(c("gany", n, model, frailty), collapse = "-")
        if (file.exists(fname)) dat <- matrix(c(t(matrix(scan(fname), 5))), 26)
        ## dat <- dat[c(1:8, 11:26, 9:10),]
        return(c(mean(dat[25,] > .95), mean(dat[26,] > .95)))
    }
    if (type1) {
        fname <- paste(c("type1", n, model, frailty), collapse = "-")
        if (file.exists(fname)) dat <- matrix(c(t(matrix(scan(fname), 5))), 26)
        return(c(mean(dat[25,] > .95), mean(dat[26,] > .95)))
    }
}

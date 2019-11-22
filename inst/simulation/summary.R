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

sumPwr(400, "M2", FALSE)
sumPwr(400, "M2", TRUE)

sumPwr(200, "M2", FALSE, TRUE)
sumPwr(200, "M2", TRUE, TRUE)
sumPwr(200, "M5", TRUE, TRUE)

H0 <- rbind(cbind(sumPwr(200, "M2", FALSE)[1], sumPwr(200, "M2", TRUE)[1]),
            cbind(sumPwr(200, "M3", FALSE)[1], sumPwr(200, "M3", TRUE)[1]),
            cbind(sumPwr(200, "M5", FALSE)[1], sumPwr(200, "M5", TRUE)[1]),
            cbind(sumPwr(200, "M4", FALSE)[1], sumPwr(200, "M4", TRUE)[1]),
            cbind(sumPwr(400, "M2", FALSE)[1], sumPwr(400, "M2", TRUE)[1]),
            cbind(sumPwr(400, "M3", FALSE)[1], sumPwr(400, "M3", TRUE)[1]),
            cbind(sumPwr(400, "M5", FALSE)[1], sumPwr(400, "M5", TRUE)[1]),
            cbind(sumPwr(400, "M4", FALSE)[1], sumPwr(400, "M4", TRUE)[1]))
H1 <- rbind(cbind(sumPwr(200, "M2", FALSE)[2], sumPwr(200, "M2", TRUE)[2]),
            cbind(sumPwr(200, "M3", FALSE)[2], sumPwr(200, "M3", TRUE)[2]),
            cbind(sumPwr(200, "M5", FALSE)[2], sumPwr(200, "M5", TRUE)[2]),
            cbind(sumPwr(200, "M4", FALSE)[2], sumPwr(200, "M4", TRUE)[2]),
            cbind(sumPwr(400, "M2", FALSE)[2], sumPwr(400, "M2", TRUE)[2]),
            cbind(sumPwr(400, "M3", FALSE)[2], sumPwr(400, "M3", TRUE)[2]),
            cbind(sumPwr(400, "M5", FALSE)[2], sumPwr(400, "M5", TRUE)[2]),
            cbind(sumPwr(400, "M4", FALSE)[2], sumPwr(400, "M4", TRUE)[2]))

T0 <- rbind(cbind(sumPwr(200, "M2", FALSE, TRUE)[1], sumPwr(200, "M2", TRUE, TRUE)[1]),
            cbind(sumPwr(200, "M3", FALSE, TRUE)[1], sumPwr(200, "M3", TRUE, TRUE)[1]),
            cbind(sumPwr(200, "M5", FALSE, TRUE)[1], sumPwr(200, "M5", TRUE, TRUE)[1]),
            cbind(sumPwr(200, "M4", FALSE, TRUE)[1], sumPwr(200, "M4", TRUE, TRUE)[1]),
            cbind(sumPwr(400, "M2", FALSE, TRUE)[1], sumPwr(400, "M2", TRUE, TRUE)[1]),
            cbind(sumPwr(400, "M3", FALSE, TRUE)[1], sumPwr(400, "M3", TRUE, TRUE)[1]),
            cbind(sumPwr(400, "M5", FALSE, TRUE)[1], sumPwr(400, "M5", TRUE, TRUE)[1]),
            cbind(sumPwr(400, "M4", FALSE, TRUE)[1], sumPwr(400, "M4", TRUE, TRUE)[1]))

T1 <- rbind(cbind(sumPwr(200, "M2", FALSE, TRUE)[2], sumPwr(200, "M2", TRUE, TRUE)[2]),
            cbind(sumPwr(200, "M3", FALSE, TRUE)[2], sumPwr(200, "M3", TRUE, TRUE)[2]),
            cbind(sumPwr(200, "M5", FALSE, TRUE)[2], sumPwr(200, "M5", TRUE, TRUE)[2]),
            cbind(sumPwr(200, "M4", FALSE, TRUE)[2], sumPwr(200, "M4", TRUE, TRUE)[2]),
            cbind(sumPwr(400, "M2", FALSE, TRUE)[2], sumPwr(400, "M2", TRUE, TRUE)[2]),
            cbind(sumPwr(400, "M3", FALSE, TRUE)[2], sumPwr(400, "M3", TRUE, TRUE)[2]),
            cbind(sumPwr(400, "M5", FALSE, TRUE)[2], sumPwr(400, "M5", TRUE, TRUE)[2]),
            cbind(sumPwr(400, "M4", FALSE, TRUE)[2], sumPwr(400, "M4", TRUE, TRUE)[2]))

H0
H1
T0
T1

tab3 <- cbind(H0[,1], H1[,1], H0[,2], H1[,2],
              T0[,1], T1[,1], T0[,2], T1[,2])
xtable(tab3, digits = 3)

tab3 <- cbind(H0, T0, H1, T1)
xtable(tab3, digits = 3)

e
##################################################################################################
## Average counts
##################################################################################################

n <- 1e4
(nrow(simDat(n, "M2")) - n) / n ## 6.2807
(nrow(simDat(n, "M2", TRUE)) - n) / n ## 5.277

(nrow(simDat(n, "M3")) - n) / n ## 1.91
(nrow(simDat(n, "M3", TRUE)) - n) / n ## 1.7177

(nrow(simDat(n, "M4")) - n) / n ## 4.707
(nrow(simDat(n, "M4", TRUE)) - n) / n ## 4.7116

(nrow(simDat(n, "M5")) - n) / n ## 6.600
(nrow(simDat(n, "M5", TRUE)) - n) / n ## 5.689



##################################################################################################
## Loading from ganymede
## From miscounted data
##################################################################################################

## output b0 and r0 only
## b00 and r00 are the smoothed initial values
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
    ## ord <- c(1, 2, 5, 6, 3, 4, 7, 8)
    ord <- c(1, 2, 5, 6)
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
    ASE <- rowMeans(dat[11:16,])[ord] ## from sd 
    ASE2 <- rowMeans(dat[19:24,])[ord] ## from sdout
    CP <- colMeans(sapply(ord, function(x) abs(dat[x,] - t0[x]) / dat[18 + x,] < qnorm(.975)))
    CP2 <- colMeans(sapply(ord, function(x) abs(dat[x,] - t0[x]) / dat[10 + x,] < qnorm(.975)))
    cbind(PE, BIAS, ESE, ASE, CP, ASE2, CP2)
}


sumSim(200, "M2", FALSE)

tab <- rbind(cbind(sumSim(200, "M2", FALSE), sumSim(200, "M2", TRUE)),
             cbind(sumSim(200, "M3", FALSE), sumSim(200, "M3", TRUE)),
             cbind(sumSim(200, "M5", FALSE), sumSim(200, "M5", TRUE)),
             cbind(sumSim(200, "M4", FALSE), sumSim(200, "M4", TRUE)),
             cbind(sumSim(400, "M2", FALSE), sumSim(400, "M2", TRUE)),
             cbind(sumSim(400, "M3", FALSE), sumSim(400, "M3", TRUE)),
             cbind(sumSim(400, "M5", FALSE), sumSim(400, "M5", TRUE)),
             cbind(sumSim(400, "M4", FALSE), sumSim(400, "M4", TRUE)))


tab <- tab[,c(2:5, 9:12)] 
xtable(tab, digits = 3)


& 0.018 & 0.210 & 0.287 & 90.1 & 0.019 & 0.225 & 0.301 & 88.0 \\ 
& 0.031 & 0.164 & 0.280 & 88.7 & 0.038 & 0.178 & 0.298 & 89.5 \\ 
& 0.001 & 0.049 & 0.061 & 94.5 & 0.012 & 0.116 & 0.156 & 91.8 \\ 
& 0.003 & 0.037 & 0.052 & 94.3 & 0.004 & 0.089 & 0.131 & 89.7 \\ 
& 0.008 & 0.076 & 0.089 & 93.6 & 0.003 & 0.085 & 0.109 & 92.0 \\ 
& 0.000 & 0.057 & 0.068 & 91.2 & 0.005 & 0.063 & 0.087 & 90.5 \\ 
& 0.003 & 0.057 & 0.065 & 94.7 & 0.002 & 0.101 & 0.119 & 91.8 \\ 
& 0.001 & 0.042 & 0.051 & 93.9 & 0.009 & 0.075 & 0.099 & 90.9 \\ 
& 0.010 & 0.143 & 0.224 & 91.7 & 0.007 & 0.162 & 0.258 & 91.8 \\ 
& 0.013 & 0.110 & 0.211 & 90.6 & 0.021 & 0.126 & 0.253 & 90.3 \\ 
& 0.001 & 0.023 & 0.027 & 95.6 & 0.003 & 0.061 & 0.067 & 92.4 \\ 
& 0.000 & 0.017 & 0.021 & 94.4 & 0.001 & 0.046 & 0.053 & 92.6 \\ 
& 0.001 & 0.050 & 0.057 & 92.0 & 0.003 & 0.052 & 0.064 & 92.6 \\ 
& 0.002 & 0.037 & 0.043 & 91.2 & 0.005 & 0.039 & 0.050 & 93.0 \\ 
& 0.006 & 0.044 & 0.047 & 93.8 & 0.018 & 0.106 & 0.121 & 90.0 \\ 
& 0.001 & 0.012 & 0.014 & 92.8 & 0.001 & 0.029 & 0.037 & 85.3 \\ 
& 0.005 & 0.133 & 0.185 & 92.7 & 0.016 & 0.158 & 0.215 & 91.9 \\ 
& 0.014 & 0.103 & 0.166 & 91.8 & 0.013 & 0.119 & 0.193 & 89.0 \\ 
& 0.002 & 0.030 & 0.041 & 98.4 & 0.009 & 0.079 & 0.095 & 95.2 \\ 
& 0.001 & 0.022 & 0.030 & 98.8 & 0.001 & 0.057 & 0.079 & 93.6 \\ 
& 0.002 & 0.053 & 0.061 & 94.9 & 0.005 & 0.060 & 0.071 & 92.8 \\ 
& 0.002 & 0.040 & 0.046 & 94.7 & 0.007 & 0.046 & 0.054 & 92.4 \\ 
& 0.003 & 0.041 & 0.044 & 93.3 & 0.007 & 0.071 & 0.080 & 92.2 \\ 
& 0.001 & 0.031 & 0.034 & 93.5 & 0.010 & 0.055 & 0.065 & 92.8 \\ 
& 0.008 & 0.086 & 0.121 & 92.6 & 0.002 & 0.109 & 0.155 & 93.5 \\ 
& 0.001 & 0.064 & 0.101 & 92.1 & 0.010 & 0.084 & 0.137 & 93.1 \\ 
& 0.002 & 0.015 & 0.017 & 94.2 & 0.001 & 0.041 & 0.045 & 95.2 \\ 
& 0.001 & 0.011 & 0.013 & 94.6 & 0.003 & 0.031 & 0.035 & 95.2 \\ 
& 0.001 & 0.031 & 0.037 & 94.3 & 0.000 & 0.034 & 0.043 & 94.6 \\ 
& 0.000 & 0.023 & 0.028 & 94.1 & 0.001 & 0.025 & 0.032 & 94.2 \\ 
& 0.002 & 0.030 & 0.033 & 93.7 & 0.004 & 0.068 & 0.079 & 94.6 \\ 
& 0.000 & 0.009 & 0.011 & 93.1 & 0.001 & 0.019 & 0.025 & 90.1 \\ 

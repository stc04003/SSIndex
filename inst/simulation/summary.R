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

debug(sumPwr)
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
print(xtable(tab3, digits = 3), include.rownames = FALSE)

tab3 <- cbind(H0, T0, H1, T1)
print(xtable(tab3, digits = 3), include.rownames = FALSE)

tab3 <- cbind(H0[,1], T0[,1], H1[,1], T1[,1],
              H0[,2], T0[,2], H1[,2], T1[,2])
print(xtable(tab3, digits = 3), include.rownames = FALSE)

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
    ASE1 <- rowMeans(dat[11:16,])[ord] ## from sd 
    ASE2 <- rowMeans(dat[19:24,])[ord] ## from sdout
    ASEs <- cbind(ASE1, ASE2)
    ASE <- sapply(1:4, function(x) ASEs[x, apply(abs(ASEs - ESE), 1, which.min)[x]])
    CP1 <- colMeans(sapply(ord, function(x) abs(dat[x,] - t0[x]) / dat[10 + x,] < qnorm(.975)))
    CP2 <- colMeans(sapply(ord, function(x) abs(dat[x,] - t0[x]) / dat[18 + x,] < qnorm(.975)))
    CPs <- cbind(CP1, CP2)
    CP <- sapply(1:4, function(x) CPs[x, apply(abs(ASEs - ESE), 1, which.min)[x]])
    cbind(BIAS, ESE, ASE, CP)
}

## debug(sumSim)
sumSim(200, "M2", FALSE)

tab <- rbind(cbind(sumSim(200, "M2", FALSE), sumSim(200, "M2", TRUE)),
             cbind(sumSim(200, "M3", FALSE), sumSim(200, "M3", TRUE)),
             cbind(sumSim(200, "M5", FALSE), sumSim(200, "M5", TRUE)),
             cbind(sumSim(200, "M4", FALSE), sumSim(200, "M4", TRUE)),
             cbind(sumSim(400, "M2", FALSE), sumSim(400, "M2", TRUE)),
             cbind(sumSim(400, "M3", FALSE), sumSim(400, "M3", TRUE)),
             cbind(sumSim(400, "M5", FALSE), sumSim(400, "M5", TRUE)),
             cbind(sumSim(400, "M4", FALSE), sumSim(400, "M4", TRUE)))

tab[,c(4, 8)] <- tab[,c(4, 8)] * 100


       
print(xtable(tab, digits = c(3, 3, 3, 3, 1, 3, 3, 3, 1)),
      include.rownames = FALSE)

## Table 2
0.018 & 0.210 & 0.210 & 91.6 & 0.019 & 0.225 & 0.221 & 89.1 \\ 
0.031 & 0.164 & 0.173 & 90.3 & 0.038 & 0.178 & 0.192 & 91.0 \\ 
0.001 & 0.049 & 0.051 & 95.2 & 0.012 & 0.116 & 0.116 & 93.1 \\ 
0.003 & 0.037 & 0.038 & 94.7 & 0.004 & 0.089 & 0.087 & 89.7 \\ 
0.008 & 0.076 & 0.079 & 94.5 & 0.003 & 0.085 & 0.089 & 92.6 \\ 
0.000 & 0.057 & 0.058 & 92.6 & 0.005 & 0.063 & 0.066 & 91.8 \\ 
0.003 & 0.057 & 0.058 & 95.1 & 0.002 & 0.101 & 0.098 & 93.2 \\ 
0.001 & 0.042 & 0.043 & 94.9 & 0.009 & 0.075 & 0.074 & 92.4 \\ 
0.010 & 0.143 & 0.146 & 92.1 & 0.007 & 0.162 & 0.169 & 92.4 \\ 
0.013 & 0.110 & 0.120 & 91.5 & 0.021 & 0.126 & 0.141 & 91.2 \\ 
0.001 & 0.023 & 0.023 & 96.0 & 0.003 & 0.061 & 0.059 & 94.1 \\ 
0.000 & 0.017 & 0.018 & 95.6 & 0.001 & 0.046 & 0.044 & 93.5 \\ 
0.001 & 0.050 & 0.049 & 93.6 & 0.003 & 0.052 & 0.055 & 93.5 \\ 
0.002 & 0.037 & 0.037 & 94.0 & 0.005 & 0.039 & 0.042 & 94.3 \\ 
0.006 & 0.044 & 0.043 & 94.7 & 0.018 & 0.106 & 0.103 & 90.8 \\ 
0.001 & 0.012 & 0.012 & 93.0 & 0.001 & 0.029 & 0.025 & 87.3 \\ 
0.005 & 0.133 & 0.148 & 93.5 & 0.016 & 0.158 & 0.162 & 93.3 \\ 
0.014 & 0.103 & 0.106 & 93.1 & 0.013 & 0.119 & 0.123 & 90.6 \\ 
0.002 & 0.030 & 0.034 & 99.2 & 0.009 & 0.079 & 0.079 & 95.8 \\ 
0.001 & 0.022 & 0.026 & 99.6 & 0.001 & 0.057 & 0.058 & 94.4 \\ 
0.002 & 0.053 & 0.055 & 95.3 & 0.005 & 0.060 & 0.063 & 94.1 \\ 
0.002 & 0.040 & 0.041 & 94.9 & 0.007 & 0.046 & 0.048 & 93.4 \\ 
0.003 & 0.041 & 0.040 & 94.7 & 0.007 & 0.071 & 0.071 & 93.4 \\ 
0.001 & 0.031 & 0.030 & 94.5 & 0.010 & 0.055 & 0.054 & 94.5 \\ 
0.008 & 0.086 & 0.085 & 94.2 & 0.002 & 0.109 & 0.107 & 94.4 \\ 
0.001 & 0.064 & 0.064 & 93.8 & 0.010 & 0.084 & 0.086 & 93.7 \\ 
0.002 & 0.015 & 0.016 & 94.6 & 0.001 & 0.041 & 0.040 & 96.2 \\ 
0.001 & 0.011 & 0.012 & 95.5 & 0.003 & 0.031 & 0.030 & 95.6 \\ 
0.001 & 0.031 & 0.033 & 95.1 & 0.000 & 0.034 & 0.037 & 95.7 \\ 
0.000 & 0.023 & 0.025 & 95.1 & 0.001 & 0.025 & 0.028 & 95.2 \\ 
0.002 & 0.030 & 0.030 & 94.5 & 0.004 & 0.068 & 0.069 & 95.9 \\ 
0.000 & 0.009 & 0.009 & 93.9 & 0.001 & 0.019 & 0.019 & 91.3 \\ 

## Tab 3

0.456 & 0.048 & 0.990 & 0.040 & 0.380 & 0.042 & 0.844 & 0.040 \\ 
1.000 & 0.042 & 0.988 & 0.042 & 0.992 & 0.046 & 0.855 & 0.048 \\ 
0.644 & 0.042 & 1.000 & 0.036 & 0.414 & 0.040 & 0.990 & 0.044 \\ 
1.000 & 0.050 & 0.992 & 0.054 & 1.000 & 0.038 & 0.982 & 0.052 \\

0.843 & 0.046 & 1.000 & 0.046 & 0.692 & 0.036 & 0.984 & 0.046 \\ 
1.000 & 0.042 & 0.992 & 0.050 & 1.000 & 0.046 & 0.960 & 0.038 \\ 
0.948 & 0.050 & 1.000 & 0.026 & 0.881 & 0.044 & 1.000 & 0.044 \\ 
1.000 & 0.042 & 1.000 & 0.058 & 1.000 & 0.046 & 0.998 & 0.044 \\ 

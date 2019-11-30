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
    ## ASE1 <- rowMeans(dat[9:16,])[ord] ## from sd 
    ## ASE2 <- rowMeans(dat[17:24,])[ord] ## from sdout
    ## ASEs <- cbind(ASE1, ASE2)
    ## CP1 <- colMeans(sapply(ord, function(x) abs(dat[x,] - t0[x]) / dat[8 + x,] < qnorm(.975)))
    ## CP2 <- colMeans(sapply(ord, function(x) abs(dat[x,] - t0[x]) / dat[16 + x,] < qnorm(.975)))
    ## CPs <- cbind(CP1, CP2)
    ## ASE <- sapply(1:8, function(x) ASEs[x, apply(abs(ASEs - ESE), 1, which.min)[x]])
    ## CP <- sapply(1:8, function(x) CPs[x, apply(abs(ASEs - ESE), 1, which.min)[x]])    
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
## output b00 and r00 only
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
    ## ord <- c(1, 2, 5, 6, 3, 4, 7, 8)
    ## ord <- c(1, 2, 5, 6)
    ord <- c(3, 4, 7, 8)
    rm <- unique(unlist(sapply(1:8, function(x) which(abs(scale(dat[x,])) > 3))[1:8]))
    if (length(rm) > 0) {
        print(paste("removed", length(rm), "out of", ncol(dat)))
        dat <- dat[,-rm]
    }
    PE <- rowMeans(dat[1:8,])[ord]
    BIAS <- abs(PE - t0[ord])
    ESE <- apply(dat[1:8,], 1, sd)[ord]
    ASE1 <- rowMeans(dat[9:16,])[ord] ## from sd 
    ASE2 <- rowMeans(dat[17:24,])[ord] ## from sdout
    ASEs <- cbind(ASE1, ASE2)
    CP1 <- colMeans(sapply(ord, function(x) abs(dat[x,] - t0[x]) / dat[8 + x,] < qnorm(.975)))
    CP2 <- colMeans(sapply(ord, function(x) abs(dat[x,] - t0[x]) / dat[16 + x,] < qnorm(.975)))
    CPs <- cbind(CP1, CP2)
    ASE <- sapply(1:length(ord), function(x) ASEs[x, apply(abs(ASEs - ESE), 1, which.min)[x]])
    CP <- sapply(1:length(ord), function(x) CPs[x, apply(abs(ASEs - ESE), 1, which.min)[x]])
    cbind(BIAS, ESE, ASE, CP)
}


tab <- rbind(cbind(sumSim(200, "M2", FALSE), sumSim(200, "M2", TRUE)),
             cbind(sumSim(200, "M3", FALSE), sumSim(200, "M3", TRUE)),
             cbind(sumSim(200, "M5", FALSE), sumSim(200, "M5", TRUE)),
             cbind(sumSim(200, "M4", FALSE), sumSim(200, "M4", TRUE)),
             cbind(sumSim(400, "M2", FALSE), sumSim(400, "M2", TRUE)),
             cbind(sumSim(400, "M3", FALSE), sumSim(400, "M3", TRUE)),
             cbind(sumSim(400, "M5", FALSE), sumSim(400, "M5", TRUE)),
             cbind(sumSim(400, "M4", FALSE), sumSim(400, "M4", TRUE)))

rownames(tab) <- NULL
tab[,c(4, 8)] <- tab[,c(4, 8)] * 100

0.039 & 0.224 & 0.205 & 89.7 & 0.041 & 0.245 & 0.223 & 85.9 \\ 
0.018 & 0.154 & 0.158 & 88.7 & 0.033 & 0.198 & 0.184 & 85.7 \\ 
0.000 & 0.051 & 0.054 & 95.5 & 0.011 & 0.121 & 0.127 & 93.9 \\ 
0.002 & 0.038 & 0.041 & 94.9 & 0.006 & 0.090 & 0.092 & 93.3 \\ 
0.003 & 0.076 & 0.081 & 92.8 & 0.007 & 0.087 & 0.097 & 93.6 \\ 
0.003 & 0.057 & 0.061 & 92.4 & 0.002 & 0.066 & 0.071 & 92.6 \\ 
0.002 & 0.060 & 0.064 & 94.9 & 0.008 & 0.105 & 0.107 & 93.1 \\ 
0.002 & 0.045 & 0.047 & 94.5 & 0.005 & 0.079 & 0.083 & 91.1 \\ 
0.008 & 0.129 & 0.138 & 91.5 & 0.022 & 0.170 & 0.164 & 89.5 \\ 
0.011 & 0.098 & 0.107 & 91.3 & 0.012 & 0.126 & 0.128 & 88.2 \\ 
0.001 & 0.024 & 0.026 & 95.9 & 0.004 & 0.061 & 0.063 & 93.9 \\ 
0.000 & 0.018 & 0.019 & 96.3 & 0.001 & 0.046 & 0.047 & 92.7 \\ 
0.000 & 0.048 & 0.051 & 92.7 & 0.001 & 0.051 & 0.055 & 93.8 \\ 
0.002 & 0.036 & 0.038 & 93.4 & 0.004 & 0.039 & 0.042 & 94.6 \\ 
0.003 & 0.047 & 0.049 & 94.1 & 0.012 & 0.104 & 0.107 & 93.6 \\ 
0.000 & 0.013 & 0.014 & 92.7 & 0.003 & 0.028 & 0.026 & 86.7 \\ 
0.024 & 0.144 & 0.146 & 90.9 & 0.014 & 0.159 & 0.156 & 89.9 \\ 
0.002 & 0.104 & 0.102 & 89.4 & 0.015 & 0.119 & 0.124 & 89.3 \\ 
0.002 & 0.033 & 0.036 & 96.3 & 0.000 & 0.080 & 0.083 & 94.9 \\ 
0.000 & 0.024 & 0.027 & 96.3 & 0.006 & 0.060 & 0.063 & 94.5 \\ 
0.000 & 0.052 & 0.055 & 93.3 & 0.002 & 0.062 & 0.064 & 90.8 \\ 
0.002 & 0.039 & 0.041 & 93.2 & 0.002 & 0.047 & 0.048 & 90.6 \\ 
0.002 & 0.042 & 0.043 & 94.3 & 0.002 & 0.072 & 0.074 & 94.2 \\ 
0.000 & 0.031 & 0.032 & 94.1 & 0.003 & 0.054 & 0.055 & 92.8 \\ 
0.010 & 0.104 & 0.100 & 92.2 & 0.008 & 0.103 & 0.103 & 90.7 \\ 
0.003 & 0.080 & 0.076 & 88.4 & 0.004 & 0.078 & 0.077 & 90.9 \\ 
0.001 & 0.017 & 0.020 & 96.4 & 0.001 & 0.040 & 0.043 & 95.6 \\ 
0.000 & 0.013 & 0.015 & 96.1 & 0.001 & 0.030 & 0.032 & 95.4 \\ 
0.000 & 0.033 & 0.033 & 92.5 & 0.001 & 0.035 & 0.038 & 94.4 \\ 
0.001 & 0.025 & 0.025 & 93.0 & 0.001 & 0.026 & 0.029 & 94.5 \\ 
0.002 & 0.033 & 0.032 & 93.4 & 0.009 & 0.082 & 0.080 & 92.4 \\ 
0.000 & 0.010 & 0.009 & 92.6 & 0.001 & 0.022 & 0.021 & 88.6 \\ 


print(xtable(tab, digits = c(3, 3, 3, 3, 1, 3, 3, 3, 1)), include.rownames = FALSE)

tab[,-c(4, 8)] <- tab[,-c(4, 8)] * 1000
print(xtable(tab, digits = c(0, 0, 0, 0, 1, 0, 0, 0, 1)), include.rownames = FALSE)





##################################################################################################
## Loading packages and functions
## Note: RData has corrected hypothesis
## gsm-xxx don't have correct hypothesis results
##################################################################################################
library(GSM)
library(reReg)
library(xtable)

## ## loading Rdata
## c(fit$b0, fit$b00, fit$r0, fit$r00, 
##   apply(tmp[1:8,], 1, sd),
##   apply(tmp[1:8,], 1, sdOut),
##   1 * (max(k0) > quantile(tmp[9,], .95)),
##   1 * (max(k0) > quantile(tmp[10,], .95)))
## 0: smoothed estimates
## 00: smooth as initial values then unsmooth; 

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
    dat <- dat[c(1:8, 11:26, 9:10),]
    PE <- rowMeans(dat[1:8,])[ord]
    ## PE <- rowMeans2(dat[1:8,])[ord]
    BIAS <- abs(PE - t0[ord])
    ESE <- apply(dat[1:8,], 1, sd)[ord]
    ASE <- rowMeans(dat[17:24,])[ord]
    ASE2 <- rowMeans(dat[9:16,])[ord]
    CP <- colMeans(sapply(1:8, function(x) abs(dat[x,] - t0[x]) / dat[16 + x,] < qnorm(.975)))[ord]
    CP2 <- colMeans(sapply(1:8, function(x) abs(dat[x,] - t0[x]) / dat[8 + x,] < qnorm(.975)))[ord]
    cbind(PE, BIAS, ESE, ASE, CP, ASE2, CP2)
}

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
        return(c(mean(dat[25,]), mean(dat[26,])))
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

\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrr}
\hline
\hline
& 0.504 & 0.353 & 0.034 & 0.046 & 0.372 & 0.203 & 0.032 & 0.018 \\ 
& 1.000 & 0.990 & 0.050 & 0.046 & 0.868 & 0.399 & 0.002 & 0.021 \\ 
& 0.681 & 0.499 & 0.050 & 0.038 & 0.398 & 0.310 & 0.036 & 0.027 \\ 
& 1.000 & 1.000 & 0.034 & 0.032 & 0.996 & 0.982 & 0.072 & 0.013 \\

& 0.846 & 0.730 & 0.038 & 0.036 & 0.626 & 0.302 & 0.022 & 0.027 \\ 
& 1.000 & 1.000 & 0.042 & 0.042 & 0.972 & 0.865 & 0.000 & 0.027 \\ 
& 0.967 & 0.906 & 0.058 & 0.052 & 0.712 & 0.410 & 0.036 & 0.033 \\ 
& 1.000 & 1.000 & 0.050 & 0.036 & 0.995 & 0.993 & 0.078 & 0.046 \\ 
\hline
\end{tabular}
\end{table}

tab <- matrix(c(0.504, 0.353, 0.034, 0.046, 0.372, 0.203, 0.032, 0.018,
                1.000, 0.990, 0.050, 0.046, 0.868, 0.399, 0.002, 0.021,
                0.681, 0.499, 0.050, 0.038, 0.398, 0.310, 0.036, 0.027, 
                1.000, 1.000, 0.034, 0.032, 0.996, 0.982, 0.072, 0.013,
                0.846, 0.730, 0.038, 0.036, 0.626, 0.302, 0.022, 0.027, 
                1.000, 1.000, 0.042, 0.042, 0.972, 0.865, 0.000, 0.027, 
                0.967, 0.906, 0.058, 0.052, 0.712, 0.410, 0.036, 0.033, 
                1.000, 1.000, 0.050, 0.036, 0.995, 0.993, 0.078, 0.046), 8, byrow = TRUE)
tab
xtable(tab[,c(1,3,2,4,5,7,6,8)], digits = 3)

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



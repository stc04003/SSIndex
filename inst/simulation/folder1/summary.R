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

sumSim <- function(n, model, frailty, type1 = FALSE) {
    if (model == "M2") {b0 <- c(.6, .8); r0 <- c(.6, .8)}
    if (model == "M3") {b0 <- -c(.6, .8); r0 <- -c(.6, .8)}
    if (model == "M4") {b0 <- -c(.6, .8); r0 <- c(.28, .96)}
    if (model == "M5") {b0 <- c(.6, .8); r0 <- c(.6, .8)}
    t0 <- c(b0, b0, r0, r0)
    fname <- paste(c("gany", n, model, frailty, type1), collapse = "-")
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
    ## First 4 rows are smooth
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
tab1[,1:2 * 4] <- 100 * tab1[,1:2 * 4]
tab2[,1:2 * 4] <- 100 * tab2[,1:2 * 4]

tab1[,1:2 * 4]
tab2[,1:2 * 4]

tab1
tab2

head(tab)
head(tab1)
head(tab2)
tail(tab1)
tail(tab2)
round(tab1, 4)
round(tab2, 4)

## Non smooth; use smooth as initial value (Table 1) in paper

sdid <- c(2:3, 6:7)
## sdid <- 2:5
tab1 <-
    rbind(cbind(sumSim(200, "M2", FALSE)[5:8, sdid], sumSim(200, "M2", TRUE)[5:8, sdid]),
          cbind(sumSim(200, "M3", FALSE)[5:8, sdid], sumSim(200, "M3", TRUE)[5:8, sdid]),
          cbind(sumSim(200, "M5", FALSE)[5:8, sdid], sumSim(200, "M5", TRUE)[5:8, sdid]),
          cbind(sumSim(200, "M4", FALSE)[5:8, sdid], sumSim(200, "M4", TRUE)[5:8, sdid]),
          cbind(sumSim(400, "M2", FALSE)[5:8, sdid], sumSim(400, "M2", TRUE)[5:8, sdid]),
          cbind(sumSim(400, "M3", FALSE)[5:8, sdid], sumSim(400, "M3", TRUE)[5:8, sdid]),
          cbind(sumSim(400, "M5", FALSE)[5:8, sdid], sumSim(400, "M5", TRUE)[5:8, sdid]),
          cbind(sumSim(400, "M4", FALSE)[5:8, sdid], sumSim(400, "M4", TRUE)[5:8, sdid]))
tab1 <- tab1 * 1000
tab1[,1:2 * 4] <- tab1[,1:2 * 4] / 10
print(xtable(tab1, digits = c(0, 0, 0, 0, 1, 0, 0, 0, 1)), include.rownames = FALSE)

## nonsmooth
## 200
7 & 212 & 202 & 91.0 & 9 & 204 & 209 & 92.7 \\ 
40 & 163 & 165 & 92.0 & 34 & 150 & 165 & 94.9 \\ 
1 & 46 & 58 & 96.9 & 10 & 104 & 118 & 96.7 \\ 
3 & 35 & 45 & 96.4 & 3 & 77 & 81 & 95.5 \\ 
9 & 91 & 101 & 96.0 & 14 & 106 & 111 & 96.4 \\ 
1 & 69 & 76 & 95.0 & 0 & 78 & 84 & 95.9 \\ 
1 & 62 & 73 & 96.7 & 5 & 115 & 129 & 95.7 \\ 
3 & 48 & 61 & 96.5 & 9 & 87 & 91 & 95.6 \\
46 & 341 & 325 & 90.8 & 48 & 333 & 331 & 91.8 \\ 
71 & 266 & 263 & 92.8 & 67 & 250 & 279 & 94.0 \\ 
1 & 53 & 61 & 96.2 & 12 & 122 & 129 & 94.5 \\ 
2 & 40 & 46 & 96.4 & 6 & 90 & 98 & 93.9 \\
1 & 39 & 46 & 95.7 & 1 & 43 & 52 & 95.9 \\ 
0 & 29 & 34 & 95.8 & 1 & 32 & 40 & 95.6 \\ 
2 & 39 & 48 & 96.9 & 7 & 98 & 108 & 96.3 \\ 
1 & 11 & 15 & 96.4 & 3 & 28 & 33 & 94.5 \\ 
## 400
2 & 153 & 150 & 92.6 & 9 & 160 & 162 & 93.3 \\ 
22 & 116 & 118 & 92.7 & 19 & 122 & 129 & 94.3 \\ 
0 & 31 & 37 & 96.1 & 0 & 80 & 86 & 95.3 \\ 
1 & 24 & 31 & 95.6 & 7 & 60 & 67 & 95.9 \\
3 & 66 & 69 & 95.2 & 4 & 72 & 81 & 95.3 \\ 
2 & 49 & 52 & 95.2 & 2 & 54 & 62 & 95.5 \\ 
2 & 47 & 54 & 96.1 & 6 & 81 & 91 & 95.7 \\ 
0 & 35 & 41 & 96.0 & 2 & 61 & 70 & 95.2 \\ 
23 & 245 & 247 & 91.6 & 10 & 246 & 249 & 92.7 \\ 
45 & 191 & 206 & 93.0 & 57 & 198 & 214 & 94.1 \\ 
0 & 37 & 42 & 95.7 & 0 & 84 & 89 & 95.4 \\ 
1 & 28 & 32 & 95.0 & 7 & 63 & 70 & 95.2 \\ 
0 & 27 & 31 & 95.1 & 1 & 28 & 34 & 95.3 \\ 
1 & 20 & 23 & 95.6 & 1 & 21 & 26 & 95.1 \\ 
2 & 28 & 33 & 95.6 & 3 & 67 & 76 & 96.3 \\ 
1 & 8 & 10 & 95.7 & 2 & 19 & 23 & 95.0 \\

## smooth
tab2 <-
    rbind(cbind(sumSim(200, "M2", FALSE)[1:4, sdid], sumSim(200, "M2", TRUE)[1:4, sdid]),
          cbind(sumSim(200, "M3", FALSE)[1:4, sdid], sumSim(200, "M3", TRUE)[1:4, sdid]),
          cbind(sumSim(200, "M5", FALSE)[1:4, sdid], sumSim(200, "M5", TRUE)[1:4, sdid]),
          cbind(sumSim(200, "M4", FALSE)[1:4, sdid], sumSim(200, "M4", TRUE)[1:4, sdid]),
          cbind(sumSim(400, "M2", FALSE)[1:4, sdid], sumSim(400, "M2", TRUE)[1:4, sdid]),
          cbind(sumSim(400, "M3", FALSE)[1:4, sdid], sumSim(400, "M3", TRUE)[1:4, sdid]),
          cbind(sumSim(400, "M5", FALSE)[1:4, sdid], sumSim(400, "M5", TRUE)[1:4, sdid]),
          cbind(sumSim(400, "M4", FALSE)[1:4, sdid], sumSim(400, "M4", TRUE)[1:4, sdid]))
tab2 <- tab2 * 1000
tab2[,1:2 * 4] <- tab2[,1:2 * 4] / 10

print(xtable(tab2, digits = c(0, 0, 0, 0, 1, 0, 0, 0, 1)), include.rownames = FALSE)

## 200 
6 & 211 & 202 & 90.9 & 8 & 204 & 219 & 91.8 \\ 
41 & 162 & 166 & 90.3 & 35 & 150 & 187 & 94.8 \\ 
1 & 42 & 51 & 96.6 & 7 & 101 & 122 & 96.0 \\ 
3 & 32 & 40 & 96.6 & 5 & 75 & 88 & 95.4 \\
10 & 91 & 100 & 96.3 & 14 & 106 & 120 & 96.2 \\ 
1 & 68 & 76 & 95.0 & 0 & 79 & 93 & 95.6 \\ 
1 & 59 & 71 & 96.5 & 5 & 110 & 124 & 95.3 \\ 
4 & 45 & 55 & 97.2 & 8 & 83 & 97 & 94.7 \\ 
55 & 342 & 322 & 90.9 & 47 & 334 & 327 & 91.0 \\ 
82 & 266 & 261 & 92.5 & 79 & 250 & 277 & 94.3 \\ 
1 & 49 & 55 & 95.9 & 8 & 115 & 125 & 95.8 \\ 
1 & 37 & 41 & 95.6 & 7 & 86 & 95 & 94.1 \\
1 & 38 & 45 & 95.9 & 1 & 42 & 52 & 97.1 \\ 
1 & 29 & 34 & 95.8 & 1 & 31 & 39 & 97.2 \\ 
1 & 35 & 39 & 96.7 & 5 & 91 & 99 & 95.3 \\ 
1 & 10 & 12 & 96.2 & 3 & 26 & 31 & 93.4 \\

1 & 151 & 149 & 91.4 & 7 & 158 & 162 & 92.5 \\ 
22 & 116 & 117 & 91.4 & 20 & 121 & 129 & 94.9 \\ 
0 & 29 & 36 & 97.0 & 1 & 78 & 84 & 95.1 \\ 
1 & 22 & 28 & 96.4 & 6 & 59 & 67 & 95.2 \\ 
4 & 66 & 69 & 95.0 & 5 & 72 & 81 & 96.6 \\
2 & 49 & 52 & 94.6 & 2 & 54 & 61 & 96.4 \\ 
2 & 44 & 50 & 96.0 & 4 & 79 & 88 & 96.5 \\ 
0 & 33 & 37 & 96.6 & 3 & 59 & 68 & 95.8 \\ 
22 & 244 & 245 & 91.7 & 9 & 245 & 248 & 91.1 \\ 
44 & 190 & 206 & 93.1 & 56 & 198 & 214 & 94.8 \\ 
0 & 35 & 38 & 95.8 & 2 & 82 & 89 & 95.2 \\ 
1 & 26 & 29 & 95.6 & 5 & 62 & 68 & 95.5 \\ 
0 & 26 & 30 & 96.4 & 1 & 28 & 34 & 98.7 \\ 
1 & 20 & 23 & 96.2 & 1 & 21 & 26 & 98.4 \\ 
2 & 25 & 27 & 94.7 & 2 & 63 & 69 & 95.8 \\ 
1 & 7 & 8 & 94.5 & 2 & 18 & 21 & 94.2 \\ 

e

sumPwr <- function(n, model, frailty, type1 = FALSE) {
    fname <- paste(c("gany", n, model, frailty, type1), collapse = "-")
    if (file.exists(fname)) dat <- matrix(c(t(matrix(scan(fname), 5))), 26)
    ## dat <- dat[c(1:8, 11:26, 9:10),]
    return(c(mean(dat[25,] > .95), mean(dat[26,] > .95)))
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

## Paper focus on M4 only

sumPwr(200, "M2", FALSE, TRUE)
sumPwr(200, "M2", TRUE, TRUE)
sumPwr(200, "M3", FALSE, TRUE)
sumPwr(200, "M3", TRUE, TRUE)
sumPwr(200, "M4", FALSE, TRUE)
sumPwr(200, "M4", TRUE, TRUE)
sumPwr(200, "M5", FALSE, TRUE)
sumPwr(200, "M5", TRUE, TRUE)

sumPwr(400, "M2", FALSE, TRUE)
sumPwr(400, "M2", TRUE, TRUE)
sumPwr(400, "M3", FALSE, TRUE)
sumPwr(400, "M3", TRUE, TRUE)

sumPwr(400, "M4", FALSE, TRUE)
sumPwr(400, "M4", TRUE, TRUE)
sumPwr(400, "M5", FALSE, TRUE)
sumPwr(400, "M5", TRUE, TRUE)

sumPwr <- function(n, model, frailty, type1 = FALSE) {
    fname <- paste(c("gany", n, model, frailty, type1), collapse = "-")
    if (file.exists(fname)) dat <- matrix(c(t(matrix(scan(fname), 5))), 26)
    return(c(mean(dat[25,] >= .95),
             mean(dat[26,] >= .95),
             mean(colSums(dat[25:26,] >= .975 ) >= 1)))
}

sumPwr(200, "M4", FALSE, FALSE)
sumPwr(200, "M4", FALSE, TRUE)
sumPwr(200, "M4", TRUE, FALSE)
sumPwr(200, "M4", TRUE, TRUE)

sumPwr(400, "M4", TRUE, TRUE)

c(sumPwr(200, "M4", FALSE, TRUE), sumPwr(400, "M4", FALSE, TRUE))
c(sumPwr(200, "M4", TRUE, TRUE), sumPwr(400, "M4", TRUE, TRUE))

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

## Summary counts
sumCount <- function(n, model, frailty) {
    dat <- simDat(n, model, frailty)
    print(quantile(unlist(lapply(split(dat$id, dat$id), length)) - 1, prob = 0:10 / 10))
    print(summary(unlist(lapply(split(dat$id, dat$id), length)) - 1))
    print(mean(unlist(lapply(split(dat$Time, dat$id), max)) == max(dat$Time)))
    ## print(sum(dat$status) / n)
}

sumCount(n, "M2", FALSE)
sumCount(n, "M2", TRUE)

sumCount(n, "M3", FALSE)
sumCount(n, "M3", TRUE)

sumCount(n, "M5", FALSE)
sumCount(n, "M5", TRUE)

sumCount(n, "M4", FALSE)
sumCount(n, "M4", TRUE)
 
library(reReg)
dat <- simDat(200, "M2", FALSE)
dat <- simDat(200, "M3", FALSE)
dat <- simDat(200, "M4", FALSE)
dat <- simDat(200, "M5", FALSE)
plot(with(dat, Recur(Time, id, event, status)))
    



library(GSM)
library(reReg)
## example

sumSim <- function(n, model, indB0 = "test") {
    fname <- paste(c("results", n, model, B), collapse = "-")
    if (file.exists(fname)) dat <- read.table(fname)
    else stop("file name does not exist")
    ## dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, nrow(dat) / 50))))
    dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, 13))))
    if (is.logical(indB0)) {
        pwr <- mean(dat[,3])
        if (indB0) {
            PE.b0 <- c(0, 0)
            ESE.b0 <- c(0, 0)
            ASE.b0 <- c(0, 0)
            PE.r0 <- eval(parse(text = paste("rowMeans(r0", n, model, ")", sep = "")))
            ESE.r0 <- eval(parse(text = paste("apply(r0", n, model, ", 1, sd)", sep = "")))
            ASE.r0 <- colMeans(dat[,12:13])
        } else {
            PE.b0 <- colMeans(dat[,1:2])
            ESE.b0 <- apply(dat[,1:2], 2, sd)
            ASE.b0 <- colMeans(dat[,6:7])
            PE.r0 <- colMeans(dat[,4:5])
            ESE.r0 <- apply(dat[,4:5], 2, sd)
            ASE.r0 <- colMeans(dat[,8:9])
        }
    } else {
        PE.b0 <- colMeans(dat[,1:2])
        ESE.b0 <- apply(dat[,1:2], 2, sd)
        ASE.b0 <- colMeans(dat[,6:7])
        PE.r0 <- colMeans(dat[,4:5])
        ESE.r0 <- apply(dat[,4:5], 2, sd) * (1 - mean(dat[,3])) +
            apply(eval(parse(text = paste("t(r0", n, model, ")", sep = ""))), 2, sd) * mean(dat[,3]) 
        ASE.r0 <- colMeans(dat[,3] * dat[,8:9] + (1 - dat[,3]) * dat[,12:13])
    }
    cbind(PE.b0, ESE.b0, ASE.b0, PE.r0, ESE.r0, ASE.r0)
}

sumPwr <- function(n, model) {
    fname <- paste(c("results", n, model, B), collapse = "-")
    if (file.exists(fname)) dat <- read.table(fname)
    else stop("file name does not exist")
    ## dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, nrow(dat) / 50))))
    dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, 13))))
    mean(dat[,3])
}

tabPwr <- rbind(c(sumPwr(50, "M1"), sumPwr(100, "M1"), sumPwr(200, "M1"), sumPwr(500, "M1")),
                c(sumPwr(50, "M2"), sumPwr(100, "M2"), sumPwr(200, "M2"), sumPwr(500, "M2")),
                c(sumPwr(50, "M3"), sumPwr(100, "M3"), sumPwr(200, "M3"), sumPwr(500, "M3")),
                c(sumPwr(50, "M4"), sumPwr(100, "M4"), sumPwr(200, "M4"), sumPwr(500, "M4")))

makeTab <- function(n) {
    cbind(rbind(sumSim(n, "M1", FALSE)[,1:3], sumSim(n, "M1", FALSE)[,4:6]),
          rbind(sumSim(n, "M1", TRUE)[,1:3], sumSim(n, "M1", TRUE)[,4:6]),
          rbind(sumSim(n, "M2", FALSE)[,1:3], sumSim(n, "M2", FALSE)[,4:6]),
          rbind(sumSim(n, "M3", FALSE)[,1:3], sumSim(n, "M3", FALSE)[,4:6]),
          rbind(sumSim(n, "M4", FALSE)[,1:3], sumSim(n, "M4", FALSE)[,4:6]))
}

tab <- rbind(makeTab(50), makeTab(100), makeTab(200), makeTab(500))

makeTab2 <- function(n) {
    rbind(rbind(sumSim(n, "M1", FALSE)[,1:3], sumSim(n, "M1", FALSE)[,4:6]),
          rbind(sumSim(n, "M1", TRUE)[,1:3], sumSim(n, "M1", TRUE)[,4:6]),
          rbind(sumSim(n, "M2", FALSE)[,1:3], sumSim(n, "M2", FALSE)[,4:6]),
          rbind(sumSim(n, "M3", FALSE)[,1:3], sumSim(n, "M3", FALSE)[,4:6]),
          rbind(sumSim(n, "M4", FALSE)[,1:3], sumSim(n, "M4", FALSE)[,4:6]))
}

tab <- cbind(makeTab2(50), makeTab2(100), makeTab2(200), makeTab2(500))

library(xtable)
print(xtable(tab, digits = 3), math.style.negative = TRUE, include.rownames=FALSE)

sumSim(100, "M3")
sumSim(500, "M3")
sumSim(1000, "M3")
sumSim(2000, "M3")


sumSim2 <- function(n, model, indB0 = "test") {
    fname <- paste(c("results", n, model, B), collapse = "-")
    if (file.exists(fname)) dat <- read.table(fname)
    else stop("file name does not exist")
    ## dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, nrow(dat) / 50))))
    dat <- t(do.call(cbind, lapply(dat, function(x) matrix(x, 13))))
    if (is.logical(indB0)) {
        pwr <- mean(dat[,3])
        if (indB0) {
            PE.b0 <- colMeans(dat[,1:2])
            ESE.b0 <- apply(dat[,1:2], 2, sd)
            ASE.b0 <- colMeans(dat[,6:7])
            PE.r0 <- colMeans(dat[,4:5])
            ESE.r0 <- apply(dat[,4:5], 2, sd)
            ASE.r0 <- colMeans(dat[,12:13])
        } else {
            PE.b0 <- colMeans(dat[,1:2])
            ESE.b0 <- apply(dat[,1:2], 2, sd)
            ASE.b0 <- colMeans(dat[,6:7])
            PE.r0 <- colMeans(dat[,4:5])
            ESE.r0 <- apply(dat[,4:5], 2, sd)
            ASE.r0 <- colMeans(dat[,8:9])
        }
    } else {
        PE.b0 <- colMeans(dat[,1:2])
        ESE.b0 <- apply(dat[,1:2], 2, sd)
        ASE.b0 <- colMeans(dat[,6:7])
        PE.r0 <- colMeans(dat[,4:5])
        ESE.r0 <- apply(dat[,4:5], 2, sd)
        ASE.r0 <- colMeans(dat[,3] * dat[,8:9] + (1 - dat[,3]) * dat[,12:13])
    }
    cbind(PE.b0, ESE.b0, ASE.b0, PE.r0, ESE.r0, ASE.r0)
}




####################################################################################3
## Parallel computing on 1000 replications
## March 10
####################################################################################3
library(parallel)
library(xtable)

do <- function(n, model, frailty = FALSE) {
    B <- 200
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    fit.indep <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
                     data = dat, shp.ind = TRUE)
    bi <- seq(0, 2 * pi, length = 100)
    k0 <- sapply(bi, function(x) getk0(dat, c(cos(x), sin(x))))
    mm <- aggregate(event ~ id, dat, sum)[, 2]
    n <- length(unique(dat$id))
    getBootk <- function(dat) {
        ind <- sample(1:n, replace = TRUE)
        dat0 <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        dat0$id <- rep(1:n, mm[ind] + 1)
        fit0 <- gsm(reSurv(time1 = t, id = id, event = event, status = status) ~ x1 + x2,
                    data = dat0, shp.ind = FALSE)
        fit1 <- gsm(reSurv(time1 = t, id = id, event = event, status = status) ~ x1 + x2,
                    data = dat0, shp.ind = TRUE)
        c(max(sapply(1:length(bi), function(x) getk0(dat0, c(cos(bi[x]), sin(bi[x]))) - k0[x])),
          fit0$b0, fit0$r0, fit1$b0, fit1$r0)
    }
    tmp <- replicate(B, getBootk(dat))
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta, (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0,
      1 * (max(k0) > quantile(tmp[1,], .95)), fit$r0, fit.indep$r0,
      sqrt(diag(var(t(tmp[2:3,])))),
      sqrt(diag(var(t(tmp[4:5,])))),
      ## sqrt(diag(var(t(tmp[6:7,])))),
      sqrt(diag(var(t(tmp[8:9,])))))
}

do2 <- function(n, model, B = 200, frailty = FALSE) {
    seed <- sample(1:1e7, 1)
    set.seed(seed)
    dat <- simDat(n, model, frailty)
    fit <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
               data = dat, shp.ind = FALSE)
    fit.indep <- gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
                     data = dat, shp.ind = TRUE)
    ## outputs are (1:2) \hat\beta, (3) reject H_0:\beta = 0? (1 = reject),
    ## (4:5) \hat\gamma,
    ## (6:7) \hat\gamma under independence,
    ## (8:9) bootstrap sd for \hat\beta, (10:11) bootstrap sd for \hat\gamma
    ## (12:13) bootstrap sd for \hat\gamma; these assumes indep.
    c(fit$b0, fit$r0, fit.indep$r0)
}

system.time(print(do(100, "M1", TRUE)))
do(50, "M4", FALSE)
do2(100, "M1", FALSE)
do2(50, "M4", FALSE)

cl <- makePSOCKcluster(8)
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, c('do', 'do2')))
invisible(clusterEvalQ(NULL, library(GSM)))
invisible(clusterEvalQ(NULL, library(reReg)))

runPara <- function(model, n, frailty) {
    obj <- paste(c(model, n, frailty), collapse = "")
    eval(parse(text = paste(obj, " <- NULL")))
    toRun <- paste(obj, " <- t(matrix(unlist(parLapply(NULL, 1:500, function(z) do(",
                   n, ",'", model, "',", frailty, "))), 13))", sep = "")
    ptm <- proc.time()
    eval(parse(text = toRun))
    print(proc.time() - ptm)
    fname <- paste(obj, ".RData", sep = "")
    eval(parse(text = paste("save(", obj, ",file = '", fname, "')", sep = "")))
}

runPara("M1", 50, TRUE)
runPara("M1", 100, TRUE)
runPara("M1", 200, TRUE)

runPara("M2", 50, TRUE)
runPara("M2", 100, TRUE)
runPara("M2", 200, TRUE)

runPara("M3", 50, TRUE)
runPara("M3", 100, TRUE)
runPara("M3", 200, TRUE)

runPara("M4", 50, TRUE)
runPara("M4", 100, TRUE)
runPara("M4", 200, TRUE)

runPara("M1", 50, FALSE)
runPara("M1", 100, FALSE)
runPara("M1", 200, FALSE)

runPara("M2", 50, FALSE)
runPara("M2", 100, FALSE)
runPara("M2", 200, FALSE)

runPara("M3", 50, FALSE)
runPara("M3", 100, FALSE)
runPara("M3", 200, FALSE)

runPara("M4", 50, FALSE)
runPara("M4", 100, FALSE)
runPara("M4", 200, FALSE)

stopCluster(cl)

makeTab2 <- function(dat) {
    cbind(colMeans(dat), apply(dat, 2, sd))
}

tab <- rbind(cbind(makeTab2(sim150), makeTab2(sim1100), makeTab2(sim1200)),
             cbind(makeTab2(sim250), makeTab2(sim2100), makeTab2(sim2200)),
             cbind(makeTab2(sim350), makeTab2(sim3100), makeTab2(sim3200)),
             cbind(makeTab2(sim450), makeTab2(sim4100), makeTab2(sim4200)),
             cbind(makeTab2(sim1502), makeTab2(sim11002), makeTab2(sim12002)),
             cbind(makeTab2(sim2502), makeTab2(sim21002), makeTab2(sim22002)),
             cbind(makeTab2(sim3502), makeTab2(sim31002), makeTab2(sim32002)),
             cbind(makeTab2(sim4502), makeTab2(sim41002), makeTab2(sim42002)))
             
print(xtable(tab, digits = 3), math.style.negative = TRUE, include.rownames = FALSE)

## X ~ truncated normal
## Without frailty
## True (0, 0), (0.28, 0.96)
0.005 & 0.692 & $-$0.038 & 0.681 & $-$0.023 & 0.694 \\ 
0.023 & 0.722 & $-$0.005 & 0.732 & $-$0.013 & 0.720 \\ 
0.253 & 0.130 & 0.262 & 0.087 & 0.265 & 0.060 \\ 
0.958 & 0.039 & 0.960 & 0.040 & 0.962 & 0.015 \\ 
0.251 & 0.132 & 0.261 & 0.078 & 0.264 & 0.053 \\ 
0.958 & 0.033 & 0.962 & 0.021 & 0.963 & 0.014 \\
## True (0.6, 0.8), (0.6, 0.8)
$-$0.542 & 0.346 & $-$0.580 & 0.214 & $-$0.596 & 0.122 \\ 
$-$0.707 & 0.296 & $-$0.765 & 0.179 & $-$0.788 & 0.094 \\ 
0.585 & 0.161 & 0.597 & 0.106 & 0.594 & 0.082 \\ 
0.786 & 0.119 & 0.788 & 0.104 & 0.796 & 0.083 \\ 
0.608 & 0.178 & 0.610 & 0.098 & 0.610 & 0.063 \\ 
0.759 & 0.151 & 0.782 & 0.082 & 0.788 & 0.055 \\
## True (0.6, 0.8), (0.6, 0.8)
$-$0.552 & 0.294 & $-$0.600 & 0.157 & $-$0.594 & 0.098 \\ 
$-$0.744 & 0.236 & $-$0.775 & 0.121 & $-$0.795 & 0.074 \\ 
$-$0.577 & 0.205 & $-$0.595 & 0.116 & $-$0.597 & 0.071 \\ 
$-$0.775 & 0.156 & $-$0.790 & 0.096 & $-$0.797 & 0.054 \\ 
$-$0.582 & 0.196 & $-$0.598 & 0.122 & $-$0.601 & 0.071 \\ 
$-$0.774 & 0.154 & $-$0.787 & 0.089 & $-$0.794 & 0.054 \\
## True (0.6, 0.8), (0.28, 0.96)
$-$0.535 & 0.344 & $-$0.570 & 0.223 & $-$0.592 & 0.113 \\ 
$-$0.701 & 0.322 & $-$0.771 & 0.174 & $-$0.793 & 0.090 \\ 
0.224 & 0.235 & 0.241 & 0.159 & 0.264 & 0.102 \\ 
0.942 & 0.088 & 0.956 & 0.041 & 0.958 & 0.043 \\ 
0.224 & 0.224 & 0.247 & 0.144 & 0.266 & 0.090 \\ 
0.946 & 0.071 & 0.957 & 0.046 & 0.959 & 0.025 \\

## With frailty
## True (0, 0), (0.28, 0.96)
$-$0.029 & 0.703 & $-$0.027 & 0.700 & $-$0.003 & 0.691 \\ 
$-$0.005 & 0.712 & $-$0.020 & 0.714 & $-$0.011 & 0.724 \\ 
0.253 & 0.136 & 0.261 & 0.078 & 0.265 & 0.064 \\ 
0.957 & 0.046 & 0.962 & 0.022 & 0.962 & 0.015 \\ 
0.250 & 0.124 & 0.259 & 0.081 & 0.263 & 0.054 \\ 
0.960 & 0.031 & 0.962 & 0.022 & 0.963 & 0.014 \\ 
## True (0.6, 0.8), (0.6, 0.8)
$-$0.527 & 0.348 & $-$0.587 & 0.208 & $-$0.596 & 0.125 \\ 
$-$0.716 & 0.299 & $-$0.765 & 0.166 & $-$0.787 & 0.099 \\ 
0.588 & 0.147 & 0.590 & 0.126 & 0.590 & 0.083 \\ 
0.783 & 0.139 & 0.789 & 0.115 & 0.797 & 0.095 \\ 
0.596 & 0.189 & 0.613 & 0.098 & 0.609 & 0.067 \\ 
0.768 & 0.142 & 0.779 & 0.092 & 0.787 & 0.069 \\
## True (0.6, 0.8), (0.6, 0.8)
$-$0.563 & 0.313 & $-$0.602 & 0.162 & $-$0.595 & 0.098 \\ 
$-$0.707 & 0.292 & $-$0.769 & 0.139 & $-$0.795 & 0.073 \\ 
$-$0.579 & 0.202 & $-$0.596 & 0.116 & $-$0.597 & 0.076 \\ 
$-$0.771 & 0.171 & $-$0.789 & 0.099 & $-$0.797 & 0.057 \\ 
$-$0.580 & 0.200 & $-$0.599 & 0.113 & $-$0.601 & 0.074 \\ 
$-$0.775 & 0.149 & $-$0.787 & 0.093 & $-$0.794 & 0.056 \\
## True (0.6, 0.8), (0.28, 0.96)
$-$0.540 & 0.363 & $-$0.573 & 0.219 & $-$0.591 & 0.111 \\ 
$-$0.690 & 0.318 & $-$0.771 & 0.175 & $-$0.795 & 0.084 \\ 
0.225 & 0.246 & 0.243 & 0.153 & 0.262 & 0.112 \\ 
0.936 & 0.115 & 0.957 & 0.040 & 0.957 & 0.058 \\ 
0.227 & 0.233 & 0.247 & 0.141 & 0.267 & 0.087 \\ 
0.941 & 0.088 & 0.958 & 0.035 & 0.959 & 0.023 \\ 


makeTab2(sim150)
makeTab2(sim250)
makeTab2(sim350)
makeTab2(sim450)

makeTab2(sim1100)
makeTab2(sim2100)
makeTab2(sim3100)
makeTab2(sim4100)


set.seed(2)
do2(100, "M2", FALSE)

tmp <- matrix(NA, 200, 6)
for (i in 1:200) {
    set.seed(i)
    tmp[i,] <- do2(100, "M2", FALSE)
}



debug(gsm)
undebug(gsm)

debug(getb0)
undebug(getb0)

set.seed(6822898)
dat <- simDat(100, "M2", FALSE)
gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
    data = dat, shp.ind = FALSE)


spg(par = acos(1 / sqrt(2)), fn = Cn, quiet = TRUE, control = list(trace = FALSE))
spg(par = double(p - 1), fn = Cn, quiet = TRUE, control = list(trace = FALSE))
optimize(f = Cn, interval = c(-10, 10))


plot(seq(0, 2 * pi, .01), sapply(seq(0, 2 * pi, .01), Cn), 'l')


acos(-.6) ## 2.214
abline(v = acos(-.6))
optimize(f = Cn, interval = c(4, 5))
abline(v = 4.62165)
c(cos(4.62165), sin(4.62165))


set.seed(21)
dat <- simDat(100, "M2", FALSE)

gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
    data = dat, shp.ind = FALSE)
## -0.6355681 -0.7720448


Cn(0)
Cn2(0)

Cn(2)
Cn2(2)

Cn(2.02)
Cn2(2.02)
replicate(100, Cn2(2.47236))
replicate(100, Cn2(1))

x0 <- seq(0, 2 * pi, .01)
y1 <- sapply(seq(0, 2 * pi, .01), Cn)
y2 <- sapply(seq(0, 2 * pi, .01), Cn2)

plot(x0, y1, 'l')
lines(x0, y2, 'l', col = 2)

plot(x0, y2, 'l')
lines(x0, y1, 'l', col = 2)

tab <- tibble(x = rep(x0, 2), y = c(y1, y2),
              type = c(rep("Unsmoothed", length(x0)), rep("Smoothed", length(x0))))

ggplot(data = tab, aes(x = x, y = y, color = type)) +
    geom_line(lwd = 2, alpha = 0.5) +
    ggtitle("M2, with n = 100") + xlab(quote(theta)) + ylab(quote(U[1](theta)))
ggsave("rank.pdf")

library(GSM)
library(reReg)
set.seed(21)
dat <- simDat(100, "M2", FALSE)
## debug(getb0)
gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
    data = dat, shp.ind = FALSE)
## -0.6355681 -0.7720448



abline(v = spg(par = acos(1 / sqrt(2)), fn = Cn, quiet = TRUE, control = list(trace = FALSE))$par)
abline(v = spg(par = double(p - 1), fn = Cn, quiet = TRUE, control = list(trace = FALSE))$minimum)
abline(v = optimize(f = Cn, interval = c(-10, 10))$minimum)
abline(v = optimize(f = Cn, interval = c(0, 10))$minimum)

abline(v = spg(par = acos(1 / sqrt(2)), fn = Cn2, quiet = TRUE, control = list(trace = FALSE))$par, col = 3)
abline(v = spg(par = double(p - 1), fn = Cn2, quiet = TRUE, control = list(trace = FALSE))$par, col = 3)
abline(v = optimize(f = Cn2, interval = c(-10, 10))$minimum, col = 3)
abline(v = optimize(f = Cn2, interval = c(0, 10))$minimum, col = 3)

spg(par = acos(1 / sqrt(2)), fn = Cn2, quiet = TRUE, control = list(trace = FALSE))$par %% (2 * pi)
spg(par = double(p - 1), fn = Cn2, quiet = TRUE, control = list(trace = FALSE))$par %% (2 * pi)
optimize(f = Cn2, interval = c(-10, 10))$minimum %% (2 * pi)
optimize(f = Cn2, interval = c(0, 10))$minimum %% (2 * pi)



cos(spg(par = acos(1 / sqrt(2)), fn = Cn2, quiet = TRUE, control = list(trace = FALSE))$par)
cos(spg(par = double(p - 1), fn = Cn2, quiet = TRUE, control = list(trace = FALSE))$par)
cos(optimize(f = Cn2, interval = c(-10, 10))$minimum)
cos(optimize(f = Cn2, interval = c(0, 10))$minimum)


f1 <- t(matrix(unlist(parLapply(NULL, 1:200, function(z) do2(100, "M2", FALSE))), 6))
f2 <- t(matrix(unlist(parLapply(NULL, 1:200, function(z) do2(200, "M2", FALSE))), 6))
f3 <- t(matrix(unlist(parLapply(NULL, 1:200, function(z) do2(500, "M2", FALSE))), 6))

f12 <- t(matrix(unlist(parLapply(NULL, 1:100, function(z) do2(100, "M3", FALSE))), 6))
f22 <- t(matrix(unlist(parLapply(NULL, 1:100, function(z) do2(200, "M3", FALSE))), 6))
f32 <- t(matrix(unlist(parLapply(NULL, 1:100, function(z) do2(500, "M3", FALSE))), 6))


xtable(cbind(colMeans(f1)[1:4], apply(f1, 2, sd)[1:4],
             colMeans(f2)[1:4], apply(f2, 2, sd)[1:4],
             colMeans(f3)[1:4], apply(f3, 2, sd)[1:4]), digits = 3)



xtable(cbind(colMeans(f12)[1:4], apply(f12, 2, sd)[1:4],
             colMeans(f22)[1:4], apply(f22, 2, sd)[1:4],
             colMeans(f32)[1:4], apply(f32, 2, sd)[1:4]), digits = 3)

 




set.seed(21)
dat <- simDat(100, "M2", FALSE)
gsm(reSurv(time1 = t, id = id, event = event, status =  status) ~ x1 + x2,
    data = dat, shp.ind = FALSE)
## -0.6355681 -0.7720448

debug(gsm)

library(grid)
library(gridExtra)
library(SSIndex)
library(ggplot2)

n <- 200
n <- 500
model <- "M5" ## This is M3 in the paper
frailty <- TRUE
dat <- simDat(n, model, frailty)

fit <- gsm(reSurv(time1 = Time, id = id, event = event, status =  status) ~ x1 + x2,
           data = dat, shp.ind = FALSE, B = 0)

i <- 10
ggplot(data.frame(Time = fit$xb[order(fit$xb)], Fhat = fit$Fhat0[[i]][order(fit$xb)]),
       aes(x = Time, y = Fhat)) + geom_step()
summary(fit$xb)

## True fhat
f <- function(t, a) {
    (10 * (1 + t)^(exp(a) / 5) - 10) / (10 * 11^(exp(a) / 5) - 10)
}
    
tt <- subset(dat, event == 0)$Time

i  <- 36
plot(fit$xb, fhat[i,], cex = .5, pch = 19, main = i)
points(fit$xb, f(tt[i], fit$xb), col = 2, cex = .5)

fhat <- do.call(rbind, fit$Fhat0)
for (i in 1:length(fit$Fhat0)) {
    plot(fit$xb, fhat[i,], cex = .5, pch = 19, main = i)
    points(fit$xb, f(tt[i], fit$xb), col = 2, cex = .5)
    ## plot(dat00$Time, fhat[,i], cex = .5, pch = 19)
    Sys.sleep(.5)
}

## Plot F(t, a) for different t

plot(subset(dat, event == 0)$Time, fhat[,1], cex = 0.5, pch = 19)

for (i in 1:length(fit$Fhat0)) {
    plot(subset(dat, event == 0)$Time, fhat[,i], cex = 0.5, pch = 19, main = i)
    ## plot(dat00$Time, fhat[,i], cex = .5, pch = 19)
    Sys.sleep(.5)
}

fhat[1,]

x <- subset(dat, event == 0)$Time
y <- fit$xb
fhat3d <- fhat[order(x),]
fhat3d <- fhat[,order(y)]

persp(x[order(x)], y[order(y)], fhat3d)


persp(x[order(x)][1:75], y[order(y)][1:75], fhat3d[1:75, 1:75])

persp(x[order(x)][1:75], y[order(y)][1:75], fhat3d[1:75, 1:75], theta = 30, phi = 45)

#' dx = a - xb
kh <- function(dx) {
    return(ifelse(dx <= 1 && dx >= -1, .75 * (1 - dx^2), 0))
}

tmp <- matrix(sapply(c(outer(fit$xb, fit$xb, "-")), kh), length(fit$x))


## ########################################################################################
## Debug area
## ########################################################################################

library(tidyverse)
library(SSIndex)

## set.seed(0)
n <- 1000
model <- "M5" ## This is M3 in the paper
frailty <- TRUE
dat <- simDat(n, model, frailty)
dat

m <- mm <- aggregate(event ~ id, data = dat, sum)[,2]
tij <- subset(dat, event == 1)$Time
yi <-  subset(dat, event == 0)$Time
midx <- c(0, cumsum(mm)[-length(mm)])
X <- as.matrix(cbind(x1 = dat$x1, x2 = dat$x2))[dat$event == 0,]
p <- ncol(X)
xb <- X %*% c(.6, .8)
h <- 2.78 * sd(xb) * n^-.25
## h <- 2.78 * sd(xb) * n^-.25


.C("shapeFun", 
   as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), 
   as.double(yi), as.double(xb), as.double(xb[1]), as.double(5), 
   as.double(h), result = double(1), PACKAGE = "SSIndex")$result

yy <- 2.5
ff <- unlist(mapply(FUN = function(x, y)
    .C("shapeFun", 
       as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), 
       as.double(yi), as.double(unique(xb)), as.double(x), as.double(y), 
       as.double(h), result = double(1), PACKAGE = "SSIndex")$result, 
    xb, rep(yy, n)))
ff <- exp(-ff)

f <- function(t, a) {
    (10 * (1 + t)^(exp(a) / 5) - 10) / (10 * 11^(exp(a) / 5) - 10)
}

plot(xb, ff, cex = .5, pch = 19, ylim = c(0, 1))
points(xb, f(yy, xb), col = 2, cex = .5)

makeGG <- function(yy) {
    ff <- unlist(mapply(FUN = function(x, y)
        .C("shapeFun", 
           as.integer(n), as.integer(mm), as.integer(midx), as.double(tij), 
           as.double(yi), as.double(unique(xb)), as.double(x), as.double(y), 
           as.double(h), result = double(1), PACKAGE = "SSIndex")$result, 
        xb, rep(yy, n)))
    ff <- exp(-ff)
    datmp <- data.frame(xb = rep(xb, 2), Fhat = c(ff, f(yy, xb)),
                        type = rep(c("estimated", "truth"), each = n))
    datmp
}



p1 <- ggplot(makeGG(1), aes(x = xb, y = Fhat, col = type)) + geom_line()
p2 <- ggplot(makeGG(2), aes(x = xb, y = Fhat, col = type)) + geom_line()
p3 <- ggplot(makeGG(3), aes(x = xb, y = Fhat, col = type)) + geom_line()
p4 <- ggplot(makeGG(4), aes(x = xb, y = Fhat, col = type)) + geom_line()
p5 <- ggplot(makeGG(5), aes(x = xb, y = Fhat, col = type)) + geom_line()
p6 <- ggplot(makeGG(6), aes(x = xb, y = Fhat, col = type)) + geom_line()
p7 <- ggplot(makeGG(7), aes(x = xb, y = Fhat, col = type)) + geom_line()
p8 <- ggplot(makeGG(8), aes(x = xb, y = Fhat, col = type)) + geom_line()
p9 <- ggplot(makeGG(9), aes(x = xb, y = Fhat, col = type)) + geom_line()

p1 <- p1 + theme(legend.position = "none") + labs(title = "t = 1")
p2 <- p2 + theme(legend.position = "none") + labs(title = "t = 2")
p3 <- p3 + theme(legend.position = "none") + labs(title = "t = 3")
p4 <- p4 + theme(legend.position = "none") + labs(title = "t = 4")
p5 <- p5 + theme(legend.position = "none") + labs(title = "t = 5")
p6 <- p6 + theme(legend.position = "none") + labs(title = "t = 6")
p7 <- p7 + theme(legend.position = "none") + labs(title = "t = 7")
p8 <- p8 + theme(legend.position = "none") + labs(title = "t = 8")
p9 <- p9 + theme(legend.position = "none") + labs(title = "t = 9")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)
## ggsave("fhat.pdf")

library(GSM)
library(reReg)
library(frailtyHL)
library(parallel)

data(cgd)
## ID 87 ends with status = 1
dat0 <- rbind(cgd, cgd[cgd$id == 87,][2,])
dat0$status[nrow(dat0)] <- 0
dat0 <- dat0[order(dat0$id),]
dat0$event <- dat0$status
dat0$status <- 0
dat0$m <- rep(aggregate(event ~ id, dat0, sum)[,2], aggregate(event ~ id, dat0, sum)[,2] + 1)
rownames(dat0) <- NULL
names(dat0)

dat0 <- subset(dat0, select = c(id, tstop, event, status, m,
                                treat, inherit, steroids, propylac,
                                sex, age, height, weight))
dat0$treat <- as.numeric(dat0$treat) - 1
dat0$inherit <- as.numeric(dat0$inherit) - 1
dat0$sex <- as.numeric(dat0$sex) - 1
dat0$id <- rep(1:length(unique(dat0$id)), aggregate(event ~ id, dat0, sum)[,2] + 1)
dat0$age0 <- as.numeric(scale(dat0$age))
names(dat0)[2] <- "Time"

head(dat0)
str(dat0)
                        
str(gsm(reSurv(Time, id, event, status) ~ treat + propylac, data = dat0))
str(gsm(reSurv(Time, id, event, status) ~ treat + propylac + inherit, data = dat0))
str(gsm(reSurv(Time, id, event, status) ~ treat + propylac + sex, data = dat0))


##########################################################################################
## Functions
##########################################################################################
fname <- reSurv(Time, id, event, status) ~ treat + propylac + inherit


pVal <- function(fname, B = 100, dat0 = dat0) {
    xNames <- attr(terms(fname), "term.labels")
    p <- length(attr(terms(fname), "term.labels"))
    fit <- gsm(fname, data = dat0)
    ## str(fit)
    dat1 <- dat0
    xCol <- as.numeric(sapply(xNames, function(x) which(names(dat0) == x)))
    colnames(dat1)[xCol] <- paste("x", 1:p, sep = "")
    dat1 <- dat1[,c(1:5, xCol)]
    ## head(dat1)
    bi <- as.matrix(expand.grid(rep(list(seq(0, 2 * pi, length = 100)), p - 1)))
    system.time(k0 <- sapply(1:NROW(bi), function(x)
        getk0(dat1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1))))
    system.time(k02 <- sapply(1:NROW(bi), function(x)
        getk02(dat1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1), fit$Fhat0)))
    getBootK <- function(dat) {
        n <- length(unique(dat$id))    
        mm <- aggregate(event ~ id, dat, length)[, 2]
        ind <- sample(1:n, n, TRUE)
        datB <- dat[unlist(sapply(ind, function(x) which(dat$id %in% x))),]
        datB$id <- rep(1:n, mm[ind])
        rownames(datB) <- NULL
        fitB <- gsm(fname, dat = datB)
        datB1 <- datB
        xCol <- as.numeric(sapply(xNames, function(x) which(names(datB) == x)))
        colnames(datB1)[xCol] <- paste("x", 1:p, sep = "")
        datB1 <- datB1[,c(1:5, xCol)]
        kb <- max(sapply(1:NROW(bi), function(x)
            getk0(datB1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1)) - k0[x]))
        kb2 <- max(sapply(1:NROW(bi), function(x)
            getk02(datB1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1), fitB$Fhat0) - k02[x]))
        c(max(kb), max(kb2),
          fitB$b0, fitB$b00, fitB$r0, fitB$r00)
    }
    cl <- makePSOCKcluster(8)
    ## cl <- makePSOCKcluster(16)
    setDefaultCluster(cl)
    invisible(clusterExport(cl, c("bi", "k0", "k02", "fname", "dat0", "xNames", "p", "getBootK"),
                            environment()))
    invisible(clusterEvalQ(NULL, library(GSM)))
    invisible(clusterEvalQ(NULL, library(reReg)))
    system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootK(dat0))) 
    stopCluster(cl)
    c(mean(max(k0) > tmp[1,]), mean(max(k02) > tmp[2,]))
    ## list(k0 = k0, k02 = k02, tmp = tmp)
}

foo <- pVal(reSurv(Time, id, event, status) ~ treat + propylac + inherit, dat0 = dat0)

system.time(f1 <- pVal(reSurv(Time, id, event, status) ~ treat + propylac + age0, dat0 = dat0))
system.time(f2 <- pVal(reSurv(Time, id, event, status) ~ treat + propylac + sex, dat0 = dat0))
system.time(f3 <- pVal(reSurv(Time, id, event, status) ~ treat + inherit, dat0 = dat0))
system.time(f4 <- pVal(reSurv(Time, id, event, status) ~ treat + propylac, dat0 = dat0))
system.time(f5 <- pVal(reSurv(Time, id, event, status) ~ treat + sex + inherit, dat0 = dat0))
system.time(f6 <- pVal(reSurv(Time, id, event, status) ~ treat + age0 + inherit, dat0 = dat0))
system.time(f7 <- pVal(reSurv(Time, id, event, status) ~ treat + sex, dat0 = dat0))
system.time(f8 <- pVal(reSurv(Time, id, event, status) ~ treat + age0, dat0 = dat0))



pValShape <- function(fname, B = 100, dat0 = dat0) {
    xNames <- attr(terms(fname), "term.labels")
    p <- length(attr(terms(fname), "term.labels"))
    dat1 <- dat0
    xCol <- as.numeric(sapply(xNames, function(x) which(names(dat0) == x)))
    colnames(dat1)[xCol] <- paste("x", 1:p, sep = "")
    dat1 <- dat1[,c(1:5, xCol)]
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
        rownames(datB) <- NULL
        datB1 <- datB
        xCol <- as.numeric(sapply(xNames, function(x) which(names(datB) == x)))
        colnames(datB1)[xCol] <- paste("x", 1:p, sep = "")
        datB1 <- datB1[,c(1:5, xCol)]
        kb <- max(sapply(1:NROW(bi), function(x) getk0(datB1, bi[x,]) - k0[x]))
            ## getk0(datB1, cumprod(c(1, sin(bi[x,]))) * c(cos(bi[x,]), 1)) - k0[x]))
        max(kb)
    }
    cl <- makePSOCKcluster(8)
    ## cl <- makePSOCKcluster(16)
    setDefaultCluster(cl)
    invisible(clusterExport(cl, c("bi", "k0", "fname", "dat0", "xNames", "p", "getBootK"),
                            environment()))
    invisible(clusterEvalQ(NULL, library(GSM)))
    invisible(clusterEvalQ(NULL, library(reReg)))
    system.time(tmp <- parSapply(NULL, 1:B, function(z) getBootK(dat0))) 
    stopCluster(cl)
    mean(max(k0) > tmp)
}

system.time(f1 <- pValShape(reSurv(Time, id, event, status) ~ treat + propylac + age0, dat0 = dat0)); print(f1) # 0.21
system.time(f2 <- pValShape(reSurv(Time, id, event, status) ~ treat + propylac + sex, dat0 = dat0)); print(f2) # 0.51
system.time(f3 <- pValShape(reSurv(Time, id, event, status) ~ treat + inherit, dat0 = dat0)); print(f3) # 0.45
system.time(f4 <- pValShape(reSurv(Time, id, event, status) ~ treat + propylac, dat0 = dat0)); print(f4) # 0.4
system.time(f5 <- pValShape(reSurv(Time, id, event, status) ~ treat + sex + inherit, dat0 = dat0)); print(f5) # 0.31
system.time(f6 <- pValShape(reSurv(Time, id, event, status) ~ treat + age0 + inherit, dat0 = dat0)); print(f6) # 0.27
system.time(f7 <- pValShape(reSurv(Time, id, event, status) ~ treat + sex, dat0 = dat0)); print(f7) # 0.59
system.time(f8 <- pValShape(reSurv(Time, id, event, status) ~ treat + age0, dat0 = dat0)); print(f8) # 0.35
system.time(f8 <- pValShape(reSurv(Time, id, event, status) ~ treat + height + weight, dat0 = dat0)); print(f8) # 0.06
system.time(f8 <- pValShape(reSurv(Time, id, event, status) ~ treat + propylac + age, dat0 = dat0)); print(f8) # 0.26
system.time(f8 <- pValShape(reSurv(Time, id, event, status) ~ treat + sex + propylac, dat0 = dat0)); print(f8) # 0.5
system.time(f8 <- pValShape(reSurv(Time, id, event, status) ~ treat + sex + propylac + age0, dat0 = dat0)); print(f8) # 0.

system.time(f2 <- pValShape(reSurv(Time, id, event, status) ~ treat + propylac + sex + age, dat0 = dat0)); print(f2) # 0.14
system.time(f3 <- pValShape(reSurv(Time, id, event, status) ~ treat + propylac + sex + height + weight, dat0 = dat0)); print(f3)
system.time(f4 <- pValShape(reSurv(Time, id, event, status) ~ treat + propylac + sex + height + weight + age, dat0 = dat0)); print(f4)

## (treat,inherit,age,height,weight,steroids,prophylactic,sex,hosp1-3

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
        getk02(dat1, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fit$Fhat0)))
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
            getk0(datB1, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1))) - k0[x]))
        kb2 <- max(sapply(1:NROW(bi), function(x)
            getk02(datB1, cumprod(c(1, sin(bi[x,])) * c(cos(bi[x,]), 1)), fitB$Fhat0) - k02[x]))
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
}

pVal(reSurv(Time, id, event, status) ~ treat + propylac + inherit, B = 100, dat0 = dat0)
pVal(reSurv(Time, id, event, status) ~ treat + propylac + age0, B = 100, dat0 = dat0)
pVal(reSurv(Time, id, event, status) ~ treat + propylac + sex, B = 100, dat0 = dat0)
pVal(reSurv(Time, id, event, status) ~ treat + steroids + inherit, B = 100, dat0 = dat0)

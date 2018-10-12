library(GSM)

dim(simDat(10, "M1"))
dim(simDat(10, "M2"))
dim(simDat(10, "M3"))
dim(simDat(10, "M4"))

c(7, 24) / 25

do <- function(n, model) {
    dat <- simDat(n, model)
    unlist(gsm(dat))
}

sim1 <- t(replicate(100, do(300, "M1")))
sim2 <- t(replicate(100, do(300, "M2")))
sim3 <- t(replicate(100, do(300, "M3")))
sim4 <- t(replicate(100, do(300, "M4")))

summary(sim1)
summary(sim2)
summary(sim3)
summary(sim4)

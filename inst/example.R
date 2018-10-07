library(GSM)


set.seed(1)
dat <- simDat(100, "M1")

gsm(dat)

## true values are b0 = c(.6, .8), g0 = c(7, 24) / 25 = (.28, .96)

do <- function(n, model) {
    dat <- simDat(n, model)
    unlist(gsm(dat))
}

foo <- replicate(100, do(300, "M1"))
summary(t(foo))

foo <- replicate(100, do(300, "M2"))
summary(t(foo))

foo <- replicate(100, do(300, "M3"))
summary(t(foo))

foo <- replicate(100, do(300, "M4"))
summary(t(foo))

foo <- replicate(50, do(100, "M5"))
summary(t(foo))



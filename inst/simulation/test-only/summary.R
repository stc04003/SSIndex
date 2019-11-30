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

sumPwr <- function(n, model, frailty, type1 = FALSE) {
    fname <- paste(c("gany", n, model, frailty, "FALSE"), collapse = "-")
    if (file.exists(fname)) dat <- matrix(unlist(read.table(fname)), 2)
    ## dat <- matrix(c(t(matrix(scan(fname), 5))), 2)
    return(c(mean(dat[1,] > .95),
             mean(dat[2,] > .95),
             mean(colSums(dat > .975 ) >= 1)))
}

n <- 200
model <- "M2"
frailty <- FALSE
type1 = FALSE



xtable(rbind(sumPwr(200, "M3", FALSE),
             sumPwr(400, "M3", FALSE), 
             sumPwr(200, "M4", FALSE),
             sumPwr(400, "M4", FALSE)), digits = 3)

xtable(rbind(sumPwr(200, "M3", TRUE),
             sumPwr(400, "M3", TRUE), 
             sumPwr(200, "M4", TRUE),
             sumPwr(400, "M4", TRUE)), digits = 3)


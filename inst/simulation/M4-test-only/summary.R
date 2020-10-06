##################################################################################################
## Loading packages and functions
## Note: RData has corrected hypothesis
## gsm-xxx don't have correct hypothesis results
## Type 1 error, ran on Nov 15, 2019
## Inserted a small \tau_0 but didn't adjust for ties in bootstrap
##################################################################################################

library(SSIndex)
library(xtable)
library(ggplot2)

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
## With 2 offsets
##################################################################################################

sumPwr <- function(n, frailty, offset1, offset2) {
    fname <- paste(c("gany", n, frailty, as.character(offset1), as.character(offset2)),
                   collapse = "-")
    if (file.exists(fname)) {
        dat <- matrix(unlist(read.table(fname)), 2)
        return(c(mean(dat[1,] > .95),
                 mean(dat[2,] > .95),
                 mean(colSums(dat > .975 ) >= 1)))
    } else {return(rep(NA, 3))}
}

sumPwr(200, FALSE, .1, .2)
sumPwr(200, TRUE, .1, .2)
sumPwr(200, TRUE, .2, .2)

n <- 200
frailty <- FALSE
offset1 <- 1
offset2 <- .8

sumPwr(n, frailty, offset1, offset2)

m <- expand.grid(1:10/10, 1:10/10)
## ## Independent
## m$shape <- mapply(FUN = function(x, y) sumPwr(200, FALSE, x, y)[1], m[,1], m[,2])
## m$size <- mapply(FUN = function(x, y) sumPwr(200, FALSE, x, y)[2], m[,1], m[,2])
## m$rate <- mapply(FUN = function(x, y) sumPwr(200, FALSE, x, y)[3], m[,1], m[,2])

## Informative
m$shape <- mapply(FUN = function(x, y) sumPwr(200, TRUE, x, y)[1], m[,1], m[,2])
m$size <- mapply(FUN = function(x, y) sumPwr(200, TRUE, x, y)[2], m[,1], m[,2])
m$rate <- mapply(FUN = function(x, y) sumPwr(200, TRUE, x, y)[3], m[,1], m[,2])

names(m)[1:2] <- c("beta", "gamma")
attr(m, "out.attrs") <- NULL
str(m)
head(m)

ggplot(m, aes(beta, gamma)) +
    geom_tile(aes(fill = shape), color = "white") +
    scale_fill_gradient("Rejection \nproportion", low = "gray", high = "gray50", limits = c(0,1)) +
    xlab(expression(paste("||", beta, "||"))) +
    ylab(expression(paste("||", gamma, "||"))) + ggtitle("Shape test") +
    scale_y_continuous(breaks = seq(0, 1, by = .1)) +
    scale_x_continuous(breaks = seq(0, 1, by = .1)) 
## ggsave("shape.pdf")


ggplot(m, aes(beta, gamma)) +
    geom_tile(aes(fill = size), color = "white") +
    scale_fill_gradient("Rejection \nproportion", low = "gray", high = "gray50", limits = c(0,1)) +
    xlab(expression(paste("||", beta, "||"))) +
    ylab(expression(paste("||", gamma, "||"))) + ggtitle("Size test") +
    scale_y_continuous(breaks = seq(0, 1, by = .1)) +
    scale_x_continuous(breaks = seq(0, 1, by = .1)) 
## ggsave("size.pdf")

ggplot(m, aes(beta, gamma)) +
    geom_tile(aes(fill = rate), color = "white") +
    scale_fill_gradient("Rejection \nproportion", low = "gray", high = "gray50", limits = c(0,1)) +
    xlab(expression(paste("||", beta, "||"))) +
    ylab(expression(paste("||", gamma, "||"))) + ggtitle("Rate test") +
    scale_y_continuous(breaks = seq(0, 1, by = .1)) +
    scale_x_continuous(breaks = seq(0, 1, by = .1)) 
## ggsave("rate.pdf")

p1 <- ggplot(m, aes(beta, gamma)) +
    geom_tile(aes(fill = shape), color = "white") +
    scale_fill_gradient("Rejection \nproportion", low = "gray", high = "gray50", limits = c(0,1)) +
    xlab(expression(paste("||", beta, "||"))) +
    ylab(expression(paste("||", gamma, "||"))) + ggtitle("Shape test") +
    scale_y_continuous(breaks = seq(0, 1, by = .1)) +
    scale_x_continuous(breaks = seq(0, 1, by = .1)) 
p2 <- ggplot(m, aes(beta, gamma)) +
    geom_tile(aes(fill = size), color = "white") +
    scale_fill_gradient("Rejection \nproportion", low = "gray", high = "gray50", limits = c(0,1)) +
    xlab(expression(paste("||", beta, "||"))) +
    ylab(expression(paste("||", gamma, "||"))) + ggtitle("Size test") +
    scale_y_continuous(breaks = seq(0, 1, by = .1)) +
    scale_x_continuous(breaks = seq(0, 1, by = .1)) 
p3 <- ggplot(m, aes(beta, gamma)) +
    geom_tile(aes(fill = rate), color = "white") +
    scale_fill_gradient("Rejection \nproportion", low = "gray", high = "gray50", limits = c(0,1)) +
    xlab(expression(paste("||", beta, "||"))) +
    ylab(expression(paste("||", gamma, "||"))) + ggtitle("Rate test") +
    scale_y_continuous(breaks = seq(0, 1, by = .1)) +
    scale_x_continuous(breaks = seq(0, 1, by = .1)) 

p1 + theme(legend.position= "top")
p1 + theme(legend.direction = "horizontal", legend.position= "bottom")

library(gridExtra)

grid.arrange(p1 + theme(legend.position = "none"), 
             p2 + theme(legend.position = "none") + ylab(""),
             p3 + theme(legend.position = "none") + ylab(""), nrow = 1)
## ggsave("all.pdf")

m0 <- data.frame(beta = rep(m$beta, 3),
                 gamma = rep(m$gamma, 3),
                 test = rep(c("shape", "size", "rate"), each = nrow(m)),
                 power = c(m$shape, m$size, m$rate))
head(m0)
m0$test <- ordered(m0$test, levels = c("shape", "size", "rate"))

ggplot(m0, aes(beta, gamma)) +
    geom_tile(aes(fill = power), color = "white") + facet_wrap(~test, nrow = 1) +
    scale_fill_gradient("Rejection \nproportion", low = "gray", high = "gray50", limits = c(0,1)) +
    xlab(expression(paste("|| ", beta, " ||"))) +
    ylab(expression(paste("|| ", gamma, " ||"))) +
    scale_y_continuous(breaks = seq(0, 1, by = .1)) +
    scale_x_continuous(breaks = seq(0, 1, by = .1)) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14,face="bold")) +
    theme(strip.text.x = element_text(size = 12))
## ggsave("all-independent.pdf")
ggsave("all-informative.pdf")


sumPwr(200, FALSE, 0.2, 0.2)
sumPwr(200, TRUE, 0.2, 0.2)

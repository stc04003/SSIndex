#' Function to generate simulated data
#'
#' The function \code{simDat} generates simulated recurrent event data from the following model:
#' \describe{
#'   \item{\code{type = "M1"}}{generates recurrent event data from the transformation model.}
#'   \item{\code{type = "M2"}}{generates recurrent event data from the additive model.}
#'   \item{\code{type = "M3"}}{generates recurrent event data from the additive model with \deqn{\gamma = 0}.}
#'   \item{\code{type = "M3"}}{generates recurrent event data from the two-sample model.}
#'   \item{\code{type = "M4"}}{generates recurrent event data from the likelihood ratio model.}
#' }
#' The data frame is similar to the one used in \code{reReg}.
#' The data consists of the following columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{t}{Patient's event time; this column consists of both the recurrent event time and the terminal event time}
#'   \item{m}{The observed number of recurrent events}
#'   \item{x1, x2}{Baseline covariates}
#'   \item{event}{recurrent event indicator; 1 = recurrent event, 0 = not a recurrent event}
#'   \item{status}{censoring indicator; this is only meanful for when event = 0}
#' }
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom tibble as.tibble
simDat <- function(n, model) {
    dat <- NULL
    beta0 <- c(.6, .8)
    gamma0 <- c(7, 24) / 25
    for (i in 1:n) {
        x0 <- 100
        while (abs(x0) > 1 ) {
            x <- mvrnorm(1, c(0,0), diag(2))
            x0 <-  sum(x * beta0)
        }
        y <- min(rexp(1, 1 / (10 * (1 + abs(x[1])))), 5)
        x_beta <- as.numeric(x %*% beta0)
        x_gamma <- as.numeric(x %*% gamma0)
        len <- 0  # max of recurrent times
        tmpt <- NULL # recurrent times 
        inx <- 0   # number of recurrent times
        while( len < Lam.f(y, x_gamma, x_beta, model)) {
	    tmpt <- c(tmpt, rexp(1))
	    len <- sum (tmpt) 
	    inx <- inx + 1
        }   
        m <- ifelse(inx > 0, inx - 1, 0)
        if (m > 0)  {
            tt <- invLam.f(cumsum(tmpt[1:m]), x_gamma, x_beta, model)
            tmp <- cbind(id = i, t = c(tt[order(tt)], y), m = m,
                         x1 = x[1], x2 = x[2], event = c(rep(1, length(tt)), 0), status = c(rep(0, length(tt)), 1))
        } else {
            tmp <- cbind(id = i, t = y, m = m, x1 = x[1], x2 = x[2], event = 0, status = 1)
        }
        dat <- rbind(dat, tmp)
    }
    return(as.tibble(dat))
}

#' Cumulative rate function
#' @noRd
#' @keywords internal
Lam.f <- function(t, r, b, model){ 
    if (model == "M1") return(5 * exp(-b) * t^2 * exp(2 * b) / (1 + exp(2 * b) * t^2) * exp(r))
    if (model == "M2") return((exp(t / 10) * 10 + t * b - 10) * exp(r) / 2)
    if (model == "M3") return((exp(t / 10) * 10 + t * b - 10))
    if (model == "M4") return(2 * log(1 + exp(b) * t ^1.5) * exp(r)/ 3)
    if (model == "M5") return((exp(t * b) - 1) * exp(r) / b)
    ## 5 * t / (1 + t)^2 + t * b
    ## (t / (t + 1) - 1 + t * b)
}

#' Inverse cumulative rate function
#' @noRd
#' @keywords internal
invLam.f <- function (t, r, b, model) {
    mapply(t, FUN = function(u) {  
        uf <- function (x) u - Lam.f (x, r, b, model)
        uniroot(uf, c(0, 100))$root
    })
}

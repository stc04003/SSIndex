#' Function to generate simulated data
#'
#' We formulate the proposed rate function as
#' \deqn{\mu(t|Z) = f(t, \beta_0^\top Z) g(\gamma_0^\top Z), 0\le t \le \tau}.
#' See Details.
#'
#' 
#' The function \code{simDat} generates simulated recurrent event data from our model's special cases.
#' The following model are implemented.
#' \describe{
#'   \item{\code{type = "M1"}}{generates recurrent event data from the proportional rate model.
#'     \deqn{\mu(t|Z) = \mu_0(t) e^{\gamma_0^\top Z}.}
#'      This model implies \eqn{\beta_0 = 0}.}
#'   \item{\code{type = "M2"}}{generates recurrent event data from the additive model.
#'     \deqn{\mu(t|Z) = \mu_0(t) + \alpha_0^\top Z.}
#'      This model implies \eqn{\beta_0 = \gamma_0}.}
#'   \item{\code{type = "M3"}}{generates recurrent event data from the accelerated rate model.
#'     \deqn{\mu(t|Z) = \mu_0(t e^{\alpha_0^\top Z}).}
#'      This model implies \eqn{\beta_0 = \gamma_0}.}
#'   \item{\code{type = "M4"}}{generates recurrent event data from the two-sample model.
#'     \deqn{\mu(t|Z) = \mu_0(t e^{\beta_0^\top Z})e^{\gamma_0^\top Z}.}}
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
#' For all scenarios, we set \eqn{\tau = 10}, \eqn{\beta_0 = (0.6, 0.8)} and \eqn{\gamma_0 = (7/25, 24/25)}.
#' We also set \eqn{\mu_0(t) = \frac{2}{1 + t}}, then \eqn{m_0(t) = 2\log(1 + t)}.
#' 
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom tibble as.tibble
simDat <- function(n, model, frailty = FALSE) {
    dat <- NULL
    beta0 <- c(.6, .8)
    gamma0 <- c(7, 24) / 25
    for (i in 1:n) {
        x0 <- 100
        if (frailty) zz <- rgamma(1, 1, 1)
        else zz <- 1
        if (model %in% c("M1", "M2", "M3", "M4", "M5")) {
            while (abs(x0) > 1 ) {
                x <- mvrnorm(1, c(0,0), diag(2))
                x0 <-  sum(x * beta0)
            }
        }
        if (model %in% c("M21", "M22", "M23", "M24", "M25")) 
            x <- c(sample(0:1, 1), runif(1, -1, 1))
        if (model %in% c("M31", "M32", "M33", "M34", "M35"))
            x <- c(runif(1, -1, 1), rexp(1))
        if (model %in% c("M4", "M24", "M34")) tau <- 1
        else tau <- 10
        y <- min(rexp(1, zz / (10 * (1 + abs(x[1])))), tau)
        x_beta <- as.numeric(x %*% beta0)
        x_gamma <- as.numeric(x %*% gamma0)
        len <- 0  # max of recurrent times
        tmpt <- NULL # recurrent times
        inx <- 0   # number of recurrent times
        while( len < Lam.f(y, x_gamma, x_beta, model, zz)) {
	    tmpt <- c(tmpt, rexp(1))
	    len <- sum(tmpt)
	    inx <- inx + 1
        }
        m <- ifelse(inx > 0, inx - 1, 0)
        if (m > 0)  {
            tt <- invLam.f(cumsum(tmpt[1:m]), x_gamma, x_beta, model, zz)
            tmp <- cbind(id = i, t = c(tt[order(tt)], y), m = m,
                         x1 = x[1], x2 = x[2], event = c(rep(1, length(tt)), 0),
                         status = c(rep(0, length(tt)), 1))
        } else {
            tmp <- cbind(id = i, t = y, m = m, x1 = x[1], x2 = x[2], event = 0, status = 1)
        }
        dat <- rbind(dat, tmp)
    }
    return(as.tibble(dat))
}

#' Cumulative rate function
#'
#' M1 ~ M3 assume \eqn{\tau = 10}, M4 assumes \eqn{\tau = 1}.
#' 
#' @param t is time
#' @param r is \eqn{\gamma_0^\top Z}
#' @param b is \eqn{\beta_0^\top Z}
#' @param model is the model indicator
#' 
#' @noRd
#' @keywords internal
Lam.f <- function(t, r, b, model, zz){
    if (model == "M1") return(zz * (2 * log(1 + t) * exp(r)))
    if (model == "M2") return(zz * (2 * (1 + t)^2 + t * b)/ 200)
    ## if (model == "M2") return(zz * (0.2 * t^2  + 0.2 * t * exp(b))) # works 
    ## if (model == "M2") return(zz * (0.5 * t^2  + 0.5 * t * exp(b))) # works but average 30 events / subject
    ## if (model == "M2") return(zz * (exp(t / 10) * 10 + t * exp(b) - 10)) # okay
    ## if (model == "M2") return(zz * (exp(t / 10) * 1 + t * exp(b) - 1))
    if (model == "M3") return(zz * ((1 - exp(-t * exp(b) / 2)) * 2 * exp(-b)))
    if (model == "M4") return(zz * (pbeta(t, 2, 1 + exp(b)) * exp(r)))
    if (model == "M5") return(zz * (10 * (1 + t) ^ (exp(b) / 5) - 10))
    ## Old settings, under Pico's draft
    ## if (model == "M1") return(5 * exp(-b) * t^2 * exp(2 * b) / (1 + exp(2 * b) * t^2) * exp(r))
    ## if (model == "M2") return((exp(t / 10) * 10 + t * b - 10) * exp(r) / 2)
    ## if (model == "M3") return((exp(t / 10) * 10 + t * b - 10))
    ## if (model == "M4") return(2 * log(1 + exp(b) * t ^1.5) * exp(r)/ 3)
    ## if (model == "M5") return((exp(t * b) - 1) * exp(r) / b)
    ## 5 * t / (1 + t)^2 + t * b
    ## (t / (t + 1) - 1 + t * b)
}

#' Inverse cumulative rate function
#' @noRd
#' @keywords internal
invLam.f <- function (t, r, b, model, zz) {
    mapply(t, FUN = function(u) {
        uf <- function(x) u - Lam.f(x, r, b, model, zz) ## / Lam.f(10, r, b, model)
        uniroot(uf, c(0, 100))$root
    })
}

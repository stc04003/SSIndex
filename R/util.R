#' Kernel function used in the estimation of shape function
#'
#' @noRd
#' @keywords internal
k_f <- function(u){
    k1 <- (15/16) * ((1-u^2)^2)*(abs(u)<=1)
    return((7/4)*(1-3*(u^2))*k1)
}


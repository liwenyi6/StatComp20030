# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A laplace function
#' @description A probability density function of laplace function
#' @param x the independent value
#' @return the value of the function
#' @export
f <- function(x) {
    .Call('_StatComp20030_f', PACKAGE = 'StatComp20030', x)
}

#' @title A random walk sampling using Rcpp
#' @description A random walk sampling using Rcpp
#' @param sigma the parameter of variance
#' @param x0 the initial value 
#' @param N the number of uniform distributed random samples
#' @return a random sample of size \code{n}
#' @importFrom stats runif dnorm rnorm
#' @export
rwMetropolis <- function(sigma, x0, N) {
    .Call('_StatComp20030_rwMetropolis', PACKAGE = 'StatComp20030', sigma, x0, N)
}


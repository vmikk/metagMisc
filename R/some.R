#' @title Check random lines in the data
#'
#' @description Returns several random lines from the data.frame or matrix.
#' Usage is quite similar to the \code{head} or \code{tail} functions.
#'
#' @param x A data.frame or matrix
#' @param n Number of records to show
#' @return Part of the data.frame or matrix with n lines.
#' @author Adapted from \code{some} in the car-package by John Fox.
#' @examples
#' x <- matrix(rnorm(500), ncol=5)
#' some(x)
#'
some <- function(x, n=10){
  show <- sample(x = 1:nrow(x), size = 10)
  return( x[show, ] )
}

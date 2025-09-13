
#' Pad numbers with leading zeros
#'
#' @description 
#' Converts numeric values to character strings with leading zeros to ensure 
#' consistent width. Useful for creating sortable labels or file names.
#'
#' @param x Numeric vector of integers to be padded with leading zeros
#' @param ndigits Integer specifying the total number of digits in the output. 
#'   If \code{NULL} (default), automatically determined from the maximum value in \code{x}
#'
#' @details
#' This function formats integers as character strings with leading zeros to create 
#' fixed-width representations. When \code{ndigits} is not specified, the function 
#' automatically calculates the minimum number of digits needed to represent the 
#' largest value in \code{x}.
#' 
#' The function handles edge cases appropriately:
#' \itemize{
#'   \item Zero values are handled correctly
#'   \item Negative numbers are supported (the minus sign does not count toward digit padding)
#'   \item Non-finite values (Inf, -Inf, NaN) are treated as having 0 digits for calculation purposes
#' }
#'
#' @return Character vector of the same length as \code{x}, with numbers formatted 
#'   as zero-padded strings
#'
#' @examples
#' # Automatic padding based on maximum value
#' zero_pad(c(1, 10, 100))
#' # Returns: "001" "010" "100"
#' 
#' # Manual specification of digits
#' zero_pad(c(1, 5, 12), ndigits = 4)
#' # Returns: "0001" "0005" "0012"
#'
#' @export
#' 
zero_pad <- function(x, ndigits = NULL){

  n_int_digits <- function(z) {
    result <- floor(log10(abs(z)))
    result[!is.finite(result)] <- 0
    result
  }

  ## Automatically choose the number of digits
  if(is.null(ndigits)){
    res <- sprintf(paste("%0", n_int_digits(max(x)) + 1, "d", sep = ""), x)
  } else {
  ## Use user-defined number
    res <- sprintf(paste("%0", ndigits, "d", sep = ""), x)
  }
  return(res)
}

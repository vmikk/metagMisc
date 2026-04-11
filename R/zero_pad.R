
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
#'   \item Non-finite values (Inf, -Inf, NaN) are converted to character strings
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
    ## Only consider finite values for calculating ndigits
    finite_vals <- x[is.finite(x) & !is.na(x)]
    if(length(finite_vals) == 0) {
      ## If no finite values, default to 1 digit
      ndigits <- 1
    } else {
      ndigits <- max(n_int_digits(finite_vals)) + 1
    }
  }
  
  ## Initialize result vector
  res <- character(length(x))
  
  ## Create masks for different value types
  finite_mask <- is.finite(x) & !is.na(x)
  na_mask <- is.na(x) & !is.nan(x)
  nan_mask <- is.nan(x)
  pos_inf_mask <- !is.na(x) & is.infinite(x) & x > 0
  neg_inf_mask <- !is.na(x) & is.infinite(x) & x < 0
  
  ## Format finite values (including negative numbers and zero)
  if(any(finite_mask)) {
    finite_vals <- x[finite_mask]
    
    ## For negative numbers, we need an extra digit to account for the minus sign
    ## since sprintf counts the minus sign in the total width
    positive_mask <- finite_vals >= 0
    negative_mask <- finite_vals < 0
    
    if(any(positive_mask)) {
      fmt_pos <- paste0("%0", ndigits, "d")
      res[finite_mask][positive_mask] <- sprintf(fmt_pos, finite_vals[positive_mask])
    }
    
    if(any(negative_mask)) {
      fmt_neg <- paste0("%0", ndigits + 1, "d")
      res[finite_mask][negative_mask] <- sprintf(fmt_neg, finite_vals[negative_mask])
    }
  }
  
  ## Handle non-finite values as character strings
  if(non_finits_found == TRUE){
    res[na_mask] <- "NA"
    res[nan_mask] <- "NaN"
    res[pos_inf_mask] <- "Inf"
    res[neg_inf_mask] <- "-Inf"
  }
  
  return(res)
}

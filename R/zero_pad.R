
## Pad number with leading zeros
## NB. output = character vector
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

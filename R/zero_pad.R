
## Pad number with leading zeros
## NB. output = character vector
zero_pad <- function(x){

  n_int_digits <- function(z) {
    result <- floor(log10(abs(z)))
    result[!is.finite(result)] <- 0
    result
  }

  res <- sprintf(paste("%0", n_int_digits(max(x)) + 1, "d", sep = ""), x)
  return(res)
}

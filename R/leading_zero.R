
## Function to add leading zeroes to numeric vector
leading_zero <- function(x, z = NULL){
  # x = numeric vector
  # z = resulting length of number
  if(is.null(z)){
    ## Find the number of characters
    z <- nchar(as.integer(max(x)))
  }

  zz <- paste("%0", z, "d", sep = "")
  sprintf(zz, x)
}
# leading_zero(1:10, z=3)
# leading_zero(c(1,10,100,1000,10000))

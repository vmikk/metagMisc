
## Function to add or remove leading zeros to numeric vector
leading_zero <- function(x, z = NULL, mode = "add"){
  # x = numeric vector
  # z = resulting length of number
  # mode = "add" or "remove" leading zero

  ## Add leading zeros
  if(mode == "add"){
    if(is.null(z)){
      ## Find the number of characters
      z <- nchar(as.integer(max(x)))
    }
  
    zz <- paste("%0", z, "d", sep = "")
    res <- sprintf(zz, x)
  }

  ## Remove leading zeros
  if(mode == "remove"){
    res <- gsub("(^|[^0-9])0+", "\\1", x, perl = TRUE)
  }

  return(res)
}
# leading_zero(1:10, z=3)
# leading_zero(c(1,10,100,1000,10000))
# leading_zero(c("001", "1001", "g01", "0abc"), mode = "remove")

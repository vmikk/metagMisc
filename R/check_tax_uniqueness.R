
check_tax_uniqueness <- function(x, col = "k", return_all = F, dropNA = T){
  # x = data frame
  #NB! Columns should be ordered, missing values are coded with NA

  require(plyr)

  x <- as.data.frame(x)

  ## Where is the selected rank?
  COLID <- which(colnames(x) == col)

  ## Remove NAs
  if(dropNA == TRUE){
    nn <- is.na(x[,COLID])
    if(any(nn)){ x <- x[!nn, ] }
  }

  ## Remove lower ranks and count number of unique values
  res <- ddply(.data = x[, 1:COLID], .variables = col, .fun = function(z){
    rez <- z[!duplicated(z),]    # remove duplicates
    rezz <- data.frame(UniqueCombs = nrow(rez))
    return(rezz)
    })

  ## Remove unique
  if(return_all == FALSE){
    res <- subset(res, UniqueCombs > 1)
    if(nrow(res) == 0){ cat("All ranks of ", col, " are unique.\n") }
  }
  return(res)
}

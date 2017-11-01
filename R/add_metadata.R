
add_metadata <- function(x, metad, xid, mid, drop_mid = T){

  ## Match metadata to the main data
  mm <- match(x = x[, xid], table = metad[, mid])

  ## Reorder metadata
  metad <- metad[mm, ]

  ## Remove sample ID column from metadata
  if(drop_mid == TRUE){
    metad <- metad[, -which(colnames(metad) == mid)]
  }

  ## Merge data with metadata
  res <- cbind(x, metad)

  return(res)
}

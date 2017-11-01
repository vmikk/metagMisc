
#' @title Add metadata to data.frame
#' @description This function combines two tables based on the row identifiers.
#' @param x Data frame with data of interest
#' @param metad Data frame with metadata that should be added; or phyloseq object with sample_data slot
#' @param xid Character, name of the column with row identifiers in the main data (x)
#' @param mid Character, name of the column with row identifiers in the metadata (metad)
#' @param drop_mid Logical; if TRUE (default), metadata identifiers will be removed from results
#' @return Data frame
#' @export
#'
#' @examples
#'
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

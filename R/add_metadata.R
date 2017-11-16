
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
#' # Load data
#' data(GlobalPatterns)
#'
#' # As an example, estimate richness
#' div <- estimate_richness(GlobalPatterns, measures = c("Observed", "Shannon"))
#' div <- dfRowName(div, name = "SampleID")
#'
#' # Take metadata from phyloseq object
#' div1 <- add_metadata(div, metad = GlobalPatterns, xid = "SampleID", mid = "X.SampleID")
#' head(div1)
#'
#' # Or metadata from data frame
#' metad <- as(sample_data(GlobalPatterns), "data.frame")
#' div2 <- add_metadata(div, metad = metad, xid = "SampleID", mid = "X.SampleID")
#' head(div2)
#'
add_metadata <- function(x, metad, xid = NULL, mid = NULL, drop_mid = T){

  ## Data validation
  if(is.null(xid) | is.null(mid)){ stop("Error: row identifiers should be provided.\n") }
  if(!length(xid) == 1){ stop("Error: row identifiers should be provided as a single character string, e.g. 'SampleID'.\n") }
  if(!length(mid) == 1){ stop("Error: row identifiers should be provided as a single character string, e.g. 'SampleID'.\n") }

  ## Extract metadata from phyloseq object
  if(class(metad) %in% "phyloseq"){
    if(is.null(sample_data(metad, errorIfNULL = F))){
      stop("Error: sample_data is missing from the phyloseq object 'metad'.\n")
    } else {
      metad <- as(object = sample_data(metad), Class = "data.frame")
    }
  }

  ## Data validation
  if(!xid %in% colnames(x)){ stop("Error: '", xid, "' column is missing in the main data.\n", sep="") }
  if(!mid %in% colnames(metad)){ stop("Error: '", mid, "' column is missing in the metadata.\n", sep="") }

  if(nrow(x) != length(unique(x[, xid]))){ stop("Error: Row identifiers are not unique in 'x'.\n") }
  if(nrow(metad) != length(unique(metad[, mid]))){ stop("Error: Row identifiers are not unique in 'metad'.\n") }


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

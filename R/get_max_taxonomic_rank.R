
## Get max taxonomic rank
# (Select last non-NA column)
get_max_taxonomic_rank <- function(x, return_rank_only = FALSE){
  # x = data.frame, columns are ordered tax ranks (NA is the missing rank)

  require(plyr)

  ## Find the indices of the last non-NA column for each row
  res <- adply(.data = x, .margins = 1, .fun = function(z){

    ## Test which tax ranks are not NAs
    rnk <- aaply(.data = z, .margins = 1, .fun = function(y) which(!is.na(y)) )

    ## Return last non-NA column number
    if(length(rnk) > 0){
      rez <- max(rnk)
    } else {
      rez <- 0   # if all tax ranks are NA
    }
    rez <- data.frame(RankColumn = rez)
    return(rez)
  })

  ## Table with correspondence of taxonomic ranks and their indices
  rnks <- data.frame(
    RankColumn = c(0, 1:ncol(x)),
    RankName = c(NA, colnames(x)),
    stringsAsFactors = F)

  ## Substitute column index with tax rank name
  res$RankName <- rnks$RankName[match(x = res$RankColumn, table = rnks$RankColumn)]

  if(return_rank_only == TRUE){
    res <- res$RankName
  }

  return(res)
}

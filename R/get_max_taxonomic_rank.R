
#' @title Determine the lowest level of taxonomic classification
#' @details This function will find the last non-NA column in the taxonomy table and return
#' @param x Either a phyloseq object, or a data frame with columns as taxonomic ranks and rows as entries (e.g., OTUs). Columns in the data frame should be ordered from the highest level of classification (e.g., Kingdom) to the lowest level (e.g., Species), missing data are coded as NA
#' @param return_rank_only Logical, if TRUE only name of the taxonomic rank will be returned
#'
#' @return Data frame with taxonomy and additional column containing the name of the lowest level of taxonomic classification. Alternatively, if 'return_rank_only = TRUE', a vector of the lowest taxonomic ranks for each OTU.
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#'
#' # phyloseq-class as input
#' taxx <- get_max_taxonomic_rank(GlobalPatterns)
#' summary(taxx$RankName)
#'
#' # data frame as input
#' taxtbl <- as.data.frame(tax_table(GlobalPatterns))
#' taxx <- get_max_taxonomic_rank(taxtbl)
#' summary(taxx$RankName)
#'
get_max_taxonomic_rank <- function(x, return_rank_only = FALSE){

  # require(plyr)

  ## Input data
  inp_class <- class(x)

  ## If input is of class 'phyloseq'
  if("phyloseq" %in% inp_class || "taxonomyTable" %in% inp_class){

    if(is.null(phyloseq::tax_table(x, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }

    otu_names <- phyloseq::taxa_names(x)
    x <- as.data.frame(phyloseq::tax_table(x), stringsAsFactors = F)
  }

  ## If there are factors in the input data frame this can cause an "Error: cannot allocate vector of size ... Gb"
  ## Convert all factors to character
  if("factor" %in% plyr::laply(.data = x, .fun = class)){
    x[] <- lapply(x, as.character)
  }

  ## Display progress bar?
  if(nrow(x) < 200){
    progr <- "none"
  } else {
    progr <- "text"
  }

  ## Find the indices of the last non-NA column for each row
  res <- plyr::adply(.data = x, .margins = 1, .fun = function(z){

    ## Test which tax ranks are not NAs
    rnk <- plyr::aaply(.data = z, .margins = 1, .fun = function(y) which(!is.na(y)) )

    ## Return last non-NA column number
    if(length(rnk) > 0){
      rez <- max(rnk)
    } else {
      rez <- 0   # if all tax ranks are NA
    }
    rez <- data.frame(RankColumn = rez)
    return(rez)
  }, .progress = progr)

  ## Table with correspondence of taxonomic ranks and their indices
  rnks <- data.frame(
    RankColumn = c(0, 1:ncol(x)),
    RankName = c(NA, colnames(x)),
    stringsAsFactors = F)

  ## Substitute column index with tax rank name
  res$RankName <- rnks$RankName[match(x = res$RankColumn, table = rnks$RankColumn)]

  ## Remove rank column
  res$RankColumn <- NULL

  ## Reorder taxonomic ranks
  res$RankName <- factor(res$RankName, levels = colnames(x))

  ## Add OTU name
  if("phyloseq" %in% inp_class || "taxonomyTable" %in% inp_class){
    res <- data.frame(TaxaName = otu_names, res, stringsAsFactors = F)
  }

  if(return_rank_only == TRUE){
    res <- res$RankName
  }

  return(res)
}

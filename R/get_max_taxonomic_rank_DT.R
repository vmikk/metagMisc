
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
get_max_taxonomic_rank_DT <- function(x, return_rank_only = FALSE){

  ## Input data
  inp_class <- class(x)

  ## If input is of class 'phyloseq'
  if("phyloseq" %in% inp_class || "taxonomyTable" %in% inp_class){

    if(is.null(phyloseq::tax_table(x, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }

    ## Extract taxonomy table
    otu_names <- phyloseq::taxa_names(x)
    x <- as.data.frame(phyloseq::tax_table(x), stringsAsFactors = F)
  }

  ## Convert to data.table object (rownames would be lost!)
  otu_names <- rownames(x)
  setDT(x)

  ## Tax rank names
  tax_cols <- copy(colnames(x))
  # tax_cols <-  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  ## If there are factors in the input data frame this can cause an "Error: cannot allocate vector of size ... Gb"
  ## Convert all factors to character
  if("factor" %in% plyr::laply(.data = x, .fun = class)){
    x <- x[, lapply(.SD, as.character), .SDcols = tax_cols]
  }

  ## Function to extract lowest taxonomy rank observed
  get_max_rank <- function(z){ 
    rnk <- which(!is.na(z))
    if(length(rnk) > 0){
      return( tax_cols[ max(rnk) ]  )
    } else {
      return(NA)
    }
  }

  ## Function to extract lowest taxon
  get_max_taxon <- function(z){
    z <- z[ !is.na(z) ]
    if(length(z)==0){ 
      return(NA) 
    } else { 
      return(z[length(z)]) 
    }
  }

  ## Find the lowest taxa for each row
  x[, `:=` (
      LowestTaxon = apply(.SD, 1, get_max_taxon),
      RankName = apply(.SD, 1, get_max_rank) ),
    .SDcols = tax_cols
    ]

  x$RankName <- factor(x = x$RankName, levels = tax_cols)

  ## Return only ranks
  if(return_rank_only == TRUE){
    res <- res$RankName
  } else {
  ## Or full table
    res <- data.frame(TaxaName = otu_names, x, stringsAsFactors = F)
  }

  return(res)
}


#' Impute missing taxonomy assignments in phyloseq objects
#'
#' @description 
#' Replaces missing (NA) taxonomic assignments with informative placeholders 
#' based on the nearest available higher-level taxonomic information. Handles 
#' species and higher taxonomic ranks differently.
#'
#' @param phys A \code{phyloseq} object or \code{taxonomyTable} containing 
#'   taxonomic assignments with potential missing values
#' @param unknown_taxon Character string to append to higher taxonomic ranks 
#'   when the assignment is missing (default: "_unidentified")
#' @param unknown_sp Character string to append to genus names when species 
#'   assignment is missing (default: " sp")
#' @param make_unique Logical indicating whether to make species names unique 
#'   by appending numerical suffixes (default: \code{FALSE})
#' @param addmaxrank Logical indicating whether to add a column showing the 
#'   lowest resolved taxonomic rank for each OTU/ASV (default: \code{FALSE})
#'
#' @details
#' This function addresses common issues in metagenomic/phylogenetic datasets 
#' where taxonomic assignments are incomplete. Missing values are replaced using 
#' a hierarchical approach:
#' 
#' \strong{Higher taxonomic ranks}: Missing values are replaced with the previous 
#' (higher) rank name plus the \code{unknown_taxon} suffix. For example, if an 
#' OTU is classified to family "Enterobacteriaceae" but genus is missing, the 
#' genus becomes "Enterobacteriaceae_unidentified".
#' 
#' \strong{Species level}: Missing species assignments are replaced with the 
#' genus name plus the \code{unknown_sp} suffix (e.g., "Escherichia sp").
#' 
#' \strong{Cleanup operations}:
#' \itemize{
#'   \item Removes duplicate unidentified strings (e.g., "_unidentified_unidentified")
#'   \item Cleans up species-level artifacts (e.g., "Genus_unidentified sp" -> "Genus sp")
#'   \item Optionally ensures unique species names using \code{make.unique()}
#' }
#'
#' @return A \code{phyloseq} object with imputed taxonomy table. Missing 
#'   taxonomic assignments are replaced with informative placeholders derived 
#'   from higher-level classifications.
#'
#' @examples
#' library(phyloseq)
#' 
#' # Basic imputation with default settings
#' ps_imputed <- phyloseq_taxonomy_imputation(physeq)
#' 
#' # Custom placeholders for missing assignments
#' ps_custom <- phyloseq_taxonomy_imputation(
#'   physeq,
#'   unknown_taxon = "_unknown",
#'   unknown_sp = "_sp"
#' )
#' 
#' # Make species names unique and add rank information
#' ps_enhanced <- phyloseq_taxonomy_imputation(
#'   physeq,
#'   make_unique = TRUE,
#'   addmaxrank = TRUE
#' )
#'
#' @seealso \code{\link{get_max_taxonomic_rank}} for identifying the lowest 
#'   resolved taxonomic rank
#'
#' @export
#'
phyloseq_taxonomy_imputation <- function(phys,
  unknown_taxon = "_unidentified", unknown_sp = " sp",
  make_unique = FALSE, addmaxrank = FALSE){

  ## If input is of class 'phyloseq'
  inp_class <- class(phys)
  if("phyloseq" %in% inp_class || "taxonomyTable" %in% inp_class){

    if(is.null(phyloseq::tax_table(phys, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }

    x <- as.data.frame(phyloseq::tax_table(phys), stringsAsFactors = F)
  } else {
    stop("Error: input data should be of class `phyloseq` or `taxonomyTable`")
  }

  ## Function to replace NAs with higher taxa names + unident string
  replace_col <- function(x, col_num = 2, unident = "unidentified"){
    nas <- is.na(x[, col_num])
    if(any(nas) == FALSE){
      return(x)              # no missing values
    } else {
      ## which values to replace
      to_repl <- which(nas)

      if(col_num == 1){
      	x[to_repl, col_num] <- unident
      } else {
        x[to_repl, col_num] <- paste(x[to_repl, col_num - 1], unident, sep = "")
      }
      return(x)
    }
  }

  clz <- colnames(x)
  sp_in_ranks <- clz %in% c("Species", "species", "sp")
  if(any(sp_in_ranks)){
    sp_id <- which(sp_in_ranks)
    non_sp_id <- (1:length(clz))[-sp_id]
  } else {
    sp_id <- NA
    non_sp_id <- 1:length(clz)
  }

  ## Replace higher taxa
  for(i in non_sp_id){
    x <- replace_col(x, col_num = i, unident = unknown_taxon)
  }

  ## Replace species names
  if(!is.na(sp_id)){
    x <- replace_col(x, col_num = sp_id, unident = unknown_sp)
  }

  ## Function to remove multiple unident strings (e.g., "_unidentified_unidentified")
  replace_unidents <- function(tt, strr = "_unidentified", spp = " sp"){
    ## tt = character vector
  
    ## Prepare regex for paterns
    # Multiple string occurrences
    multpatt <- paste("(", strr, ")(\\1+)", sep = "")
    unsp <- paste(strr, spp, sep = "")

    ## Replace "_unidetified_unidentified" with single occurrence
    # gsub(x = tt, pattern = "(_unidentified)(\\1+)", replacement = "_unidentified", perl = T)
    rez <- gsub(x = tt, pattern = multpatt, replacement = strr, perl = T)

    ## Replace "_unidetified sp" with " sp"
    rez <- gsub(x = rez, pattern = unsp, replacement = spp, perl = T)    

    return(rez)
  }

  ## Remove multiple unident strings from all tax columns
  x <- sapply(x, replace_unidents)

  ## Make species names unique
  if(make_unique == TRUE){
    x[, ncol(x)] <- base::make.unique(names = x[, ncol(x)], sep = ".")
  }

  ## Add the OTU classification at the lowest annotated taxonomic rank
  if(addmaxrank == TRUE){
    LowestTaxRank <- as.character( get_max_taxonomic_rank(phys, return_rank_only = TRUE) )
    x <- cbind(x, LowestTaxRank = LowestTaxRank)
  }

  ## Add taxa names
  rownames(x) <- phyloseq::taxa_names(phys)

  ## Replace tax_table with the modified one
  phyloseq::tax_table(phys) <- phyloseq::tax_table(as.matrix(x))
  return(phys)
}



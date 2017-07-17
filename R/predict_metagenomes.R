
#' @title Produce the actual metagenome functional predictions.
#'
#' @param otu_tab Data frame with OTU abundances (rows = OTUs, columns = Samples, first column = OTU names)
#' @param func_tab Data frame with precalculated function predictions on per OTU basis (rows = OTUs, columns = feature counts, first column = OTU names)
#' @param NSTI_present Logical; idnicating weather NSTI values are present in the last column of func_tab
#' @param rel_abund Logical; if TRUE, OTU counts will be transformed to relative abundances
#' @param add_sub_tabs Logical; if TRUE, subsetted OTU and functional tables will be added to the results
#' @details
#' This function is analogous to the 'predict_metagenomes.py' from PICRUSt.
#' It will produce the actual metagenome functional predictions for a given OTU table and table with feature counts,
#' e.g. count number of gene copies for each feature (e.g., KEGG) in each sample.
#' Feature counts for each OTU will be multiplied by the abundance of that OTU in each each sample and summed across all OTUs.
#'
#' @return Data frame with samples as rows and features as columns.
#' @export
#' @references https://picrust.github.io/picrust/scripts/predict_metagenomes.html
#'
#' @examples
#' ## Create dummy data
#' set.seed(111)
#' NSAMP=10    # number of samples
#' NSPEC=30    # number of species/OTUs
#' NGENES=8   # number of features
#'
#' dummy_name <- function(len = 5){ paste(sample(letters, size = len, replace = T), collapse = "") }
#'
#' # Table with precalculated number of gene copies per OTU
#' func_tab <- data.frame(
#'   OTU_ID = replicate(n = 100, expr = dummy_name()),
#'   matrix(data = sample(1:100, size = 100*NGENES, replace = T), nrow = 100),
#'   stringsAsFactors = F)
#' colnames(func_tab)[-1] <- paste("Gene", 1:NGENES, sep="")
#'
#' # Table with OTU abundance
#' otu_tab <- data.frame(
#'   OTU = sample(func_tab$OTU_ID, size = NSPEC),
#'   matrix(data = sample(0:200, size = NSPEC*NSAMP, replace = T), nrow = NSPEC),
#'   stringsAsFactors = F)
#' colnames(otu_tab)[-1] <- paste("Sample", 1:NSAMP, sep="")
#'
#' ## Predict metagenomes
#' predict_metagenomes(otu_tab, func_tab, NSTI_present = FALSE, rel_abund = FALSE, add_sub_tabs = FALSE)
#'
predict_metagenomes <- function(otu_tab, func_tab, NSTI_present = TRUE, rel_abund = TRUE, add_sub_tabs = FALSE){

  # Subset OTUs for the features present in the functional table
  otu_vs_func <- otu_tab[, 1] %in% func_tab[, 1]
  if(any(!otu_vs_func)){
    cat("Warning: OTUs that are not present in the table with functional features will be removed.\n")
    otu_tab <- otu_tab[otu_vs_func, ]
  }
  otus <- as.matrix( otu_tab[, -1] )

  # Extract OTUs from functional table
  match_otus <- match(x = otu_tab[, 1], table = func_tab[,1])
  if(any(is.na(match_otus))){ match_otus <- na.omit(match_otus) }
  funcs <- as.matrix( func_tab[match_otus, -1] )

  # Remove NSTI scores if present
  if(NSTI_present == TRUE){ funcs <- funcs[, -ncol(funcs)] }

  # Standardize OTU count to relative abundance
  if(rel_abund == TRUE){ otus <- decostand(otus, method = "total", MARGIN = 2) }

  # Multiply gene counts for each OTU by the abundance of that OTU in each each sample
  # and sum across all OTUs
  res <- crossprod(otus, funcs)

  # Add subsetted tables to the results
  if(add_sub_tabs == TRUE){
    attr(res, which = "OTUs") <- otu_tab
    attr(res, which = "Funcs") <- func_tab[match_otus, ]
  }

  return(res)
}

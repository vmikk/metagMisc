
#' @title Remove unused taxonomy levels from phyloseq-object.
#' @details This function is particularly useful after \code{\link{tax_glom}}.
#' @param physeq  A phyloseq-class object.
#' @details 'physeq_rm_na_tax' will remove columns filled with NA values from the taxonomy table of phyloseq object.
#' @return Phyloseq object.
#' @export
#' @seealso \code{\link{tax_table}}
#'
#' @examples
#' data(GlobalPatterns)
#' gp <- tax_glom(GlobalPatterns, taxrank="Family")  # 7 taxonomic ranks
#' head(tax_table(gp))   # last columns are filled with NAs
#'
#' # Remove taxonomic levels filled with NAs
#' gp2 <- physeq_rm_na_tax(gp)  # 5 taxonomic ranks
#' head(tax_table(gp2))
#'
physeq_rm_na_tax <- function(physeq){
  # rm_all <- function(x) { Filter(function(x)!all(is.na(x)), df) }
  rm_all <- function(df) { df[, !apply(is.na(df), 2, all)] }
  tax_table(physeq) <- rm_all( tax_table(physeq) )
  return(physeq)
}

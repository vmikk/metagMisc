
#' @title Standardize OTU abundance table
#' @description Standardize phyloseq OTU table with with methods from \code{\link[vegan]{decostand}} from vegan package.
#' @param physeq A phyloseq-class object
#' @param method Standardization method (see \code{\link[vegan]{decostand}} for available options)
#' @param ... Additional parameters may be passed to vegan \code{\link[vegan]{decostand}} function
#' @return phyloseq object with standardized OTU table.
#' @seealso \code{\link[vegan]{decostand}}, \code{\link{phyloseq_transform_css}}, \code{\link{phyloseq_transform_vst_blind}}, \code{\link{phyloseq_transform_rlog_blind}}, \code{\link{physeq_transform_anderson_log}}
#' @export
#'
#' @examples
#'
phyloseq_standardize_otu_abundance <- function(physeq, method="total", ...){

  ## Check the orientation of the OTU table
  trows <- taxa_are_rows(physeq)
  if(trows == TRUE){ marg <- 2 } else { marg <- 1 }

  ## Extact OTU table
  comm <- as(object = otu_table(physeq), Class = "matrix")

  ## Standardize community table
  comm_std <- decostand(comm, method, MARGIN = marg, ...)

  ## Replace old otu_table with the new one
  otu_table(physeq) <- otu_table(comm_std, taxa_are_rows = trows)

  return(physeq)
}

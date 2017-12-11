
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
#' # Load data
#' data("esophagus")
#'
#' # Total-sum scaling (TSS) normalization
#' phyloseq_standardize_otu_abundance(esophagus, method = "total")
#' # the same as
#' transform_sample_counts(esophagus, function(OTU) OTU/sum(OTU) )
#' identical(
#'   phyloseq_standardize_otu_abundance(esophagus, method = "total"),
#'   transform_sample_counts(esophagus, function(OTU) OTU/sum(OTU)) )
#'
#' # Presence-absence scaling (0/1)
#' phyloseq_standardize_otu_abundance(esophagus, method = "pa")
#'
#' # Logarithmic transformation as in Anderson et al., 2006
#' phyloseq_standardize_otu_abundance(esophagus, method = "log")
#'
#' # Hellinger standardization
#' phyloseq_standardize_otu_abundance(esophagus, method = "hellinger")
#'
phyloseq_standardize_otu_abundance <- function(physeq, method="total", ...){

  ## Check the orientation of the OTU table
  trows <- phyloseq::taxa_are_rows(physeq)
  if(trows == TRUE){ marg <- 2 } else { marg <- 1 }

  ## Extact OTU table
  comm <- as(object = phyloseq::otu_table(physeq), Class = "matrix")

  ## Standardize community table
  comm_std <- vegan::decostand(comm, method, MARGIN = marg, ...)

  ## Replace old otu_table with the new one
  phyloseq::otu_table(physeq) <- phyloseq::otu_table(comm_std, taxa_are_rows = trows)

  return(physeq)
}

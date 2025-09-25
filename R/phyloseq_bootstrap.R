
#' Bootstrap samples in a phyloseq object
#'
#' Resamples the samples (columns) of a `phyloseq` object with replacement to
#' generate bootstrap replicates. Each replicate preserves the `otu_table` and
#' `sample_data`, and can optionally copy additional slots such as
#' `tax_table`, `refseq`, and `phy_tree`.
#'
#' @param phys A `phyloseq` object to bootstrap
#' @param n Integer number of bootstrap replicates to generate (default, 100)
#' @param seed Optional integer random seed for reproducibility. If `NULL`, the
#'   current RNG state is used
#' @param other_slots Character vector of additional slots to copy into each
#'   replicate. Allowed values are any of "tax_table", "refseq",
#'   "phy_tree". Use "none" or `NULL` to copy none. Default:
#'   `c("tax_table", "refseq", "phy_tree")`
#' @param progress Progress bar type passed to `plyr::llply` (e.g., "text" or "none")
#'
#' @return A list of length `n` of `phyloseq` objects containing
#'   bootstrap-resampled samples.
#'
#' @export
#' 
phyloseq_bootstrap <- function(phys, n = 100, seed = NULL, other_slots = c("tax_table", "refseq", "phy_tree"), progress = "text"){

  ## Preserve no other slots except "otu_table" and "sample_data"
  if(is.null(other_slots)) { other_slots <- "" }
  if(other_slots == "none"){ other_slots <- "" }

  ## Function to subset samples by IDs
  sample_boot <- function(ps, samp_ids, slots = other_slots){

    if(taxa_are_rows(ps) == TRUE) { oo <- phyloseq::otu_table(ps)[, samp_ids] }
    if(taxa_are_rows(ps) == FALSE){ oo <- phyloseq::otu_table(ps)[samp_ids, ] }

    if(!is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE))){
      dd <- phyloseq::sample_data(ps)[samp_ids, ]
      phyloseq::sample_names(oo) <- phyloseq::sample_names(dd)
      newps <- phyloseq::phyloseq(oo, dd)
    } else {
      newps <- phyloseq::phyloseq(oo)
    }

    ## Add other slots
    if("tax_table" %in% slots && !is.null(phyloseq::tax_table(ps, errorIfNULL = F))){ phyloseq::tax_table(newps) <- phyloseq::tax_table(ps) }
    if("phy_tree" %in% slots && !is.null(phyloseq::phy_tree(ps, errorIfNULL = F))){ phyloseq::phy_tree(newps) <- phyloseq::phy_tree(ps) }
    if("refseq" %in% slots && !is.null(phyloseq::refseq(ps, errorIfNULL = F))){ phyloseq::refseq(newps) <- phyloseq::refseq(ps) }

    return(newps)
  }

  ## Set random seed
  if(!is.null(seed)){ set.seed(seed) }

  ## Generate bootstrapped IDs
  IDS <- plyr::rlply(
    .n = n,
    .expr = sample(x = 1:nsamples(phys), size = phyloseq::nsamples(phys), replace = TRUE))

  ## Subset p
  res <- plyr::llply(
    .data = IDS,
    .fun = function(z, ...){ sample_boot(phys, samp_ids = z, ...) },
    .progress = progress,
    slots = other_slots
    )

  return(res)
}

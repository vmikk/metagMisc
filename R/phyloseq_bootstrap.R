
## Bootstrap samples in phyloseq object
phyloseq_bootstrap <- function(phys, n = 100, seed = NULL, other_slots = c("tax_table", "refseq", "phy_tree"), progress = "text"){

  ## Preserve no other slots except "otu_table" and "sample_data"
  if(is.null(other_slots)) { other_slots <- "" }
  if(other_slots == "none"){ other_slots <- "" }

  ## Function to subset samples by IDs
  sample_boot <- function(ps, samp_ids, slots = other_slots){

    if(taxa_are_rows(ps) == TRUE) { oo <- phyloseq::otu_table(ps)[, samp_ids] }
    if(taxa_are_rows(ps) == FALSE){ oo <- phyloseq::otu_table(ps)[samp_ids, ] }

    if(!is.null(phyloseq::sample_data(ps, errorIfNULL=TRUE))){
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
# gg <- phyloseq_bootstrap(phys = GlobalPatterns, n = 10, other_slots = "")

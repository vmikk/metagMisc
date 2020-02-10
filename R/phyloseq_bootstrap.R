
## Bootstrap samples in phyloseq object
phyloseq_bootstrap <- function(phys, n = 100, seed = NULL, other_slots = c("tax_table", "refseq", "phy_tree"), progress = "text"){

  ## Preserve no other slots except "otu_table" and "sample_data"
  if(is.null(other_slots)) { other_slots <- "" }
  if(other_slots == "none"){ other_slots <- "" }

  ## Function to subset samples by IDs
  sample_boot <- function(ps, samp_ids, slots = other_slots){

    if(taxa_are_rows(ps) == TRUE) { oo <- otu_table(ps)[, samp_ids] }
    if(taxa_are_rows(ps) == FALSE){ oo <- otu_table(ps)[samp_ids, ] }

    if(!is.null(sample_data(ps, errorIfNULL=TRUE))){
      dd <- sample_data(ps)[samp_ids, ]
      sample_names(oo) <- sample_names(dd)
      newps <- phyloseq(oo, dd)
    } else {
      newps <- phyloseq(oo)
    }

    ## Add other slots
    if("tax_table" %in% slots && !is.null(tax_table(ps, errorIfNULL = F))){ tax_table(newps) <- tax_table(ps) }
    if("phy_tree" %in% slots && !is.null(phy_tree(ps, errorIfNULL = F))){ phy_tree(newps) <- phy_tree(ps) }
    if("refseq" %in% slots && !is.null(refseq(ps, errorIfNULL = F))){ refseq(newps) <- refseq(ps) }

    return(newps)
  }

  ## Set random seed
  if(!is.null(seed)){ set.seed(seed) }

  ## Generate bootstrapped IDs
  IDS <- plyr::rlply(
    .n = n,
    .expr = sample(x = 1:nsamples(phys), size = nsamples(phys), replace = TRUE))

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

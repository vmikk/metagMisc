
## CLR transform function
phyloseq_clr <- function(physeq, pseudocount = 0.5){

  ## CLR function with special handling of zeros
  # https://github.com/zdk123/SpiecEasi/blob/e4d7c0af49eae1ead40a89c1f569ee589de94776/R/normalization.R#L53
  clr_func <- function(x, base = exp(1), tol=.Machine$double.eps, ...) {
      nzero <- (x >= tol)
      LOG <- log(ifelse(nzero, x, 1), base)
      ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
  }

  ## Add pseudocount if reqiured
  if(!is.null(pseudocount)){
    otu_table(physeq) <- otu_table(physeq) + pseudocount
  }

  ## Convert to relative abundances
  physeq <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))

  ## CLR-transform OTU table
  otu_table(physeq) <- phyloseq::transform_sample_counts(otu_table(physeq), clr_func)

  return(physeq)
}


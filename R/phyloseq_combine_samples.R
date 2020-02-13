## Function to combine samples (sum OTU abundance) into a single sample
phyloseq_combine_samples <- function(phys){

  ## Remove sample metadata
  if(!is.null(phyloseq::sample_data(phys, errorIfNULL=F))){ phys@sam_data <- NULL }
  
  ## Estimate OTU total abundance 
  tss <- phyloseq::taxa_sums(phys)

  ## Replace OTU table
  ott <- matrix(data = tss, ncol = 1, dimnames = list(names(tss), "Total"))
  phyloseq::otu_table(phys) <- phyloseq::otu_table(ott, taxa_are_rows = TRUE)
  
  return(phys)
}

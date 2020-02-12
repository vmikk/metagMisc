## Function to combine samples (sum OTU abundance) into a single sample
phyloseq_combine_samples <- function(phys){

  if(!is.null(sample_data(phys, errorIfNULL=F))){ phys@sam_data <- NULL }
  tss <- taxa_sums(phys)
  ott <- matrix(data = tss, ncol = 1, dimnames = list(names(tss), "Total"))
  otu_table(phys) <- otu_table(ott, taxa_are_rows = TRUE)
  return(phys)
}

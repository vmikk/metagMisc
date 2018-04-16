
# Estimate the observed abundance-based sample coverage for phyloseq object
phyloseq_coverage <- function(physeq, correct_singletons = FALSE){
  
  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        as.data.frame(otu_table(physeq, taxa_are_rows = T)),
        correct_singletons = correct_singletons)

  ## Estimate sample coverages
  res <- ldply(.data = x, .fun = function(z){ iNEXT:::Chat.Ind(z, sum(z)) })
  colnames(res) <- c("SampleID", "SampleCoverage")

  return(res)
}

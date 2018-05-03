
## Compute OTU diversity for a particular sample size or coverage
phyloseq_effort_div <- function(physeq, base = "size", level = NULL, conf = 0.95, correct_singletons = FALSE){
  # conf=0.95 (could be NULL)
  # If base="size" and level=NULL,
  #   then this function computes the diversity estimates for the minimum sample size among all sites.
  # If base="coverage" and level=NULL,
  #   then this function computes the diversity estimates for the minimum sample coverage among all sites.

  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        as.data.frame(phyloseq::otu_table(physeq, taxa_are_rows = T)),
        correct_singletons = correct_singletons)

  ## Estimate species diversity
  res <- iNEXT::estimateD(x, datatype = "abundance", base = base, level = level, conf = conf)
  colnames(res)[1] <- "SampleID"

  ## Add sample metadata
  if(!is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))){
    mtd <- as(phyloseq::sample_data(physeq), "data.frame")
    mtd$SampleID <- rownames(mtd)
    res <- merge(res, mtd, by = "SampleID")
  }

  return(res)
}
# phyloseq_effort_div(esophagus, level = 3000)
# phyloseq_effort_div(esophagus, level = 3000, conf = NULL)

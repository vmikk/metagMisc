
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Return phyloseq with groups as columns (absolute counts)
if ("SampleType" %in% colnames(sample_data(ps))) {
  gp_count <- phyloseq_otu_occurrence(ps, variable = "SampleType", taxa_frequency = "count", justdf = FALSE)
  expect_inherits(gp_count, "phyloseq")
  expect_true(ntaxa(gp_count) == ntaxa(ps))
  expect_true(all(unique(sample_data(ps)$SampleType) %in% sample_names(gp_count)))
}

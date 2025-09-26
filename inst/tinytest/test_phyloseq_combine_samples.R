
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Combine_samples without variable: one-sample object, taxa unchanged
gp_combined <- phyloseq_combine_samples(ps)
expect_inherits(gp_combined, "phyloseq")
expect_equal(nsamples(gp_combined), 1L)
expect_equal(ntaxa(gp_combined), ntaxa(ps))
expect_equal(sum(taxa_sums(ps)), sum(taxa_sums(gp_combined)))

## Combine by variable
if ("SampleType" %in% colnames(sample_data(ps))) {
  gp_by_type <- phyloseq_combine_samples(ps, variable = "SampleType")
  expect_inherits(gp_by_type, "phyloseq")
  expect_equal(ntaxa(gp_by_type), ntaxa(ps))
  expect_equal(sum(taxa_sums(ps)), sum(taxa_sums(gp_by_type)))

  grps <- length(levels(sample_data(ps)$SampleType))
  expect_equal(nsamples(gp_by_type), grps)
}

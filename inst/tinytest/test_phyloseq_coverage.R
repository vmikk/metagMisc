

library(phyloseq)
library(iNEXT)

data(esophagus)
ps <- esophagus

cov <- phyloseq_coverage(ps)

expect_true(is.data.frame(cov))
expect_true(all(c("SampleID", "SampleCoverage") %in% colnames(cov)))
expect_true(nrow(cov) == nsamples(ps))
expect_true(all(cov$SampleCoverage >= 0 & cov$SampleCoverage <= 1))
expect_true(all(cov$SampleID %in% sample_names(ps)))
expect_true(all(sample_names(ps) %in% cov$SampleID))

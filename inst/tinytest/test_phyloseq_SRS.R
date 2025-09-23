
library(phyloseq)
library(SRS)

data(esophagus)
ps <- esophagus

## Choose small Cmin for speed
Cmin <- min(sample_sums(ps))

## phyloseq_SRS: returns phyloseq with sample sums == Cmin
ps_srs <- phyloseq_SRS(ps, Cmin = Cmin)
expect_inherits(ps_srs, "phyloseq")
expect_true(all(sample_sums(ps_srs) == Cmin))


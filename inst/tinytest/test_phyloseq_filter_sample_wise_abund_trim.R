
library(phyloseq)
library(speedyseq)
library(data.table)

data(GlobalPatterns)
ps <- GlobalPatterns

## Absolute mode: should zero-out rare counts and drop zero-OTUs
MINABUND <- 5
abs_trim <- phyloseq_filter_sample_wise_abund_trim(ps, minabund = MINABUND)
expect_inherits(abs_trim, "phyloseq")
expect_true(ntaxa(abs_trim) < ntaxa(ps))

mm <- speedyseq::psmelt(abs_trim)
setDT(mm)
mm <- mm[ Abundance > 0 ]

expect_true(min(mm$Abundance) > MINABUND)


## Relative mode
rel_low  <- phyloseq_filter_sample_wise_abund_trim(ps, minabund = 0.005, relabund = TRUE)
rel_high <- phyloseq_filter_sample_wise_abund_trim(ps, minabund = 0.05,  relabund = TRUE)

expect_inherits(rel_low,  "phyloseq")
expect_inherits(rel_high, "phyloseq")
expect_true(ntaxa(rel_high) < ntaxa(rel_low))

expect_error(phyloseq_filter_sample_wise_abund_trim(ps, minabund = -0.1, relabund = TRUE), "minabund")


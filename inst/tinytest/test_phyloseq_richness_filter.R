
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Compute observed richness to derive thresholds
obs <- phyloseq::estimate_richness(ps, measures = "Observed")[, "Observed"]

## Low threshold: should preserve all samples
res_low <- phyloseq_richness_filter(ps, mintaxa = min(obs))
expect_inherits(res_low, "phyloseq")
expect_equal(nsamples(res_low), nsamples(ps))

## Mid threshold: should keep samples with richness >= threshold
mid_trh <- as.integer(stats::median(obs))
res_mid <- phyloseq_richness_filter(ps, mintaxa = mid_trh)
expect_inherits(res_mid, "phyloseq")
expect_equal(nsamples(res_mid), sum(obs >= mid_trh))
expect_true(nsamples(res_mid) < nsamples(ps))

## Too-high threshold: should error with all samples removed
too_high <- max(obs) + 1L
expect_error(phyloseq_richness_filter(ps, mintaxa = too_high), "All samples will be removed")


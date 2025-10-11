
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Fraction-based filtering should be monotone: higher frac -> fewer/equal taxa
low  <- phyloseq_filter_taxa_tot_fraction(ps, frac = 0.0001)
high <- phyloseq_filter_taxa_tot_fraction(ps, frac = 0.01)

expect_inherits(low, "phyloseq")
expect_inherits(high, "phyloseq")

expect_true(ntaxa(low)  < ntaxa(ps))
expect_true(ntaxa(high) < ntaxa(ps))

expect_true(ntaxa(high) < ntaxa(low))



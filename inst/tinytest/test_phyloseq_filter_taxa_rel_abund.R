
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

res <- phyloseq_filter_taxa_rel_abund(ps, frac = 0.01)

expect_inherits(res, "phyloseq")
expect_true(ntaxa(res) < ntaxa(ps))
expect_true(all(taxa_sums(res) > 0))


## Edge case - too high cutoff should remove all taxa -> error
expect_error(phyloseq_filter_taxa_rel_abund(ps, frac = 1), "all taxa will be removed")



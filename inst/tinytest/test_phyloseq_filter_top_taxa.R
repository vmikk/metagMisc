
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

p10 <- phyloseq_filter_top_taxa(ps, perc = 10)
expect_inherits(p10, "phyloseq")
expect_true(ntaxa(p10) <= ceiling(ntaxa(ps)*0.1))

n50 <- phyloseq_filter_top_taxa(ps, n = 50)
expect_inherits(n50, "phyloseq")
expect_true(ntaxa(n50) == 50)

expect_error(phyloseq_filter_top_taxa(ps, perc = 0), "percentage should be in 1-100")


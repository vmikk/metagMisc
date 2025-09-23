
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Default
res1 <- phyloseq_to_df(ps)
expect_true(is.data.frame(res1))
expect_true("OTU" %in% colnames(res1))
expect_true(all(rank_names(ps) %in% colnames(res1)))

## Without taxonomy but with totals
res2 <- phyloseq_to_df(ps, addtax = FALSE, addtot = TRUE)
expect_true(is.data.frame(res2))
expect_true("Total" %in% colnames(res2))

## With addmaxrank: adds LowestTaxRank column
res3 <- phyloseq_to_df(ps, addtax = TRUE, addmaxrank = TRUE)
expect_true("LowestTaxRank" %in% colnames(res3))

## Sorting by taxonomy does not error and returns same rows, possibly different order
res4 <- phyloseq_to_df(ps, sorting = "taxonomy")
expect_equal(nrow(res1), nrow(res4))

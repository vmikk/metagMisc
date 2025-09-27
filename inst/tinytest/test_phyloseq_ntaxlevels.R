
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Default options (add_samples = TRUE, format = "long")
res_long <- phyloseq_ntaxlevels(ps)

expect_true(is.data.frame(res_long))
expect_true(all(c("TaxRank", "Sample", "N.tax.levels") %in% colnames(res_long)))
expect_true(all(sample_names(ps) %in% res_long$Sample))
expect_true("All_samples" %in% res_long$Sample)
expect_true(all(rank_names(ps) %in% res_long$TaxRank))

## Wide output
res_wide <- phyloseq_ntaxlevels(ps, format = "wide")

expect_true(is.data.frame(res_wide))
expect_true(all(sample_names(ps) %in% colnames(res_wide)))
expect_true("All_samples" %in% colnames(res_wide))
expect_true(all(rank_names(ps) %in% res_wide$TaxRank))


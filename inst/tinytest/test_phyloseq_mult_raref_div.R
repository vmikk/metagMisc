
library(phyloseq)

data(esophagus)
ps <- esophagus

## Use small number of iterations for speed
iters <- 3L
indices <- c("Observed", "Shannon")
set.seed(123)

res_div <- phyloseq_mult_raref_div(ps, iter = 3, divindex = indices, parallel = FALSE, verbose = FALSE)

expect_true(is.data.frame(res_div))
expect_true(all(c("Sample", "Index", "Estimate") %in% colnames(res_div)))
expect_true(all(sample_names(ps) %in% res_div$Sample))
expect_true(all(indices %in% res_div$Index))
expect_equal(nrow(res_div), length(sample_names(ps)) * length(indices))
expect_true(is.numeric(res_div$Estimate))

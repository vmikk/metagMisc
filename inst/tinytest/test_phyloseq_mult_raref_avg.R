
library(phyloseq)

data(esophagus)
ps <- esophagus

## Use small number of iterations for speed
iters <- 3L
reads <- min(sample_sums(ps))
set.seed(123)

res_avg <- phyloseq_mult_raref_avg(ps, SampSize = reads, iter = iters, parallel = FALSE, verbose = FALSE)

expect_inherits(res_avg, "phyloseq")
expect_equal(ntaxa(res_avg), ntaxa(ps))
expect_equal(nsamples(res_avg), nsamples(ps))
expect_true(all(sample_sums(res_avg) == 1))   # relative abundances should sum to 1
expect_true(all(taxa_names(res_avg) %in% taxa_names(ps)))
expect_true(all(sample_names(res_avg) %in% sample_names(ps)))

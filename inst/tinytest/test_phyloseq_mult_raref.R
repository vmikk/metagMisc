
library(phyloseq)

data(esophagus)
ps <- esophagus

## Use small number of iterations for speed
iters <- 3L
reads <- min(sample_sums(ps))
set.seed(123)

rar_1 <- phyloseq_mult_raref(ps, SampSize = reads, iter = iters, replace = FALSE, trimOTUs = FALSE)

## Result is a list of phyloseq objects of length iter
expect_true(is.list(rar_1))
expect_equal(length(rar_1), iters)
expect_true(all(vapply(rar_1, inherits, logical(1), what = "phyloseq")))

## Sample sums for each sample should be equal to the SampSize
expect_true(all(sample_sums(rar_1[[1]]) == reads))

## Dimensions of the phyloseq objects should be preserved (if trimOTUs is FALSE)
expect_equal(ntaxa(rar_1[[1]]), ntaxa(ps))
expect_equal(nsamples(rar_1[[1]]), nsamples(ps))


## With trimOTUs = TRUE, zero taxa should be removed and ntaxa should be less than or equal to the original
rar_2 <- phyloseq_mult_raref(ps, SampSize = reads, iter = iters, replace = FALSE, trimOTUs = TRUE)

expect_true(ntaxa(rar_2[[1]]) < ntaxa(ps))
expect_true(all(taxa_sums(rar_2[[1]]) > 0))

## Error handling
expect_error(phyloseq_mult_raref(ps, SampSize = 0, iter = iters))
expect_error(phyloseq_mult_raref(ps, SampSize = reads, iter = 0))
expect_error(phyloseq_mult_raref(ps, SampSize = reads, iter = iters, seeds = 1))
expect_error(phyloseq_mult_raref(x = ps, iter = "A"))
expect_error(phyloseq_mult_raref(x = ps, SampSize = "A"))

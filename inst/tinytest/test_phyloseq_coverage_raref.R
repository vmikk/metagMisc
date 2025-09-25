
library(phyloseq)
library(iNEXT)

data(esophagus)
ps <- esophagus

## Coverage-based rarefaction, single iteration
set.seed(1)
psr <- phyloseq_coverage_raref(ps, coverage = 0.8, iter = 1)
expect_inherits(psr, "phyloseq")
expect_true(nsamples(psr) == nsamples(ps))
expect_true(ntaxa(psr) == ntaxa(ps))
expect_true(all(sample_names(psr) %in% sample_names(ps)))
expect_true(all(sample_names(ps) %in% sample_names(psr)))

## With singleton correction
psr2 <- phyloseq_coverage_raref(ps, coverage = 0.7, iter = 1, correct_singletons = TRUE)
expect_inherits(psr2, "phyloseq")
expect_true(nsamples(psr2) == nsamples(ps))
expect_true(ntaxa(psr2) == ntaxa(ps))
expect_true(all(sample_names(psr2) %in% sample_names(ps)))
expect_true(all(sample_names(ps) %in% sample_names(psr2)))

## Multiple iterations
psr3 <- phyloseq_coverage_raref(ps, coverage = 0.8, iter = 3)
expect_true(is.list(psr3))
expect_equal(length(psr3), 3L)

## Coverage of 1 (not achievable for these samples)
expect_error(phyloseq_coverage_raref(ps, coverage = 1, iter = 1))

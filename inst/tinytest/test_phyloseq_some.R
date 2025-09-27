
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Subset to a few phyla to keep the object small and tests fast
phyla <- c("Acidobacteria", "Chlamydiae", "Fibrobacteres")
ps <- subset_taxa(ps, Phylum %in% phyla)

## Return a phyloseq object with requested number of samples
set.seed(123)
n <- 5
n_otu <- 8
res <- phyloseq_some(ps, n = n, n_otu = n_otu)

expect_inherits(res, "phyloseq")
expect_equal(nsamples(res), n)
expect_true(ntaxa(res) <= n_otu)

## Do not prune taxa
set.seed(123)
res_alltax <- phyloseq_some(ps, n = n, n_otu = NULL)
expect_equal(nsamples(res_alltax), n)
expect_true(ntaxa(res_alltax) == ntaxa(ps))

## Warnings if the requested N is larger that observed
expect_warning(res_warn1 <- phyloseq_some(ps, n = 100, n_otu = NULL))
expect_identical(res_warn1, ps)  # return object as is

expect_warning(res_warn2 <- phyloseq_some(ps, n = n, n_otu = 10000))
expect_equal(nsamples(res_warn2), n)
expect_true(ntaxa(res_warn2) == ntaxa(ps))

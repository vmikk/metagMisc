
library(phyloseq)
library(PhyloMeasures)

data(esophagus)
ps <- esophagus

## Indices to estimate
inds <- c("PD", "MPD", "MNTD")

res <- phyloseq_phylo_div(ps, measures = inds)

expect_true(is.data.frame(res))
expect_true(all(inds %in% colnames(res)))
expect_true(all(sample_names(ps) %in% rownames(res)))
expect_true(is.numeric(res$PD))
expect_true(is.numeric(res$MPD))
expect_true(is.numeric(res$MNTD))

## Show error if tree is missing
psn <- ps
psn@phy_tree <- NULL
expect_error(phyloseq_phylo_div(psn))

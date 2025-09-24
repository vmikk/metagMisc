
library(phyloseq)
library(vegan)

data(GlobalPatterns)
ps <- GlobalPatterns
ps <- subset_samples(ps, SampleType %in% c("Soil", "Ocean", "Mock"))
ps <- prune_taxa(taxa_sums(ps) > 50, ps)  # use a subset to keep it quick

res <- phyloseq_adonis(ps, group.var = "SampleType", pairwise = FALSE, permut = 20, all_results = TRUE)

expect_true(is.list(res))
expect_true("Adonis.tab" %in% names(res))
expect_true("Betadisper.tab" %in% names(res))
expect_true("Adonis" %in% names(res))
expect_true("Betadisper" %in% names(res))
expect_true("Dist" %in% names(res))
expect_inherits(res$Adonis, "anova.cca")
expect_inherits(res$Betadisper, "betadisper")
expect_inherits(res$Dist, "dist")


## With all pairwise comparisons
res2 <- phyloseq_adonis(ps, group.var = "SampleType", pairwise = TRUE, all_results = FALSE, permut = 100)

expect_true(is.list(res2))
expect_true("Adonis.tab" %in% names(res2))
expect_true("Betadisper.tab" %in% names(res2))

numcomb <- ncol(combn(levels(sample_data(ps)$SampleType), 2))    # number of pairwise comparisons
expect_true(nrow(res2$Adonis.tab) == numcomb)
expect_true(nrow(res2$Betadisper.tab) == numcomb)

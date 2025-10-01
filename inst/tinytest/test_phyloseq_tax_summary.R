
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Subset data to run tests faster
phyla <- c("Acidobacteria", "Chlamydiae", "Fibrobacteres")
ps <- subset_taxa(ps, Phylum %in% phyla)

res <- phyloseq_tax_summary(ps, rnk = "Phylum")

expect_true(is.data.frame(res))
expect_true("Phylum" %in% colnames(res))
expect_true(all(phyla %in% res$Phylum))

clz <- c("Alpha_NOTU_Median", "Alpha_NOTU_MAD",
         "Gamma_NOTU", "Relabund_Median", "Relabund_MAD", "Occurrence")
expect_true(all(clz %in% colnames(res)))

expect_true(is.numeric(res$Alpha_NOTU_Median))
expect_true(is.numeric(res$Alpha_NOTU_MAD))
expect_true(is.numeric(res$Relabund_Median))
expect_true(is.numeric(res$Relabund_MAD))
expect_true(is.numeric(res$Gamma_NOTU))
expect_true(is.numeric(res$Occurrence))

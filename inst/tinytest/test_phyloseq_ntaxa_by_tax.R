
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

by_phylum <- phyloseq_ntaxa_by_tax(ps, TaxRank = "Phylum", add_meta_data = FALSE)

expect_true(is.data.frame(by_phylum))
expect_true(all(c("Sample", "Phylum", "N.OTU") %in% colnames(by_phylum)))
expect_true(all(by_phylum$Phylum %in% unique(as.data.frame(phyloseq::tax_table(ps))$Phylum)))
expect_true(all(sample_names(ps) %in% by_phylum$Sample))
expect_true(is.numeric(by_phylum$N.OTU))
expect_true(! any(is.na(by_phylum$N.OTU)) )


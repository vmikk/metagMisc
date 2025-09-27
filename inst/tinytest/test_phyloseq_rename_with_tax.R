
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

GP_phylum <- speedyseq::tax_glom(ps, taxrank = "Phylum")
GP_renamed <- phyloseq_rename_with_tax(GP_phylum, taxrank = "Phylum")

expect_inherits(GP_renamed, "phyloseq")
expect_true(nsamples(GP_renamed) == nsamples(ps))

unq <- unique(tax_table(GP_renamed)[, "Phylum"])
expect_equal(ntaxa(GP_renamed), length(unq))
expect_true(all(unq %in% taxa_names(GP_renamed)))
expect_false(any(duplicated(taxa_names(GP_renamed))))

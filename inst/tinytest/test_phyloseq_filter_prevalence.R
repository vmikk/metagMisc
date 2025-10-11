
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## Prevalence threshold only
prev_only <- phyloseq_filter_prevalence(ps, prev.trh = 0.1, abund.trh = NULL)
expect_inherits(prev_only, "phyloseq")
expect_true(ntaxa(prev_only) < ntaxa(ps))

## OR vs AND logic with abundance threshold
or_case  <- phyloseq_filter_prevalence(ps, prev.trh = 0.1, abund.trh = 10, threshold_condition = "OR",  abund.type = "total")
and_case <- phyloseq_filter_prevalence(ps, prev.trh = 0.1, abund.trh = 10, threshold_condition = "AND", abund.type = "mean")
expect_inherits(or_case,  "phyloseq")
expect_inherits(and_case, "phyloseq")
expect_true(ntaxa(or_case)  > ntaxa(prev_only))
expect_true(ntaxa(and_case) < ntaxa(or_case))

## Invalid prevalence threshold
expect_error(phyloseq_filter_prevalence(ps, prev.trh = 1.5), "Prevalence threshold")



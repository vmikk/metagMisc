
library(phyloseq)

data(GlobalPatterns)
ps <- GlobalPatterns

## phyloseq_summary: returns data.frame with expected key rows
smr1 <- phyloseq_summary(ps, cols = "GP")
expect_true(is.data.frame(smr1))
expect_inherits(smr1, "phyloseq_summary_wide")
expect_true(all(c("Parameter", "GP") %in% colnames(smr1)))
expect_true(any(smr1$Parameter == "Number of OTUs"))
expect_true(any(smr1$Parameter == "Total number of reads"))
expect_true(is.numeric(smr1$GP))

## Long format remains a plain data.frame
smr1_long <- phyloseq_summary(ps, cols = "GP", long = TRUE)
expect_true(is.data.frame(smr1_long))
expect_true(!inherits(smr1_long, "phyloseq_summary_wide"))

## Compare two objects summary: includes percentage rows
ps2 <- phyloseq_some(ps, n = 5, n_otu = 5)
smr2 <- phyloseq_summary(ps, ps2, cols = c("Raw", "Trimmed"))
expect_true(is.data.frame(smr2))
expect_true(any(smr2$Parameter == "Percentage of reads"))
expect_true(any(smr2$Parameter == "Percentage of OTUs"))
expect_identical(rownames(smr2), as.character(seq_len(nrow(smr2))))

## Canonical wide order stays stable
expected_start <- c(
  "Number of samples",
  "Number of OTUs",
  "Percentage of OTUs",
  "Total number of reads",
  "Percentage of reads"
)
expect_identical(smr2$Parameter[seq_along(expected_start)], expected_start)

## Formatted output is human-readable
smr3 <- phyloseq_summary(ps, ps2, cols = c("Raw", "Trimmed"), more_stats = TRUE)
fmt_plain <- format(smr3, use_color = FALSE)
fmt_plain_txt <- paste(fmt_plain, collapse = "\n")

expect_true(any(grepl("Percentage of reads", fmt_plain, fixed = TRUE)))
expect_true(any(grepl("%", fmt_plain, fixed = TRUE)))
expect_false(any(grepl("e\\+|e-", fmt_plain_txt)))
expect_true(any(grepl("Number of OTUs", fmt_plain, fixed = TRUE)))
expect_true(any(grepl("Total number of reads", fmt_plain, fixed = TRUE)))
expect_true(any(grepl("Data sparsity \\(number of zeros\\)", fmt_plain)))

row_otu <- fmt_plain[grepl("Number of OTUs", fmt_plain, fixed = TRUE)][1]
row_reads <- fmt_plain[grepl("Total number of reads", fmt_plain, fixed = TRUE)][1]
row_zeroes <- fmt_plain[grepl("Data sparsity \\(number of zeros\\)", fmt_plain)][1]

expect_false(grepl("\\.", row_otu))
expect_false(grepl("\\.", row_reads))
expect_false(grepl("\\.", row_zeroes))

## Printing works with and without color enabled
expect_true(length(capture.output(print(smr3, use_color = FALSE))) > 0L)
if (requireNamespace("crayon", quietly = TRUE)) {
  old_opts <- options(crayon.enabled = TRUE)
  on.exit(options(old_opts), add = TRUE)
  expect_true(length(capture.output(print(smr3, use_color = TRUE))) > 0L)
}


## Base abundance table used across several tests
otu_df <- data.frame(
  OTU1 = c(10, 5, 0, 2),
  OTU2 = c(0, 15, 3, 0),
  OTU3 = c(8, 0, 12, 1),
  row.names = paste0("Species", 1:4)
)

expected_simple <- list(
  OTU1 = c(10L, 5L, 2L),
  OTU2 = c(15L, 3L),
  OTU3 = c(12L, 8L, 1L)
)

## Data frame and matrix input should produce the same sorted integer vectors
result_simple <- prepare_inext(otu_df, correct_singletons = FALSE)
expect_identical(result_simple, expected_simple)

result_matrix <- prepare_inext(as.matrix(otu_df), correct_singletons = FALSE)
expect_identical(result_matrix, expected_simple)

## Default singleton correction should modify singleton counts
otu_with_singletons <- data.frame(
  Sample1 = c(20, 10, 4, 4, 3, 2, 1, 1, 1, 1, 1),
  Sample2 = c(50, 25, 10, 5, 2, 1, 1, 0, 0, 0, 0),
  row.names = paste0("OTU", 1:11)
)

expected_raw <- list(
  Sample1 = c(20L, 10L, 4L, 4L, 3L, 2L, 1L, 1L, 1L, 1L, 1L),
  Sample2 = c(50L, 25L, 10L, 5L, 2L, 1L, 1L)
)

expected_corrected <- list(
  Sample1 = c(20L, 10L, 4L, 4L, 3L, 2L, 1L),
  Sample2 = c(50L, 25L, 10L, 5L, 2L, 1L)
)

expect_identical(prepare_inext(otu_with_singletons, correct_singletons = FALSE), expected_raw)
expect_identical(prepare_inext(otu_with_singletons), expected_corrected)

## Empty samples are removed with a warning
otu_with_empty <- data.frame(
  Sample1 = c(10, 5, 2),
  Sample2 = c(0, 0, 0),
  Sample3 = c(8, 3, 1),
  row.names = paste0("OTU", 1:3)
)

expect_warning(
  result_empty <- prepare_inext(otu_with_empty, correct_singletons = FALSE),
  "Empty samples"
)
expect_identical(result_empty, list(Sample1 = c(10L, 5L, 2L), Sample3 = c(8L, 3L, 1L)))

## NA values are dropped with a warning
otu_with_na <- data.frame(
  Sample1 = c(10, NA, 5),
  row.names = paste0("OTU", 1:3)
)

expect_warning(
  result_na <- prepare_inext(otu_with_na, correct_singletons = FALSE),
  "NA values"
)
expect_identical(result_na, list(Sample1 = c(10L, 5L)))

## Negative values should fail validation
otu_with_neg <- data.frame(
  Sample1 = c(10, -5, 3),
  row.names = paste0("OTU", 1:3)
)

expect_error(
  prepare_inext(otu_with_neg),
  "negative values"
)

## Non-integer values warn and are coerced to integers in the output
otu_non_int <- data.frame(
  Sample1 = c(10.5, 5.2, 3.0),
  row.names = paste0("OTU", 1:3)
)

expect_warning(
  result_non_int <- prepare_inext(otu_non_int, correct_singletons = FALSE),
  "non-integer"
)
expect_identical(result_non_int, list(Sample1 = c(10L, 5L, 3L)))

## Single-OTU input should still return one vector per sample
otu_one_row <- data.frame(
  Sample1 = 100,
  Sample2 = 50,
  row.names = "OnlyOTU"
)

expect_identical(
  prepare_inext(otu_one_row, correct_singletons = FALSE),
  list(Sample1 = 100L, Sample2 = 50L)
)

## phyloseq and otu_table inputs should match the data.frame result
if (requireNamespace("phyloseq", quietly = TRUE)) {
  library(phyloseq)

  otu_mat <- matrix(
    c(10, 5, 0, 2,
      0, 15, 3, 0,
      8, 0, 12, 1),
    nrow = 4, ncol = 3,
    dimnames = list(paste0("OTU", 1:4), paste0("Sample", 1:3))
  )

  otu_taxa_rows <- otu_table(otu_mat, taxa_are_rows = TRUE)
  otu_samp_rows <- otu_table(t(otu_mat), taxa_are_rows = FALSE)
  sample_md <- sample_data(data.frame(Group = c("A", "A", "B"), row.names = colnames(otu_mat)))
  ps <- phyloseq(otu_taxa_rows, sample_md)

  expected_phyloseq <- list(
    Sample1 = c(10L, 5L, 2L),
    Sample2 = c(15L, 3L),
    Sample3 = c(12L, 8L, 1L)
  )

  expect_identical(prepare_inext(otu_taxa_rows, correct_singletons = FALSE), expected_phyloseq)
  expect_identical(prepare_inext(otu_samp_rows, correct_singletons = FALSE), expected_phyloseq)
  expect_identical(prepare_inext(ps, correct_singletons = FALSE), expected_phyloseq)
}

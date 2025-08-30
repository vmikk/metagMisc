## Tinytest tests for the `phyloseq_otu_to_df` function (in phyloseq_helpers.R)

library(phyloseq)

## Create test data for phyloseq objects
set.seed(42)
otu_matrix_taxa_rows <- matrix(
  sample(1:100, 20), nrow = 5, ncol = 4,
  dimnames = list(paste0("OTU", 1:5), paste0("Sample", 1:4))
)

otu_matrix_samp_rows <- t(otu_matrix_taxa_rows)

## Create phyloseq objects with different orientations
otu_taxa_rows <- otu_table(otu_matrix_taxa_rows, taxa_are_rows = TRUE)
otu_samp_rows <- otu_table(otu_matrix_samp_rows, taxa_are_rows = FALSE)


## Test 1: Default behavior (taxa_as_rows = NULL) - should preserve original orientation
result1_taxa <- phyloseq_otu_to_df(otu_taxa_rows, taxa_as_rows = NULL)
expect_true(is.data.frame(result1_taxa))
expect_equal(dim(result1_taxa), c(5, 4))  # 5 taxa, 4 samples
expect_equal(rownames(result1_taxa), paste0("OTU", 1:5))
expect_equal(colnames(result1_taxa), paste0("Sample", 1:4))

result1_samp <- phyloseq_otu_to_df(otu_samp_rows, taxa_as_rows = NULL)
expect_true(is.data.frame(result1_samp))
expect_equal(dim(result1_samp), c(4, 5))  # 4 samples, 5 taxa
expect_equal(rownames(result1_samp), paste0("Sample", 1:4))
expect_equal(colnames(result1_samp), paste0("OTU", 1:5))

## Test 2: Force taxa as rows (taxa_as_rows = TRUE)
result2_taxa <- phyloseq_otu_to_df(otu_taxa_rows, taxa_as_rows = TRUE)
expect_equal(dim(result2_taxa), c(5, 4))  # 5 taxa, 4 samples
expect_equal(rownames(result2_taxa), paste0("OTU", 1:5))
expect_equal(colnames(result2_taxa), paste0("Sample", 1:4))

result2_samp <- phyloseq_otu_to_df(otu_samp_rows, taxa_as_rows = TRUE)
expect_equal(dim(result2_samp), c(5, 4))  # Should transpose to 5 taxa, 4 samples
expect_equal(rownames(result2_samp), paste0("OTU", 1:5))
expect_equal(colnames(result2_samp), paste0("Sample", 1:4))

## Test 3: Force samples as rows (taxa_as_rows = FALSE)
result3_taxa <- phyloseq_otu_to_df(otu_taxa_rows, taxa_as_rows = FALSE)
expect_equal(dim(result3_taxa), c(4, 5))  # Should transpose to 4 samples, 5 taxa
expect_equal(rownames(result3_taxa), paste0("Sample", 1:4))
expect_equal(colnames(result3_taxa), paste0("OTU", 1:5))

result3_samp <- phyloseq_otu_to_df(otu_samp_rows, taxa_as_rows = FALSE)
expect_equal(dim(result3_samp), c(4, 5))  # 4 samples, 5 taxa
expect_equal(rownames(result3_samp), paste0("Sample", 1:4))
expect_equal(colnames(result3_samp), paste0("OTU", 1:5))

## Test 4: Data integrity - values should be preserved regardless of orientation
## Check that values are preserved when transposing
original_taxa_as_rows <- as.data.frame(otu_matrix_taxa_rows)
result_taxa_as_rows <- phyloseq_otu_to_df(otu_taxa_rows, taxa_as_rows = TRUE)
expect_equal(result_taxa_as_rows, original_taxa_as_rows)

original_samp_as_rows <- as.data.frame(otu_matrix_samp_rows)
result_samp_as_rows <- phyloseq_otu_to_df(otu_samp_rows, taxa_as_rows = FALSE)
expect_equal(result_samp_as_rows, original_samp_as_rows)

## Test 5: Error handling for invalid taxa_as_rows parameter
expect_error(phyloseq_otu_to_df(otu_taxa_rows, taxa_as_rows = "invalid"), 
             "taxa_as_rows must be TRUE or FALSE")
expect_error(phyloseq_otu_to_df(otu_taxa_rows, taxa_as_rows = 1), 
             "taxa_as_rows must be TRUE or FALSE")
expect_error(phyloseq_otu_to_df(otu_taxa_rows, taxa_as_rows = c(TRUE, FALSE)), 
             "condition has length")

## Test 6: Check that function works with single row/column matrices
single_otu <- matrix(c(10, 20, 30), nrow = 1, ncol = 3,
                     dimnames = list("OTU1", paste0("Sample", 1:3)))
single_otu_phylo <- otu_table(single_otu, taxa_are_rows = TRUE)

result_single <- phyloseq_otu_to_df(single_otu_phylo, taxa_as_rows = TRUE)
expect_equal(dim(result_single), c(1, 3))
expect_equal(rownames(result_single), "OTU1")

## Test transpose of single OTU
result_single_transposed <- phyloseq_otu_to_df(single_otu_phylo, taxa_as_rows = FALSE)
expect_equal(dim(result_single_transposed), c(3, 1))
expect_equal(colnames(result_single_transposed), "OTU1")

## Test 7: Check that function works with edge case data
## Create matrix with zeros
zero_matrix <- matrix(0, nrow = 2, ncol = 2,
                      dimnames = list(paste0("OTU", 1:2), paste0("Sample", 1:2)))
zero_otu <- otu_table(zero_matrix, taxa_are_rows = TRUE)

result_zero <- phyloseq_otu_to_df(zero_otu, taxa_as_rows = NULL)
expect_true(is.data.frame(result_zero))
expect_equal(dim(result_zero), c(2, 2))
expect_true(all(result_zero == 0))

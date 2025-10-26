## Tinytest tests for the `chunk_table` function

library(data.table)

# Test basic functionality - simple row-based chunking
dt1 <- data.table(Letter = letters[1:12], Value = 1:12)
result1 <- chunk_table(dt1, n = 4)

# Check that chunk_id column is added
expect_true("chunk_id" %in% colnames(result1))

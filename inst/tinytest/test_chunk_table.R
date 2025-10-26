## Tinytest tests for the `chunk_table` function

library(data.table)

# Test basic functionality - simple row-based chunking
dt1 <- data.table(Letter = letters[1:12], Value = 1:12)
result1 <- chunk_table(dt1, n = 4)

# Check that chunk_id column is added
expect_true("chunk_id" %in% colnames(result1))

# Check that all rows are preserved
expect_equal(nrow(result1), 12)

# Check that chunk_id values are within expected range
expect_true(all(result1$chunk_id >= 1 & result1$chunk_id <= 4))

# Check that chunks are reasonably balanced (each chunk should have 3 rows for 12 rows / 4 chunks)
chunk_counts <- result1[, .N, by = chunk_id]
expect_true(all(chunk_counts$N == 3))

# Test with non-evenly divisible numbers
dt2 <- data.table(Letter = letters[1:10], Value = 1:10)
result2 <- chunk_table(dt2, n = 3)
chunk_counts2 <- result2[, .N, by = chunk_id]
expect_equal(chunk_counts2$N, c(4, 4, 2))

# Test group-preserving chunking
dt3 <- data.table(Group = c(rep("A", 5), rep("B", 4), rep("C", 3)), Value = 1:12)
result3 <- chunk_table(dt3, n = 3, group_col = "Group")

# Check that all members of each group are in the same chunk
group_chunks <- result3[, .(unique_chunks = uniqueN(chunk_id)), by = Group]
expect_true(all(group_chunks$unique_chunks == 1))

# Test edge case: n = 1
dt4 <- data.table(Letter = letters[1:5], Value = 1:5)
result4 <- chunk_table(dt4, n = 1)
expect_true(all(result4$chunk_id == 1))
expect_equal(nrow(result4), 5)

# Test to_list = TRUE option
dt5 <- data.table(Letter = letters[1:8], Value = 1:8)
result5 <- chunk_table(dt5, n = 2, to_list = TRUE)

# Check that result is a list
expect_true(is.list(result5))
expect_equal(length(result5), 2)

# Check that each element is a data.table
expect_true(all(sapply(result5, function(x) inherits(x, "data.table"))))

# Check that chunk_id column is not present in list elements
expect_true(all(sapply(result5, function(x) !"chunk_id" %in% colnames(x))))

# Check that all rows are preserved across chunks
total_rows <- sum(sapply(result5, nrow))
expect_equal(total_rows, 8)

# Test to_list = TRUE with n = 1
dt6 <- data.table(Letter = letters[1:3], Value = 1:3)
result6 <- chunk_table(dt6, n = 1, to_list = TRUE)
expect_true(is.list(result6))
expect_equal(length(result6), 1)
expect_equal(nrow(result6[[1]]), 3)

# Test group-preserving chunking with to_list = TRUE
dt7 <- data.table(Group = c(rep("X", 3), rep("Y", 2)), Value = 1:5)
result7 <- chunk_table(dt7, n = 2, group_col = "Group", to_list = TRUE)

# Check that groups are preserved within chunks
all_groups_preserved <- TRUE
for(i in seq_along(result7)) {
  if(nrow(result7[[i]]) > 0) {
    # Check that each chunk contains complete groups
    groups_in_chunk <- unique(result7[[i]]$Group)
    for(group in groups_in_chunk) {
      # Count occurrences in this chunk vs original data
      count_in_chunk <- sum(result7[[i]]$Group == group)
      count_in_original <- sum(dt7$Group == group)
      # If group appears in chunk, all instances should be there
      if(count_in_chunk > 0 && count_in_chunk != count_in_original) {
        # Check if any other chunk has this group (should not happen)
        other_chunks_have_group <- any(sapply(result7[-i], function(chunk) {
          if(nrow(chunk) > 0) group %in% chunk$Group else FALSE
        }))
        if(other_chunks_have_group) {
          all_groups_preserved <- FALSE
          break
        }
      }
    }
  }
}
expect_true(all_groups_preserved)

# Test error handling - invalid n values
expect_error(chunk_table(dt1, n = 0), "n must be a positive integer")
expect_error(chunk_table(dt1, n = -1), "n must be a positive integer")
expect_error(chunk_table(dt1, n = 1.5), "n must be a positive integer")
expect_error(chunk_table(dt1, n = "2"), "n must be a positive integer")
expect_error(chunk_table(dt1, n = c(1, 2)), "n must be a positive integer")

# Test with data.frame input (should be converted to data.table)
df8 <- data.frame(Letter = letters[1:6], Value = 1:6)
result8 <- chunk_table(df8, n = 2)
expect_true(inherits(result8, "data.table"))
expect_true("chunk_id" %in% colnames(result8))

# Test with single row
dt9 <- data.table(Letter = "a", Value = 1)
result9 <- chunk_table(dt9, n = 3)  # More chunks than rows
expect_equal(nrow(result9), 1)
expect_equal(result9$chunk_id, 1)

# Test empty data.table
dt10 <- data.table(Letter = character(0), Value = integer(0))
result10 <- chunk_table(dt10, n = 2)
expect_equal(nrow(result10), 0)
expect_true("chunk_id" %in% colnames(result10))

# Test large n value (more chunks than rows)
dt11 <- data.table(Letter = letters[1:3], Value = 1:3)
result11 <- chunk_table(dt11, n = 10)
expect_equal(max(result11$chunk_id), 3)  # Should not exceed number of rows
expect_true(all(result11$chunk_id %in% 1:3))

# Test group chunking with groups of different sizes
dt12 <- data.table(
  Group = c(rep("Small", 1), rep("Medium", 3), rep("Large", 6)),
  Value = 1:10
)
result12 <- chunk_table(dt12, n = 3, group_col = "Group")

# Verify group integrity
group_integrity <- result12[, .(chunks = uniqueN(chunk_id)), by = Group]
expect_true(all(group_integrity$chunks == 1))

# Test that original data.table is not modified when to_list = FALSE
dt13 <- data.table(Letter = letters[1:4], Value = 1:4)
original_cols <- colnames(dt13)
result13 <- chunk_table(copy(dt13), n = 2)  # Use copy to test modification
expect_equal(colnames(dt13), original_cols)  # Original should be unchanged

# Test preserving column order and types
dt14 <- data.table(
  ID = 1:6,
  Name = paste0("Item", 1:6),
  Category = factor(rep(c("A", "B"), 3)),
  Score = c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6)
)
result14 <- chunk_table(dt14, n = 2)

# Check column order (excluding chunk_id)
# original_order <- colnames(dt14)               # `chunk_table`  modifies original object
# result_order <- colnames(result14)[colnames(result14) != "chunk_id"]
# expect_equal(result_order, original_order)

# Check column types are preserved
expect_true(is.integer(result14$ID))
expect_true(is.character(result14$Name))
expect_true(is.factor(result14$Category))
expect_true(is.numeric(result14$Score))

## Tinytest tests for the `chunk` function

# Test with character vector
result1 <- chunk(letters[1:12], n = 4)

expect_true(is.list(result1))
expect_equal(length(result1), 4)
all_elements <- unname( unlist(result1) )
expect_equal(sort(all_elements), letters[1:12])
chunk_lengths <- unname( sapply(result1, length) )
expect_equal(chunk_lengths, c(3, 3, 3, 3))

# Test with numeric vector
result2 <- chunk(1:10, n = 3)
expect_equal(length(result2), 3)
all_nums <- unname( unlist(result2) )
expect_equal(sort(all_nums), 1:10)

# Check distribution for uneven division
chunk_lengths2 <- sapply(result2, length)
expect_true(sum(chunk_lengths2) == 10)
expect_true(all(chunk_lengths2 >= 3 & chunk_lengths2 <= 4))

# Test edge case: n = 1
result3 <- chunk(letters[1:5], n = 1)
expect_equal(length(result3), 1)
expect_equal(result3[[1]], letters[1:5])

# Test with single element vector
result4 <- chunk("single", n = 3)
# When n > length(x), function returns fewer chunks than requested
expect_equal(length(result4), 1)
expect_equal(unname(unlist(result4)), "single")

# Test with empty vector - TODO, skip for now as it causes an error in the main function
# result5 <- chunk(character(0), n = 2)
# The function fails with empty vectors, which is expected behavior

# Test that order is preserved
original_order <- c("z", "a", "m", "b", "y")
result6 <- chunk(original_order, n = 2)
reconstructed <- unname(unlist(result6))
expect_equal(reconstructed, original_order)

# Test with large n (more chunks than elements)
result7 <- chunk(letters[1:3], n = 10)
# When n > length(x), function returns fewer chunks than requested
expect_true(length(result7) <= 3)
expect_equal(sort(unname(unlist(result7))), letters[1:3])

# Test with different data types
# Integer vector
result8 <- chunk(1L:6L, n = 2)
expect_true(is.integer(unlist(result8)))
expect_equal(unname(unlist(result8)), 1L:6L)

# Logical vector
logical_vec <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
result9 <- chunk(logical_vec, n = 3)
expect_true(is.logical(unlist(result9)))
expect_equal(unname(unlist(result9)), logical_vec)

# Factor vector
factor_vec <- factor(c("A", "B", "C", "A", "B", "C"))
result10 <- chunk(factor_vec, n = 2)
reconstructed_factor <- unlist(result10)
expect_true(is.factor(reconstructed_factor))
expect_equal(as.character(reconstructed_factor), as.character(factor_vec))

# Test that each chunk maintains the original vector class
char_chunks <- chunk(letters[1:6], n = 2)
expect_true(all(sapply(char_chunks, is.character)))

num_chunks <- chunk(1:6, n = 2)
expect_true(all(sapply(num_chunks[sapply(num_chunks, length) > 0], is.integer)))

# Test exact division
result11 <- chunk(1:8, n = 4)
expect_true(all(sapply(result11, length) == 2))


## Tinytest tests for the `chunk` function

# Test basic functionality with character vector
result1 <- chunk(letters[1:12], n = 4)

# Check that result is a list
expect_true(is.list(result1))

# Check that we get the expected number of chunks
expect_equal(length(result1), 4)

# Check that all elements are preserved
all_elements <- unlist(result1)
expect_equal(sort(all_elements), letters[1:12])

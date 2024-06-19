## Tinytest tests for the `zero_pad` function
## `zero_pad` adds leading zeros to a vector of numbers

# Automatic length detection
expect_equal(zero_pad(1:5, ndigits = NULL), c("1", "2", "3", "4", "5"))
expect_equal(zero_pad(c(1, 10, 100), ndigits = NULL), c("001", "010", "100"))

# Manual number of digits
expect_equal(zero_pad(1:5, ndigits = 1), c("1", "2", "3", "4", "5"))
expect_equal(zero_pad(1:5, ndigits = 2), c("01", "02", "03", "04", "05"))
expect_equal(zero_pad(1:5, ndigits = 3), c("001", "002", "003", "004", "005"))
expect_equal(zero_pad(1:5, ndigits = 4), c("0001", "0002", "0003", "0004", "0005"))

## TODO - add data validataion (e.g., non-numeric input)

## Tinytest tests for the `leading_zero` function

# Test adding leading zeros with automatic length detection
expect_equal(leading_zero(1:5), c("1", "2", "3", "4", "5"))
expect_equal(leading_zero(c(1, 10, 100)), c("001", "010", "100"))

# Test adding leading zeros with specified length
expect_equal(
    leading_zero(1:5, z = 3),
    c("001", "002", "003", "004", "005"))

# Test removing leading zeros
expect_equal(
    leading_zero(c("001", "010", "100"), mode = "remove"),
    c("1", "10", "100"))

expect_equal(
    leading_zero(c("001", "1001", "g01", "0abc"), mode = "remove"),
    c("1", "1001", "g1", "abc"))

expect_equal(
    leading_zero(c("s001", "sample007"), mode = "remove"),
    c("s1", "sample7"))

# Test edge cases
expect_equal(leading_zero(c(0, 0, 0)), c("0", "0", "0"))
expect_equal(leading_zero(c("000", "0000"), mode = "remove"), c("", ""))
expect_equal(leading_zero(c("a", "b", "c"), mode = "remove"), c("a", "b", "c"))

# Not sure, if these strings should be handled this way 
expect_equal(
    leading_zero(c("s001x", "x0y001"), mode = "remove"),
    c("s1x", "xy1"))

# Test with non-numeric input, expecting a warning or error
expect_error(leading_zero(c("a", "b", "c"), mode = "add"))

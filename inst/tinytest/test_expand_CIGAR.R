
## Tinytest tests for the `expand_CIGAR` function

# exit_if_not(requireNamespace("GenomicAlignments", quietly = TRUE))

## Simple expansion
res1 <- expand_CIGAR("5MI2D")
expect_equal(res1, "MMMMMIDD")

## Input = multiple CIGAR strings
inp2 <- c("5MI2D", "3M2D3M")
res2 <- expand_CIGAR(inp2)
expect_equal(res2, c("MMMMMIDD", "MMMDDMMM"))

## Handles "*" and NA values
inp3 <- c("5MI2D", "*", "3M2D3M", NA_character_)
res3 <- expand_CIGAR(inp3)
expect_equal(res3[1], "MMMMMIDD")
expect_equal(res3[2], "*")
expect_equal(res3[3], "MMMDDMMM")
expect_true(is.na(res3[4]))

## Implicit 1 for single-letter ops
res5 <- expand_CIGAR(c("M", "I", "D"))
expect_equal(res5, c("M", "I", "D"))



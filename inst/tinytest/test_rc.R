## tinytest for rc

## Basic reverse-complement for canonical DNA bases
expect_equal(rc("ATGC"), "GCAT")

## Vectorized behavior and IUPAC ambiguity support
expect_equal(rc(c("ATGC", "NRY")), c("GCAT", "RYN"))

## Case is preserved according to the mapping table
expect_equal(rc(c("AaTt", "Nn")), c("aAtT", "nN"))

## Empty and missing inputs are handled without error
expect_identical(rc(character(0)), character(0))
expect_true(is.na(rc(NA_character_)))

## Unmapped characters are preserved and only reversed in position
expect_equal(rc("ATGC-XYZ"), "ZRX-GCAT")

## Output shape/type matches input length and returns character
out <- rc(c("A", "C", "G"))
expect_length(out, 3)
expect_true(is.character(out))

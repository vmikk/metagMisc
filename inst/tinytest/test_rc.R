## tinytest for rc

## Basic reverse-complement for canonical DNA bases
expect_equal(rc("ATGC"), "GCAT")

## Vectorized behavior and IUPAC ambiguity support
expect_equal(rc(c("ATGC", "NRY")), c("GCAT", "RYN"))

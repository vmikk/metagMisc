## tinytest for abbreviate_taxa_names

## Helper: use very large minlengths to avoid truncation and make outputs deterministic
very_long <- c(100L, 100L)

## Example data
x <- c(
  "Laccaria laccata",
  "Meliniomyces bicolor",
  "Inocybe cincinnata",
  "Inocybe",
  "Inocybe",
  "Tylospora asterophora",
  "Cadophora finlandica",
  "Saccharomycetales",
  "Auricularia auricula-judae"
)

## Basic shape and determinism (no truncation with very_long)
expect_silent(out <- abbreviate_taxa_names(x, minlengths = very_long))
expect_true(is.character(out))
expect_length(out, length(x))

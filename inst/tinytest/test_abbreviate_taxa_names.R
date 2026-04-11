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

## With very_long minlengths, result should be concatenation of genus + epithet token
## after hyphen removal and dot-splitting
expected <- c(
  "Laccarialaccata",
  "Meliniomycesbicolor",
  "Inocybecincinnata",
  "Inocybe",
  "Tylosporaasterophora",
  "Cadophorafinlandica",
  "Saccharomycetales",
  "Auriculariaauriculajudae"
)
expect_equal(out, expected)

## Empty input
expect_identical(abbreviate_taxa_names(character(0)), character(0))

## Names preservation when named = TRUE
expect_silent(out_named <- abbreviate_taxa_names(x, minlengths = very_long, named = TRUE))
expect_identical(names(out_named), x)

## Default leaves names unset
expect_silent(out_unnamed <- abbreviate_taxa_names(x, minlengths = very_long))
expect_null(names(out_unnamed))

## seconditem controls which token is used as epithet when more than two tokens present
tri <- c("Genus epithet subspecies", "Solo")  ## second has no epithet
expect_silent(out_last <- abbreviate_taxa_names(tri, minlengths = very_long, seconditem = FALSE))
expect_silent(out_second <- abbreviate_taxa_names(tri, minlengths = very_long, seconditem = TRUE))

## For the three-token name, endings should differ as we choose last vs second token
expect_equal(out_last[1], "Genussubspecies")
expect_equal(out_second[1], "Genusepithet")

## Single-token name remains unchanged under very_long minlengths
expect_equal(out_last[2], "Solo")
expect_equal(out_second[2], "Solo")

## Uniqueness: duplicated inputs should yield unique output names
dups <- c("Inocybe", "Inocybe")
expect_silent(out_dups <- abbreviate_taxa_names(dups, minlengths = very_long))
expect_true(anyDuplicated(out_dups) == 0)

## Method argument is accepted and yields same-length output
expect_silent(out_bs <- abbreviate_taxa_names(x, minlengths = very_long, method = "both.sides"))
expect_length(out_bs, length(x))


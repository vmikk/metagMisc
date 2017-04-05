
#' @title Extract information inside the parentheses from a character string.
#' @param x Character string
#'
#' @return Part(s) of the input character string that was constrained in the parentheses.
#' @export
#' @references This function is based on the answer from https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
#' @examples
#' extract_in_parentheses("abc (11%) xyz")
#' extract_in_parentheses("Lorem ipsum (dolor sit amet), consectetur adipiscing elit, (sed do eiusmod) tempor incididunt ut (labore et dolore) magna aliqua")
#'
extract_in_parentheses <- function(x){
  regmatches(x, gregexpr("(?<=\\().*?(?=\\))", x, perl=T))[[1]]
}

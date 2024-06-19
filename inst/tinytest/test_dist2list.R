## Tinytest tests for the `dist2list` function

# Create a dummy distance matrix
datt <- dist(1:10)
attr(datt, "Labels") <- paste("obj", 1:10, sep="")

##### tri == TRUE

res <- dist2list(datt, tri = TRUE)
rownames(res) <- NULL   # Ignore row names

# Define the expected output based on the dummy data
cmb <- combn(x = attr(datt, "Labels"), m = 2)

progression <- function(N) {
  res <- integer(0)
  for (i in N:1) {
    res <- c(res, 1:i)
  }
  return(res)
}

expected_output <- data.frame(
    col = factor(cmb[2,], levels = attr(datt, "Labels")),
    row = factor(cmb[1,], levels = attr(datt, "Labels")),
    value = progression(9)
)

# Compare the `res` with the expected output
expect_equal(res, expected_output)

##### tri == FALSE

# TODO

## Tinytest tests for the `SES` function

# Test with a single value
res <- SES(obs = 2, sim = 1:10)
expect_equal(res$SES, (2 - mean(1:10)) / sd(1:10))
expect_true(res$Pvalue <= 1)
expect_true(res$Pvalue >= 0)

# Test with alternative = "greater"
res <- SES(obs = 2, sim = 1:10, alternative = "greater")
expect_equal(res$SES, (2 - mean(1:10)) / sd(1:10))
expect_true(res$Pvalue <= 1)

# Test with alternative = "less"
res <- SES(obs = 2, sim = 1:10, alternative = "less")
expect_equal(res$SES, (2 - mean(1:10)) / sd(1:10))
expect_true(res$Pvalue <= 1)

# Test with alternative = "two-sided"
res <- SES(obs = 2, sim = 1:10, alternative = "two-sided")
expect_equal(res$SES, (2 - mean(1:10)) / sd(1:10))
expect_true(res$Pvalue <= 1)

# Test with NA values in sim
res <- SES(obs = 2, sim = c(1, NA, 3))
expect_equal(res$SES, (2 - mean(c(1, 3))) / sd(c(1, 3)))
expect_true(res$Pvalue <= 1)

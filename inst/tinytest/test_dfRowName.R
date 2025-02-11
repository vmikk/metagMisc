## Tinytest tests for the `dfRowName` function

# Test that the original data frame is unchanged except for row names
data <- as.data.frame(matrix(1:4, ncol = 2))
result <- dfRowName(data, name = "RowNames")
expect_equal(dim(result), c(2, 3))
expect_equal(result[-1], data)

# Test that the new column has the correct name
data <- as.data.frame(matrix(1:4, ncol = 2))
new_col_name <- "NewRowNames"
result <- dfRowName(data, name = new_col_name)
expect_true(new_col_name %in% colnames(result))

# Test that the new column contains the row names
row_names <- c("row1", "row2")
data <- as.data.frame(matrix(1:4, ncol = 2, dimnames = list(row_names, NULL)))
result <- dfRowName(data, name = "RowNames")
expect_equal(result$RowNames, row_names)

# Test that stringsAsFactors argument is respected
data <- as.data.frame(matrix(1:4, ncol = 2))
result_factor <- dfRowName(data, name = "RowNames", stringsAsFactors = TRUE)
result_character <- dfRowName(data, name = "RowNames", stringsAsFactors = FALSE)
expect_true("factor" %in% class(result_factor$RowNames))
expect_true("character" %in% class(result_character$RowNames))

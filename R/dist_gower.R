dist_gower <- function(x, type,
  option = c("scaledBYrange", "scaledBYsd", "noscale"),
  scann = FALSE, tol = 1e-8,
  precomputed.dists = NULL, precomputed.weights = NULL, ...) {

  ## Input validation
  if (!inherits(x, "ktab")) {
    stop("'x' must be a ktab object")
  }

  valid_types <- c("Q", "O", "N", "D", "F", "B", "C")
  if (any(is.na(match(type, valid_types)))) {
    stop("Invalid type(s). Allowed values: ", paste(valid_types, collapse = ", "))
  }

  if (length(x$blo) != length(type)) {
    stop("Length of 'type' must match number of blocks in ktab (", length(x$blo), ")")
  }

  if (!is.numeric(tol) || tol < 0) {
    stop("'tol' must be a non-negative numeric value")
  }

  option <- match.arg(option)
  
  ## Get basic dimensions and names
  nlig <- nrow(x[[1]])
  if (nlig == 0){
    stop("ktab object has no rows")
  }

  d.names <- rownames(x[[1]])
  if (is.null(d.names)){
    d.names <- as.character(seq_len(nlig))
  }
  
  ## Interactive metric selection
  method_params <- .get_method_parameters(type, scann)
  
  ## Determine if matrix weights are needed for NA handling
  napres <- any(type == "D") || any(is.na(unlist(x[seq_along(x$blo)])))
  
  ## Process ktab blocks
  ktab_results <- .process_ktab_blocks(x, type, option, method_params, tol, nlig, d.names, napres)
  
  ## Process precomputed distances
  precomp_results <- .process_precomputed_distances(precomputed.dists, precomputed.weights, nlig, d.names)
  
  ## Combine all distance objects and weights
  all_dist_objects <- c(ktab_results$dist_objects, precomp_results$dist_objects)
  all_weight_infos <- c(ktab_results$weight_infos, precomp_results$weight_infos)
  
  if (length(all_dist_objects) == 0) {
    warning("No distances to aggregate. Returning NULL.")
    return(NULL)
  }
  
  ## Aggregate distances using Gower's formula
  result <- .aggregate_distances(all_dist_objects, all_weight_infos, nlig, d.names, napres || precomp_results$uses_matrix_weights)
  return(result)
}


#' Scale quantitative data according to option
#' @noRd
.scale_quantitative_data <- function(df, option, tol) {
  switch(option,
    "scaledBYsd" = as.data.frame(scale(df)),
    "scaledBYrange" = {
      ranges <- apply(df, 2, function(col) diff(range(col, na.rm = TRUE)))
      centers <- apply(df, 2, min, na.rm = TRUE)
      scales <- ifelse(ranges < tol, 1, ranges)
      as.data.frame(scale(df, center = centers, scale = scales))
    },
    "noscale" = df
  )
}

#' Reorder distance matrix to match row names
#' @noRd
.reorder_distance_matrix <- function(pd, d.names, index) {
  
  pd_matrix <- as.matrix(pd)
  pd_labels <- rownames(pd_matrix)
  
  if (is.null(pd_labels) && is.null(d.names)) {
    return(pd)  # No reordering needed
  }
  
  if (is.null(d.names) && !is.null(pd_labels)) {
    warning("ktab has no row names, but precomputed distance ", index, 
            " does. Using distance matrix as-is.")
    return(pd)
  }
  
  if (!is.null(d.names) && !is.null(pd_labels)) {
    if (!all(pd_labels %in% d.names)) {
      warning("Row/col names of precomputed distance ", index, 
              " do not fully match ktab row names. Attempting to reorder.")
    }
    
    # Attempt reordering
    tryCatch({
      pd_reordered <- pd_matrix[d.names, d.names, drop = FALSE]
      return(stats::as.dist(pd_reordered))
    }, error = function(e) {
      stop("Could not reorder precomputed distance ", index, 
           " to match ktab row names. Error: ", e$message)
    })
  }
  
  return(pd)
}

#' Aggregate distance matrices using Gower's formula
#'
#' @param dist_objects List of dist objects to aggregate (from ktab blocks and precomputed distances)
#' @param weight_infos List of weights for each distance object (scalars or nlig×nlig matrices)
#' @param nlig Number of observations (rows) in the distance matrices
#' @param d.names Character vector of row/item labels for the final distance object
#' @param use_matrix_weights Logical; if TRUE uses matrix-based weights (for missing data), if FALSE uses scalar weights (simple weighted average)
#' @details d_ij = sqrt(Σ(w_k * d_ijk²) / Σ(w_k)). Function handles both scalar and matrix-based weighting for missing data scenarios.
#' @noRd
.aggregate_distances <- function(dist_objects, weight_infos, nlig, d.names, use_matrix_weights) {
  
  # Compute sum of squared distances
  sum_sq_distances <- matrix(0, nlig, nlig)
  for (d_obj in dist_objects) {
    sum_sq_distances <- sum_sq_distances + as.matrix(d_obj)^2
  }
  
  if (use_matrix_weights) {
    # Matrix-based weight aggregation
    total_weights <- matrix(0, nlig, nlig)
    
    for (i in seq_along(weight_infos)) {
      w_info <- weight_infos[[i]]
      d_matrix <- as.matrix(dist_objects[[i]])
      
      if (is.matrix(w_info)) {
        total_weights <- total_weights + w_info
      } else {
        # Convert scalar to matrix, setting 0 where distance is NA
        scalar_matrix <- matrix(w_info, nlig, nlig)
        scalar_matrix[is.na(d_matrix)] <- 0
        diag(scalar_matrix) <- 0
        total_weights <- total_weights + scalar_matrix
      }
    }
    
    # Add diagonal offset (ade4::dist.ktab behavior)
    diag(total_weights) <- diag(total_weights) + 1
    
    # Compute final distances
    final_distances <- matrix(NA, nlig, nlig)
    valid_weights <- total_weights > 0
    final_distances[valid_weights] <- sqrt(sum_sq_distances[valid_weights] / 
                                           total_weights[valid_weights])
    
  } else {
    # Scalar weight aggregation
    total_weight <- sum(unlist(weight_infos))
    
    if (total_weight <= 0) {
      warning("Total weight is non-positive (", total_weight, 
              "). Distance computation may be invalid.")
    }
    
    final_distances <- sqrt(sum_sq_distances / total_weight)
  }
  
  # Ensure symmetric matrix with zero diagonal
  diag(final_distances) <- 0
  
  # Convert to dist object
  res <- stats::as.dist(final_distances)
  attr(res, "Labels") <- d.names
  
  return(res)
}


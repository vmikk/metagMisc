#' Aggregate distance matrices using Gower's formula
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
      final_distances <- sqrt(sum_sq_distances / total_weight)
    } else {
      final_distances <- sqrt(sum_sq_distances / total_weight)
    }
  }
  
  # Ensure symmetric matrix with zero diagonal
  diag(final_distances) <- 0
  
  # Convert to dist object
  res <- stats::as.dist(final_distances)
  attr(res, "Labels") <- d.names
  
  return(res)
}


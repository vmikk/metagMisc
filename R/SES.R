
#' @title Estimate standardized effect size (SES) and p-value
#'
#' @param obs Single value of the observed statistics
#' @param sim Numeric vector of simulated values (null distribution)
#' @param alternative Character string specifying the alternative hypothesis; "greater", "less", or "two-sided" (default)
#'
#' @return Data frame with the effect size estimate (SES) and p-value (Pval).
#' @export
#'
#' @details
#' Standardized effect size (SES) represents Z-score, and is calculated as
#' (observed - mean_simulated) / standard_deviation_simulated.
#' SES values more extreme than +/-1.96 belong to the outer 5% of the null distribution
#' and thus are significant at alpha=0.05 with two-sided hypothesis test.
#'
#' @examples
#' SES(2, 1:10)
#'
SES <- function(obs, sim, alternative = "two-sided"){

  ## Check fo NA values in sim
  if(any(is.na(sim))){
    warning("Warning: NA values in `sim` were removed.\n")
    sim <- na.omit(sim)
  }

  ## Estimate SES and simulation summary
  res <- data.frame(
    Obs = obs,
    SES = (obs - mean(sim)) / sd(sim),
    Sim.mean = mean(sim),
    Sim.Q1 = stats::quantile(sim, probs = 0.25),
    Sim.Q3 = stats::quantile(sim, probs = 0.75),
    Sim.Var = var(sim)
    )

  ## Estimate P-value (based on ade4::as.randtest)
  if(alternative == "greater"){
    res$Pvalue <- (sum(sim >= obs) + 1)/(length(sim) + 1)
  }
  else if(alternative == "less"){
    res$Pvalue <- (sum(sim <= obs) + 1)/(length(sim) + 1)
  }
  else if(alternative == "two-sided") {
    sim0 <- abs(sim - mean(sim))
    obs0 <- abs(obs - mean(sim))
    res$Pvalue <- (sum(sim0 >= obs0) + 1)/(length(sim) +1)
  }

  rownames(res) <- NULL
  return(res)
}


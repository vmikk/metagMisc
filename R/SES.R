
## Estimate standardized effect size (SES) and p-value
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


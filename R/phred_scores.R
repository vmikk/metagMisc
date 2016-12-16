

#' @title Estimate Phred Quality Score score from probability of an error
#'
#' @param p Probability that the base is called wrong.
#'
#' @return Phred Q-Score
#' @export
#' @seealso \code{\link{phred_inverse}}
#'
#' @examples phred_score(c(0.001, 0.01, 0.05, 0.1))
#'
phred_score <- function(p){
  Q.score <- -1 * 10 * log(p, base = 10)
  return(Q.score)
}

#' @title Estimate probability of an error based on the Phred-score.
#'
#' @param q = Phred Quality Score
#'
#' @return Probability that the base is called wrong.
#' @export
#' @seealso \code{\link{phred_score}}
#'
#' @examples phred_inverse()
#' phred_inverse(30)
#' data.frame(Q = 2:42, ProbError = phred_inverse(2:42))
#'
phred_inverse <- function(q){
  Prob <-10^(-1 * q / 10)
  return(Prob)
}


# typical range from Q2 to Q40
# Q score of 3 means P=0.5, meaning that there is a 50% chance the base is wrong
# The lowest value usually found in practice is Q=2 (P=0.63), which means the base call is more likely to be wrong than correct.



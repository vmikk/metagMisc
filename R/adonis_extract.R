
## Extract ADONIS results
adonis_extract <- function(z){
  nf <- nrow(z) - 2
  rz <- data.frame(
    Factor = rownames(z)[1:nf],
    R2 = z$R2[1:nf],
    F = z$F[1:nf],
    df = paste(z$Df[1:nf], z$Df[nf+1],sep=";"),
    p = z$Pr[1:nf])
  return(rz)
}

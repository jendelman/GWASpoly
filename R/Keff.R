#' Effective number of markers
#' 
#' Effective number of markers
#' 
#' @param r2 matrix of squared correlations
#' @param alpha significance level
#' 
#' @references Moskvina V, Schmidt KM (2008) On multiple-testing correction in genome-wide association studies. Genetic Epidemiology 32:567-573. doi:10.1002/gepi.20331
#' 
#' @keywords internal
#' @return numeric

Keff <- function(r2,alpha) {
  m <- nrow(r2)
  Q <- sqrt(r2)
  Q[upper.tri(Q,diag=T)] <- NA
  rmax <- apply(Q[-1,],1,max,na.rm=T)
  kappa <- sqrt(1-rmax^(-1.31*log10(alpha)))
  return(1+sum(kappa))
}

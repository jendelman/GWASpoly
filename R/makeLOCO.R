makeLOCO <- function(K,exclude) {
  #K is list
  #exclude is vector of indices
  n.chr <- length(K)
  tmp <- K[[1]]*0
  keep <- setdiff(1:n.chr,exclude)
  for (i in keep) {
    tmp <- tmp + K[[i]]
  }
  return(tmp/length(keep))
}
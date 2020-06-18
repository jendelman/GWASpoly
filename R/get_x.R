get_x <- function(map) {
  #takes a map with chrom and position and returns x axis values for plotting multiple chromosomes
  a <- tapply(map[,2],map[,1],max)
  n <- length(a)
  m <- tapply(map[,2],map[,1],length)
  b <- c(0,apply(array(1:(n-1)),1,function(k){sum(a[1:k])}))
  x <- map[,2] + rep(b,times=m)
  return(x)
}

#' Plot LD vs distance
#' 
#' Plot LD vs distance
#' 
#' A monotone decreasing, convex spline is fit using R package \code{scam}. 
#' 
#' @param data variable inheriting from class \code{\link{GWASpoly}}
#' @param max.pair maximum number of r2 pairs for the spline
#' @param dof degrees of freedom for the spline
#' @param max.loci maximum number of markers to use per chromosome
#' @param position "bp" or "cM"
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @import scam
#' @importFrom stats dist

LD.plot <- function(data,max.pair=1e4,dof=8,max.loci=NULL,position="bp") {
  
  chroms <- levels(data@map$Chrom)
  n.chrom <- length(chroms)
  result <- NULL
  for (i in 1:n.chrom) {
    ix <- which(data@map$Chrom==chroms[i])
    m <- length(ix)
    if (!is.null(max.loci)) {
      if (m > max.loci) {
        ix <- sample(ix,max.loci)
        m <- max.loci
      }
    }
    tmp <- expand.grid(col=1:m,row=1:m)
    tmp <- tmp[tmp$row >= tmp$col,]  #only need lower triangular
    r2 <- cor(data@geno[,ix])^2
    r2.vec <- as.vector(r2[cbind(tmp$row,tmp$col)])
    d <- as.matrix(dist(matrix(data@map$Position[ix],ncol=1))) #distance matrix
    d.vec <- as.vector(d[cbind(tmp$row,tmp$col)])
    if (position=="bp") 
      d.vec <- d.vec/1e6
    
    result <- rbind(result,data.frame(d=d.vec,r2=r2.vec))
  }
  
  max.pair <- min(max.pair,nrow(result))
  scam.ans <- scam(formula=r2~s(d,bs=c("mdcx"),k=dof),data=result[sample(1:nrow(result),max.pair),])
  dmax <- max(result$d)
  predans <- predict.scam(scam.ans,newdata=data.frame(d=seq(0,dmax,length.out = 500)))
  spline.data <- data.frame(d=seq(0,dmax,length.out = 500),r2=predans)
  p <- ggplot(data=spline.data,aes(x=d,y=r2)) +  ylab(expression(r^2)) + theme_bw() + geom_line()
  if (position=="bp") {
    p <- p + xlab("Distance (Mb)") 
  } else {
    p <- p + xlab("Distance (cM)") 
  }
  return(p)
}

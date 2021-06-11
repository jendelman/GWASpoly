#' Set covariance matrix for polygenic effect
#' 
#' Set covariance matrix for polygenic effect
#' 
#' When \code{LOCO} = TRUE, K is computed for each chromosome as $K=MM'$, where M is the centered genotype matrix (lines x markers), and scaled to have unit diagonal (the overall scaling is not important for GWAS). When \code{LOCO} = FALSE, a single K matrix is computed for all markers (this is not recommended but provided for legacy reasons). Alternatively, the user can supply their own positive semidefinite K, with row.names that match the genotype identifiers (this option cannot be used with LOCO).
#' 
#' @param data Output from \code{read.GWASpoly}
#' @param K Optional: user-supplied matrix
#' @param n.core Number of cores for parallel computing
#' @param LOCO TRUE/FALSE, whether to use leave-one-chromosome-out 
#' 
#' @return Variable of class \code{GWASpoly.K}
#' 
#' @export
#' @importFrom parallel makeCluster clusterExport stopCluster parLapply
#' @importFrom methods as
#' 
set.K <- function(data,K=NULL,n.core=1,LOCO=TRUE) {
	stopifnot(inherits(data,"GWASpoly.data"))

	if (!is.null(K)) {
		gid.K <- rownames(K)
		gid <- rownames(data@geno)
		missing <- setdiff(gid,gid.K)
		if (length(missing) > 0) {
			cat("K matrix does not include following genotypes:\n",paste(missing,collapse="\n"))
			stop()
		}
		K <- list(all=K[gid,gid])
	} else {
	  
	  if (LOCO) {
	    markers <- split(data@map$Marker,data@map$Chrom)
	  } else {
	    markers <- list(all=data@map$Marker)
	  }
	  
	  f1 <- function(marks,data) {
	    M <- scale(data@geno[,marks],center=T,scale=F)
	    K <- tcrossprod(M)
	    K <- K/mean(diag(K))
	  }
	  
	  if (n.core > 1) {
	    cl <- makeCluster(n.core)
	    clusterExport(cl=cl,varlist=NULL)
	    K <- parLapply(cl,markers,f1,data=data)
	    stopCluster(cl)
	  } else {
	    K <- lapply(markers,f1,data=data)
	  }
	}
	return(new("GWASpoly.K",map=data@map,pheno=data@pheno,fixed=data@fixed,geno=data@geno,ploidy=data@ploidy,K=K))
}
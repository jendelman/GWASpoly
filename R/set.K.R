#' Set covariance matrix for polygenic effect
#' 
#' Set covariance matrix for polygenic effect
#' 
#' By default, K is computed as $K=MM'$, where M is the centered genotype matrix (lines x markers).  For GWAS, the overall scaling of K is irrelevant.  At present, K is scaled such that the mean of its diagonal elements is 1.  Alternatively, the user can supply any positive semidefinite K (with row.names that match the genotype identifiers).
#' 
#' @param data Output from \code{read.GWASpoly}
#' @param K Optional: user-supplied matrix
#' 
#' @return Variable of class \code{GWASpoly.K}
#' 
#' @export
set.K <- function(data,K=NULL) {
	stopifnot(inherits(data,"GWASpoly.data"))
	if (!is.null(K)) {
		gid.K <- rownames(K)
		gid <- rownames(data@geno)
		missing <- setdiff(gid,gid.K)
		if (length(missing) > 0) {
			cat("K matrix does not include following genotypes:\n",paste(missing,collapse="\n"))
			stop()
		}
		K <- K[gid,gid]
	} else {
		M <- scale(data@geno,center=T,scale=F)
		K <- tcrossprod(M)
		K <- K/mean(diag(K))
	}
	return(new("GWASpoly.K",map=data@map,pheno=data@pheno,fixed=data@fixed,geno=data@geno,ploidy=data@ploidy,K=K))
}
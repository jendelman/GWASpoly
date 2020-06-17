set.K <- function(data,K=NULL) {
	stopifnot(inherits(data,"GWASpoly"))
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
		#make K using realized relationship matrix
		M <- scale(data@geno,center=T,scale=F)
		K <- tcrossprod(M)
		K <- K/mean(diag(K))
	}
	return(new("GWASpoly.K",data,K=K))
}
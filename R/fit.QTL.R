#' Test markers as QTL under backward elimination
#' 
#' Test markers as QTL under backward elimination
#' 
#' Each element of \code{qtl} is a character vector of length two with  format c("marker","model"). Each element of \code{fixed} is a character vector of length two: the first element is the name of the effect (must match column in phenotype input file) and the second element is either "factor" or "numeric". The p-value and R2 for each maker are based on the likelihood ratio test under backward elimination, comparing the deviance to the chi-squared distribution. 
#' 
#' @param data variable inheriting from class \code{\link{GWASpoly.K}}
#' @param trait name of trait
#' @param qtl list of markers to fit in the multi-QTL model (see Details)
#' @param fixed list to specify fixed effects (see Details)
#' 
#' @return data frame with partial r2 and p-values
#' 
#' @export
#' @importFrom rrBLUP mixed.solve
#' @importFrom stats model.matrix pchisq
#' 
fit.QTL <- function(data,trait,qtl,fixed=NULL) {
	stopifnot(inherits(data,"GWASpoly.K"))
  stopifnot(is.element(trait,names(data@scores)))
  markers <- sapply(qtl,"[",1)
  models <- sapply(qtl,"[",2)
	stopifnot(models %in% c("additive","general",
	                        paste(1:(data@ploidy/2),"dom-ref",sep="-"),
	                        paste(1:(data@ploidy/2),"dom-alt",sep="-")))
	not.miss <- which(!is.na(data@pheno[,trait]))
	y <- data@pheno[not.miss,trait]
	pheno.gid <- data@pheno[not.miss,1]
	geno.gid <- rownames(data@geno)
	n.gid <- length(geno.gid)
	n <- length(y)
	Z <- matrix(0,n,n.gid)
	Z[cbind(1:n,match(pheno.gid,geno.gid))] <- 1
	X <- matrix(1,n,1)
	if (!is.null(fixed)) {
	  for (i in 1:length(fixed)) {
	    if (fixed[[i]][2]=="factor") {
	      xx <- factor(data@fixed[not.miss,fixed[[i]][1]])	
	      if (length(levels(xx)) > 1) {X <- cbind(X,model.matrix(~x,data.frame(x=xx))[,-1])}
	    } else {
	      X <- cbind(X,data@fixed[not.miss,fixed[[i]][1]])	
	    }
	  }
	}
	
	n.qtl <- length(qtl)
	S <- vector("list",n.qtl)
	df <- integer(n.qtl)
	X0 <- X
	for (i in 1:n.qtl) {
	  S[[i]] <- .design.score(data@geno[,markers[i]],model=models[i],ploidy=data@ploidy,min.MAF=0,max.geno.freq=1)
	  stopifnot(!is.null(S[[i]]))
	  df[i] <- ncol(S[[i]])
    X <- cbind(X,Z%*%S[[i]])
	}
	full.model <- mixed.solve(y=y,X=.make.full(X),Z=Z,K=data@K,method = "ML")

	pval <- R2 <- numeric(n.qtl)
	for (i in 1:n.qtl) {
	  X <- X0
	  if (n.qtl > 1) {
	    for (j in setdiff(1:n.qtl,i)) {
	      X <- cbind(X,Z%*%S[[j]])
	    }
	  }
	  reduced.model <- mixed.solve(y=y,X=.make.full(X),Z=Z,K=data@K,method = "ML")
	  deviance <- 2*(full.model$LL - reduced.model$LL)
	  pval[i] <- pchisq(q=deviance,df=df[i],lower.tail=FALSE)
	  R2[i] <- 1-exp(-deviance/n)
	}
	return(data.frame(marker=markers,model=models,R2=R2,pval=pval))
}

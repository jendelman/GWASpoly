#' Test markers as QTL under backward elimination
#' 
#' Test markers as QTL under backward elimination
#' 
#' \code{qtl} is a data frame with columns "Marker" and "Model", where each row corresponds to a QTL. \code{fixed} is a data frame with columns "Effect" and "Type": the first column is the name of the effect, which must match a column in the phenotype input file, and the second column is either "factor" or "numeric". The p-value and R2 for each marker are based on the likelihood ratio test under backward elimination, comparing the deviance to the chi-squared distribution. 
#' 
#' @param data variable inheriting from class \code{\link{GWASpoly.K}}
#' @param trait name of trait
#' @param qtl data frame to specify the multi-QTL model (see Details)
#' @param fixed data frame to specify the fixed effects (see Details)
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
	stopifnot(qtl$Model %in% c("additive","general",
	                        paste(1:(data@ploidy/2),"dom-ref",sep="-"),
	                        paste(1:(data@ploidy/2),"dom-alt",sep="-"),
	                        "diplo-general", "diplo-additive"))
	stopifnot(qtl$Marker %in% data@map$Marker)
	
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
	  for (i in 1:nrow(fixed)) {
	    if (fixed$Type[i]=="factor") {
	      xx <- factor(data@fixed[not.miss,fixed$Effect[i]])	
	      if (length(levels(xx)) > 1) {X <- cbind(X,model.matrix(~x,data.frame(x=xx))[,-1])}
	    } else {
	      X <- cbind(X,data@fixed[not.miss,fixed$Effect[i]])	
	    }
	  }
	}
	chrom <- as.character(data@map$Chrom[match(qtl$Marker,data@map$Marker)])
	if (length(data@K) > 1) {
	  if (length(setdiff(levels(data@map$Chrom),chrom))==0) {
	    stop("LOCO model cannot be used because there are QTL on every chromosome. Run set.K with LOCO=FALSE")
	  }
	  K <- makeLOCO(data@K,exclude=match(chrom,levels(data@map$Chrom)))
	} else {
	  K <- data@K[[1]]
	}
	
	n.qtl <- nrow(qtl)
	S <- vector("list",n.qtl)
	df <- integer(n.qtl)
	X0 <- X
	for (i in 1:n.qtl) {
	  S[[i]] <- .design.score(data@geno[,qtl$Marker[i]],model=qtl$Model[i],ploidy=data@ploidy,min.MAF=0,max.geno.freq=1)
	  stopifnot(!is.null(S[[i]]))
	  df[i] <- ncol(S[[i]])
    X <- cbind(X,Z%*%S[[i]])
	}
	
	full.model <- mixed.solve(y=y,X=.make.full(X),Z=Z,K=K,method = "ML")

	pval <- R2 <- numeric(n.qtl)
	for (i in 1:n.qtl) {
	  X <- X0
	  if (n.qtl > 1) {
	    for (j in setdiff(1:n.qtl,i)) {
	      X <- cbind(X,Z%*%S[[j]])
	    }
	  }
	  reduced.model <- mixed.solve(y=y,X=.make.full(X),Z=Z,K=K,method = "ML")
	  deviance <- 2*(full.model$LL - reduced.model$LL)
	  pval[i] <- pchisq(q=deviance,df=df[i],lower.tail=FALSE)
	  R2[i] <- 1-exp(-deviance/n)
	}
	return(data.frame(data@map[match(qtl$Marker,data@map$Marker),1:3],Model=qtl$Model,R2=R2,pval=pval))
}

#' Compute marker significance scores
#' 
#' Compute marker significance scores
#' 
#' The following marker-effect models are available:
#' \itemize{
#'  \item "additive": Indicates the marker effect is proportional to the dosage of the alternate allele
#'  \item "X-dom": where X can be any integer between 1 and ploidy/2 and refers to the allele dosage needed for complete dominance (e.g., "1-dom" = simplex dominance, "2-dom" = duplex dominance).  The software tries both dominance patterns for a given dosage model, e.g., whether the reference or alternate allele is dominant
#'  \item "diplo-general": All heterozygotes have the same effect
#'  \item "diplo-additive": All heterozygotes have the same effect, constrained to be halfway between the homozygous effects
#'  \item "general": There are no constraints on the effects of the different dosage levels
#'  }
#' To specify additional model parameters, such as the inclusion of fixed effects (Q matrix) and the minimum minor allele frequency, use \code{set.params}
#' 
#' @param data Output from \code{set.K}
#' @param models Vector of model names
#' @param traits Vector trait names (by default, all traits)
#' @param params Optional list of params created by \code{set.params}
#' @param n.core Number of cores for parallel computing
#' @param quiet TRUE/FALSE whether to suppress output charting progress
#'
#' @return Variable of class \code{GWASpoly.fitted}
#' @export
#' @importFrom rrBLUP mixed.solve
#' @importFrom stats model.matrix
#' 
GWASpoly <- function(data,models,traits=NULL,params=NULL,n.core=1,quiet=F) {
	
stopifnot(inherits(data,"GWASpoly.K"))
if (is.null(params)) {params <- set.params()}
if (!is.null(params$fixed)) {
	stopifnot(is.element(params$fixed,colnames(data@fixed)))
}
if (length(union(grep("ref",models,fixed=T),grep("alt",models,fixed=T)))) {
  stop("Do not include `dom` or `alt` in the model name.")
}
if (is.null(traits)) {traits <- colnames(data@pheno)[-1]}
stopifnot(is.element(traits,colnames(data@pheno)[-1]))

geno.gid <- rownames(data@geno)
m <- nrow(data@map)
data@K <- data@K[geno.gid,geno.gid]
n.gid <- length(geno.gid)

#params$models <- models
n.model <- length(models)
dom.models <- grep("dom",models,fixed=T)
if (length(dom.models)>0) {
	dom.orders <- as.integer(substr(models[dom.models],1,1))
	if (max(dom.orders) > data@ploidy/2) {
		stop("Maximum dominance model is ploidy/2")
	}
	dom.models2 <- sort(c(paste(models[dom.models],"ref",sep="-"),paste(models[dom.models],"alt",sep="-")))
} else {
	dom.models2 <- character(0)
}
other.models <- setdiff(1:n.model,dom.models)
if (!all(is.element(models[other.models],c("additive","general","diplo-general","diplo-additive")))) {
	stop("Invalid model")
}
models <- c(models[other.models],dom.models2)
n.model <- length(models)
n.trait <- length(traits)
if (params$n.PC > 0) {
	eig.vec <- eigen(data@K)$vectors
}
all.scores <- vector("list",n.trait)
names(all.scores) <- traits
all.effects <- all.scores

for (j in 1:n.trait) {
trait <- traits[j]
if (!quiet) {cat(paste("Analyzing trait:",trait,"\n"))}
not.miss <- which(!is.na(data@pheno[,trait]))
y <- data@pheno[not.miss,trait]
pheno.gid <- data@pheno[not.miss,1]
n <- length(y)
Z <- matrix(0,n,n.gid)
Z[cbind(1:n,match(pheno.gid,geno.gid))] <- 1
X <- matrix(1,n,1)
if (!is.null(params$fixed)) {
	for (i in 1:length(params$fixed)) {
		if (params$fixed.type[i]=="factor") {
			xx <- factor(data@fixed[not.miss,params$fixed[i]])	
			if (length(levels(xx)) > 1) {X <- cbind(X,model.matrix(~x,data.frame(x=xx))[,-1])}
		} else {
			X <- cbind(X,data@fixed[not.miss,params$fixed[i]])	
		}
	}
}
if (params$n.PC > 0) {
	X <- cbind(X,Z%*%eig.vec[,1:params$n.PC])
}
X2 <- .make.full(X) 

if (params$P3D) {
	if (!quiet) {cat("P3D approach: Estimating variance components...")}
	Hinv <- mixed.solve(y=y,X=X2,Z=Z,K=data@K,return.Hinv=TRUE)$Hinv
	if (!quiet) {cat("Completed \n")}
} else {
	Hinv <- NULL
}
scores <- matrix(NA,m,n.model)
colnames(scores) <- models
rownames(scores) <- colnames(data@geno)
betas <- scores
for (k in 1:n.model) {
	if (!quiet) {cat(paste("Testing markers for model:",models[k],"\n"))}
	if ((n.core > 1)& requireNamespace("parallel",quietly=TRUE)) {
    	it <- split(1:m,factor(cut(1:m,n.core,labels=FALSE)))
	    score.list <- parallel::mclapply(it,.score.calc,y,Z,X2,data@K,data@geno,Hinv,data@ploidy,models[k],params$min.MAF,params$max.geno.freq,mc.cores=n.core)
	   	scores[,k] <- unlist(lapply(score.list,function(el){el$score}))
	    betas[,k] <- unlist(lapply(score.list,function(el){el$beta}))
	} else {
    	ans <- .score.calc(1:m,y,Z,X2,data@K,data@geno,Hinv,data@ploidy,models[k],params$min.MAF,params$max.geno.freq)
    	scores[,k] <- ans$score
    	betas[,k] <- ans$beta
	}
}
all.scores[[j]] <- data.frame(scores,check.names=F)
all.effects[[j]] <- data.frame(betas,check.names=F)
}

return(new("GWASpoly.fitted",map=data@map,pheno=data@pheno,fixed=data@fixed,geno=data@geno,ploidy=data@ploidy,K=data@K,scores=all.scores,effects=all.effects,params=params))  
}

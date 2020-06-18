#' Set the significance threshold
#' 
#' Set the significance threshold
#' 
#' The FDR method is based on version 1.30.0 of the qvalue package
#' 
#' @param data Variable of class \code{GWASpoly.fitted}
#' @param method One of the following: "Bonferroni","FDR","permute"
#' @param level Genome-wide false positive rate for the Bonferroni or permutation methods; false discovery rate for method FDR
#' @param n.permute Number of permutations for method "permute"
#' @param n.core Number of cores to use for multicore processing
#' 
#' @return Variable of class \code{GWASpoly.thresh}
#' 
#' @export
#' 
set.threshold <- function(data,method,level=0.05,n.permute=1000,n.core=1) {

	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)
	n.trait <- length(traits)
	models <- colnames(data@scores[[1]])
	n.model <- length(models)
	methods <- c("Bonferroni","FDR","permute")
	stopifnot(is.element(method,methods))
	threshold <- matrix(NA,n.trait,n.model)
	colnames(threshold) <- models
	rownames(threshold) <- traits
	for (i in 1:n.trait) {
		trait <- traits[i]
		if (method=="permute") {
			print(paste("Trait:",trait),quote=F)
			y <- data@pheno[,trait]
			ix <- which(!is.na(y))
			max.scores <- matrix(NA,n.permute,n.model)
			colnames(max.scores) <- models
			for (q in 1:n.permute) {
				print(paste("Permutation",q),quote=F)
				data2 <- data
				data2@pheno[ix,trait] <- sample(y[ix])
				data2 <- GWASpoly(data2,models=data@params$models,traits=trait,params=data@params,quiet=T,n.core=n.core)
				for (j in 1:n.model) {max.scores[q,j] <- max(data2@scores[[trait]][,models[j]],na.rm=T)}				
			}
		}
		for (j in 1:n.model) {
			model <-  models[j]
			scores <- as.vector(na.omit(data@scores[[trait]][,model]))
			m <- length(scores)
			if (method=="Bonferroni") {threshold[i,j] <- -log10(level/m)}
			if (method=="FDR") {
				tmp <- cbind(10^(-scores),.qvalue(10^(-scores)))
				tmp <- tmp[order(tmp[,2]),]
				if (tmp[1,2] > level) {
					threshold[i,j] <- -log10(tmp[1,1])*1.2
				} else {
					k <- max(which(tmp[,2] < level))
					threshold[i,j] <- -log10(mean(tmp[k:(k+1),1]))
				}
			}
			if (method=="permute") {
				threshold[i,j] <- sort(max.scores[,model],decreasing=TRUE)[floor(level*n.permute)]
			}	
		}
	}
	return(new("GWASpoly.thresh",map=data@map,pheno=data@pheno,fixed=data@fixed,geno=data@geno,ploidy=data@ploidy,K=data@K,scores=data@scores,effects=data@effects,params=data@params,threshold=threshold))
}
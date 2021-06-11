#' Set the significance threshold
#' 
#' Set the significance threshold
#' 
#' The default method, "M.eff", is a Bonferroni-type correction but using an effective number of markers that accounts for LD between markers (Moskvina and Schmidt, 2008). The FDR method is based on version 1.30.0 of the qvalue package. 
#' 
#' @param data Variable of class \code{GWASpoly.fitted}
#' @param method One of the following: "M.eff","Bonferroni","FDR","permute"
#' @param level Genome-wide false positive or false discovery rate (depending on \code{method}). 
#' @param n.permute Number of permutations for method "permute"
#' @param n.core Number of cores to use for multicore processing
#' 
#' @references Moskvina V, Schmidt KM (2008) On multiple-testing correction in genome-wide association studies. Genetic Epidemiology 32:567-573. doi:10.1002/gepi.20331
#' 
#' @return Variable of class \code{GWASpoly.thresh}
#' 
#' @export
#' @importFrom stats cor
#' @import Matrix
#' @importFrom methods as
#' 
set.threshold <- function(data,method="M.eff",level=0.05,n.permute=1000,n.core=1) {

	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)
	n.trait <- length(traits)
	models <- colnames(data@scores[[1]])
	n.model <- length(models)
	methods <- c("M.eff","Bonferroni","FDR","permute")
	stopifnot(is.element(method,methods))
	threshold <- matrix(NA,n.trait,n.model)
	colnames(threshold) <- models
	rownames(threshold) <- traits
	
	if (method=="M.eff") {
	  chrom <- levels(data@map[,2])
	  n.chrom <- length(chrom)
	  r2 <- vector("list",n.chrom)
	  names(r2) <- chrom
	  for (i in chrom) {
	    ix <- which(data@map[,2]==i)
	    r2[[i]] <- as(cor(data@geno[,ix])^2,"dspMatrix")
	  }
	}
	
	for (i in 1:n.trait) {
		trait <- traits[i]
		if (method=="permute") {
			print(paste("Trait:",trait),quote=F)
			#y <- data@pheno[,trait]
			#ix <- which(!is.na(y))
			max.scores <- matrix(NA,n.permute,n.model)
			colnames(max.scores) <- models
			for (q in 1:n.permute) {
				print(paste("Permutation",q),quote=F)
				data2 <- data
				data2@pheno[,1] <- sample(data@pheno[,1]) #permute id
				data2 <- GWASpoly(data2,models=models,traits=trait,params=data@params,quiet=T,n.core=n.core)
				for (j in 1:n.model) {max.scores[q,j] <- max(data2@scores[[trait]][,models[j]],na.rm=T)}				
			}
		}
		for (j in 1:n.model) {
			model <-  models[j]
			iv <- which(!is.na(data@scores[[trait]][,model]))
			scores <- as.vector(data@scores[[trait]][iv,model])
			m <- length(scores)
			if (method=="Bonferroni") {threshold[i,j] <- -log10(level/m)}
			if (method=="M.eff") {
			  me <- 0
			  for (chr in chrom) {
			    ix <- data@map[intersect(iv,which(data@map[,2]==chr)),1]
			    if (length(ix)>1) {
			      me <- me + Keff(r2=r2[[chr]][ix,ix],alpha=level)
			    } else {
			      me <- me + 1
			    }
			  }
			  threshold[i,j] <- -log10(level/me)
			}
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
				threshold[i,j] <- sort(max.scores[,model],decreasing=TRUE)[max(floor(level*n.permute),1)]
			}	
		}
	}
	cat("Thresholds\n")
	print(round(threshold,2))
	return(new("GWASpoly.thresh",map=data@map,pheno=data@pheno,fixed=data@fixed,geno=data@geno,ploidy=data@ploidy,K=data@K,scores=data@scores,effects=data@effects,params=data@params,threshold=threshold))
}
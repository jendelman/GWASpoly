#' Extract significant QTL
#' 
#' Output a table with significant markers
#' 
#' To return all significant markers (original behavior of the function), use \code{bp.window=NULL}. Assumes input map position in bp. 
#' 
#' @param data Output from \code{set.threshold}
#' @param traits Vector of trait names (by default, all traits)
#' @param models Vector of model names (by default, all models)
#' @param bp.window prune output to return only the most significant marker within this window size
#' 
#' @return Data frame with results. Score = -log10(p). Effect = marker effect (not available for the general and diplo-general models because there are multiple effects).
#' 
#' @export
#' @importFrom stats dist
#' 
get.QTL <- function(data,traits=NULL,models=NULL,bp.window=1e6) {
	stopifnot(inherits(data,"GWASpoly.thresh"))
	if (is.null(traits)) {
		traits <- names(data@scores)
	} else {
		stopifnot(is.element(traits,names(data@scores)))
	}
  all.models <- colnames(data@scores[[1]])
  if (is.null(models)) {
    models <- all.models
  } else {
    models <- unlist(lapply(as.list(models),function(x){all.models[grep(x,all.models,fixed=T)]}))
    stopifnot(all(is.element(models,all.models)))
  }
	n.model <- length(models)
	n.trait <- length(traits)
	output <- data.frame(NULL)
	for (i in 1:n.trait) {
		for (j in 1:n.model) {
			ix <- which(data@scores[[traits[i]]][,models[j]] > data@threshold[traits[i],models[j]])
			n.ix <- length(ix)
			if (n.ix > 0) {
			  tmp <- data.frame(Trait=rep(traits[i],n.ix),Model=rep(models[j],n.ix),
			                    Threshold=round(rep(data@threshold[traits[i],models[j]],n.ix),2),
			                    data@map[ix,],Score=round(data@scores[[traits[i]]][ix,models[j]],2),
			                    Effect=data@effects[[traits[i]]][ix,models[j]],stringsAsFactors=F,check.names=F)
			  if (!is.null(bp.window)) {
			    chrom <- unique(tmp$Chrom)
			    for (k in chrom) {
			      ix <- which(tmp$Chrom==k)
			      tmp2 <- tmp[ix,]
			      tmp2 <- tmp2[order(tmp2$Score,decreasing=T),]
			      d <- as.matrix(dist(tmp2$Position))
			      d[upper.tri(d,diag=F)] <- NA
			      diag(d) <- bp.window + 1
			      min.d <- apply(d,1,min,na.rm=T)
			      tmp3 <- tmp2[which(min.d > bp.window),]
			      output <- rbind(output,tmp3[order(tmp3$Position),])
			    }
			  } else {
			    output <- rbind(output,tmp)
			  }
			}
		}
	}
	return(output)
}
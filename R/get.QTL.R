#' Extract significant QTL
#' 
#' Output a table with all significant markers
#' 
#' Score = -log10(p). Effect = marker effect (not available for the general and diplo-general models).
#' 
#' @param data Output from \code{set.threshold}
#' @param traits Vector of trait names (by default, all traits)
#' @param models Vector of model names (by default, all models)
#' 
#' @return Data frame with results for significant markers
#' 
#' @export
#' 
get.QTL <- function(data,traits=NULL,models=NULL) {
	stopifnot(inherits(data,"GWASpoly.thresh"))
	if (is.null(traits)) {
		traits <- names(data@scores)
	} else {
		stopifnot(is.element(traits,names(data@scores)))
	}
	if (is.null(models)) {
		models <- colnames(data@scores[[1]])
	} else {
		stopifnot(is.element(models,colnames(data@scores[[1]])))
	}
	n.model <- length(models)
	n.trait <- length(traits)
	output <- data.frame(NULL)
	for (i in 1:n.trait) {
		for (j in 1:n.model) {
			ix <- which(data@scores[[traits[i]]][,models[j]] > data@threshold[traits[i],models[j]])
			n.ix <- length(ix)
			if (n.ix > 0) {
			  output <- rbind(output,data.frame(Trait=rep(traits[i],n.ix),Model=rep(models[j],n.ix),Threshold=round(rep(data@threshold[traits[i],models[j]],n.ix),2),data@map[ix,],Score=round(data@scores[[traits[i]]][ix,models[j]],2),Effect=formatC(data@effects[[traits[i]]][ix,models[j]],2,format="e"),stringsAsFactors=F,check.names=F))
			}
		}
	}
	return(output)
}
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
			output <- rbind(output,data.frame(Trait=rep(traits[i],n.ix),Model=rep(models[j],n.ix),Threshold=round(rep(data@threshold[traits[i],models[j]],n.ix),2),data@map[ix,],Score=round(data@scores[[traits[i]]][ix,models[j]],2),Effect=round(data@effects[[traits[i]]][ix,models[j]],2),stringsAsFactors=F,check.names=F))
		}
	}
	return(output)
}
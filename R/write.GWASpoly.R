#' Write results to file
#' 
#' Write results to file
#' 
#' Score = -log10(p). Effect = marker effect (not available for the general and diplo-general models).
#' 
#' @param data Variable of class \code{GWASpoly.fitted}
#' @param trait Trait name
#' @param filename Filename
#' @param what Either "scores" or "effects"
#' @param delim Delimiter to use in the output file (default is comma)
#' 
#' @export
#' @importFrom utils write.table
#' 
write.GWASpoly <- function(data,trait,filename,what="scores",delim=",") {
	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)
	stopifnot(is.element(trait,traits))
	stopifnot(is.element(what,c("scores","effects")))
	if (what=="scores") {
		write.table(cbind(data@map,data@scores[[trait]]),file=filename,quote=F,sep=delim,row.names=F)
	} else {
		write.table(cbind(data@map,data@effects[[trait]]),file=filename,quote=F,sep=delim,row.names=F)
	}
	return(NULL)
}
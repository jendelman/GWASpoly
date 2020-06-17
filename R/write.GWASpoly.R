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
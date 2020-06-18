#' Quantile-Quantile (QQ) Plot
#' 
#' Inspect p-value inflation using a QQ plot
#' 
#' @param data Variable of class \code{GWASpoly.fitted}
#' @param traits Vector of trait names (by default, all traits plotted)
#' @param models Vector of model names (by default, all models plotted)
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom stats ppoints


qq.plot <- function(data,traits=NULL,models=NULL) {
	stopifnot(inherits(data,"GWASpoly.fitted"))
  if (is.null(traits)) {
    traits <- names(data@scores)
  } else {
	  stopifnot(all(is.element(traits,names(data@scores))))
  }
  n.trait <- length(traits)
  all.models <- colnames(data@scores[[1]])
  if (is.null(models)) {
    models <- all.models
  } else {
    stopifnot(all(is.element(models,all.models)))
  }

  plotme <- NULL
	for (k in 1:n.trait) {
	  scores <- as.data.frame(data@scores[[traits[k]]][,models])
	  colnames(scores) <- models
	  tmp <- pivot_longer(data=scores,colnames(scores),names_to="model",values_to="y",values_drop_na=TRUE)
	  tmp$model <- factor(tmp$model,levels=models,ordered=T)
	  tmp <- tmp[order(tmp$model,tmp$y,decreasing = c(FALSE,TRUE)),]
	  tmp$x <- unlist(
	    tapply(tmp$y,factor(tmp$model,levels=models,ordered=T),function(x) {
	     n <- length(x)
	     unif.p <- -log10(ppoints(n))
	     }))
	  tmp$trait <- traits[k]
	  plotme <- rbind(plotme,tmp)
	}
  plotme$trait <- factor(plotme$trait,ordered=T,levels=traits)
	 
  p <- ggplot(data=plotme,aes(x=x,y=y,colour=model)) + facet_wrap(~trait) + geom_point() + theme_bw() + xlab(expression(paste("Expected -log"[10],"(p)",sep=""))) + ylab(expression(paste("Observed -log"[10],"(p)",sep=""))) + scale_colour_brewer(palette="Set1") + geom_abline(slope=1,intercept=0,linetype=2) + theme(text = element_text(size=20))
  return(p)	  
}
  
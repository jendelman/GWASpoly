#' Quantile-Quantile (QQ) Plot
#' 
#' Inspect p-value inflation using a QQ plot
#' 
#' @param data Variable of class \code{GWASpoly.fitted}
#' @param trait Trait name
#' @param models Vector of model names (by default, all models plotted)
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom stats ppoints
#' @importFrom rlang .data

qq.plot <- function(data,trait,models=NULL) {
	stopifnot(inherits(data,"GWASpoly.fitted"))
  stopifnot(is.element(trait,names(data@scores)))

  all.models <- colnames(data@scores[[trait]])
  if (is.null(models)) {
    models <- all.models
  } else {
    models <- unlist(lapply(as.list(models),function(x){all.models[grep(x,all.models,fixed=T)]}))
    stopifnot(all(is.element(models,all.models)))
  }

	scores <- as.data.frame(data@scores[[trait]][,models])
	colnames(scores) <- models
  scores$Chrom <- data@map$Chrom
	tmp <- pivot_longer(data=scores,cols=1:length(models),names_to="model",values_to="y",values_drop_na=TRUE)
  tmp <- as.data.frame(tmp)
  tmp$model <- factor(tmp$model,levels=models,ordered=T)
  tmp <- tmp[order(tmp$model,tmp$Chrom,tmp$y,decreasing = c(FALSE,FALSE,TRUE)),]
  tmp2 <- tapply(tmp$y,list(tmp$Chrom,tmp$model),function(x) {
    n <- length(x)
    unif.p <- -log10(ppoints(n))
	  })
  tmp$x <- unlist(tmp2)
	 
  p <- ggplot(data=tmp,aes(x=.data$x,y=.data$y,colour=.data$model)) + facet_wrap(~Chrom) + geom_point() + theme_bw() + xlab(expression(paste("Expected -log"[10],"(p)",sep=""))) + ylab(expression(paste("Observed -log"[10],"(p)",sep=""))) + scale_colour_brewer(palette="Set1") + geom_abline(slope=1,intercept=0,linetype=2) + theme(text = element_text(size=15))
  return(p)	  
}
  
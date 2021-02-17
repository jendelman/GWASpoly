#' Create Manhattan plot
#' 
#' Create Manhattan plot
#' 
#' Results for the ref and alt versions of the dominance model are combined. If \code{data} is the output from \code{\link{set.threshold}}, then the threshold is displayed as a horizontal dashed line when \code{models} contains a single model. Because the threshold varies between models, it is not drawn when multiple models are included. Although the ref and alt versions of each dominance model are slightly different (as seen with \code{\link{qq.plot}}), they are treated as a single model for the Manhattan plot, and the average threshold is shown.
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

manhattan.plot <- function(data,traits=NULL,models=NULL) {
		
	stopifnot(inherits(data,"GWASpoly.fitted"))
  if (is.null(traits)) {
    traits <- names(data@scores)
  } else {
    stopifnot(all(is.element(traits,names(data@scores))))
  }
  n.trait <- length(traits)
  n.model <- length(models)
  all.models <- colnames(data@scores[[1]])
  if (is.null(models)) {
    models <- all.models
  } else {
    models <- unlist(lapply(as.list(models),function(x){all.models[grep(x,all.models,fixed=T)]}))
    stopifnot(all(is.element(models,all.models)))
  }

  plotme <- thresh.data <- NULL
  x <- get_x(data@map[,2:3])
  for (k in 1:n.trait) {
    scores <- as.data.frame(data@scores[[traits[k]]][,models])
    colnames(scores) <- models
    scores$x <- x
    scores$color <- factor(ifelse(as.integer(factor(data@map[,2]))%%2==1,1,0))
    tmp <- pivot_longer(data=scores,cols=match(models,colnames(scores)),names_to="model",values_to="y",values_drop_na=TRUE)
    tmp$trait <- traits[k]
    if (n.model==1 & inherits(data,"GWASpoly.thresh")) {
      thresh.data <- rbind(thresh.data,data.frame(y=mean(data@threshold[traits[k],models]),trait=traits[k]))
    }
    plotme <- rbind(plotme,tmp)
  }
  plotme$trait <- factor(plotme$trait)
  plotme$model <- gsub(pattern="-ref",replacement="",x=plotme$model)
  plotme$model <- gsub(pattern="-alt",replacement="",x=plotme$model)
  plotme$model <- factor(plotme$model)
  
  allchr <- unique(data@map[,2])
  breaks <- (tapply(x,data@map[,2],max) + tapply(x,data@map[,2],min))/2
  p <- ggplot(data=plotme,aes(x=x,y=y,colour=color,shape=model)) +
    ylab(expression(paste("-log"[10],"(p)"))) +
    theme_bw() + theme(text = element_text(size=15),panel.grid = element_blank()) +
    scale_x_continuous(name="Chromosome",breaks=breaks,labels=allchr) +
    guides(colour="none") + geom_point() + scale_shape(solid=FALSE) + scale_colour_manual(values = c("cornflowerblue","darkblue")) + facet_wrap(~trait)
  if (n.model==1 & inherits(data,"GWASpoly.thresh")) {
    p <- p + geom_hline(data=thresh.data,mapping=aes(yintercept=y),linetype=2,colour="grey50") + guides(shape=FALSE)
  }
	return(p)	
}

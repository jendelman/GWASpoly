manhattan.plot <- function(data,trait,model,cex=1,y.max=NULL,filename=NULL) {
		
	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)
	stopifnot(is.element(trait,traits))
	models <- colnames(data@scores[[trait]])
	stopifnot(is.element(model,models))
	input <- data.frame(data@map[,1:3],scores=data@scores[[trait]][,model])

	chroms <- levels(data@map$Chrom)	
	n.chrom <- length(chroms)
	chrom.start <- rep(0,n.chrom)
	chrom.mid <- rep(0,n.chrom)
	if (is.null(y.max)) {
		y.max <- max(input[,4],na.rm=T)+1
		if (class(data)=="GWASpoly.thresh") {
			y.max <- max(y.max,data@threshold[trait,model]+1)
		}
	}

	if (!is.null(filename)) {postscript(file=filename,horizontal=FALSE)}
	par(pty="s")
	if (n.chrom > 1) {
		for (i in 1:(n.chrom-1)) {chrom.start[i+1] <- chrom.start[i]+max(input[which(input[,2]==chroms[i]),3])+1}
	}
	x.max <- chrom.start[n.chrom]+max(input[which(input[,2]==chroms[n.chrom]),3])
	plot(0,0,type="n",xlim=c(0,x.max),ylim=c(0,y.max),ylab=expression(paste("-log"[10],"(p)",sep="")),xlab="Chromosome",xaxt="n",main=paste(trait," (",model,") ",sep=""))

	for (i in seq(1,n.chrom,by=2)) {
		ix <- which(input[,2]==chroms[i])
		chrom.mid[i] <- chrom.start[i]+max(input[ix,3])/2
		points(chrom.start[i]+input[ix,3],input[ix,4],col="dark blue",pch=16,cex=cex)
	}	

	if (n.chrom > 1){
	for (i in seq(2,n.chrom,by=2)) {
		ix <- which(input[,2]==chroms[i])
		chrom.mid[i] <- chrom.start[i]+max(input[ix,3])/2
		points(chrom.start[i]+input[ix,3],input[ix,4],col="cornflowerblue",pch=16,cex=cex)
	}	
	}
	
	if (class(data)=="GWASpoly.thresh") {
		abline(h=data@threshold[trait,model],lty=2)
	}
	axis(side=1,at=chrom.mid,labels=chroms)
	if (!is.null(filename)) {dev.off()}
	return(NULL)	
}
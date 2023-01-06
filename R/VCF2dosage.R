#' Convert VCF to dosage file
#' 
#' Convert VCF to dosage file
#' 
#' Only bi-allelic variants supported. The "GT" option for \code{geno.code} is the posterior maximum genotype (e.g., 0/0/1/1). "DS" represents the posterior mean dosage of the alternate allele. VCF file must conform to 4.1 or later.  
#'
#' @param VCF.file VCF filename (can be gzipped)
#' @param dosage.file CSV filename to output with allele dosage
#' @param geno.code genotype code in the FORMAT field: "GT" of "DS"
#' @param ploidy ploidy
#' @param samples optional vector of sample names, to export subset of the population
#' @param min.DP minimum per sample depth (DP) to export genotype. Default is 1, for no filtering.
#' @param max.missing threshold for missing data per marker, as a proportion.
#' @param min.minor minimum number of samples with the minor allele. Default is 5.
#' 
#' @export

VCF2dosage <- function(VCF.file, dosage.file, geno.code, ploidy, samples=NULL,
                       min.DP=1, max.missing, min.minor=5) {

  stopifnot(geno.code %in% c("GT","DS"))
  stopifnot(min.minor > 0)
  
  con <- file(VCF.file,"r")
  temp <- readLines(con,1)
  comment.line <- 0
  while(substr(temp,1,2)=="##") {
  	temp <- readLines(con,1)
  	comment.line <- comment.line+1
  }
  header <- strsplit(temp,split="\t",fixed=TRUE)[[1]]
  sample.names <- header[-(1:9)]
  if (is.null(samples)) {samples <- sample.names}
  n.sample <- length(samples)
  sample.ix <- match(samples,sample.names)
  iw <- which(is.na(sample.ix))
  if (length(iw)>0) {
    print("The following samples are not in the VCF file:")
    for (k in 1:length(iw)) {print(samples[iw])}
    close(con)
    stop()
  }
  out1 <- file(dosage.file,"w")
  delim=","
  writeLines(con=out1,text=paste(c("Marker","Chrom","Position","REF","ALT",samples),collapse=delim))
  
  FORMAT.pos <- match("FORMAT",header,nomatch = 0)
  if (FORMAT.pos==0) {
    close(con)
    close(out1)
    stop("FORMAT information not present in the VCF file")
  }
  temp <- readLines(con,1)
  temp2 <- strsplit(temp,split="\t",fixed=TRUE)[[1]]
  FORMAT <- strsplit(temp2[FORMAT.pos],split=":",fixed=TRUE)[[1]]  #FORMAT field
  geno.pos <- match(geno.code,FORMAT,nomatch = 0)
  if (geno.pos==0) {
    close(con)
    close(out1)
    stop("geno code not present in FORMAT")
  }
  DP.pos <- match("DP",FORMAT,nomatch = 0)
  if (DP.pos==0 & min.DP > 1) {
    close(con)
    close(out1)
    stop("DP not present in FORMAT")
  }
  
  m1 <- m <- 0
  n.minor <- 0
  while(length(temp) > 0) {
    m1 <- m1 + 1
    if (m1%%1e4==0) {cat(sub("X",m1,"Number of variants processed...X\n"))}
    
  	temp2 <- strsplit(temp,split="\t",fixed=TRUE)[[1]]  
  	marker <- paste(temp2[1],temp2[2],sep="_")

  	REF <- temp2[4]
  	ALT <- temp2[5]
  	if (max(regexpr(",",c(REF,ALT),fixed=T)) < 0) {
		  #only include bi-allelic variants
		  
		  FORMAT <- strsplit(temp2[FORMAT.pos],split=":",fixed=TRUE)[[1]]  #FORMAT field
		  geno.pos <- match(geno.code,FORMAT,nomatch = 0)
		  if (geno.pos==0) {
		    close(con)
		    close(out1)
		    stop("geno code not present in FORMAT")
		  }
		  DP.pos <- match("DP",FORMAT,nomatch = 0)
		  if (DP.pos==0 & min.DP > 1) {
		    close(con)
		    close(out1)
		    stop("DP not present in FORMAT")
		  } 
  
		  data <- strsplit(temp2[-(1:9)],split=":",fixed=T)[sample.ix]  #holds string of data for each sample, as a list
		  
		  convert <- function(w,geno.code,geno.pos) {
		    if(length(w) < geno.pos || w[geno.pos]==".")
		      return(as.numeric(NA))  
		    if (geno.code=="GT") {
		      z <- as.integer(strsplit(w[geno.pos],split="/",fixed=T)[[1]])
		      return(sum(z))
		    }
		    if (geno.code=="DS"){
		      return(as.numeric(w[geno.pos]))
		    }
		  }
		  geno <- sapply(data,convert,geno.code=geno.code,geno.pos=geno.pos)
		  
		  if (min.DP > 1) {
		    DPvec <- sapply(data,function(w){
		      if(length(w) < DP.pos || w[DP.pos]==".") {
		        return(0)  
		      } else {
		        return(as.integer(w[DP.pos]))
		      }
		    })  
		    geno[DPvec < min.DP] <- NA
		  }
      
		  frac.miss <- sum(is.na(geno))/length(geno)
		  if (ploidy > 0 & frac.miss < 1) {
		    p <- mean(geno,na.rm=T)/ploidy
		    tab <- table(factor(round(geno),levels=0:ploidy))
		    if (p > 0.5) {
		      n.minor <- sum(tab) - tab[ploidy+1]
		    } else {
		      n.minor <- sum(tab) - tab[1]
		    }
		  }
		  
		  geno <- round(geno,1)
		  if (frac.miss <= max.missing & n.minor >= min.minor) {
  		  writeLines(con=out1,text=paste(c(marker,temp2[1],temp2[2],REF,ALT,geno),collapse=delim))
  		  m <- m + 1
		  }
  	}
  	temp <- readLines(con,1)
	}
  close(con)
  close(out1)
  print(paste("Number of markers in output:",m))
}

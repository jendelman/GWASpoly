set.params <- function(fixed=NULL,fixed.type=NULL,n.PC=0,MAF=0.05,geno.freq=0.95,P3D=T) {
stopifnot(MAF > 0)
stopifnot(MAF < 0.5)
stopifnot(geno.freq > 0)
stopifnot(geno.freq < 1)
stopifnot(length(fixed)==length(fixed.type))
stopifnot(is.element(fixed.type,c("numeric","factor")))
return(list(fixed=fixed,fixed.type=fixed.type,n.PC=n.PC,min.MAF=MAF,max.geno.freq=geno.freq,P3D=P3D))
}
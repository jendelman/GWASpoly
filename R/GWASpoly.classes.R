setClass("GWASpoly",slots=c(map="data.frame",pheno="data.frame",geno="matrix",fixed="data.frame",ploidy="numeric"))
setClass("GWASpoly.K",slots=c(K="matrix"),contains="GWASpoly")
setClass("GWASpoly.fitted",slots=c(scores="list",effects="list",params="list"),contains="GWASpoly.K")
setClass("GWASpoly.thresh",slots=c(threshold="matrix"),contains="GWASpoly.fitted")

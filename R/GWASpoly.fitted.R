#' S4 class with results from genome-wide scan
#' 
#' @slot map data frame with columns Marker,Chrom,Position,Ref,Alt
#' @slot pheno data frame of phenotypes
#' @slot geno matrix with allele dosages
#' @slot fixed data frame of fixed effects
#' @slot ploidy ploidy
#' @slot K list of covariance matrices 
#' @slot scores -log10(p) results 
#' @slot effects estimated marker effects
#' @slot params parameters used for the analysis
#' 
#' @include GWASpoly.K.R
GWASpoly.fitted <- setClass("GWASpoly.fitted",slots=c(scores="list",effects="list",params="list"),contains="GWASpoly.K")

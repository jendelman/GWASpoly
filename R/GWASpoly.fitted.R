#' S4 class with results from genome-wide scan
#' 
#' @slot map data frame with marker,chrom,and position (either bp or cM)
#' @slot pheno data frame of phenotypes
#' @slot geno matrix with allele dosages
#' @slot fixed data frame of fixed effects
#' @slot ploidy ploidy
#' @slot K covariance matrix for polygenic effect
#' @slot scores -log10(p) results 
#' @slot effects estimated marker effects
#' @slot params parameters used for the analysis
#' 
#' @include GWASpoly.K.R
GWASpoly.fitted <- setClass("GWASpoly.fitted",slots=c(scores="list",effects="list",params="list"),contains="GWASpoly.K")

#' S4 class with results from genome-wide scan and detection threshold
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
#' @slot threshold thresholds for significance
#' 
#' @include GWASpoly.fitted.R
GWASpoly.thresh <- setClass("GWASpoly.thresh",slots=c(threshold="matrix"),contains="GWASpoly.fitted")

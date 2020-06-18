#' S4 class with genotypes, phenotypes, and polygenic covariance
#' 
#' @slot map data frame with marker,chrom,and position (either bp or cM)
#' @slot pheno data frame of phenotypes
#' @slot geno matrix with allele dosages
#' @slot fixed data frame of fixed effects
#' @slot ploidy ploidy
#' @slot K covariance matrix for polygenic effect
#' 
#' @export
#' @include GWASpoly.data.R
GWASpoly.K <- setClass("GWASpoly.K",slots=c(K="matrix"),contains="GWASpoly.data")

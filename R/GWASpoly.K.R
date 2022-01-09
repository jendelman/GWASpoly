#' S4 class with genotypes, phenotypes, and polygenic covariance
#' 
#' @slot map data frame with columns Marker,Chrom,Position,Ref,Alt
#' @slot pheno data frame of phenotypes
#' @slot geno matrix with allele dosages
#' @slot fixed data frame of fixed effects
#' @slot ploidy ploidy
#' @slot K list of covariance matrices (one for each chromosome)
#' 
#' @include GWASpoly.data.R
GWASpoly.K <- setClass("GWASpoly.K",slots=c(K="list"),contains="GWASpoly.data")

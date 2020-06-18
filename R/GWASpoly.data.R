#' S4 class with genotype and phenotype data
#' 
#' @slot map data frame with marker,chrom,and position (either bp or cM)
#' @slot pheno data frame of phenotypes
#' @slot geno matrix with allele dosages
#' @slot fixed data frame of fixed effects
#' @slot ploidy ploidy
#' 
#' @export
GWASpoly.data <- setClass("GWASpoly.data",slots=c(map="data.frame",pheno="data.frame",geno="matrix",fixed="data.frame",ploidy="numeric"))
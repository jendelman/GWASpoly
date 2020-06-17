\name{read.GWASpoly}
\alias{read.GWASpoly}
\title{
Read in data for GWASpoly
}
\description{
Reads marker and phenotype data
}
\usage{
read.GWASpoly(ploidy, pheno.file, geno.file, format, n.traits, delim = ",")
}
\arguments{
  \item{ploidy}{
Integer specifying the ploidy level (e.g., 2 for diploid, 4 for tetraploid)
}
  \item{pheno.file}{
Name of the phenotype file.
}
  \item{geno.file}{
Name of the genotype file.
}
  \item{format}{
Format for the marker data, either "numeric", "AB", or "ACGT" (see details).
}
  \item{n.traits}{
Integer specifying the number of traits in the phenotype file.  
}
  \item{delim}{
Character to indicate the delimiter in the data files (e.g., "," for csv, "\\t" for tab-delimited).
}
}
\details{
The first column of the phenotype file contains the genotype identifier, columns 2 through (n.traits + 1) contain trait values, and subsequent columns contain the levels (for factors) or numeric values (for covariates) of any fixed effects.  

The first three columns of the genotype file are (1) marker name, (2) chromosome, and (3) position.  Subsequent columns contain the marker data for each individual in the population.  

Marker data can be coded in one of three formats: 
\itemize{
\item "numeric" = markers are coded based on the dosage of the alternate allele, taking on integer values between 0 and ploidy (fractional values not allowed)
\item "AB" format (e.g., AAAB, ABBB for tetraploids)
\item "ACGT" format (e.g., AAAT, GGCC for tetraploids)
}
Only bi-allelic markers are allowed.  Missing marker data will be imputed with the population mode (most frequent value) for each marker.
}
\value{
An object of class \code{\link[=GWASpoly-class]{GWASpoly}}.
}
\references{
Rosyara et al. (2016) Software for genome-wide association studies in autopolyploids and its application to potato.  Plant Genome 9.
}
\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}
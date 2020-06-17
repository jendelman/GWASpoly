\name{GWASpoly.K-class}
\Rdversion{1.1}
\docType{class}
\alias{GWASpoly.K-class}

\title{Class \code{GWASpoly.K}}
\description{
Extends the class \code{\link[=GWASpoly-class]{GWASpoly}} by the addition of
a kinship matrix, using the function \code{\link{set.K}}.
}
\section{Slots}{
  \describe{
    \item{\code{ploidy}:}{Ploidy level (e.g., 2 for diploid, 4 for tetraploid)}
    \item{\code{map}:}{Data frame with five columns: Marker, Chrom, Position, Ref, Alt}
    \item{\code{pheno}:}{Data frame of phenotypes, with genotype identifiers as rownames.}  
    \item{\code{geno}:}{Numeric matrix (lines x markers) with marker data coded as the dosage of the alternate allele.}
    \item{\code{fixed}:}{Data frame of factors and covariates to model as fixed effects.}
    \item{\code{K}:}{Positive semidefinite kinship matrix.}
	}
}
\examples{
showClass("GWASpoly.K")
}

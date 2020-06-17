\name{GWASpoly.fitted-class}
\Rdversion{1.1}
\docType{class}
\alias{GWASpoly.fitted-class}

\title{Class \code{GWASpoly.fitted}}
\description{
Extends the class \code{\link[=GWASpoly.K-class]{GWASpoly.K}} by the addition of
GWAS results, using the function \code{\link{GWASpoly}}.
}
\section{Slots}{
  \describe{
    \item{\code{ploidy}:}{Ploidy level (e.g., 2 for diploid, 4 for tetraploid)}
    \item{\code{map}:}{Data frame with five columns: Marker, Chrom, Position, Ref, Alt}
    \item{\code{pheno}:}{Data frame of phenotypes, with genotype identifiers as rownames.}  
    \item{\code{geno}:}{Numeric matrix (lines x markers) with marker data coded as the dosage of the alternate allele.}
    \item{\code{fixed}:}{Data frame of factors and covariates to model as fixed effects.}
    \item{\code{K}:}{Positive semidefinite covariance (kinship) matrix.}
    \item{\code{scores}:}{List of data.frames, one for each trait. The data frames contain the -log10(p) results, organized as markers x models.}
    \item{\code{effects}:}{List of data.frames, one for each trait. The data frames contain the marker effect estimates, organized as markers x models. Effect estimates are not provided for the diplo-general or general models. Effect estimates for the additive model are based on marker coding {0,1,2,..ploidy}. For the diplo-additive model, the coding is {0,1,2}. For dominance models, the coding is {0,1}.}
    \item{\code{params}:}{List with the following model parameters:
    \itemize{
    \item fixed = Names of fixed effects
    \item fixed.type = Vector of fixed effect types: either "numeric" (covariate) or "factor"
    \item n.PC = Number of principal components included as covariates
    \item min.MAF = Minimum minor allele frequency
    \item max.geno.freq = Maximum genotype frequency (after applying dominance relations)
    \item P3D = TRUE/FALSE whether the P3D approximation was used
    \item models = Names of marker-effect models 
    }
    }    
	}
}
\examples{
showClass("GWASpoly.fitted")
}

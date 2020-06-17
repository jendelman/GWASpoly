\name{set.K}
\alias{set.K}
\title{
Creates kinship (relationship) matrix 
}
\description{
Creates the kinship matrix for GWAS
}
\usage{
set.K(data, K = NULL)
}

\arguments{
  \item{data}{
Object of class \code{\link[=GWASpoly-class]{GWASpoly}}.
}
  \item{K}{
Optional: User-supplied kinship matrix
}
}
\details{
By default, K is computed as \eqn{MM^T}, where M is the centered genotype matrix (lines x markers).  For GWAS, the overall scaling of K is irrelevant.  At present, K is scaled such that the mean of its diagonal elements is 1.  Alternatively, the user can supply any positive semidefinite K (with row.names that match the genotype identifiers in data@geno).
}
\value{
Object of class \code{\link[=GWASpoly.K-class]{GWASpoly.K}}.
}
\references{
Rosyara et al. (in press) Software for genome-wide association studies in autopolyploids and its application to potato.  Plant Genome
}
\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}
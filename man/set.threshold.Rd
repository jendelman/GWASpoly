\name{set.threshold}
\alias{set.threshold}
\title{
Sets the significance threshold.
}
\description{
Sets the significance threshold.
}
\usage{
set.threshold(data, method, level = 0.05, n.permute = 1000, n.core = 1)
}
\arguments{
  \item{data}{
Object of class \code{\link[=GWASpoly.fitted-class]{GWASpoly.fitted}}.
}
  \item{method}{
One of the following: "Bonferroni","FDR","permute"
}
  \item{level}{
Genome-wide false positive rate for the Bonferroni or permutation methods; false discovery rate for method FDR.
}
  \item{n.permute}{
Number of permutations for method "permute"
}
  \item{n.core}{
Number of cores to use for multicore processing (only Linux and MacOS).
}
}
\details{
The FDR method is based on version 1.30.0 of the qvalue package.
}
\value{
Object of class \code{\link[=GWASpoly.thresh-class]{GWASpoly.thresh}}.}
\references{
Rosyara et al. (2016) Software for genome-wide association studies in autopolyploids and its application to potato.  Plant Genome 9.
}
\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}
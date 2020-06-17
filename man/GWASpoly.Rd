\name{GWASpoly}
\alias{GWASpoly}
\title{
GWAS
}
\description{
Genome-wide association studies
}
\usage{
GWASpoly(data, models, traits = NULL, params = NULL, n.core = 1, quiet = F)
}
\arguments{
  \item{data}{
Object of class \code{\link[=GWASpoly.K-class]{GWASpoly.K}}.
}
  \item{models}{
Vector of model names for the marker effect.  See Details for options.
}
  \item{traits}{
Vector of trait names. If NULL, all traits are analyzed.
}
 \item{params}{
List of model parameters, created using \code{\link{set.params}}.  If NULL, default values from \code{\link{set.params}} are used.
 }
   \item{n.core}{
Number of cores to use for multicore processing (only Linux and MacOS).
}
  \item{quiet}{
TRUE/FALSE whether to suppress output charting the progress of the analysis
}
}
\details{
The following marker-effect models are available:
\itemize{
\item "additive" = Indicates the marker effect is proportional to the dosage of the alternate allele.
\item "X-dom" = where X can be any integer between 1 and ploidy/2 and refers to the allele dosage needed for complete dominance (e.g., "1-dom" = simplex dominance, "2-dom" = duplex dominance).  The software tries both dominance patterns for a given dosage model, e.g., whether the reference or alternate allele is dominant.
\item "diplo-general" = All heterozygotes have the same effect.
\item "diplo-additive" = All heterozygotes have the same effect, constrained to be halfway between the homozygous effects.
\item "general" = There are no constraints on the effects of the different dosage levels.
}

To specify additional model parameters, such as the inclusion of fixed effects (Q matrix) and the minimum minor allele frequency, use \code{\link{set.params}}.
}
\value{
Object of class \code{\link[=GWASpoly.fitted-class]{GWASpoly.fitted}}.
}
\references{
Rosyara et al. (in press) Software for genome-wide association studies in autopolyploids and its application to potato.  Plant Genome
}
\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}
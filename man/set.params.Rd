\name{set.params}
\alias{set.params}
\title{
Set parameters for GWAS
}
\description{
Set parameters for GWAS
}
\usage{
set.params(fixed=NULL,fixed.type=NULL,n.PC=0,MAF=0.05,geno.freq=0.95,P3D=T)
}
\arguments{
  \item{fixed}{
Vector of names of fixed effects
}
  \item{fixed.type}{
Vector of effect types ("numeric" or "factor"), corresponding to the effects listed in "fixed". 
}
  \item{n.PC}{
Number of principal components to include as covariates.
}
  \item{MAF}{
Minimum minor allele frequency
}
  \item{geno.freq}{
Maximum genotype frequency (after applying dominance relations)
}
  \item{P3D}{
TRUE/FALSE whether to use the P3D approximation (variance components not re-estimated for every marker).
}
}
\details{
The list returned by the function should be passed to \code{\link{GWASpoly}}.
}
\value{
A list with the following components 
\item{fixed}{Names of fixed effects}
\item{fixed.type}{Types of fixed effects}
\item{n.PC}{Number of principal components to include as covariates}
\item{min.MAF}{Minimum minor allele frequency}
\item{max.geno.freq}{Maximum genotype frequency (after applying dominance relations)}
\item{P3D}{TRUE/FALSE whether to use the P3D approximation}
}
\references{
Rosyara et al. (in press) Software for genome-wide association studies in autopolyploids and its application to potato.  Plant Genome
}
\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}
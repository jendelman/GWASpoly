\name{get.QTL}
\alias{get.QTL}
\title{
Extracts significant QTL
}
\description{
Extracts significant QTL
}
\usage{
get.QTL(data, traits = NULL, models = NULL)
}
\arguments{
  \item{data}{
Object of class \code{\link[=GWASpoly.thresh-class]{GWASpoly.thresh}}.
}
  \item{traits}{
Vector of trait names.
}
  \item{models}{
Vector of model names.
}
}
\details{
The function returns which markers exceed the significance threshold.
}
\value{
Data frame with the following columns: Trait, Model, Threshold, Marker, Chrom, Position, Ref, Alt, Score, Effect. 
}
\references{
Rosyara et al. (in press) Software for genome-wide association studies in autopolyploids and its application to potato.  Plant Genome
}
\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}
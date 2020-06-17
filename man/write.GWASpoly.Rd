\name{write.GWASpoly}
\alias{write.GWASpoly}
\title{
Write GWASpoly results to file
}
\description{
Write GWASpoly results to file
}
\usage{
write.GWASpoly(data, trait, filename, what = "scores", delim = ",")
}
\arguments{
  \item{data}{
Object inheriting from class \code{\link[=GWASpoly.fitted-class]{GWASpoly.fitted}}.
}
  \item{trait}{
Name of trait to output
}
  \item{filename}{
Name of file
}
  \item{what}{
Two valid options: "scores" or "effects"
}
  \item{delim}{
Delimiter to use in the output file (default is comma).
}
}
\details{
By default the "scores", which are the -log10(p) results, are written to a comma-separated values file.  If what = "effects", the estimated effects for the markers are written to file instead of the scores.  (Effects are not returned for the general and diplo-general models.  See \code{\link[=GWASpoly.fitted-class]{GWASpoly.fitted}}.)
}
\value{
NULL
}
\references{
Rosyara et al. (in press) Software for genome-wide association studies in autopolyploids and its application to potato.  Plant Genome
}
\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}
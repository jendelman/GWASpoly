\name{manhattan.plot}
\alias{manhattan.plot}
\title{
Manhattan plot
}
\description{
Manhattan plot
}
\usage{
manhattan.plot(data, trait, model, cex = 1, y.max = NULL, filename = NULL)
}
\arguments{
  \item{data}{
Object inheriting from class \code{\link[=GWASpoly.fitted-class]{GWASpoly.fitted}}.
}
  \item{trait}{
Name of trait to plot
}
  \item{model}{
Name of model to plot
}
  \item{cex}{
Controls size of points
}
  \item{y.max}{
Controls the upper limit for y = -log10(p).  If NULL, y.max is determined from the data.
}
  \item{filename}{
If NULL, plot appears on screen.  Providing a file name will generate a postscript file.
}
}
\details{
If a threshold has been set with \code{\link{set.threshold}}, it is displayed as well.
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
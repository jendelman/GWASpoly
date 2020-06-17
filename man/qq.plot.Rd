\name{qq.plot}
\alias{qq.plot}
\title{
QQ plot
}
\description{
Generates a quantile-quantile plot of the observed vs. expected -log10(p) values.
}
\usage{
qq.plot(data, trait, model, cex = 1, filename = NULL)
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
  \item{filename}{
If NULL, plot appears on screen.  Providing a file name will generate a postscript file.
}
}
\details{
QQ plots are used to check the inflation of p-values in GWAS.
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
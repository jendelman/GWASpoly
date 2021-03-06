% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.QTL.R
\name{get.QTL}
\alias{get.QTL}
\title{Extract significant QTL}
\usage{
get.QTL(data, traits = NULL, models = NULL, bp.window = NULL)
}
\arguments{
\item{data}{Output from \code{set.threshold}}

\item{traits}{Vector of trait names (by default, all traits)}

\item{models}{Vector of model names (by default, all models)}

\item{bp.window}{prune output to return only the most significant marker within this window size}
}
\value{
Data frame with results. Score = -log10(p). Effect = marker effect (not available for the general and diplo-general models because there are multiple effects).
}
\description{
Output a table with significant markers
}
\details{
To return all significant markers (original behavior of the function), use \code{bp.window=NULL}. Assumes input map position in bp.
}

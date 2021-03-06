% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qq.plot.R
\name{qq.plot}
\alias{qq.plot}
\title{Quantile-Quantile (QQ) Plot}
\usage{
qq.plot(data, trait, models = NULL)
}
\arguments{
\item{data}{Variable of class \code{GWASpoly.fitted}}

\item{trait}{Trait name}

\item{models}{Vector of model names (by default, all models plotted)}
}
\value{
ggplot2 object
}
\description{
Inspect p-value inflation using a QQ plot
}

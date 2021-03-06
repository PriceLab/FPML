% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_functions.R
\name{mattcc.curve}
\alias{mattcc.curve}
\title{Create Matthews Correlation Coefficient Data for Plotting}
\usage{
mattcc.curve(true.values, pred.probs, num.points = 101)
}
\arguments{
\item{true.values}{The true classifications for the data}

\item{pred.probs}{The predicted probabilities for the data}

\item{num.points}{The number of different thresholds to set between 0 and 1 (default = 101)}
}
\value{
A numeric vector of MCC values over the supplied number of points
}
\description{
Using the formula for Matthews Correlation Coefficient and a given number of points ranging
from 0 to 100, create data for MCC over a set of thresholds between 0 and 1. These data can then
be used to plot MCC curves
}

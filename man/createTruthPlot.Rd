% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildModels.R
\name{createTruthPlot}
\alias{createTruthPlot}
\title{Create the "truth plot"}
\usage{
createTruthPlot(model, modelType, seed, threshold, X_test, y_test, plotFile)
}
\arguments{
\item{model}{A model trained on the annotated dataset, either a linear or gradient-boosted model}

\item{modelType}{Either "boosted" or "linear" to indicate the type of model provided}

\item{seed}{One of "16", "20", or "both", denoting which seed is being used}

\item{threshold}{A numeric from 0-1 indicating the cutoff for what constitutes a "positive"
probability. This is an INCLUSIVE value; positives are greater than or equal to the threshold}

\item{X_test}{The matrix of test data for predictors, generally created using "prepModelData"}

\item{y_test}{The matrix of test data for response, generally created using "prepModelData"}

\item{plotFile}{A location to use for creating the png plot file}
}
\value{
The "truth data frame", containing the Wellington/HINT scores for each point in the
test set, the predicted probability of a ChIPseq hit, the binary prediction based on that
probability and the supplied threshold, the true ChIPseq hit value, and the corresponding
label of TP/TN/FP/FN.

It also constructs and saves a plot of HINT vs Wellington scores, where both have been
transformed via asinh() and the Wellington scores have been converted via absolute value.
Points are coded by confusion matrix outcome.
}
\description{
Create a "truth" plot, which plots HINT v. Wellington scores given a model and a threshold
to use for determining a "positive"
}

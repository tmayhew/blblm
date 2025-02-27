% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{confint.blblm}
\alias{confint.blblm}
\title{blblm Confidence Interval}
\usage{
\method{confint}{blblm}(object, parm = NULL, level = 0.95, ...)
}
\arguments{
\item{object}{model fit produced by blblm(formula, data, ...)}

\item{parm}{character vector, parameters to give confidence interval for. Defaults to all predictors.}

\item{level}{numeric, confidence interval level.}

\item{...}{other arguments.}
}
\description{
Constructs a confidence interval for parameters of the model
}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/model.R
\name{trajectorymixture}
\alias{trajectorymixture}
\title{Fit a mixture of mean trajectories.}
\usage{
trajectorymixture(trajectories, K, bfn, maxiter = 100, eps = 1e-05)
}
\arguments{
\item{trajectories}{Data frame containing trajectory ids in the first column, observation times in the second, and observation values in the third.}

\item{K}{The number of components in the mixture.}

\item{bfn}{The basis function mapping times to a design matrix.}

\item{maxiter}{The maximum number of EM iterations.}

\item{eps}{The smallest relative change in log likelihood needed to stop.}
}
\description{
Fit a mixture of mean trajectories.
}
\author{
pschulam
}


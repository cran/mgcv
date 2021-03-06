% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mgcvExports.R
\name{mini.roots}
\alias{mini.roots}
\title{Obtain square roots of penalty matrices}
\usage{
mini.roots(S, off, np, rank = NULL)
}
\arguments{
\item{S}{a list of penalty matrices, in packed form.}

\item{off}{a vector where the i-th element is the offset for the i-th matrix. 
The elements in columns \code{1:off[i]} of \code{B[[i]]} will be equal to zero.}

\item{np}{total number of parameters.}

\item{rank}{here \code{rank[i]} is optional supplied rank of \code{S[[i]]}. Set \code{rank[i] < 1}, or 
\code{rank=NULL} to estimate.}
}
\value{
A list of matrix square roots such that \code{S[[i]]=B[[i]]\%*\%t(B[[i]])}.
}
\description{
INTERNAL function to obtain square roots, \code{B[[i]]}, of the penalty matrices \code{S[[i]]}'s having as few
             columns as possible.
}
\author{
Simon N. Wood <simon.wood@r-project.org>.
}

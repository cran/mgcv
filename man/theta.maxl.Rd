%
\name{theta.maxl}
\alias{theta.maxl}
\title{
Estimate theta of the Negative Binomial by Maximum Likelihood
}
\description{
Given the estimated mean vector, estimate \code{theta} of the Negative Binomial
Distribution. Based on \code{theta.ml} in \code{library(MASS)}, but modified to avoid \code{theta} overflow problems. }
\usage{
theta.maxl(y, mu, n, limit=20, eps=.Machine$double.eps^0.25, trace=FALSE)
}
\arguments{
\item{y}{
Vector of observed values from the Negative Binomial.
}
\item{mu}{
Estimated mean vector.
}
\item{n}{
Residual degrees of freedom (assuming \code{theta} known).
}
\item{limit}{
Limit on the number of iterations.
}
\item{eps}{
Tolerance to determine convergence.
}
\item{trace}{
logical: should iteration progress be printed?
}}
\value{
The required estimate of \code{theta}, as a scalar.
}
\seealso{
\code{\link{gam}}, \code{\link{gam.nbut}}
}
\keyword{models}

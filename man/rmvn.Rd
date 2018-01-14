\name{rmvn}
\alias{rmvn}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generate multivariate normal deviates}
\description{ Generates multivariate normal random deviates.
}

\usage{
rmvn(n,mu,V)
}
\arguments{
\item{n}{number of simulated vectors required.}

\item{mu}{the mean of the vectors: either a single vector of length \code{p=ncol(V)} or an \code{n} by \code{p} matrix.}

\item{V}{A positive semi definite covariance matrix.}

}

\value{ An \code{n} row matrix, with each row being a draw from a multivariate normal density with covariance matrix \code{V} and mean vector \code{mu}. Alternatively each row may have a different mean vector if \code{mu} is a vector.

}


\details{Uses a `square root' of \code{V} to transform standard normal deviates to multivariate normal with the correct covariance matrix. 
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}


\seealso{\code{\link{ldTweedie}}, \code{\link{Tweedie}}}

\examples{
library(mgcv)
V <- matrix(c(2,1,1,2),2,2) 
mu <- c(1,3)
n <- 1000
z <- rmvn(n,mu,V)
crossprod(sweep(z,2,colMeans(z)))/n ## observed covariance matrix
colMeans(z) ## observed mu 
}
\keyword{models} \keyword{regression}%-- one or more ..



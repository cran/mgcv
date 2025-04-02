\name{dpnorm}
\alias{dpnorm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stable evaluation of difference between normal c.d.f.s}
\description{Evaluates the difference between two \eqn{N(0,1)}{N(0,1)} cumulative distribution functions avoiding cancellation error.  
}
\usage{
dpnorm(x0,x1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x0}{vector of lower values at which to evaluate standard normal distribution function.}
 \item{x1}{vector of upper values at which to evaluate standard normal distribution function.}
}

\details{ Equivalent to \code{pnorm(x1)-pnorm(x0)}, but stable when \code{x0} and \code{x1} values are very close, or in the upper tail of the standard normal.}
 

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
require(mgcv)
x <- seq(-10,10,length=10000)
eps <- 1e-10
y0 <- pnorm(x+eps)-pnorm(x) ## cancellation prone
y1 <- dpnorm(x,x+eps)       ## stable
## illustrate stable computation in black, and
## cancellation prone in red...
par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(log(y1),log(y0),type="l")
lines(log(y1[x>0]),log(y0[x>0]),col=2)
plot(x,log(y1),type="l")
lines(x,log(y0),col=2)

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










\name{FFdes}
\alias{FFdes}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Level 5 fractional factorial designs}
\description{Computes level 5 fractional factorial designs for up to 120 factors
using the agorithm of Sanchez and Sanchez (2005), and optionally central composite designs.
}
\usage{
FFdes(size=5,ccd=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{size}{number of factors up to 120.}
\item{ccd}{if \code{TRUE}, adds points along each axis at the same distance from the origin as the points in the
fractional factorial design, to create the outer points of a central composite design. Add central points to complete.}
}

\details{Basically a translation of the code provided in the appendix of Sanchez and Sanchez (2005). 
}

\references{
Sanchez, S. M. & Sanchez, P. J. (2005) Very large fractional factorial and central composite designs.
ACM Transactions on Modeling and Computer Simulation. 15: 362-377
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
  require(mgcv)
  plot(rbind(0,FFdes(2,TRUE)),xlab="x",ylab="y",
       col=c(2,1,1,1,1,4,4,4,4),pch=19,main="CCD")
  FFdes(5)
  FFdes(5,TRUE)
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










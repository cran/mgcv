\name{cSplineDes}
\alias{cSplineDes}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Evaluate cyclic B spline basis}
\description{ Uses \code{splineDesign} to set up the model matrix for a cyclic B-spline basis.
}

\usage{
cSplineDes(x, knots, ord = 4)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ covariate values for smooth.}

\item{knots}{The knot locations: the range of these must include all the data.}

\item{ord}{ order of the basis. 4 is a cubic spline basis. Must be >1.}

}

\details{ The routine is a wrapper that sets up a B-spline basis, where the basis functions wrap at the first and 
last knot locations.} 

\value{ A matrix with \code{length(x)} rows and \code{length(knots)-1} columns.

}
 

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\seealso{\code{\link{cyclic.p.spline}}}

\examples{
 x <- 0:100/100;k<- 0:5/5
 X <- cSplineDes(x,k)
 plot(x,X[,1],type="l"); for (i in 2:5) lines(x,X[,i],col=i)
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










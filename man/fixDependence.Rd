\name{fixDependence}
\alias{fixDependence}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Detect linear dependencies of one matrix on another}
\description{Identifies columns of a matrix \code{X2} which are linearly
dependent on columns of a matrix \code{X1}. Primarily of use in setting up 
identifiability constraints for nested GAMs. 
}
\usage{
fixDependence(X1,X2,tol=.Machine$double.eps^.5)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{X1}{ A matrix.}
  \item{X2}{ A matrix, the columns of which may be partially linearly
  dependent on the columns of \code{X1}.}
  \item{tol}{The tolernce to use when assessing linear dependence.}
}
\details{ The algorithm uses a simple approach based on QR decomposition: see
code for details.

}
\value{ An array of the columns of \code{X2} which are linearly dependent on
columns of \code{X1}. \code{NULL} if the two matrices are independent.
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}



\examples{
n<-20;c1<-4;c2<-7
X1<-matrix(runif(n*c1),n,c1)
X2<-matrix(runif(n*c2),n,c2)
X2[,3]<-X1[,2]+X2[,4]*.1
X2[,5]<-X1[,1]*.2+X1[,2]*.04
fixDependence(X1,X2)
}
\keyword{models} \keyword{regression}%-- one or more ..



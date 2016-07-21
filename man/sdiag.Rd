\name{sdiag}
\alias{sdiag}
\alias{sdiag<-}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract or modify diagonals of a matrix}
\description{ Extracts or modifies sub- or super- diagonals of a matrix.
}
\usage{
sdiag(A,k=0)
sdiag(A,k=0) <- value
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{A}{a matrix}
\item{k}{sub- (negative) or super- (positive) diagonal of a matrix. 0 is the leading diagonal.}
\item{value}{single value, or vector of the same length as the diagonal.}

}

\value{A vector containing the requested diagonal, or a matrix with the requested diagonal replaced by \code{value}.

}



\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\examples{
require(mgcv)
A <- matrix(1:35,7,5)
A
sdiag(A,1) ## first super diagonal
sdiag(A,-1) ## first sub diagonal

sdiag(A) <- 1 ## leading diagonal set to 1
sdiag(A,3) <- c(-1,-2) ## set 3rd super diagonal 

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










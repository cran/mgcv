\name{mroot}
\alias{mroot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Smallest square root of matrix}
\description{ Find a square root of a positive semi-definite matrix, 
having as few columns as possible. Uses either pivoted choleski 
decomposition or singular value decomposition to do this.  
}
\usage{
mroot(A,rank=NULL,method="chol")
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{A}{ The positive semi-definite matrix, a square root of which is 
to be found.}
\item{rank}{if the rank of the matrix \code{A} is known then it should 
be supplied.}
\item{method}{ \code{"chol"} to use pivoted choloeski decompositon, 
which is fast but tends to over-estimate rank. \code{"svd"} to use 
singular value decomposition, which is slow, but is the most accurate way 
to estimate rank.}
}

\details{ The routine uses an LAPACK SVD routine, or the LINPACK pivoted 
Choleski routine. It is primarily of use for turning penalized regression 
problems into ordinary regression problems.} 
\value{ A matrix, \eqn{ {\bf B}}{B} with as many columns as the rank of 
\eqn{ {\bf A}}{A}, and such that \eqn{ {\bf A} = {\bf BB}^\prime}{A=BB'}.}
 

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










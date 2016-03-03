\name{trichol}
\alias{trichol}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Choleski decomposition of a tri-diagonal matrix}
\description{ 
Computes Choleski decomposition of a (symmetric positive definite) tri-diagonal matrix stored as a leading diagonal and sub/super diagonal.
}
\usage{
trichol(ld,sd)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{ld}{leading diagonal of matrix}
\item{sd}{sub-super diagonal of matrix}

}

\value{ A list with elements \code{ld} and \code{sd}. \code{ld} is the leading diagonal and \code{sd} is the super diagonal of bidiagonal matrix \eqn{\bf B}{B} where \eqn{{\bf B}^T{\bf B} = {\bf T}}{B'B=T} and \eqn{\bf T}{T} is the original tridiagonal matrix. 

}

\details{Calls \code{dpttrf} from \code{LAPACK}. The point of this is that it has \eqn{O(n)}{O(n)} computational cost, rather than the \eqn{O(n^3)}{O(n^3)} required by dense matrix methods.

}
 
\seealso{\code{\link{bandchol}}}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\references{
Anderson, E., Bai, Z., Bischof, C., Blackford, S., Dongarra, J., Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A. and Sorensen, D., 1999. LAPACK Users' guide (Vol. 9). Siam.
}

\examples{
require(mgcv)
## simulate some diagonals...
set.seed(19); k <- 7
ld <- runif(k)+1
sd <- runif(k-1) -.5

## get diagonals of chol factor...
trichol(ld,sd)

## compare to dense matrix result...
A <- diag(ld);for (i in 1:(k-1)) A[i,i+1] <- A[i+1,i] <- sd[i]
R <- chol(A)
diag(R);diag(R[,-1])

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










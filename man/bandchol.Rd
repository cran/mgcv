\name{bandchol}
\alias{bandchol}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Choleski decomposition of a band diagonal matrix}
\description{ 
Computes Choleski decomposition of a (symmetric positive definite) band-diagonal matrix, \code{A}. 
}
\usage{
bandchol(B)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{B}{An n by k matrix containing the diagonals of the matrix \code{A} to be decomposed. First row is leading diagonal, next is first sub-diagonal, etc. sub-diagonals are zero padded at the end. Alternatively gives \code{A} directly, i.e. a square matrix with 2k-1 non zero diagonals (those from the lower triangle are not accessed).}
}

\value{Let \code{R} be the factor such that \code{t(R)\%*\%R = A}. \code{R} is upper triangular and if the rows of \code{B} contained the diagonals of \code{A} on entry, then what is returned is an n by k matrix containing the diagonals of \code{R}, packed as \code{B} was packed on entry. If \code{B} was square on entry, then \code{R} is returned directly. See examples.  
}

\details{Calls \code{dpbtrf} from \code{LAPACK}. The point of this is that it has \eqn{O(k^2n)}{O(k^2n)} computational cost, rather than the \eqn{O(n^3)}{O(n^3)} required by dense matrix methods.

}
 

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\references{
Anderson, E., Bai, Z., Bischof, C., Blackford, S., Dongarra, J., Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A. and Sorensen, D., 1999. LAPACK Users' guide (Vol. 9). Siam.
}

\examples{
require(mgcv)
## simulate a banded diagonal matrix
n <- 7;set.seed(8)
A <- matrix(0,n,n)
sdiag(A) <- runif(n);sdiag(A,1) <- runif(n-1)
sdiag(A,2) <- runif(n-2)
A <- crossprod(A) 

## full matrix form...
bandchol(A)
chol(A) ## for comparison

## compact storage form...
B <- matrix(0,3,n)
B[1,] <- sdiag(A);B[2,1:(n-1)] <- sdiag(A,1)
B[3,1:(n-2)] <- sdiag(A,2)
bandchol(B)

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










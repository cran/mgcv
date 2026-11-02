\name{bandchol}
\alias{bandchol}
\alias{bandmult}
\alias{bandsolve}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Cholesky decomposition of a band diagonal matrix}
\description{ 
Computes Cholesky decomposition of a (symmetric positive definite) band-diagonal matrix, \code{A}, or optionally a partial Cholesky factor if the matrix consists only of diagonals of a PD matrix without necessarily being PD. Solves or multiplies with the corresponding banded Cholesky factor \code{R}.
}
\usage{
bandchol(B,partial=FALSE)
bandmult(B,x,trans=FALSE)
bandsolve(B,x,trans=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{B}{An n by k matrix containing the diagonals of the matrix \code{A} to be decomposed, or factor \code{R}. First row is leading diagonal, next is first super-diagonal, etc. Super-diagonals are zero padded at the end. Alternatively, for \code{bandchol} only, gives \code{A} directly, i.e. a square matrix with 2k-1 non zero diagonals (those from the lower triangle are not accessed).}
\item{partial}{compute partial factor if \code{B} contains diagonals of a PD matrix, but the matrix with just these bands is not guaranteed PD. See details.}
\item{x}{an n-vector or n row matrix to multiply by R (\code{bandmult}) or its inverse (\code{bandsolve}).}
\item{trans}{if \code{TRUE} then transpose \code{R} before multiplication or solving.}
}

\value{For \code{bandchol}, let \code{R} be the factor such that \code{t(R)\%*\%R = A}. \code{R} is upper triangular and if the rows of \code{B} contained the diagonals of \code{A} on entry, then what is returned is an n by k matrix containing the diagonals of \code{R}, packed as \code{B} was packed on entry. If \code{B} was square on entry, then \code{R} is returned directly. Returns \code{NULL} if matrix defined by \code{B} is not positive definite and \code{partial=FALSE}. See examples.

For \code{bandmult} or \code{bandsolve} a vector or matrix of the same dimension as \code{x} containing the required result. 
}

\details{For \code{bandchol} calls \code{dpbtrf} from \code{LAPACK}. The point of this is that it has \eqn{O(k^2n)}{O(k^2n)} computational cost, rather than the \eqn{O(n^3)}{O(n^3)} required by dense matrix methods. \code{bandmult} calls BLAS routine \code{dtbmv} (repeatedly if \code{x} has several columns). \code{bandsolve} calls \code{dtbtrs} from \code{LAPACK}.

If \code{partial=TRUE} then \code{bandchol} performs a partial decomposition (calls compiled C code, not LAPACK). The diagonals in \code{B} should come from a PD matrix, but the matrix consisting of just those diagonals need not be PD. The leading diagonal must be positive. The factor computed is of a block diagonal approximation to the matrix whose bands are in \code{B}. The approximation's elements either match the target matrix exactly, or are set to zero. The construction runs the  Cholesky decomposition until an imaginary diagonal factor element is encountered, at which point a new diagonal block is started, and all already computed elements beyond the previous block are set to zero. This is repeated until a full factor is produced. See examples. 

}
 

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\references{
Anderson, E., Bai, Z., Bischof, C., Blackford, S., Dongarra, J., Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A. and Sorensen, D., 1999. LAPACK Users' guide (Vol. 9). SIAM.
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
R <- chol(A); R ## for comparison

## compact storage form...
B <- matrix(0,3,n)
B[1,] <- sdiag(A);B[2,1:(n-1)] <- sdiag(A,1)
B[3,1:(n-2)] <- sdiag(A,2)
RB <- bandchol(B); RB

## bandmult and bandsolve with a vector...
x <- runif(n)

drop(R\%*\%x) ## direct
bandmult(RB,x) ## banded
drop(t(R)\%*\%x)
bandmult(RB,x,trans=TRUE)

backsolve(R,x)
bandsolve(RB,x)
forwardsolve(t(R),x)
bandsolve(RB,x,trans=TRUE)

## bandmult and bandsolve with a matrix...
x <- matrix(runif(2*n),n,2)

drop(R\%*\%x)
bandmult(RB,x)
drop(t(R)\%*\%x)
bandmult(RB,x,trans=TRUE)

backsolve(R,x)
bandsolve(RB,x)
forwardsolve(t(R),x)
bandsolve(RB,x,trans=TRUE)

## Partial factor example
set.seed(0); n <- 7
B0 <- crossprod(matrix(runif(n*n),n,n)) ## PD matrix
b <- 3 ## super-diagonals to use
B <- matrix(0,b+1,n) ## packed storage for b+1 diags
B[1,] <- diag(B0)
for (i in 1:b) B[i+1,1:(n-i)] <- sdiag(B0,i)
bandchol(B) ## not PD
R <- bandchol(B,partial=TRUE) ## get partial factor
## pack back into full matrix for demo...
R0 <- diag(R[1,]); for (i in 1:b) sdiag(R0,i) <- R[i+1,1:(n-i)]
BR <- crossprod(R0)
## c.f. non zero crossprod(R0) and corresponding B0  
BR;B0*as.numeric(BR!=0) 
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










\name{choldrop}
\alias{choldrop}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Down date a Cholesky factor on dropping a row/col}
\description{Given a Cholesky factor, \code{R}, of a matrix, \code{A}, finds the Cholesky factor of \code{A[-k,-k]},
where \code{k} is an integer.  
}
\usage{
choldrop(R,k)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{R}{Cholesky factor of a matrix, \code{A}.}
\item{k}{row and column of \code{A} to drop.}
}

\details{If \code{R} is upper triangular then \code{t(R[,-k])\%*\%R[,-k] == A[-k,-k]}, but \code{R[,-k]} has elements on the first sub-diagonal, from its kth column onwards. To get from this to a triangular Cholesky factor of \code{A[-k,-k]} we can apply a sequence of Given rotations from the left to eliminate the sub-diagonal elements. The routine does this. If \code{R} is a lower triangular factor then Givens rotations from the right are needed to remove the extra elements. If \code{n} is the dimension of \code{R} then the update has O(n^2) computational cost.

Note that the update is vector oriented, and is hence not susceptible to speed up by use of an optimized BLAS. The update is set up to be relatively Cache friendly, in that in the upper triangular case successive Givens rotations are stored for sequential application columnwise, rather than being applied rowwise as soon as they are computed. Even so, the upper triangular update is slightly slower than the lower triangular update. 

}
 

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
  require(mgcv)
  set.seed(0)
  n <- 6
  A <- crossprod(matrix(runif(n*n),n,n))
  R0 <- chol(A)
  k <- 3
  Rd <- choldrop(R0,k)
  range(Rd-chol(A[-k,-k]))
  Rd;chol(A[-k,-k])
  
  ## same but using lower triangular factor A = LL'
  L <- t(R0)
  Ld <- choldrop(L,k)
  range(Ld-t(chol(A[-k,-k])))
  Ld;t(chol(A[-k,-k]))
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










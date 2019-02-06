\name{choldrop}
\alias{choldrop}
\alias{cholup}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Deletion and rank one Cholesky factor update}
\description{Given a Cholesky factor, \code{R}, of a matrix, \code{A}, \code{choldrop} finds the Cholesky factor of \code{A[-k,-k]},
where \code{k} is an integer. \code{cholup} finds the factor of \eqn{A + uu^T}{A+uu'} (update) or \eqn{A - uu^T}{A-uu'} (downdate).
}
\usage{
choldrop(R,k)
cholup(R,u,up)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{R}{Cholesky factor of a matrix, \code{A}.}
\item{k}{row and column of \code{A} to drop.}
\item{u}{vector defining rank one update.}
\item{up}{if \code{TRUE} compute update, otherwise downdate.}
}

\details{First consider \code{choldrop}. If \code{R} is upper triangular then \code{t(R[,-k])\%*\%R[,-k] == A[-k,-k]}, but \code{R[,-k]} has elements on the first sub-diagonal, from its kth column onwards. To get from this to a triangular Cholesky factor of \code{A[-k,-k]} we can apply a sequence of Givens rotations from the left to eliminate the sub-diagonal elements. The routine does this. If \code{R} is a lower triangular factor then Givens rotations from the right are needed to remove the extra elements. If \code{n} is the dimension of \code{R} then the update has \eqn{O(n^2)}{O(n^2)} computational cost.

\code{cholup} (which assumes \code{R} is upper triangular) updates based on the observation that \eqn{ R^TR + uu^T = [u,R^T][u,R^T]^T =  [u,R^T]Q^TQ[u,R^T]^T}{R'R + uu' = [u,R'][u,R']' =  [u,R']Q'Q[u,R']'}, and therefore we can construct \eqn{Q}{Q} so that \eqn{Q[u,R^T]^T=[0,R_1^T]^T}{Q[u,R']'=[0,R1']'}, where \eqn{R_1}{R1} is the modified factor. \eqn{Q}{Q} is constructed from a sequence of Givens rotations in order to zero the elements of \eqn{u}{u}. Downdating is similar except that hyperbolic rotations have to be used in place of Givens rotations --- see Golub and van Loan (2013, section 6.5.4) for details. Downdating only works if  \eqn{A - uu^T}{A-uu'} is positive definite. Again the computational cost is \eqn{O(n^2)}{O(n^2)}.

Note that the updates are vector oriented, and are hence not susceptible to speed up by use of an optimized BLAS. The updates are set up to be relatively Cache friendly, in that in the upper triangular case successive Givens rotations are stored for sequential application column-wise, rather than being applied row-wise as soon as they are computed. Even so, the upper triangular update is slightly slower than the lower triangular update. 
}

\references{
Golub GH and CF Van Loan (2013) Matrix Computations (4th edition) Johns Hopkins
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

  ## Rank one update example
  u <- runif(n)
  R <- cholup(R0,u,TRUE)
  Ru <- chol(A+u \%*\% t(u)) ## direct for comparison
  R;Ru
  range(R-Ru)

  ## Downdate - just going back from R to R0
  Rd <-  cholup(R,u,FALSE)
  R0;Rd
  range(R-Ru)
  
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










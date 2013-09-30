\name{Rrank}
\alias{Rrank}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Find rank of upper triangular matrix}
\description{ 

Finds rank of upper triangular matrix R, by estimating condition
number of upper \code{rank} by \code{rank} block, and reducing \code{rank} 
until this is acceptably low. Assumes R has been computed by a method that uses 
pivoting, usually pivoted QR or Choleski.
}
\usage{
Rrank(R,tol=.Machine$double.eps^.9)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{R}{An upper triangular matrix, obtained by pivoted QR or pivoted Choleski.}
\item{tol}{the tolerance to use for judging rank.}

}

\details{ The method is based on Cline et al. (1979) as described in Golub and van Loan (1996).
}
 

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\references{
Cline, A.K., C.B. Moler, G.W. Stewart and J.H. Wilkinson (1979) 
An estimate for the condition number of a matrix. 
SIAM J. Num. Anal. 16, 368-375

Golub, G.H, and C.F. van Loan (1996) 
Matrix Computations 3rd ed. 
Johns Hopkins University Press, Baltimore.
}

\examples{
  set.seed(0)
  n <- 10;p <- 5
  X <- matrix(runif(n*(p-1)),n,p)
  qrx <- qr(X,LAPACK=TRUE)
  Rrank(qr.R(qrx))
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










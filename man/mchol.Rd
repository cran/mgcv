\name{mchol}
\alias{mchol}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Sparse \code{chol} function}
\description{A wrapper for the \code{Cholesky} function from the \code{Matrix} package that produces the sparse pivoted Cholesky factor of a positive definite matrix, and returns the pivot vector for this as an attribute. The \code{Matrix} package \code{chol} function no longer retuns the pivots. 
}

\usage{
mchol(A)
}
\arguments{
\item{A}{A sparse matrix suitable for passing to \code{Cholesky} function from the \code{Matrix} package.}
}

\value{
If \code{A} is positive definite then the upper triangular Cholesky factor matrix, with \code{"pivot"} attribute and \code{"rank"} attribute which is the diminsion of \code{A}. Otherwise \code{-1} with \code{"rank"} attribute \code{-1}.
}

\details{The \code{Matrix} version of \code{chol} performs sparse Cholesky decomposition with sparsity maintaining pivoting, but no longer returns the pivot vector, rendering the returned factor useless for many purposes. This wrapper function simply fixes this. It also ensures that there are no numerical zeroes below the leading diagonal (not the default behaviour of \code{expand1}, which can put some numerical zeroes there, in place of structural ones, at least at version 1.6-5). Calls \code{chol(A,pivot=TRUE)} if \code{A} inherits from class \code{matrix}. 
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood  \email{simon.wood@r-project.org}
}


\examples{
library(mgcv)
## A sparse +ve def matrix
u <- sample(1:100,10)*.001
x <- i <- j <- 1:10
ii <- sample(1:10,10,replace=TRUE);
jj <- sample(1:10,10,replace=TRUE)
x <- c(x,u,u);
i <- c(i,ii,jj)
j <- c(j,jj,ii)
A <- Matrix::sparseMatrix(i=i,j=j,x=x)
R <- mchol(A)
piv <- attr(R,"pivot")
range(crossprod(R)-A[piv,piv])

i <- sample(1:5,10,replace=TRUE)
j <- sample(1:10,10,replace=TRUE)
u <- sample(1:100,10)*.001
A <- crossprod(Matrix::sparseMatrix(i=i,j=j,x=u))
mchol(A)

}
\keyword{models} \keyword{regression}%-- one or more ..



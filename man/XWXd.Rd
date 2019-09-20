\name{XWXd}
\alias{XWXd}
\alias{XWyd}
\alias{Xbd}
\alias{diagXVXd}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Internal functions for discretized model matrix handling}
\description{Routines for computing with discretized model matrices as described in Wood et al. (2017) and Li and Wood (2019).
}
\usage{
XWXd(X,w,k,ks,ts,dt,v,qc,nthreads=1,drop=NULL,ar.stop=-1,ar.row=-1,ar.w=-1,
     lt=NULL,rt=NULL)
XWyd(X,w,y,k,ks,ts,dt,v,qc,drop=NULL,ar.stop=-1,ar.row=-1,ar.w=-1,lt=NULL)
Xbd(X,beta,k,ks,ts,dt,v,qc,drop=NULL,lt=NULL)
diagXVXd(X,V,k,ks,ts,dt,v,qc,drop=NULL,nthreads=1,lt=NULL,rt=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{X}{A list of the matrices containing the unique rows of model matrices for terms of a full model matrix, or the model matrices of the terms margins.
 if term subsetting arguments \code{lt} and \code{rt} are non-NULL then this requires an \code{"lpip"} attribute: see details.}
 \item{w}{An n-vector of weights}
 \item{y}{n-vector of data.}
 \item{beta}{coefficient vector.}
 \item{k}{A matrix whose columns are index n-vectors each selecting the rows of an X[[i]] required to create the full matrix.}
\item{ks}{The ith term has index vectors \code{ks[i,1]:(ks[i,2]-1)}. The corresponing full model matrices are summed over.} 
\item{ts}{The element of \code{X} at which each model term starts.}
\item{dt}{How many elements of \code{X} contribute to each term.}
\item{v}{\code{v[[i]]} is Householder vector for ith term, if \code{qc[i]>0}.}
\item{qc}{if \code{qc[i]>0} then term has a constraint.}
\item{nthreads}{number of threads to use}
\item{drop}{list of columns of model matrix/parameters to drop}
\item{ar.stop}{Negative to ignore. Otherwise sum rows \code{(ar.stop[i-1]+1):ar.stop[i]} of the rows selected by \code{ar.row} and weighted by \code{ar.w} to get ith row of model matrix to use.}
\item{ar.row}{extract these rows...}
\item{ar.w}{weight by these weights, and sum up according to \code{ar.stop}. Used to implement AR models.}
\item{lt}{use only columns of X corresponding to these model matrix terms (for left hand \code{X} in \code{XWXd}).}
\item{rt}{as \code{lt} for right hand \code{X}.}
\item{V}{Coefficient covariance matrix.}
}

\details{These functions are really intended to be internal, but are exported so that they can be used in the initialization code of families without problem. They are primarily used by \code{\link{bam}} to implement the methods given in the references. \code{XWXd} produces \eqn{X^TWX}{X'WX}, \code{XWy} produces \eqn{X^TWy}{X'Wy}, \code{Xbd} produces \eqn{X\beta}{Xb} and \eqn{diagXVXd} produces the diagonal of \eqn{XVX^T}{XVX'}.

The \code{"lpip"} attribute of \code{X} is a list of the coefficient indices for each term. Required if subsetting via \code{lt} and \code{rt}.
}
 

\references{
Wood, S.N., Li, Z., Shaddick, G. & Augustin N.H. (2017) Generalized additive models for gigadata: modelling the UK black smoke network daily data. Journal of the American Statistical Association. 112(519):1199-1210
\url{http://dx.doi.org/10.1080/01621459.2016.1195744}

Li, Z & S.N. Wood (2019) Faster model matrix crossproducts for large generalized linear models with discretized covariates. Statistics and Computing.
\url{https://doi.org/10.1007/s11222-019-09864-2}

}


\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










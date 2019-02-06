\name{sp.vcov}
\alias{sp.vcov}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract smoothing parameter estimator covariance matrix from (RE)ML GAM fit}
\description{ Extracts the estimated covariance matrix for the log smoothing parameter
estimates from a (RE)ML estimated \code{gam} object, provided the fit was with a method 
that evaluated the required Hessian.
}
\usage{
sp.vcov(x,edge.correct=TRUE,reg=1e-3)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x}{ a fitted model object of class \code{gam} as produced by \code{gam()}.}
\item{edge.correct}{ if the model was fitted with \code{edge.correct=TRUE} (see \code{\link{gam.control}}), then thereturned covariance matrix will be for the edge corrected log smoothing parameters.}
\item{reg}{regularizer for Hessian - default is equivalent to prior variance of 1000 on log smoothing parameters.}
}
\details{ Just extracts the inverse of the hessian matrix of the negative (restricted) log likelihood w.r.t
the log smoothing parameters, if this has been obtained as part of fitting. 
}

\value{ A matrix corresponding to the estimated covariance matrix of the log smoothing parameter estimators,
if this can be extracted, otherwise \code{NULL}. If the scale parameter has been (RE)ML estimated (i.e. if the method was \code{"ML"} or \code{"REML"} and the scale parameter was unknown) then the 
last row and column relate to the log scale parameter. If \code{edge.correct=TRUE} and this was used in fitting then the edge corrected smoothing parameters are in attribute \code{lsp} of the returned matrix.  
}

\author{Simon N. Wood \email{simon.wood@r-project.org}
}

\references{Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models (with discussion).
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}
}


\seealso{ \code{\link{gam}}, \code{\link{gam.vcomp}}}


\examples{ 
require(mgcv)
n <- 100
x <- runif(n);z <- runif(n)
y <- sin(x*2*pi) + rnorm(n)*.2
mod <- gam(y~s(x,bs="cc",k=10)+s(z),knots=list(x=seq(0,1,length=10)),
           method="REML")
sp.vcov(mod)
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







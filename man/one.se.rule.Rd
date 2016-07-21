\name{one.se.rule}
\alias{one.se.rule}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{The one standard error rule for smoother models}
\description{ The `one standard error rule' (see e.g. Hastie, Tibshirani and Friedman, 2009) is a way of producing smoother models than those directly estimated by automatic smoothing parameter selection methods. In the single smoothing parameter case, we select the largest smoothing parameter within one standard error of the optimum of the smoothing parameter selection criterion. This approach can be generalized to multiple smoothing parameters estimated by REML or ML.}

\details{Under REML or ML smoothing parameter selection an asyptotic distributional approximation is available for the log smoothing parameters. Let \eqn{\rho}{r} denote the log smoothing parameters that we want to increase to obtain a smoother model. The large sample distribution of the estimator of \eqn{\rho}{r} is \eqn{N(\rho,V)}{N(r,V)} where \eqn{V}{V} is the matrix returned by \code{\link{sp.vcov}}. Drop any elements of \eqn{\rho}{r} that are already at `effective infinity', along with the corresponding rows and columns of \eqn{V}{V}. The standard errors of the log smoothing parameters can be obtained from the leading diagonal of \eqn{V}{V}. Let the vector of these be \eqn{d}{d}. Now suppose that we want to increase the estimated log smoothing parameters by an amount \eqn{\alpha d}{a*d}. We choose \eqn{\alpha}{a} so that \eqn{\alpha d^T V^{-1}d = \sqrt{2p}}{a d'V^{-1}d = (2p)^0.5}, where p is the dimension of d and 2p the variance of a chi-squared r.v. with p degrees of freedom.

The idea is that we increase the log smoothing parameters in proportion to their standard deviation, until the RE/ML is increased by 1 standard deviation according to its asypmtotic distribution. }


\author{Simon N. Wood \email{simon.wood@r-project.org}
}

\references{Hastie, T, R. Tibshirani and J. Friedman (2009) The Elements of Statistical Learning 2nd ed. Springer.}


\seealso{ \code{\link{gam}}}


\examples{ 
require(mgcv)
set.seed(2) ## simulate some data...
dat <- gamSim(1,n=400,dist="normal",scale=2)
b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat,method="REML")
b
## only the first 3 smoothing parameters are candidates for
## increasing here...
V <- sp.vcov(b)[1:3,1:3] ## the approx cov matrix of sps
d <- diag(V)^.5          ## sp se.
## compute the log smoothing parameter step...
d <- sqrt(2*length(d))/d%*%solve(V,d)*d
sp <- b$sp ## extract original sp estimates
sp[1:3] <- sp[1:3]*exp(d) ## apply the step
## refit with the increased smoothing parameters...
b1 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat,method="REML",sp=sp)
b;b1 ## compare fits
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







\name{ziP}
\alias{ziP}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM zero-inflated Poisson regression family}
\description{Family for use with \code{\link{gam}}, implementing regression for zero inflated Poisson data
when the complinetart log log of the zero probability is linearly dependent on the log of the Poisson parameter. Use with care, 
noting that simply having many zero response observations is not an indication of zero inflation: the question is whether you have too many zeroes given the specified model.  
}

\usage{
ziP(theta = NULL, link = "identity")
}
\arguments{
\item{theta}{the 2 parameters controlling the slope and intercept of the 
linear transform of the mean controlling the zero inflation rate. If supplied (and second element is positive) 
then treated as fixed parameters (\eqn{\theta_1}{theta_1} and \eqn{\exp(\theta_2)}{exp(theta_2)}), otherwise estimated. If supplied and second element is negative then treated as starting values (with sign of second element reversed).}

\item{link}{The link function: only the \code{"identity"} is currently supported.}

}
\value{
An object of class \code{extended.family}.
}

\details{The probability of a zero count is given by \eqn{1-p}{1- p}, whereas the probability of
count \eqn{y>0}{y>0} is given by the truncated Poisson probability function \eqn{p\mu^y/((\exp(\mu)-1)y!)}{(pmu^y/((exp(mu)-1)y!)}. The linear predictor 
gives \eqn{\log \mu}{log(mu)}, while \eqn{\eta = \log(-\log(1-p)) }{eta=log(-log(1-p))} and \eqn{\eta = \theta_1 + \exp(\theta_2) \log \mu }{eta = theta_1 + exp(theta_2) log(mu)}. The \code{theta} parameters are estimated alongside the smoothing parameters.

Note that the fitted values for this model are the log of the Poisson parameter. Use the \code{predict} function with \code{type=="response"} to get the predicted expected response.

}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}

\section{WARNINGS }{

Zero inflated models are often over-used. Having lots of zeroes in the data does not in itself imply zero inflation. Having too many zeroes *given the model mean* may imply zero inflation. 
}


\seealso{\code{\link{ziplss}}}

\examples{

rzip <- function(gamma,theta= c(-2,.3)) {
## generate zero inflated Poisson random variables, where 
## lambda = exp(gamma), eta = theta[1] + exp(theta[2])*gamma
## and 1-p = exp(-exp(eta)).
   y <- gamma; n <- length(y)
   lambda <- exp(gamma)
   eta <- theta[1] + exp(theta[2])*gamma
   p <- 1- exp(-exp(eta))
   ind <- p > runif(n)
   y[!ind] <- 0
   np <- sum(ind)
   ## generate from zero truncated Poisson, given presence...
   y[ind] <- qpois(runif(np,dpois(0,lambda[ind]),1),lambda[ind])
   y
} 

library(mgcv)
## Simulate some beta data...
set.seed(1);n<-400
dat <- gamSim(1,n=n)
dat$y <- rzip(dat$f/4-1)

b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ziP(),data=dat)
b$outer.info ## check convergence!!
b
plot(b,pages=1)
## plot deviance residuals against expected response
plot(predict(b,type="response"),residuals(b))

}
\keyword{models} \keyword{regression}%-- one or more ..



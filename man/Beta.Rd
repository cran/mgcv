\name{betar}
\alias{betar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM beta regression family}
\description{Family for use with \code{\link{gam}} or \code{\link{bam}}, implementing regression for beta distributed data on (0,1).
A linear predictor controls the mean, \eqn{\mu}{mu} of the beta distribution, while the variance is then
\eqn{\mu(1-\mu)/(1+\phi)}{mu(1-mu)/(1+phi)}, with parameter \eqn{\phi}{phi} being estimated during 
fitting, alongside the smoothing parameters.
}

\usage{
betar(theta = NULL, link = "logit",eps=.Machine$double.eps*100)
}
\arguments{
\item{theta}{the extra parameter (\eqn{\phi}{phi} above). }

\item{link}{The link function: one of \code{"logit"}, \code{"probit"}, \code{"cloglog"} and \code{"cauchit"}.}

\item{eps}{the response variable will be truncated to the interval \code{[eps,1-eps]} if there are values outside this range.
This truncation is not entirely benign, but too small a value of \code{eps} will cause stability problems if there are 
zeroes or ones in the response.}
}
\value{
An object of class \code{extended.family}.
}

\details{These models are useful for proportions data which can not be modelled as binomial. Note the assumption that data are in 
(0,1), despite the fact that for some parameter values 0 and 1 are perfectly legitimate observations. The restriction is needed to 
keep the log likelihood bounded for all parameter values. Any data exactly at 0 or 1 are reset to be just above 0 or just below 1 using the \code{eps} argument (in fact any observation \code{<eps} is reset to \code{eps} and any observation \code{>1-eps} is reset to \code{1-eps}). Note the effect of this resetting. If \eqn{\mu\phi>1}{mu*phi>1} then impossible 0s are replaced with highly improbable \code{eps} values. If the inequality is reversed then 0s with infinite probability density are replaced with \code{eps} values having high finite probability density. The equivalent condition for 1s is \eqn{(1-\mu)\phi>1}{(1-mu)*phi>1}. Clearly all types of resetting are somewhat unsatisfactory, and care is needed if data contain 0s or 1s (often it makes sense to manually reset the 0s and 1s in a manner that somehow reflects the sampling setup).  
}

%- maybe also `usage' for other objects documented here.

\author{ Natalya Pya (nat.pya@gmail.com) and Simon Wood (s.wood@r-project.org)
}

\section{WARNINGS}{
Do read the details section if your data contain 0s and or 1s. 
}

\examples{
library(mgcv)
## Simulate some beta data...
set.seed(3);n<-400
dat <- gamSim(1,n=n)
mu <- binomial()$linkinv(dat$f/4-2)
phi <- .5
a <- mu*phi;b <- phi - a;
dat$y <- rbeta(n,a,b) 

bm <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=betar(link="logit"),data=dat)

bm
plot(bm,pages=1)
}
\keyword{models} \keyword{regression}%-- one or more ..



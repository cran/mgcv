\name{betar}
\alias{betar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM beta regression family}
\description{Family for use with \code{\link{gam}}, implementing regression for beta distributed data on (0,1).
A linear predictor controls the mean, \eqn{\mu}{mu} of the beta distribution, while the variance is then
\eqn{\mu(1-\mu)/(1+\phi)}{mu(1-mu)/(1+phi)}, with parameter \eqn{\phi}{phi} being estimated during 
fitting, alongside the smoothing parameters.
}

\usage{
betar(theta = NULL, link = "logit")
}
\arguments{
\item{theta}{the extra parameter (\eqn{\phi}{phi} above). }

\item{link}{The link function: one of \code{"logit"}, \code{"probit"}, \code{"cloglog"} and \code{"cauchit"}.}

}
\value{
An object of class \code{extended.family}.
}

\details{These models are useful for proportions data which can not be modelled as binomial. Note the assumption that data are in 
(0,1), despite the fact that for some parameter values 0 and 1 are perfectly legitimate observations. The restriction is needed to 
keep the log likelihood bounded for all parameter values. Any data exactly at 0 or 1 
are reset to be just above 0 or just below 1. 
}

%- maybe also `usage' for other objects documented here.

\author{ Natalya Pya (nyp20@bath.ac.uk) and Simon Wood (s.wood@r-project.org)
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



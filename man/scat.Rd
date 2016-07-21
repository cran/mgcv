\name{scat}
\alias{scat}
\alias{t.scaled}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM scaled t family for heavy tailed data}
\description{Family for use with \code{\link{gam}}, implementing regression for the heavy tailed response
variables, y, using a scaled t model. The idea is that \eqn{(y-\mu)/\sigma \sim t_\nu }{(y - mu)/sig ~ t_nu} where 
\eqn{mu}{mu} is determined by a linear predictor, while \eqn{\sigma}{sig} and \eqn{\nu}{nu} are parameters 
to be estimated alongside the smoothing parameters.
}

\usage{
scat(theta = NULL, link = "identity")
}
\arguments{
\item{theta}{the parameters to be estimated \eqn{\nu = 2 + \exp(\theta_1)}{nu = 2 + exp(theta_1) } and 
\eqn{\sigma = \exp(\theta_2)}{sig = exp(theta_2)}. If supplied and positive, then taken to be fixed values of 
\eqn{\nu}{nu} and \eqn{\sigma}{sig}. If any negative, then absolute values taken as starting values. }

\item{link}{The link function: one of \code{"identity"}, \code{"log"} or \code{"inverse"}.}

}
\value{
An object of class \code{extended.family}.
}

\details{Useful in place of Gaussian, when data are heavy tailed.
}

%- maybe also `usage' for other objects documented here.

\author{ Natalya Pya (nyp20@bath.ac.uk)
}
\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association.
\url{http://arxiv.org/abs/1511.03864}
}


\examples{
library(mgcv)
## Simulate some t data...
set.seed(3);n<-400
dat <- gamSim(1,n=n)
dat$y <- dat$f + rt(n,df=3)*2

b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=scat(link="identity"),data=dat)

b
plot(b,pages=1)

}
\keyword{models} \keyword{regression}%-- one or more ..



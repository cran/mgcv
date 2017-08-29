\name{scat}
\alias{scat}
\alias{t.scaled}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM scaled t family for heavy tailed data}
\description{Family for use with \code{\link{gam}} or \code{\link{bam}}, implementing regression for the heavy tailed response
variables, y, using a scaled t model. The idea is that \eqn{(y-\mu)/\sigma \sim t_\nu }{(y - mu)/sig ~ t_nu} where 
\eqn{mu}{mu} is determined by a linear predictor, while \eqn{\sigma}{sig} and \eqn{\nu}{nu} are parameters 
to be estimated alongside the smoothing parameters.
}

\usage{
scat(theta = NULL, link = "identity",min.df=3)
}
\arguments{
\item{theta}{the parameters to be estimated \eqn{\nu = b + \exp(\theta_1)}{nu = b + exp(theta_1) } (where `b' is \code{min.df}) and 
\eqn{\sigma = \exp(\theta_2)}{sig = exp(theta_2)}. If supplied and both positive, then taken to be fixed values of 
\eqn{\nu}{nu} and \eqn{\sigma}{sig}. If any negative, then absolute values taken as starting values. }

\item{link}{The link function: one of \code{"identity"}, \code{"log"} or \code{"inverse"}.}

\item{min.df}{minimum degrees of freedom. Should not be set to 2 or less as this implies infinite response variance.}
}
\value{
An object of class \code{extended.family}.
}

\details{Useful in place of Gaussian, when data are heavy tailed. \code{min.df} can be modified, but lower values can occasionally
lead to convergence problems in smoothing parameter estimation. In any case \code{min.df} should be >2, since only then does a t
random variable have finite variance.   
}

%- maybe also `usage' for other objects documented here.

\author{ Natalya Pya (nat.pya@gmail.com)
}
\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}

}


\examples{
library(mgcv)
## Simulate some t data...
set.seed(3);n<-400
dat <- gamSim(1,n=n)
dat$y <- dat$f + rt(n,df=4)*2

b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=scat(link="identity"),data=dat)

b
plot(b,pages=1)

}
\keyword{models} \keyword{regression}%-- one or more ..



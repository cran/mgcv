\name{Tweedie}
\alias{Tweedie}
\alias{tw}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM Tweedie families}
\description{Tweedie families, designed for use with \code{\link{gam}} from the \code{mgcv} library.
Restricted to variance function powers between 1 and 2. A useful alternative to \code{\link{quasi}} when a
full likelihood is desirable. \code{Tweedie} is for use with fixed \code{p}. \code{tw} is for use when \code{p}
is to be estimated during fitting. For fixed \code{p} between 1 and 2 the Tweedie is an exponential family 
distribution with variance given by the mean to the power \code{p}.

\code{tw} is only useable with \code{\link{gam}} and \code{\link{bam}} but not \code{gamm}. \code{Tweedie} works with all three.
}

\usage{
Tweedie(p=1, link = power(0))
tw(theta = NULL, link = "log",a=1.01,b=1.99)
}
\arguments{
\item{p}{the variance of an observation is proportional to its mean to the power \code{p}. \code{p} must
be greater than 1 and less than or equal to 2. 1 would be Poisson, 2 is gamma. }

\item{link}{The link function: one of \code{"log"}, \code{"identity"}, \code{"inverse"}, \code{"sqrt"}, or a 
\code{\link{power}} link (\code{Tweedie} only).}

\item{theta}{Related to the Tweedie power parameter by \eqn{p=\exp(a+b \exp(\theta))/(1+\exp(\theta))}{p=exp(a+b*exp(theta))/(1+exp(theta))}. If this is supplied as a positive value then it is taken as the fixed value for \code{p}.
If it is a negative values then its absolute value is taken as the initial value for \code{p}.}

\item{a}{lower limit on \code{p} for optimization.}

\item{b}{upper limit on \code{p} for optimization.}

}
\value{
 For \code{Tweedie}, an object inheriting from class \code{family}, with additional elements
\item{dvar}{the function giving the first derivative of the variance function w.r.t. \code{mu}.}
\item{d2var}{the function giving the second derivative of the variance function w.r.t. \code{mu}.}
\item{ls}{A function returning a 3 element array: the saturated log likelihood followed by its first 2 derivatives
w.r.t. the scale parameter.}

For \code{tw}, an object of class \code{extended.family}.
}

\details{ A Tweedie random variable with 1<p<2 is a sum of \code{N} gamma random variables 
where \code{N} has a Poisson distribution. The p=1 case is a generalization of a Poisson distribution and is a discrete 
distribution supported on integer multiples of the scale parameter. For 1<p<2 the distribution is supported on the 
positive reals with a point mass at zero. p=2 is a gamma distribution. As p gets very close to 1 the continuous 
distribution begins to converge on the discretely supported limit at p=1, and is therefore highly multimodal. 
See \code{\link{ldTweedie}} for more on this behaviour.

\code{Tweedie} is based partly on the \code{\link{poisson}} family, and partly on \code{tweedie} from the 
\code{statmod} package. It includes extra components to work with all \code{mgcv} GAM fitting methods as well as 
an \code{aic} function. 

The Tweedie density involves a normalizing constant with no closed form, so this is evaluated using the series 
evaluation method of Dunn and Smyth (2005), with extensions to also compute the derivatives w.r.t. \code{p} and the scale parameter. 
Without restricting \code{p} to (1,2) the calculation of Tweedie densities is more difficult, and there does not 
currently seem to be an implementation which offers any benefit over \code{\link{quasi}}. If you need  this 
case then the \code{tweedie} package is the place to start.
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}.
}

\references{
Dunn, P.K. and G.K. Smyth (2005) Series evaluation of Tweedie exponential dispersion model densities. 
Statistics and Computing 15:267-280

Tweedie, M. C. K. (1984). An index which distinguishes between
     some important exponential families. Statistics: Applications and
     New Directions. Proceedings of the Indian Statistical Institute
     Golden Jubilee International Conference (Eds. J. K. Ghosh and J.
     Roy), pp. 579-604. Calcutta: Indian Statistical Institute.

Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}

}

\seealso{\code{\link{ldTweedie}}, \code{\link{rTweedie}}}

\examples{
library(mgcv)
set.seed(3)
n<-400
## Simulate data...
dat <- gamSim(1,n=n,dist="poisson",scale=.2)
dat$y <- rTweedie(exp(dat$f),p=1.3,phi=.5) ## Tweedie response

## Fit a fixed p Tweedie, with wrong link ...
b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=Tweedie(1.25,power(.1)),
         data=dat)
plot(b,pages=1)
print(b)

## Same by approximate REML...
b1 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=Tweedie(1.25,power(.1)),
          data=dat,method="REML")
plot(b1,pages=1)
print(b1)

## estimate p as part of fitting

b2 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=tw(),
          data=dat,method="REML")
plot(b2,pages=1)
print(b2)

rm(dat)
}
\keyword{models} \keyword{regression}%-- one or more ..



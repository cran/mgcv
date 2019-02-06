\name{twlss}
\alias{twlss}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Tweedie location scale family}
\description{Tweedie family in which the mean, power and scale parameters can all depend on smooth linear predictors. Restricted to estimation via the extended Fellner Schall method of Wood and Fasiolo (2017). Only usable with \code{\link{gam}}. Tweedie distributions are exponential family with variance given by \eqn{\phi \mu^p}{s*m^p} where \eqn{\phi}{s} is a scale parameter, \eqn{p}{p} a parameter (here between 1 and 2) and \eqn{\mu}{m} is the mean. 
}

\usage{
twlss(link=list("log","identity","identity"),a=1.01,b=1.99)
}
\arguments{

\item{link}{The link function list: currently no choise.}

\item{a}{lower limit on the power parameter relating variance to mean.}

\item{b}{upper limit on power parameter.}

}
\value{
An object inheriting from class \code{general.family}.
}

\details{ A Tweedie random variable with 1<p<2 is a sum of \code{N} gamma random variables 
where \code{N} has a Poisson distribution. The p=1 case is a generalization of a Poisson distribution and is a discrete 
distribution supported on integer multiples of the scale parameter. For 1<p<2 the distribution is supported on the 
positive reals with a point mass at zero. p=2 is a gamma distribution. As p gets very close to 1 the continuous 
distribution begins to converge on the discretely supported limit at p=1, and is therefore highly multimodal. 
See \code{\link{ldTweedie}} for more on this behaviour.

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

Wood, S.N. and Fasiolo, M., (2017). A generalized Fellner-Schall method for smoothing
parameter optimization with application to Tweedie location, scale and shape models. Biometrics, 73(4), pp.1071-1081.
\url{https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12666}

Wood, S.N., N. Pya and B. Saefken (2016). Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}

}

\seealso{\code{\link{Tweedie}}, \code{\link{ldTweedie}}, \code{\link{rTweedie}}}

\examples{
library(mgcv)
set.seed(3)
n<-400
## Simulate data...
dat <- gamSim(1,n=n,dist="poisson",scale=.2)
dat$y <- rTweedie(exp(dat$f),p=1.3,phi=.5) ## Tweedie response

## Fit a fixed p Tweedie, with wrong link ...
b <- gam(list(y~s(x0)+s(x1)+s(x2)+s(x3),~1,~1),family=twlss(),
         data=dat)
plot(b,pages=1)
print(b)

rm(dat)
}
\keyword{models} \keyword{regression}%-- one or more ..



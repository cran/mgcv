\name{gam.neg.bin}
\alias{gam.neg.bin}
\alias{mgcv.find.theta}
\alias{mgcv.get.scale}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAMs with the negative binomial distribution}
\description{The \code{gam} modelling function is designed to be able to use 
the \code{negative.binomial} and \code{neg.bin} families from the MASS library, 
with or without a known \eqn{\theta}{theta} parameter. A value for \code{theta} 
must always be passed to these families, but if \eqn{\theta}{theta} is to be
estimated then the passed value is treated as a starting value for estimation. 


If the \code{scale} argument passed to \code{gam} is positive, then it is used 
as the scale parameter \code{theta} is treated as a fixed known parameter and 
any smoothing parameters are chosen by UBRE. If \code{scale} is not positive then
\eqn{\theta}{theta} is estimated. The method of estimation is to choose \eqn{\hat \theta}{theta} 
so that the GCV (Pearson) estimate of the scale parameter is one (since the scale parameter 
is one for the negative binomial). 

\eqn{\theta}{theta} estimation is nested within the IRLS loop used for GAM fitting. After
each call to fit an iteratively weighted additive model to the IRLS pseudodata, the \eqn{\theta}{theta}
estimate is updated. This is done by conditioning on all components of the current GCV/Pearson 
estimator of the scale parameter except \eqn{\theta}{theta} and then searching for the 
\eqn{\hat \theta}{theta} which equates this conditional  estimator to one. The search is 
a simple bisection search after an initial crude line search to bracket one. The search will 
terminate at the upper boundary of the search region is a Poisson fit would have yielded an estimated 
scale parameter <1. Search limits can be set in \code{gam.control}.

Note that 
\code{neg.bin} only allows a log link, while \code{negative.binomial} also allows \code{"sqrt"} and
\code{"identity"}. In addition the \code{negative.binomial} family results in a more 
informative \code{gam} summary.  

The negative binomial families can not yet be used with `outer' estimation of
smoothing parameters (see \code{\link{gam.method}}).

}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
library(MASS) # required for negative binomial families
set.seed(3)
n<-400
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
pi <- asin(1) * 2
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f + 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10 - 1.396
g<-exp(f/5)
# negative binomial data  
y<-rnbinom(g,size=3,mu=g)
# unknown theta ...
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=negative.binomial(1))
plot(b,pages=1)
print(b)
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=neg.bin(1)) # unknown theta
plot(b,pages=1)
print(b)
# known theta example ...
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=negative.binomial(3),scale=1)
plot(b,pages=1)
print(b)
# Now use "sqrt" link available in negative.binomial (but not neg.bin)
set.seed(1)
f<-f-min(f);g<-f^2
y<-rnbinom(g,size=3,mu=g)
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=negative.binomial(1,link="sqrt")) 
plot(b,pages=1)
print(b)
}
\keyword{models} \keyword{regression}%-- one or more ..



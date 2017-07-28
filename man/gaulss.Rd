\name{gaulss}
\alias{gaulss}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Gaussian location-scale model family}
\description{The \code{gaulss} family implements Gaussian location scale additive models in which 
the mean and the logb of the standard deviation (see details) can depend on additive smooth predictors. Useable 
only with \code{\link{gam}}, the linear predictors are specified via a list of formulae.
}

\usage{
gaulss(link=list("identity","logb"),b=0.01)
}
\arguments{
\item{link}{two item list specifying the link for the mean and the standard deviation. See details.}
\item{b}{The minumum standard deviation, for the \code{"logb"} link.}
}
\value{
 An object inheriting from class \code{general.family}.
}

\details{Used with \code{\link{gam}} to fit Gaussian location - scale models. \code{gam} is called with 
a list containing 2 formulae, the first specifies the response on the left hand side and the structure of the linear predictor for the mean on the right hand side. The second is one sided, specifying the linear predictor for the standard deviation on the right hand side. 

Link functions \code{"identity"}, \code{"inverse"}, \code{"log"} and \code{"sqrt"} are available for the mean. For the standard deviation only the \code{"logb"} link is implemented: \eqn{\eta = \log(\sigma - b)}{eta = log(sigma-b)} and \eqn{\sigma = b + \exp(\eta)}{sigma = b + exp(eta)}. This link is designed to avoid singularities in the likelihood caused by the standard deviation tending to zero. Note that internally the family is parameterized in terms of the \eqn{\tau=\sigma^{-1}}{tau=1/sigma}, i.e. the standard deviation of the precision, so the link and inverse link are coded to reflect this, however the reltaionships between the linear predictor and the standard deviation are as given above. 

The fitted values for this family will be a two column matrix. The first column is the mean, and the second column is the inverse of the standard deviation. Predictions using \code{\link{predict.gam}} will also produce 2 column matrices for \code{type} \code{"link"} and \code{"response"}. The second column when \code{type="response"} is again on the reciprocal standard deviation scale (i.e. the square root precision scale). The second column when \code{type="link"} is  \eqn{\log(\sigma - b)}{log(sigma-b)}. Also \code{\link{plot.gam}} will plot smooths relating to \eqn{\sigma}{sigma} on the \eqn{\log(\sigma - b)}{log(sigma-b)} scale (so high values correspond to high standard deviation and low values to low standard deviation). Similarly the smoothing penalties are applied on the (log) standard deviation scale, not the log precision scale.

The null deviance reported for this family is the sum of squares of the difference between the response and the mean response divided by the standard deviation of the response according to the model. The deviance is the sum of squares of residuals divided by model standard deviations. 

}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}
}




\examples{
library(mgcv);library(MASS)
b <- gam(list(accel~s(times,k=20,bs="ad"),~s(times)),
            data=mcycle,family=gaulss())
summary(b) 
plot(b,pages=1,scale=0)
}
\keyword{models} \keyword{regression}%-- one or more ..



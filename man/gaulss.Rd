\name{gaulss}
\alias{gaulss}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Gaussian location-scale model family}
\description{The \code{gaulss} family implements Gaussian location scale additive models in which 
the mean and the inverse of the standard deviation can depend on additive smooth predictors. Useable 
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

Link functions \code{"identity"}, \code{"inverse"}, \code{"log"} and \code{"sqrt"} are available for the mean. For the standard deviation only the \code{"logb"} link is implemented: \eqn{\eta = \log(\sigma - b)}{eta = log(sigma-b)} and \eqn{\sigma = b + \exp(\eta)}{sigma = b + exp(eta)}. This link is designed to avoid singularities in the likelihood caused by the standard deviation tending to zero.

The fitted values for this family will be a two column matrix. The first column is the mean, and the second column is the inverse of the 
standard deviation. Predictions using \code{\link{predict.gam}} will also produce 2 column matrices for \code{type} 
\code{"link"} and \code{"response"}.

The null deviance reported for this family is the sum of squares of the difference between the response and the mean response divided by the standard deviation of the response according to the model. The deviance is the sum of squares of residuals divided by model standard deviations. 

}




\examples{
library(mgcv);library(MASS)
b <- gam(list(accel~s(times,k=20,bs="ad"),~s(times)),
            data=mcycle,family=gaulss())
summary(b) 
plot(b,pages=1,scale=0)
}
\keyword{models} \keyword{regression}%-- one or more ..



\name{gammals}
\alias{gammals}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Gamma location-scale model family}
\description{The \code{gammals} family implements gamma location scale additive models in which 
the log of the mean, \eqn{\mu}{m}, and the log of the scale parameter, \eqn{\phi}{s} (see details) can depend on additive smooth predictors.The parameterization is of the usual GLM type where the variance of the response is given by \eqn{\phi\mu^2}{sm^2}. Useable only with \code{\link{gam}}, the linear predictors are specified via a list of formulae.
}

\usage{
gammals(link=list("identity","log"),b=-7)
}
\arguments{
\item{link}{two item list specifying the link for the mean and the standard deviation. See details for meaning which may not be intuitive.}
\item{b}{The minumum log scale parameter.}
}
\value{
 An object inheriting from class \code{general.family}.
}

\details{Used with \code{\link{gam}} to fit gamma location - scale models parameterized in terms of the log mean and the log scale parameter (the response variance is the squared mean multiplied by the scale parameter). Note that \code{identity} links mean that the linear predictors give the log mean and log scale directly. By default the \code{log} link for the scale parameter simply forces the log scale parameter to have a lower limit given by argument \code{b}: if \eqn{\eta}{l} is the linear predictor for the log scale parameter, \eqn{\phi}{s}, then \eqn{\log \phi = b + \log(1+e^\eta)}{log(s) = b + log(1+e^l)}. 


\code{gam} is called with 
a list containing 2 formulae, the first specifies the response on the left hand side and the structure of the linear predictor for the log mean on the right hand side. The second is one sided, specifying the linear predictor for the log scale on the right hand side. 


The fitted values for this family will be a two column matrix. The first column is the mean (on original, not log, scale), and the second column is the log scale. Predictions using \code{\link{predict.gam}} will also produce 2 column matrices for \code{type} \code{"link"} and \code{"response"}. The first column is on the original data scale when \code{type="response"} and on the log mean scale of the linear predictor when \code{type="link"}. The second column when \code{type="response"} is again the log scale parameter, but is on the linear predictor when \code{type="link"}.

The null deviance reported for this family computed by setting the fitted values to the mean response, but using the model estimated scale. 

}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\doi{10.1080/01621459.2016.1180986}
}




\examples{
library(mgcv)
## simulate some data
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
f3 <- function(x) 0 * x
n <- 400;set.seed(9)
x0 <- runif(n);x1 <- runif(n);
x2 <- runif(n);x3 <- runif(n);
mu <- exp((f0(x0)+f2(x2))/5)
th <- exp(f1(x1)/2-2)
y <- rgamma(n,shape=1/th,scale=mu*th)

b1 <- gam(list(y~s(x0)+s(x2),~s(x1)+s(x3)),family=gammals)
plot(b1,pages=1)
summary(b1)
gam.check(b1)
plot(mu,fitted(b1)[,1]);abline(0,1,col=2)
plot(log(th),fitted(b1)[,2]);abline(0,1,col=2)

}
\keyword{models} \keyword{regression}%-- one or more ..



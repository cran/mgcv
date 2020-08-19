\name{gumbls}
\alias{gumbls}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Gumbel location-scale model family}
\description{The \code{gumbls} family implements Gumbel location scale additive models in which the location and scale parameters (see details) can depend on additive smooth predictors. Useable only with \code{\link{gam}}, the linear predictors are specified via a list of formulae.
}

\usage{
gumbls(link=list("identity","log"),b=-7)
}
\arguments{
\item{link}{two item list specifying the link for the location \eqn{\mu}{m} and log scale parameter \eqn{\beta}{B}. See details for meaning, which may not be intuitive.}
\item{b}{The minumum log scale parameter.}
}
\value{
 An object inheriting from class \code{general.family}.
}

\details{
Let \eqn{z = (y-\mu) e^{-\beta}}{z = (y - m)exp(-B)}, then the log Gumbel density is \eqn{l = -\beta - z - e^{-z}}{l = -B - z - exp(-z)}. The expected value of a Gumbel r.v. is \eqn{\mu + \gamma e^{\beta}}{m + g exp(B)} where \eqn{\gamma}{g} is Eulers constant (about 0.57721566). The corresponding variance is \eqn{\pi^2 e^{2\beta}/6}{pi^2 exp(2B)/6}.

\code{gumbls} is used with \code{\link{gam}} to fit Gumbel location - scale models parameterized in terms of scale parameter \eqn{\mu}{m} and the log scale parameter \eqn{\beta}{B}. Note that \code{identity} link for the scale parameter means that the corresponding linear predictor gives \eqn{\beta}{B} directly. By default the \code{log} link for the scale parameter simply forces the log scale parameter to have a lower limit given by argument \code{b}: if \eqn{\eta}{l} is the linear predictor for the log scale parameter, \eqn{\beta}{B}, then \eqn{\beta = b + \log(1+e^\eta)}{B = b + log(1+e^l)}.

\code{gam} is called with 
a list containing 2 formulae, the first specifies the response on the left hand side and the structure of the linear predictor for location parameter, \eqn{\mu}{m}, on the right hand side. The second is one sided, specifying the linear predictor for the lg scale, \eqn{\beta}{B}, on the right hand side. 


The fitted values for this family will be a two column matrix. The first column is the mean, and the second column is the log scale parameter, \eqn{\beta}{B}. Predictions using \code{\link{predict.gam}} will also produce 2 column matrices for \code{type} \code{"link"} and \code{"response"}. The first column is on the original data scale when \code{type="response"} and on the log mean scale of the linear predictor when \code{type="link"}. The second column when \code{type="response"} is again the log scale parameter, but is on the linear predictor when \code{type="link"}.



}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{https://dx.doi.org/10.1080/01621459.2016.1180986}
}




\examples{
library(mgcv)
## simulate some data
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
n <- 400;set.seed(9)
x0 <- runif(n);x1 <- runif(n);
x2 <- runif(n);x3 <- runif(n);
mu <- f0(x0)+f1(x1)
beta <- exp(f2(x2)/5)
y <- mu - beta*log(-log(runif(n))) ## Gumbel quantile function

b <- gam(list(y~s(x0)+s(x1),~s(x2)+s(x3)),family=gumbls)
plot(b,pages=1,scale=0)
summary(b)
gam.check(b)

}
\keyword{models} \keyword{regression}%-- one or more ..



\name{ziplss}
\alias{ziplss}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Zero inflated Poisson location-scale model family}
\description{The \code{ziplss} family implements a zero inflated Poisson model in which one linear predictor 
controls the probability of presence  and the other controls the mean given  presence.
Useable only with \code{\link{gam}}, the linear predictors are specified via a list of formulae. 
Should be used with care: simply having a large number of zeroes is not an indication of zero inflation. 

Requires integer count data.
}

\usage{
ziplss(link=list("identity","identity"))
}
\arguments{
\item{link}{two item list specifying the link - currently only identity links are possible, as 
parameterization is directly in terms of log of Poisson response and logit of probability of presence.}
}
\value{
 An object inheriting from class \code{general.family}.
}

\details{Used with \code{\link{gam}} to fit 2 stage zero inflated Poisson models. \code{gam} is called with 
a list containing 2 formulae, the first specifies the response on the left hand side and the structure of the linear predictor for the Poisson parameter on the right hand side. The second is one sided, specifying the linear predictor for the probability of presence on the right hand side. 

The fitted values for this family will be a two column matrix. The first column is the log of the Poisson parameter, 
and the second column is the complimentary log log of probability of presnece.. 
Predictions using \code{\link{predict.gam}} will also produce 2 column matrices for \code{type} 
\code{"link"} and \code{"response"}.

The null deviance computed for this model assumes that a single probability of presence and a single Poisson parameter 
are estimated. 

For data with large areas of covariate space over which the response is zero it may be advisable to use low order penalties to 
avoid problems. For 1D smooths uses e.g. \code{s(x,m=1)} and for isotropic smooths use \code{\link{Duchon.spline}}s in place of thin plaste terms with order 1 penalties, e.g \code{s(x,z,m=c(1,.5))} --- such smooths penalize towards constants, thereby avoiding extreme estimates when the data are uninformative.

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}
}
\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}
}

\section{WARNINGS }{
Zero inflated models are often over-used. Having lots of zeroes in the data does not in itself imply zero inflation. Having too many zeroes *given the model mean* may imply zero inflation. 
}

\examples{
library(mgcv)
## simulate some data...
f0 <- function(x) 2 * sin(pi * x); f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
n <- 500;set.seed(5)
x0 <- runif(n); x1 <- runif(n)
x2 <- runif(n); x3 <- runif(n)

## Simulate probability of potential presence...
eta1 <- f0(x0) + f1(x1) - 3
p <- binomial()$linkinv(eta1) 
y <- as.numeric(runif(n)<p) ## 1 for presence, 0 for absence

## Simulate y given potentially present (not exactly model fitted!)...
ind <- y>0
eta2 <- f2(x2[ind])/3
y[ind] <- rpois(exp(eta2),exp(eta2))

## Fit ZIP model... 
b <- gam(list(y~s(x2)+s(x3),~s(x0)+s(x1)),family=ziplss())
b$outer.info ## check convergence

summary(b) 
plot(b,pages=1)
}
\keyword{models} \keyword{regression}%-- one or more ..



\name{gevlss}
\alias{gevlss}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Extreme Value location-scale model family}
\description{The \code{gevlss} family implements Generalized Extreme Value location scale additive models in which the location, scale and shape parameters depend on additive smooth predictors. Usable 
only with \code{\link{gam}}, the linear predictors are specified via a list of formulae.
}

\usage{
gevlss(link=list("identity","identity","logit"))
}
\arguments{
\item{link}{three item list specifying the link for the location scale and shape parameters. See details.}
}
\value{
 An object inheriting from class \code{general.family}.
}

\details{Used with \code{\link{gam}} to fit Generalized Extreme Value location scale and shape models. \code{gam} is called with a list containing 3 formulae: the first specifies the response on the left hand side and the structure of the linear predictor for the location parameter on the right hand side. The second is one sided, specifying the linear predictor for the log scale parameter on the right hand side. The third is one sided specifying the linear predictor for the shape parameter.

Link functions \code{"identity"} and \code{"log"} are available for the location (mu) parameter. There is no choice of link for the log scale parameter (\eqn{\rho = \log \sigma}{rho = log sigma}). The shape parameter (xi) defaults to a modified logit link restricting its range to (-1,.5), the upper limit is required to ensure existence of an MLE.

The fitted values for this family will be a three column matrix. The first column is the location parameter, the second column is the log scale parameter, the third column is the shape parameter.

This family does not produce a null deviance. Note that the distribution for \eqn{\xi=0}{xi=0} is approximated by setting \eqn{\xi}{xi} to a small number. 

The derivative system code for this family is mostly auto-generated, and the family is still somewhat experimental.  
}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association.
\url{http://arxiv.org/abs/1511.03864}
}




\examples{
library(mgcv)
Fi.gev <- function(z,mu,sigma,xi) {
## GEV inverse cdf.
  xi[abs(xi)<1e-8] <- 1e-8 ## approximate xi=0, by small xi
  x <- mu + ((-log(z))^-xi-1)*sigma/xi
}

## simulate test data...
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
set.seed(1)
n <- 500
x0 <- runif(n);x1 <- runif(n);x2 <- runif(n)
mu <- f2(x2)
rho <- f0(x0)
xi <- (f1(x1)-4)/9
y <- Fi.gev(runif(n),mu,exp(rho),xi)
dat <- data.frame(y,x0,x1,x2);pairs(dat)

## fit model....
b <- gam(list(y~s(x2),~s(x0),~s(x1)),family=gevlss,data=dat)

## plot and look at residuals...
plot(b,pages=1,scale=0)
summary(b)

par(mfrow=c(2,2))
mu <- fitted(b)[,1];rho <- fitted(b)[,2]
xi <- fitted(b)[,3]
## Get the predicted expected response... 
fv <- mu + exp(rho)*(gamma(1-xi)-1)/xi
rsd <- residuals(b)
plot(fv,rsd);qqnorm(rsd)
plot(fv,residuals(b,"pearson"))
plot(fv,residuals(b,"response"))

}
\keyword{models} \keyword{regression}%-- one or more ..



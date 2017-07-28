\name{mvn}
\alias{mvn}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Multivariate normal additive models}
\description{Family for use with \code{\link{gam}} implementing smooth multivariate Gaussian regression. 
The means for each dimension are given by a separate linear predictor, which may contain smooth components. Extra linear predictors may also be specified giving terms which are shared between components (see \code{\link{formula.gam}}). The Choleski factor of the response precision matrix is estimated as part of fitting.
}

\usage{
mvn(d=2)
}
\arguments{
\item{d}{The dimension of the response (>1).}
}

\value{
An object of class \code{general.family}.
}

\details{The response is \code{d} dimensional multivariate normal, where the covariance matrix is estimated, 
and the means for each dimension have sperate linear predictors. Model sepcification is via a list of gam like 
formulae - one for each dimension. See example.

Currently the family ignores any prior weights, and is implemented using first derivative information sufficient for BFGS estimation of smoothing parameters. \code{"response"} residuals give raw residuals, while \code{"deviance"} 
residuals are standardized to be approximately independent standard normal if all is well.
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood  \email{simon.wood@r-project.org}
}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}

}

\seealso{\code{\link{gaussian}}}

\examples{
library(mgcv)
## simulate some data...
V <- matrix(c(2,1,1,2),2,2)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
n <- 300
x0 <- runif(n);x1 <- runif(n);
x2 <- runif(n);x3 <- runif(n)
y <- matrix(0,n,2)
for (i in 1:n) {
  mu <- c(f0(x0[i])+f1(x1[i]),f2(x2[i]))
  y[i,] <- rmvn(1,mu,V)
}
dat <- data.frame(y0=y[,1],y1=y[,2],x0=x0,x1=x1,x2=x2,x3=x3)

## fit model...

b <- gam(list(y0~s(x0)+s(x1),y1~s(x2)+s(x3)),family=mvn(d=2),data=dat)
b
summary(b)
plot(b,pages=1)
solve(crossprod(b$family$data$R)) ## estimated cov matrix

}
\keyword{models} \keyword{regression}%-- one or more ..



\name{ocat}
\alias{ocat}
\alias{ordered.categorical}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM ordered categorical family}
\description{Family for use with \code{\link{gam}} or \code{\link{bam}}, implementing regression for ordered categorical data.
A linear predictor provides the expected value of a latent variable following a logistic distribution. The 
probability of this latent variable lying between certain cut-points provides the probability of the ordered 
categorical variable being of the corresponding category. The cut-points are estimated along side the model 
smoothing parameters (using the same criterion). The observed categories are coded 1, 2, 3, ... up to the 
number of categories. 
}

\usage{
ocat(theta=NULL,link="identity",R=NULL)
}
\arguments{
\item{theta}{cut point parameter vector (dimension \code{R-2}). If supplied and all positive, then taken to be the cut point increments (first cut point is fixed at -1). If any are negative then absolute values are taken as starting values for cutpoint increments. }

\item{link}{The link function: only \code{"identity"} allowed at present (possibly for ever).}

\item{R}{the number of catergories.}

}
\value{
An object of class \code{extended.family}.
}

\details{Such cumulative threshold models are only identifiable up to an intercept, or one of the cut points. 
Rather than remove the intercept, \code{ocat} simply sets the first cut point to -1. Use \code{\link{predict.gam}} with 
\code{type="response"} to get the predicted probabilities in each category.
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}
}


\examples{
library(mgcv)
## Simulate some ordered categorical data...
set.seed(3);n<-400
dat <- gamSim(1,n=n)
dat$f <- dat$f - mean(dat$f)

alpha <- c(-Inf,-1,0,5,Inf)
R <- length(alpha)-1
y <- dat$f
u <- runif(n)
u <- dat$f + log(u/(1-u)) 
for (i in 1:R) {
  y[u > alpha[i]&u <= alpha[i+1]] <- i
}
dat$y <- y

## plot the data...
par(mfrow=c(2,2))
with(dat,plot(x0,y));with(dat,plot(x1,y))
with(dat,plot(x2,y));with(dat,plot(x3,y))

## fit ocat model to data...
b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ocat(R=R),data=dat)
b
plot(b,pages=1)
gam.check(b)
summary(b)
b$family$getTheta(TRUE) ## the estimated cut points

## predict probabilities of being in each category
predict(b,dat[1:2,],type="response",se=TRUE)
}
\keyword{models} \keyword{regression}%-- one or more ..



\name{ziP}
\alias{ziP}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM zero-inflated Poisson regression family}
\description{Family for use with \code{\link{gam}} or \code{\link{bam}}, implementing regression for zero inflated Poisson data
when the complimentary log log of the zero probability is linearly dependent on the log of the Poisson parameter. Use with great care, noting that simply having many zero response observations is not an indication of zero inflation: the question is whether you have too many zeroes given the specified model.

This sort of model is really only appropriate when none of your covariates help to explain the zeroes in your data. If your covariates predict which observations are likely to have zero mean then adding a zero inflated model on top of this is likely to lead to identifiability problems. Identifiability problems may lead to fit failures, or absurd values for the linear predictor or predicted values. 
}

\usage{
ziP(theta = NULL, link = "identity",b=0)
}
\arguments{
\item{theta}{the 2 parameters controlling the slope and intercept of the 
linear transform of the mean controlling the zero inflation rate. If supplied
then treated as fixed parameters (\eqn{\theta_1}{theta_1} and \eqn{\theta_2}{theta_2}), otherwise estimated.}

\item{link}{The link function: only the \code{"identity"} is currently supported.}

\item{b}{a non-negative constant, specifying the minimum dependence of the zero inflation rate on the linear predictor.}
}
\value{
An object of class \code{extended.family}.
}

\details{The probability of a zero count is given by \eqn{1-p}{1- p}, whereas the probability of
count \eqn{y>0}{y>0} is given by the truncated Poisson probability function \eqn{p\mu^y/((\exp(\mu)-1)y!)}{(pmu^y/((exp(mu)-1)y!)}. The linear predictor 
gives \eqn{\log \mu}{log(mu)}, while \eqn{\eta = \log(-\log(1-p)) }{eta=log(-log(1-p))} and \eqn{\eta = \theta_1 + \{b+\exp(\theta_2)\} \log \mu }{eta = theta_1 + (b+exp(theta_2)) log(mu)}. The \code{theta} parameters are estimated alongside the smoothing parameters. Increasing the \code{b} parameter from zero can greatly reduce identifiability problems, particularly when there are very few non-zero data.  

The fitted values for this model are the log of the Poisson parameter. Use the \code{predict} function with \code{type=="response"} to get the predicted expected response. Note that the theta parameters reported in model summaries are \eqn{\theta_1}{theta_1} and \eqn{b + \exp(\theta_2)}{b + exp(theta_2)}.

These models should be subject to very careful checking, especially if fitting has not converged. It is quite easy to set up models with identifiability problems, particularly if the data are not really zero inflated, but simply have many zeroes because the mean is very low in some parts of the covariate space. See example for some obvious checks. Take convergence warnings seriously.

}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}

\section{WARNINGS }{

Zero inflated models are often over-used. Having lots of zeroes in the data does not in itself imply zero inflation. Having too many zeroes *given the model mean* may imply zero inflation. 
}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}

}

\seealso{\code{\link{ziplss}}}

\examples{

rzip <- function(gamma,theta= c(-2,.3)) {
## generate zero inflated Poisson random variables, where 
## lambda = exp(gamma), eta = theta[1] + exp(theta[2])*gamma
## and 1-p = exp(-exp(eta)).
   y <- gamma; n <- length(y)
   lambda <- exp(gamma)
   eta <- theta[1] + exp(theta[2])*gamma
   p <- 1- exp(-exp(eta))
   ind <- p > runif(n)
   y[!ind] <- 0
   np <- sum(ind)
   ## generate from zero truncated Poisson, given presence...
   y[ind] <- qpois(runif(np,dpois(0,lambda[ind]),1),lambda[ind])
   y
} 

library(mgcv)
## Simulate some ziP data...
set.seed(1);n<-400
dat <- gamSim(1,n=n)
dat$y <- rzip(dat$f/4-1)

b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ziP(),data=dat)

b$outer.info ## check convergence!!
b
plot(b,pages=1)
plot(b,pages=1,unconditional=TRUE) ## add s.p. uncertainty 
gam.check(b)
## more checking...
## 1. If the zero inflation rate becomes decoupled from the linear predictor, 
## it is possible for the linear predictor to be almost unbounded in regions
## containing many zeroes. So examine if the range of predicted values 
## is sane for the zero cases? 
range(predict(b,type="response")[b$y==0])

## 2. Further plots...
par(mfrow=c(2,2))
plot(predict(b,type="response"),residuals(b))
plot(predict(b,type="response"),b$y);abline(0,1,col=2)
plot(b$linear.predictors,b$y)
qq.gam(b,rep=20,level=1)

## 3. Refit fixing the theta parameters at their estimated values, to check we 
## get essentially the same fit...
thb <- b$family$getTheta()
b0 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ziP(theta=thb),data=dat)
b;b0

## Example fit forcing minimum linkage of prob present and
## linear predictor. Can fix some identifiability problems.
b2 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ziP(b=.3),data=dat)

}
\keyword{models} \keyword{regression}%-- one or more ..



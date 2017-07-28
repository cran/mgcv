\name{multinom}
\alias{multinom}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM multinomial logistic regression}
\description{Family for use with \code{\link{gam}}, implementing regression for categorical response data. Categories must be coded 0 to K, where K is a positive integer. \code{\link{gam}} should be called with a list of K formulae, one for each category except category zero (extra formulae for shared terms may also be supplied: see \code{\link{formula.gam}}). The first formula also specifies the response variable.      
}

\usage{
multinom(K=1)
}
\arguments{
\item{K}{There are K+1 categories and K linear predictors. }
}
\value{
An object of class \code{general.family}.
}

\details{
The model has K linear predictors, \eqn{\eta_j}{h_j}, each dependent on smooth functions of predictor variables, in the usual way. If response variable, y, contains the class labels 0,...,K then the likelihood for y>0 is \eqn{\exp(\eta_y)/\{1+\sum_j \exp(\eta_j) \}}{exp(h_y)/(1 + sum_j exp(h_j) )}. If y=0 the likelihood is \eqn{1/\{1+\sum_j \exp(\eta_j) \}}{1/(1 + sum_j exp(h_j) )}. In the two class case this is just a binary logistic regression model. The implementation uses the approach to GAMLSS models described in Wood, Pya and Saefken (2016).

The residuals returned for this model are simply the square root of -2 times the deviance for each observation, with a positive sign if the observed y is the most probable class for this observation, and a negative sign otherwise. 

Use \code{predict} with \code{type="response"} to get the predicted probabilities in each category.

Note that the model is not completely invariant to category relabelling, even if all linear predictors have the same form. Realistically this model is unlikely to be suitable for problems with large numbers of categories. Missing categories are not supported.   
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

\seealso{\code{\link{ocat}}}

\examples{
library(mgcv)
set.seed(6)
## simulate some data from a three class model
n <- 1000
f1 <- function(x) sin(3*pi*x)*exp(-x)
f2 <- function(x) x^3
f3 <- function(x) .5*exp(-x^2)-.2
f4 <- function(x) 1
x1 <- runif(n);x2 <- runif(n)
eta1 <- 2*(f1(x1) + f2(x2))-.5
eta2 <- 2*(f3(x1) + f4(x2))-1
p <- exp(cbind(0,eta1,eta2))
p <- p/rowSums(p) ## prob. of each category 
cp <- t(apply(p,1,cumsum)) ## cumulative prob.
## simulate multinomial response with these probabilities
## see also ?rmultinom
y <- apply(cp,1,function(x) min(which(x>runif(1))))-1
## plot simulated data...
plot(x1,x2,col=y+3)

## now fit the model...
b <- gam(list(y~s(x1)+s(x2),~s(x1)+s(x2)),family=multinom(K=2))
plot(b,pages=1)
gam.check(b)

## now a simple classification plot...
expand.grid(x1=seq(0,1,length=40),x2=seq(0,1,length=40)) -> gr
pp <- predict(b,newdata=gr,type="response")
pc <- apply(pp,1,function(x) which(max(x)==x)[1])-1
plot(gr,col=pc+3,pch=19)

}
\keyword{models} \keyword{regression}%-- one or more ..



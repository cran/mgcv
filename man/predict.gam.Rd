\name{predict.gam}
\alias{predict.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction from fitted GAM model}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} 
and produces predictions given a new set of values for the model covariates 
or the original values used for the model fit. Predictions can be accompanied
by standard errors, based on the posterior distribution of the model
coefficients. The routine can optionally return the matrix by which the model
coefficients must be pre-multiplied in order to yield the values of the linear predictor at
the supplied covariate values: this is useful for obtaining credible regions
for quantities derived from the model.}

\usage{
predict.gam(object,newdata,type="link",se.fit=FALSE,terms=NULL,
            block.size=1000,newdata.guaranteed=FALSE,na.action=na.pass,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
  
 \item{object}{ a fitted \code{gam} object as produced by \code{gam()}.
                }
 \item{newdata}{ A data frame containing the values of the model covariates at which predictions
  are required. If this is not provided then predictions corresponding to the
                original data are returned. If \code{newdata} is provided then
                it should contain all the variables needed for prediction: a
                warning is generated if not. } 
 
\item{type}{ When this has the value \code{"link"} (default) the linear predictor (possibly with
associated standard errors) is returned. When \code{type="terms"} each component of the 
linear predictor is returned seperately (possibly with standard errors): this includes 
parametric model components, followed by each smooth component, but excludes
any offset and any intercept. When \code{type="response"} predictions 
on the scale of the response are returned (possibly with approximate
standard errors). When \code{type="lpmatrix"} then a matrix is returned
which yields the values of the linear predictor (minus any offset) when
postmultiplied by the
parameter vector (in this case \code{se.fit} is ignored). The latter
option is most useful for getting variance estimates for quantities derived from
the model: for example integrated quantities, or derivatives of smooths. }

\item{se.fit}{ when this is TRUE (not default) standard error estimates are returned for each prediction.}

\item{terms}{if \code{type=="terms"} then only results for the terms given in this array
will be returned.}

\item{block.size}{maximum number of predictions to process per call to underlying
code: larger is quicker, but more memory intensive. Set to < 1 to use total number
of predictions as this.}

\item{newdata.guaranteed}{Set to \code{TRUE} to turn off all checking of
\code{newdata} except for sanity of factor levels: this can speed things up
for large prediction tasks, but \code{newdata} must be complete, with no
\code{NA} values for predictors required in the model. }

\item{na.action}{what to do about \code{NA} values in \code{newdata}. With the
default \code{na.pass}, any row of \code{newdata} containing \code{NA} values
for required predictors, gives rise to \code{NA} predictions (even if the term concerned has no
\code{NA} predictors). \code{na.exclude} or \code{na.omit} result in the
dropping of \code{newdata} rows, if they contain any \code{NA} values for
required predictors. If \code{newdata} is missing then \code{NA} handling is 
determined from \code{object$na.action}.}

\item{...}{ other arguments.}

}

\value{ If \code{type=="lpmatrix"} then a matrix is returned which will
  give a vector of linear predictor values (minus any offest) at the supplied covariate
  values, when applied to the model coefficient vector. 
Otherwise, if \code{se.fit} is \code{TRUE} then a 2 item list is returned with items (both arrays) \code{fit}
and \code{se.fit} containing predictions and associated standard error estimates, otherwise an 
array of predictions is returned. The dimensions of the returned arrays depends on whether 
\code{type} is \code{"terms"} or not: if it is then the array is 2 dimensional with each 
term in the linear predictor separate, otherwise the array is 1 dimensional and contains the 
linear predictor/predicted values (or corresponding s.e.s). The linear predictor returned termwise will 
not include the offset or the intercept.

\code{newdata} can be a data frame, list or model.frame: if it's a model frame
then all variables must be supplied.

}

\details{The standard errors produced by \code{predict.gam} are based on the
Bayesian posterior covariance matrix of the parameters \code{Vp} in the fitted
gam object. 

To facilitate plotting with \code{\link{termplot}}, if \code{object} possesses
an attribute \code{"para.only"} and \code{type=="terms"} then only parametric
terms of order 1 are returned (i.e. those that \code{termplot} can handle).

Note that, in common with other prediction functions, any offset supplied to
\code{\link{gam}} as an argument is always ignored when predicting, unlike
offsets specified in the gam model formula.

See the examples for how to use the \code{lpmatrix} for obtaining credible
regions for quantities derived from the model. 
}

\references{

Chambers and Hastie (1993) Statistical Models in S. Chapman & Hall.

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428


Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

Wood, S.N. (2004) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass. 99:637-686

\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org} 

The design is inspired by the S function of the same name described in
Chambers and Hastie (1993) (but is not a clone).

}

\section{WARNING }{ Note that the behaviour of this function is not identical to 
\code{predict.gam()} in Splus.

\code{type=="terms"} does not exactly match what \code{predict.lm} does for
parametric model components.
} 

\seealso{  \code{\link{gam}}, \code{\link{gamm}}, \code{\link{plot.gam}}}

\examples{
library(mgcv)
n<-200
sig <- 2
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
y <- 2 * sin(pi * x0)
y <- y + exp(2 * x1) 
y <- y + 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10
y <- y + x3
e <- rnorm(n, 0, sig)
y <- y + e
b<-gam(y~s(x0)+s(I(x1^2))+s(x2)+offset(x3))
rm(y,x0,x1,x2,x3)
newd <- data.frame(x0=(0:30)/30,x1=(0:30)/30,x2=(0:30)/30,x3=(0:30)/30)
pred <- predict.gam(b,newd)

## now get variance of sum of predictions using lpmatrix

Xp <- predict(b,newd,type="lpmatrix") 

## Xp %*% coef(b) yields vector of predictions

a <- rep(1,31)
Xs <- t(a) %*% Xp ## Xs %*% coef(b) gives sum of predictions
var.sum <- Xs %*% b$Vp %*% t(Xs)

## Now get the variance of non-linear function of predictions
## by simulation from posterior distribution of the params

library(MASS)
br<-mvrnorm(1000,coef(b),b$Vp) ## 1000 replicate param. vectors
res <- rep(0,1000)
for (i in 1:1000)
{ pr <- Xp %*% br[i,] ## replicate predictions
  res[i] <- sum(log(abs(pr))) ## example non-linear function
}
mean(res);var(res)

## loop is replace-able by following .... 

res <- colSums(log(abs(Xp \%*\% t(br))))
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..


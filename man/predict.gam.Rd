\name{predict.gam}
\alias{predict.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction from fitted GAM model}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} 
and produces predictions given a new set of values for the model covariates 
or the original values used for the model fit.}
}
\usage{
predict.gam(object,newdata,type="link",se.fit=FALSE,plot.call=FALSE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
  
 \item{object}{ a fitted \code{gam} object as produced by \code{gam()}.
                }
 \item{newdata}{ A data frame containing the values of the model covariates at which predictions
  are required. If this is not provided then predictions corresponding to the original data are returned.  } 
 
 
\item{type}{ When this has the value \code{"link"} (default) the linear predictor (possibly with
associated standard errors) is returned. When \code{type="terms"} each component of the 
linear predictor is returned seperately (possibly with standard errors): this includes 
parametric model components, followed by each smooth component. When \code{type="response"} predictions 
on the scale of the response are returned (possibly with approximate standard errors).}

 \item{se.fit}{ when this is TRUE (not default) standard error estimates are returned for each prediction.}
 \item{plot.call}{used to signal a call from plot.gam, to avoid having to construct a data frame containing all
 the parametric model components.}
 \item{\dots}{further arguments to be passed to or from methods.}
}

\value{ If \code{se.fit} is \code{TRUE} then a 2 item list is returned with items (both arrays) \code{fit}
and \code{se.fit} containing predictions and associated standard error estimates, otherwise an 
array of predictions is returned. The dimensions of the returned arrays depends on whether 
\code{type} is \code{"terms"} or not: if it is then the array is 2 dimensional with each 
term in the linear predictor separate, otherwise the array is 1 dimensional and contains the 
linear predictor/predicted values (or corresponding s.e.s). 


}

\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428


\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}}

\section{WARNING }{ Note that the behaviour of this function is not identical to 
\code{predict.gam()} in Splus.} 

\seealso{  \code{\link{gam}} \code{\link{plot.gam}}}

\examples{

 library(mgcv)
n<-200
sig2<-4
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
pi <- asin(1) * 2
y <- 2 * sin(pi * x0)
y <- y + exp(2 * x1) - 3.75887
y <- y + 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10 - 1.396
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- y + e
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3))
newd<-data.frame(x0=(0:30)/30,x1=(0:30)/30,x2=(0:30)/30,x3=(0:30)/30)
pred<-predict.gam(b,newd)
 
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
















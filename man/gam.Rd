\name{gam}
\alias{gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Models using penalized regression splines and 
GCV}
\description{ Fits the specified  generalized additive model to data.
The GAM is represented  using one dimensional penalized regression splines
 with smoothing parameters selected by GCV.
}
}
\usage{
 gam(formula,data,weights,family,scale)
}
%- maybe also `usage' for other objects documented here.
\details{
 Each smooth model terms is represented using a cubic penalized
     regression spline. Knots  of the spline are placed evenly
     throughout the covariate values to which the term refers:  For
     example, if fitting 101 data with a 10 knot spline of \code{x} then
     there would be a knot at every 10th (ordered)  \code{x} value. The
use
     of penalized regression splines turns  the gam fitting problem
     into  a penalized glm fitting problem, which can be fitted using a
 slight modification of  \code{glm.fit} : \code{gam.fit}.  The penalized
glm
     approach also allows smoothing parameters for all smooth terms to
     be  selected  simultaneously by GCV or UBRE. This is achieved as
     part of fitting by calling \code{mgcv}  within \code{gam.fit}.
Details are
     given in Wood (2000).
 
}

\arguments{ 
  
 \item{formula}{A GAM formula. This is exactly like the formula for a
glm exept that smooth terms can be added to the right hand side of the
formula, and the l.h.s. must contain only the names of a variable, and
not some transformation function applied to a named variable.
Smooth terms are specified by expressions of the form:
\code{s(var,knots)} where \code{var} is the covariate which the smooth
is a function of and \code{var} is the number of knots for the spline
representing this smooth.}
 \item{data}{ A data frame containing the model covariates required by the
formula. If this is missing the searth list is used to try and find the
variables needed.} 
\item{weights}{ 
prior weights on the data.
}
 \item{family}{
This is a family object specifying the distribution and link to use on
fitting etc. See \code{\link{glm}} and \code{\link{family}} for more
details. 
} 
\item{scale}{ If this is zero then GCV is used for all distributions
except Poisson and binomial where UBRE is used with scale parameter
assumed to be 1. If this is greater than 1 it is assumed to be the scale
parameter/variance and UBRE is used, otherwise GCV is used.} }

\value{ 
 The function returns an object of class \code{"gam"}. 

\item{coefficients}{ the coefficients of the fitted model. parametric
          coefficients are  first followed  by coefficients for each
          spline term in turn.}

\item{residuals}{ the deviance residuals for the fitted model.}

\item{fitted.values}{ fitted model predictions of expected value for each
          datum.}
\item{family}{family object specifying distribution and link used.
}
\item{linear.predictor}{fitted model prediction of link function of
expected
          value for  each datum.}
\item{deviance}{(unpenalized)}

\item{null.deviance:} 

\item{df.null}{} 

\item{iter}{number of iterations of IRLS taken to get convergence.}

\item{weights}{final weights used in IRLS iteration.}

\item{prior.weights}{} 

\item{y}{response data.}

\item{converged:}{} 

\item{sig2}{estimated or supplied variance/scale parameter.}

\item{edf}{estimated degrees of freedom for each smooth.}

\item{boundary}{} 
\item{sp}{ smoothing parameter for each smooth.}

\item{df}{ number of knots for each smooth (one more than maximum
          degrees of freedom).}

\item{nsdf}{number of parametric, non-smooth, model terms excluding the
          intercept.}

\item{Vp}{ estimated covariance matrix for parameters.}

\item{xp}{knot locations for each smooth. xp[i,] are the locations for
          the ith smooth.}

\item{formula}{ the model formula.}

\item{x}{parametric design matrix columns (excluding intercept)
          followed by the data that form arguments of the smooth.}

}

\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2)


\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}}

\section{WARNING }{The code does not check for rank defficiency of the
model matrix -
      it will likely just fail instead!
 } 

\seealso{  \code{\link{predict.gam}} \code{\link{plot.gam}}}

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
plot(b) 
}
\keyword{GAM, GCV, AIC, UBRE, penalized regression spline }%-- one or more ...

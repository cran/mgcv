\name{gam.check}
\alias{gam.check}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Some diagnostics for a fitted gam model.}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and produces some diagnostic information
about the fitting procedure and results. 
}
\usage{
gam.check(b)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{b}{ a fitted \code{gam} object as produced by \code{gam()}.}
}
\details{ This function plots 4 standard diagnostic plots, and some other convergence diagnostics. 

The first plot shows the GCV or UBRE score against model 
degrees of freedom given the final estimates of the relative smoothing parameters for the model. This is a slice through the 
GCV/UBRE score function that passes through the minimum found during fitting. Although not conclusive (except in the single 
smoothing parameter case), a lack of multiple local minima on this plot is suggestive of a lack of multiple local minima in the 
GCV/UBRE function and is therefore a good thing. Multiple local minima on this plot indicates that the GCV/UBRE function 
may have multiple local minima, but in a multiple smoothing parameter case this is not conclusive - multiple local minima on one slice 
through a function do not necessarily imply that the function has multiple local minima. A "good" plot here is a smooth curve with 
only one local minimum (which is therefore its global minimum). In any case, note that possible problems with multiple local minima are 
substantially reduced by the global line search with respect to the overall model smoothing parameter that is performed at 
each iteration of the minimizer (and produces the data used in this plot).

The other 3 plots are two residual plots and plot of fitted values against original data.

The function also prints out information about the convergence of the GCV minimization algorithm, indicating how 
many iterations  were required to minimise the GCV/UBRE score. A message is printed if the minimization terminated by 
failing to improve the score with a steepest descent step: otherwise minimization terminated by meeting convergence criteria.
The mean absolute gradient of the GCV/UBRE function at the minimum is given (compare this to the range of the score shown on the GCV/UBRE 
plot). An indication of whether or not the Hessian of the GCV/UBRE function is positive definite is given. If some smoothing parameters 
are not well defined (effectively zero, or infinite) then it may not be, although this is not usually a problem. A message is printed 
if the second guess smoothing parameters did not improve on the first guess - this is primarily there for the developers. 

The ideal results from this function have a smooth, single minima GCV/UBRE plot, good residual plots, and convergence to small 
gradients with a positive definite Hessian. However, failure to meet some of these criteria 
is often  acceptable, 
and the information provided is primarily of use in diagnosing suspected problems. 
High gradients at convergence  are a clear indication of problems, however.

Fuller data can be extracted from \code{mgcv.conv} part of the \code{gam} object.  

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


\seealso{  \code{\link{gam}} \code{\link{mgcv}}}

\examples{
library(mgcv)
set.seed(0)
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
plot(b,pages=1)
gam.check(b)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







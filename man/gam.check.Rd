\name{gam.check}
\alias{gam.check}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Some diagnostics for a fitted gam model}
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
\details{ This function plots 4 standard diagnostic plots, and some other
convergence diagnostics. Output differs depending on whether the underlying
fitting method was \code{\link{mgcv}} or another method (see \code{\link{gam.method}}).

For \code{\link{mgcv}} based fits , the first plot shows the GCV or UBRE score against model 
degrees of freedom, given the final estimates of the relative smoothing
parameters for the model. This is a slice through the 
GCV/UBRE score function that passes through the minimum found during fitting. Although not conclusive (except in the single 
smoothing parameter case), a lack of multiple local minima on this plot is
suggestive of a lack of multiple local minima 
in the GCV/UBRE function and is therefore a good thing. Multiple local minima on this plot indicates that the GCV/UBRE function 
may have multiple local minima, but in a multiple smoothing parameter case this is not conclusive - multiple local minima on one slice 
through a function do not necessarily imply that the function has multiple local minima. A `good' plot here is a smooth curve with 
only one local minimum (which is therefore its global minimum).

The location of the minimum used for the fitted model is also marked on the first plot. Sometimes this location may be a local minimum
that is not the global minimum on the plot. There is a legitimate reason for this to happen, and it does not always indicate problems.
Smoothing parameter selection is based on applying GCV/UBRE to the approximating linear model produced by the GLM IRLS fitting method
employed in \code{gam.fit()}. It is sometimes possible for these approximating models to develop `phantom' minima in their GCV/UBRE scores. These 
minima usually imply a big change in model parameters, and have the characteristic that the minimia will not be present in the GCV/UBRE score 
of the approximating model that would result from actually applying this parameter change. In other words, these are spurious minima in regions 
of parameter space well beyond those for which the weighted least squares problem can be expected to represent the real underlying likelihood well.
Such minima can lead to convergence problems. To help ensure convergence even in the presence of phantom minima, 
\code{gam.fit} switches to a cautious optimization mode after a user controlled number of iterations of the IRLS algorithm (see \link{gam.control}).
In the presence of local minima in the GCV/UBRE score, this method selects the minimum that leads to the smallest change in 
model estimated degrees of freedom. This approach is usually sufficient to deal with phantom minima. Setting \code{trace} to \code{TRUE} in 
\code{gam.control} will allow you to check exactly what is happening. 

If the location of the point indicating the minimum is not on the curve showing the GCV/UBRE function then there are numerical problems with the 
estimation of the effective degrees of freedom: this usually reflects problems with  the relative scaling of covariates that are arguments of a single smooth. 
In this circumstance reported estimated degrees of freedom can not be trusted, although the fitted model and term estimates are likely to be quite 
acceptable. 

If the fit method is based on \code{\link{magic}} or \code{\link{gam.fit2}} then there is no global search and the problems with phantom local minima are much reduced. The first plot in this case will simply be a normal QQ plot of the standardized residuals.

The other 3 plots are two residual plots and plot of fitted values against original data.

The function also prints out information about the convergence of the GCV minimization algorithm, indicating how 
many iterations  were required to minimise the GCV/UBRE score. A message is printed if the minimization terminated by 
failing to improve the score with a steepest descent step: otherwise minimization terminated by meeting convergence criteria.
The mean absolute gradient or RMS gradient of the GCV/UBRE function at the minimum is given. An indication of whether or not the Hessian of the GCV/UBRE function is positive definite is given. If some smoothing parameters 
are not well defined (effectively zero, or infinite) then it may not be, although this is not usually a problem. If the fit method is \code{\link{mgcv}}, a message is printed 
if the second guess smoothing parameters did not improve on the first guess - this is primarily there for the developer. 
For other fit methods the estimated rank of the model is printed.


The ideal results from this function have a smooth, single minima GCV/UBRE
plot (where appropriate), 
good residual plots, and convergence to small 
gradients with a positive definite Hessian. However, failure to meet some of these criteria 
is often  acceptable,  and the information provided is primarily of use in diagnosing suspected problems. 
High gradients at convergence  are a clear indication of problems, however.

Fuller data can be extracted from \code{mgcv.conv} part of the \code{gam} object.  

}

\references{

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

\url{http://www.stats.gla.ac.uk/~simon/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  \code{\link{gam}}, \code{\link{mgcv}}, \code{\link{magic}}}

\examples{
library(mgcv)
set.seed(0)
n<-200
sig<-2
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
y <- 2 * sin(pi * x0)
y <- y + exp(2 * x1) - 3.75887
y <- y+0.2*x2^11*(10*(1-x2))^6+10*(10*x2)^3*(1-x2)^10-1.396
e <- rnorm(n, 0, sig)
y <- y + e
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3))
plot(b,pages=1)
gam.check(b)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







\name{gam.check}
\alias{gam.check}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Some diagnostics for a fitted gam model}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and produces some diagnostic information
about the fitting procedure and results. The default is to produce 4 residual
plots, and some information about the convergence of the smoothness selection optimization.
}
\usage{
gam.check(b)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{b}{ a fitted \code{gam} object as produced by \code{gam()}.}
}
\details{ This function plots 4 standard diagnostic plots, and some other
convergence diagnostics. Usually the 4 plots are various residual plots. The
printed information relates to the optimization used to select smoothing
parameters. For the default optimization methods the information is summarized in a
readable way, but for other optimization methods, whatever is returned by way of
convergence diagnostics is simply printed.

For \code{\link{mgcv}} based fits (not the default), the first plot shows the GCV or UBRE score against model 
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

If the fit method is based on \code{\link{magic}} or \code{\link{gam.fit2}} or
\code{\link{gam.fit3}}then there is no global search and the problems with
phantom local minima are much reduced. These more recent methods are also much
more robust than the \code{\link{mgcv}} based methods.

}

\references{

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman
and Hall/CRC Press.


\url{http://www.maths.bath.ac.uk/~sw283/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{ \code{\link{choose.k}},  \code{\link{gam}}, \code{\link{mgcv}}, \code{\link{magic}}}

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







\name{plot.gam}
\alias{plot.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Default GAM plotting}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and plots the 
   component smooth functions that make it up, on the scale of the linear predictor.}
}
\usage{
plot.gam(object,rug=TRUE,se=TRUE,pages=0,scale=-1,n)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
  
 \item{object}{ a fitted \code{gam} object as produced by \code{gam()}.
                }
 \item{rug}{ when TRUE (default) then the covariate to which the plot applies is displayed as a rug plot
 at the foot of each plot.} 
  
\item{se}{ when TRUE (default) upper and lower lines are added to the plot at 2 standard errors
        above and below the estimate of the smooth being plotted.}
 \item{pages}{ (default 0) the number of pages over which to spread the output. For example, 
if \code{pages=1} then all terms will be plotted on one page with the layout performed automatically. 
Set to 0 to have the routine leave all graphics settings as they are. }
\item{scale}{ set to -1 (default) to have the same y-axis scale for each plot, and to 0 for a 
          different y axis for each plot.  }
\item{n}{ number of points used for each plot - for a nice smooth plot this needs to be several times the estimated 
    degrees of freedom for the smooth. Default value 100.}
}
\details{ Produces default plot showing the smooth components of a fitted GAM. The x axis of each plot is labelled 
   with the covariate name, while the y axis is labelled \code{s(cov,edf) } where \code{cov}
   is the covariate name, and \code{edf} the estimated (or user defined for regression splines) degrees of freedom of the smooth.
   Within the function, the data for the plots is obtained by a call to \code{predict.gam()} with a 
   suitably constructed data frame. A rug plot along the bottom of each frame shows the covariate values.}
   }

\value{ The function simply generates plots.
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
\code{plot.gam()} in Splus.} 

\seealso{  \code{\link{gam}} \code{\link{predict.gam}}}

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
plot(b,pages=1)

 
}
\keyword{GAM, GCV, AIC, UBRE, penalized regression spline }%-- one or more ...

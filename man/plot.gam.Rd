\name{plot.gam}
\alias{plot.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Default GAM plotting}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and plots the 
   component smooth functions that make it up, on the scale of the linear predictor.}
}
\usage{
plot.gam(x,rug=TRUE,se=TRUE,pages=0,select=NULL,scale=-1,n=100,n2=40,pers=FALSE,theta=30,phi=30,jit=FALSE,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
  
 \item{x}{ a fitted \code{gam} object as produced by \code{gam()}.
                }
 \item{rug}{ when TRUE (default) then the covariate to which the plot applies is displayed as a rug plot
 at the foot of each plot of a 1-d smooth, and the loactions of the
 covariates are plotted as points on the contour plot representing a 2-d
 smooth.} 
  
\item{se}{ when TRUE (default) upper and lower lines are added to the
  1-d plots at 2 standard errors
        above and below the estimate of the smooth being plotted. For
	2-d plots, surfaces at +1 and -1 standard errors are contoured
	and overlayed on the contour plot for the estimate.}
 \item{pages}{ (default 0) the number of pages over which to spread the output. For example, 
if \code{pages=1} then all terms will be plotted on one page with the layout performed automatically. 
Set to 0 to have the routine leave all graphics settings as they are. }
\item{select}{Allows the  plot for a single model term to be selected for printing. e.g. if you just want the plot for the second smooth term set \code{select=2}. }
\item{scale}{ set to -1 (default) to have the same y-axis scale for each plot, and to 0 for a 
          different y axis for each plot.  }
\item{n}{ number of points used for each 1-d plot - for a nice smooth plot this needs to be several times the estimated 
  degrees of freedom for the smooth. Default value 100.}
\item{n2}{Square root of number of points used to grid estimates of 2-d
  functions for contouring.}
\item{pers}{Set to \code{TRUE} if you want perspective plots for 2-d
  terms.}
\item{theta}{One of the perspective plot angles.}
\item{phi}{The other perspective plot angle.}
\item{jit}{Set to TRUE if you want rug plots for 1-d terms to be jittered.}
\item{...}{ other arguments.}
}
\details{ Produces default plot showing the smooth components of a
  fitted GAM.

  For plots of 1-d smooths, the x axis of each plot is labelled 
   with the covariate name, while the y axis is labelled \code{s(cov,edf) } where \code{cov}
   is the covariate name, and \code{edf} the estimated (or user defined for regression splines) degrees of freedom of the smooth.

Contour plots are produced for 2-d smooths with the x-axes labelled with the first covariate
name and the y axis with the second covariate name. The main title of
the plot is something like \code{s(var1,var2,edf)}, indicating the
variables of which the term is a function, and the estimated degrees of
freedom for the term. When \code{se=TRUE}, estimator variability is shown by overlaying
contour plots at plus and minus 1 s.e. relative to the main
estimate. Contour levels are chosen to try and ensure reasonable
separation of the contours of the different plots, but this is not
always easy to achieve.  
   
   Within the function, the data for the plots is obtained by a call to \code{predict.gam()} with a 
   suitably constructed data frame.

   Smooths of more than 2 variables are not currently dealt with, but
   simply generate a warning.
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
  \code{plot.gam()} in Splus.

  The function can not deal with smooths of more than 2 variables!
} 

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
# example with 2-d plots
b1<-gam(y~s(x0,x1)+s(x2)+s(x3))
op<-par(mfrow=c(2,2))
plot(b1)
par(op) 
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







\name{gam.check}
\alias{gam.check}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Some diagnostics for a fitted gam model}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and produces some diagnostic information
about the fitting procedure and results. The default is to produce 4 residual
plots, and some information about the convergence of the smoothness selection optimization.
}
\usage{
gam.check(b,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{b}{ a fitted \code{gam} object as produced by \code{gam()}.}
\item{...}{extra graphics parameters to pass to plotting functions.}
}
\details{ This function plots 4 standard diagnostic plots, and some other
convergence diagnostics. Usually the 4 plots are various residual plots. The
printed information relates to the optimization used to select smoothing
parameters. For the default optimization methods the information is summarized in a
readable way, but for other optimization methods, whatever is returned by way of
convergence diagnostics is simply printed.
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
dat <- gamSim(1,n=200)
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
plot(b,pages=1)
gam.check(b)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







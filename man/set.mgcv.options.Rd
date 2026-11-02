\name{set.mgcv.options}
\alias{set.mgcv.options}
\alias{mgcv.options}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Set mgcv options}
\description{Allows convenient setting and clearing of mgcv options. Options are set using \code{\link{options}} after having \code{mgcv.} appended to the front of the option name. \code{NULL} leads to function defined defaults being used. 
}

\usage{
set.mgcv.options(vc.logrange=25,s.bs=NULL,xt=NULL,te.bs=NULL,te.np=NULL,
ti.mc=NULL,gam.method=NULL,gam.optimizer=NULL,bam.discrete=NULL,plot.scheme=NULL)
}
\arguments{
\item{vc.logrange}{used by \code{\link{notLog}}.}
\item{s.bs}{default \code{\link{s}} basis.}
\item{xt}{default \code{xt} argument to \code{\link{s}}.}
\item{te.bs}{default \code{\link{te}} basis.}
\item{te.np}{logical indicating whether \code{\link{te}} marginals should th re-parameterized.  }
\item{ti.mc}{default \code{mc} argument for \code{\link{te}}.}
\item{gam.method}{default \code{method} argument for \code{\link{gam}}.}
\item{gam.optimizer}{default \code{optimizer} argument for \code{\link{gam}}.}
\item{bam.discrete}{default \code{dscrete} argument for \code{\link{bam}}.}
\item{plot.scheme}{default plot \code{scheme} for \code{\link{plot.gam}}.}
}

\value{
Value returned by \code{\link{options}}.
}

\details{Some \code{mgcv} defaults can usefully be reset: for example the \code{gam} default method is GCV for back-compatibulity and REML is usually better. This routine facilitiates this. 
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood  \email{simon.wood@r-project.org}
}


\examples{
library(mgcv)
set.seed(2) ## simulate some data... 
dat <- gamSim(1,n=400,dist="normal",scale=2)
set.mgcv.options(gam.method="REML",plot.scheme=2)
b <- gam(y~s(x0)+s(x1)+te(x2,x3),data=dat)
b
plot(b,pages=1) 
set.mgcv.options() ## reset to 'factory fresh'
b <- gam(y~s(x0)+s(x1)+te(x2,x3),data=dat)
b
plot(b,pages=1) 
}
\keyword{models} \keyword{regression}%-- one or more ..



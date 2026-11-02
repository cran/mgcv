\name{options.mgcv}
\alias{options.mgcv}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Options in mgcv}
\description{Some \code{mgcv} options can be set using \code{\link{options}}. These are described here.}


\details{
The following options control smoother behaviour. 
\itemize{
\item \code{mgcv.s.bs} give the default basis used for \code{\link{s}} terms. Defaults to \code{c("tp","tp")}, the first being the default for univariate smoothing and the second for mutivariate.
\item \code{mgcv.te.bs} is the default basis for \code{\link{te}} and \code{ti} smooths. Defaults to \code{"cr"}.
\item \code{mgcv.te.np} is the default \code{np} parameter for  \code{\link{te}} and \code{ti}. Defaults to \code{TRUE}.
\item \code{mgcv.xt} is the default \code{xt} parameter for  \code{\link{s}}, \code{\link{te}} and \code{ti}. \code{NULL} by default.
\item \code{mgcv.ti.mc} is the default \code{mc} parameter for \code{\link{ti}}. \code{NULL} by default.}

The following options control \code{\link{gam}} and \code{\link{bam}} methods.
\itemize{
\item \code{mgcv.gam.method} gives the default \code{method} argument for \code{\link{gam}}. \code{"GCV.Cp"} by default (historical, \code{"REML"} usually better).
\item \code{mgcv.gam.optimizer} gives the default \code{optimizer} argument for \code{\link{gam}}. \code{c("outer","newton")} by default.
\item \code{mgcv.bam.discrete} gives the default \code{discrete} argument for \code{\link{bam}}. \code{FALSE} by default (historical, \code{TRUE} usually better).
}

The following controls plotting defaults.
\itemize{
\item \code{mgcv.plot.scheme} controls the plotting \code{scheme} argument to \code{\link{plot.gam}}. Defaults to 0. 
}

}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood (s.wood@r-project.org).
}



\keyword{models} \keyword{regression}%-- one or more ..



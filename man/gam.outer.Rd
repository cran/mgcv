\name{gam.outer}
\alias{gam.outer}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Minimize GCV or UBRE score of a GAM using `outer' iteration}
\description{Estimation of GAM smoothing parameters is most stable if
optimization of the UBRE or GCV score is outer to the penalized iteratively
re-weighted least squares scheme used to estimate the model given smoothing 
parameters.

This routine optimizes a GCV or UBRE score in this way. Basically the GCV or
UBRE score is evaluated for each trial set of smoothing parameters by
estimating the GAM for those smoothing parameters. The score is minimized
w.r.t. the parameters numerically, using \code{optim} or \code{nlm}. Exact
derivatives of the score can be used by fitting with \code{\link{gam.fit2}},
which improves efficiency and reliability relative to relying solely on finite
difference derivatives.  

Note that there is a choise between basing GCV/UBRE scores on the deviance or
the Pearson statistic: see \code{\link{gam.method}}.

Not normally called directly, but rather a service routine for \code{\link{gam}}.
}
\usage{
gam.outer(lsp,fscale,family,control,method,gamma,G)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{lsp}{The log smoothing parameters.}

\item{fscale}{Typical scale of the GCV or UBRE score.}

\item{family}{the model family.}

\item{control}{control argument to pass to \code{\link{gam.fit}} if pure
finite differencing is being used.}

\item{method}{method list reurned from \code{\link{gam.method}}. This defines
the optimization method to use.}

\item{gamma}{ The degree of freedom inflation factor for the GCV/UBRE score.}

\item{G}{List produced by \code{\link{gam.setup}}, containing most of what's
needed to actually fit GAM.}

}
\details{ 
Estimation of smoothing parameters by optimizing GCV scores obtained at
convergence of the P-IRLS iteration was proposed by O'Sullivan et al. (1986),
and is here termed `outer' iteration. 
}

\references{

O 'Sullivan, Yandall & Raynor (1986) Automatic smoothing of regression
functions in generalized linear models. J. Amer. Statist. Assoc. 81:96-103.

\url{http://www.maths.bath.ac.uk/~sw283/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{\code{\link{gam.fit2}},  \code{\link{gam}}, \code{\link{mgcv}}, \code{\link{magic}}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...


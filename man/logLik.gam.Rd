\name{logLik.gam}
\alias{logLik.gam}
\alias{AIC.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Log likelihood for a fitted GAM, for AIC}
\description{ Function to extract the log-likelihood for a fitted \code{gam}
model (note that the models are usually fitted by penalized likelihood maximization). 
Used by \code{\link{AIC}}.
}
\usage{
\method{logLik}{gam}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{ fitted model objects of class \code{gam} as produced by \code{gam()}.}
\item{...}{un-used in this case} 

}
\details{ Modification of \code{logLik.glm} which corrects the degrees of
freedom for use with \code{gam} objects. 

The function is provided so that \code{\link{AIC}} functions correctly with
\code{gam} objects, and uses the appropriate degrees of freedom (accounting
for penalization). Note,  when using \code{AIC} for penalized models, that the
degrees of freedom are the effective degrees of freedom and not the number of
parameters, and the model maximizes the penalized likelihood, not the actual
likelihood. (See e.g. Hastie and Tibshirani, 1990, section 6.8.3 and also Wood 2008), 

By default this routine uses a definition of the effective degrees of freedom that 
includes smoothing parameter uncertainty, if this is available (i.e. if smoothing 
parameter selection is by some variety of marginal likelihood).
  
}

\value{ Standard \code{logLik} object: see \code{\link{logLik}}.
}


\references{Hastie and Tibshirani, 1990, Generalized Additive Models.

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for
generalized additive models. J.R.Statist. Soc. B 70(3):495-518



}

\author{ Simon N. Wood \email{simon.wood@r-project.org} based directly on \code{logLik.glm}}


\seealso{  \code{\link{AIC}}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







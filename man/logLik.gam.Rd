\name{logLik.gam}
\alias{logLik.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract the log likelihood for a fitted GAM}
\description{ Function to extract the log-likelihood for a fitted \code{gam}
model (note that the models are usually fitted by penalized likelihood maximization). 
}
\usage{
logLik.gam(object,...)
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
parameters, and the model maximizes the penalized likelihood, not the actual likelihood!

}

\value{ Standard \code{logLik} object: see \code{\link{logLik}}.
}


\author{ Simon N. Wood \email{simon.wood@r-project.org} based directly on \code{logLik.glm}}

\section{WARNING }{ Note that the likelihood is evaluated at the maximum of
the penalized likelihood for most GAMs, rather than at the MLE.
} 

\seealso{  \code{\link{AIC}}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







\name{residuals.gam}
\alias{residuals.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Model residuals}
\description{Returns residuals for a fitted \code{gam} model
  object. Pearson, deviance, working and response residuals are
  available. 
}

\usage{
residuals.gam(object, type = c("deviance", "pearson", "working", "response"),...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ a \code{gam} fitted model object. }
  \item{type}{the type of residuals wanted. }
   \item{...}{other arguments.}
}
\details{Response residuals are the raw residuals (data minus fitted
  values). Pearson residuals are raw residuals divided by the standard
  deviation of the data according to the model mean variance
  relationship (\eqn{(y-\mu)/\sqrt{V(\mu)}}{(y-m)/V(m)^0.5}, where
  \eqn{y}{y} is data \eqn{\mu}{m} is model fitted value and \eqn{V}{V}
  is model mean-variance relationship.). Deviance residuals simply
  return the deviance residuals defined by the model family. Working
  residuals are the residuals returned from model fitting at convergence.

  There is  a special function for \code{gam} objects because of a bug
  in the calculation of Pearson residuals in some versions of
  residual.glm().
 
}
\value{ An array of residuals.
}


\author{ Simon N. Wood \email{snw@st-and.ac.uk}}



\seealso{   \code{\link{gam}} \code{\link{mgcv}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..













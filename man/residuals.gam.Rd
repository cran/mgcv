\name{residuals.gam}
\alias{residuals.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Model residuals}
\description{Returns residuals for a fitted \code{gam} model
  object. Pearson, deviance, working and response residuals are
  available. 
}

\usage{
\method{residuals}{gam}(object, type = "deviance",...)
}
\arguments{
  \item{object}{ a \code{gam} fitted model object. }
  \item{type}{the type of residuals wanted. Usually one of 
             \code{"deviance"}, \code{"pearson"},\code{"scaled.pearson"}, 
             \code{"working"}, or \code{"response"}.  }
   \item{...}{other arguments.}
}
\details{Response residuals are the raw residuals (data minus fitted
  values). Scaled Pearson residuals are raw residuals divided by the standard
  deviation of the data according to the model mean variance
  relationship and estimated scale parameter. Pearson residuals are the same, but multiplied by the square root 
  of the scale parameter (so they are independent of the scale parameter):
  (\eqn{(y-\mu)/\sqrt{V(\mu)}}{(y-m)/V(m)^0.5}, where  \eqn{y}{y} is data \eqn{\mu}{m} 
  is model fitted value and \eqn{V}{V} is model mean-variance relationship.). Both are provided since not all texts 
  agree on the definition of Pearson residuals. Deviance residuals simply
  return the deviance residuals defined by the model family. Working
  residuals are the residuals returned from model fitting at convergence.

  Families can supply their own residual function, which is used in place of the standard 
  function if present, (e.g. \code{\link{cox.ph}}).

}
\value{ A vector of residuals.
}


\author{ Simon N. Wood \email{simon.wood@r-project.org}}



\seealso{   \code{\link{gam}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..













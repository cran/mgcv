\name{fix.family.link}
\alias{fix.family.link}
\alias{fix.family.var}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Modify families for use in GAM fitting}
\description{ Generalized Additive Model fitting by `outer' iteration,
requires extra derivatives of the variance and link functions to be 
added to family objects. The functions add what is needed.
}
\usage{
fix.family.link(fam)
fix.family.var(fam)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{fam}{ A \code{family}.}
}
\details{ Outer iteration GAM estimation requires derivatives of the GCV/UBRE
score, which are obtained by differentiating the P-IRLS GAM fitting iteration
w.r.t. the model smoothing parameters. The expressions for the derivatives
require the second derivative of the link w.r.t. the mean, and the first
derivative of the variance function w.r.t. the mean. These functions add
functions evaluating these quantities to a family.
}
\value{ A family object with extra component functions \code{dvar} and \code{d2link}.
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\keyword{models} \keyword{regression}%-- one or more ..



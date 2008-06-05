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
\item{fam}{A \code{family}.}
}
\details{ Outer iteration GAM estimation requires derivatives of the GCV/UBRE
score, which are obtained by differentiating the P-IRLS GAM fitting iteration
w.r.t. the model smoothing parameters. The expressions for the derivatives
require the second and third derivatives of the link w.r.t. the mean, and the first
and second derivatives of the variance function w.r.t. the mean. These functions add
functions evaluating these quantities to a family. 

If the family already has functions \code{dvar}, \code{d2var}, \code{d2link}
and \code{d3link} then these
functions simply return the family unmodified: this allows non-standard links
to be used with \code{\link{gam}} when using outer iteration (performance
iteration operates with unmodified families).

The \code{dvar} function is a function of a mean vector, \code{mu}, and returns
a vector of corresponding first derivatives of the family variance
function. The \code{d2link} function is also a function of a vector of mean
values, \code{mu}: it returns a vector of second derivatives of the link,
evaluated at \code{mu}. Higher derivatives are defined similarly.

If modifying your own family, note that you can often get away with supplying
only a \code{dvar} and \code{d2var}, function if your family only requires links that occur in
one of the standard families.

}
\value{ A family object with extra component functions \code{dvar},
\code{d2var}, \code{d2link} and \code{d3link}.
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\seealso{  
\code{\link{gam.fit3}}} 


\keyword{models} \keyword{regression}%-- one or more ..



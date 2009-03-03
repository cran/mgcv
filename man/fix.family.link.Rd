\name{fix.family.link}
\alias{fix.family.link}
\alias{fix.family.var}
\alias{fix.family.ls}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Modify families for use in GAM fitting}
\description{ Generalized Additive Model fitting by `outer' iteration,
requires extra derivatives of the variance and link functions to be 
added to family objects. The functions add what is needed.
}
\usage{
fix.family.link(fam)
fix.family.var(fam)
fix.family.ls(fam)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{fam}{A \code{family}.}
}
\details{ Outer iteration GAM estimation requires derivatives of the GCV, UBRE/gAIC,
GACV, REML or ML score, which are obtained by finding the derivatives of the model
coefficients w.r.t. the log smoothing parameters, using the implicit function theorem. 
 The expressions for the derivatives require the second and third derivatives of the link 
w.r.t. the mean (and the 4th derivatives if Fisher scoring is not used). Also required are the first
and second derivatives of the variance function w.r.t. the mean (plus the third 
derivative if Fisher scoring is not used). Finally REML or ML estimation of smoothing parameters
requires the log saturated likelihood and its first two derivatives w.r.t. the scale parameter.
These functions add functions evaluating these quantities to a family. 

If the family already has functions \code{dvar}, \code{d2var}, \code{d3var}, \code{d2link},
\code{d3link}, \code{d4link} and for RE/ML \code{ls}, then these functions simply 
return the family unmodified: this allows non-standard links
to be used with \code{\link{gam}} when using outer iteration (performance
iteration operates with unmodified families). Note that if you only need Fisher scoring then 
\code{d4link} and \code{d3var} can be dummy, as they are ignored. similalry \code{ls} is only needed for 
RE/ML.

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
\code{d2var}, \code{d2link}, \code{d3link}, \code{d4link} and \code{ls}. 
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\seealso{  
\code{\link{gam.fit3}}} 


\keyword{models} \keyword{regression}%-- one or more ..



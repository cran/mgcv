\name{gam.parser}
\alias{gam.parser}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Model fitting using penalized regression splines and 
GCV}
\description{ This is an internal function of package \code{mgcv}. It is a service routine for
 \code{gam} which splits of the strictly parametric part of the model formula and returns it as a formula and interprets the smooth parts 
of the model formula.  
 The processing of 
 the smooth terms results in the return of the basis dimension for each smooth term, whether 
 or not the smooth term is to be treated as a pure regression spline
 (i.e. one without a penalty) and the type of basis used to represent
 the spline. It also returns an array of covariate names which are arguments of smooths, and the name of the response variable.  
 For further information on usuage see code for \code{gam}.
}
}



\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}


\seealso{   \code{\link{gam}} \code{\link{mgcv}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..




\name{print.gam}
\alias{print.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Model default print statement}
\description{ This is the default print statement for a GAM object. If you need a list of everything
that is part of a gam object see \code{\link{gam}}, or use \code{names()}. The family (including link), 
model formula, and estimated degrees of freedom for each model term (plus total) are printed, as 
well as the minimized GCV or UBRE score, depending on which was used.
}


\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org} }



\seealso{   \code{\link{gam}},\code{\link{mgcv}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










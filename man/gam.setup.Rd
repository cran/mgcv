\name{gam.setup}
\alias{gam.setup}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Model set up.}
\description{ This is an internal function of package \code{mgcv}. It is called 
by \code{gam} to obtain the design matrix and penalty matrices for a GAM
set up using penalized regression splines. This is done by calling a mixture of
R routines and compiled C code.  For further information on usuage see code for \code{gam}.
}
}



\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428

Wood (2003) Thin plate regression splines. JRSSB in press

\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}}



\seealso{   \code{\link{gam}}, \code{\link{mgcv}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









\name{gam.fit}
\alias{gam.fit}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Models fitting using penalized regression splines and 
GCV}
\description{ This is an internal function of package \code{mgcv}. It is a modification
  of the function \code{glm.fit}, designed to be called from \code{gam}. The major
  modification is that rather than solving a weighted least squares problem at each IRLS step, 
  a weighted, penalized least squares problem
  is solved at each IRLS step with smoothing parameters associated with each penalty chosen by GCV or UBRE,
  using routine \code{mgcv}. For further information on usage see code for \code{gam}. Some regularization of the 
IRLS weights is also permitted as a way of addressing identifiability related problems (see 
\code{\link{gam.control}}). Negative binomial parameter estimation is supported.
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

\section{WARNING }{The code does not check for rank defficiency of the
model matrix -
      it will likely just fail instead!
 } 

\seealso{   \code{\link{gam}} \code{\link{mgcv}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..

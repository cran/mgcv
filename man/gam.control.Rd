\name{gam.control}
\alias{gam.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Setting Generalized Additive Models fitting defaults}
\description{ This is an internal function of package \code{mgcv}. It is a modification
  of the function \code{glm.control}. It enables the user to set defaults for convergence 
  tolerance and maximum number of iterations when using
  \code{gam}. Argument \code{mgcv.tol} controls the tolerence used to judge
  multiple smoothing parameter convergence (default 1e-6). \code{mgcv.max.half} is the
  maximum number of times to halve the step length in an iterative update of the
  smoothing parameters if the step is failing to decrease the GCV/UBRE
  score (default 15). See \code{\link{glm.control}} 
  for more information on other arguments. 
}
}



\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428


\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}}


\seealso{   \code{\link{gam}} \code{\link{gam.fit}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..



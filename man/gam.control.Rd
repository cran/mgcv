\name{gam.control}
\alias{gam.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Setting Generalized Additive Models fitting defaults}
\description{ This is an internal function of package \code{mgcv}. It is a modification
  of the function \code{glm.control}. It enables the user to set defaults for convergence 
  tolerance and maximum number of iterations when using
  \code{gam}. Argument \code{mgcv.tol} controls the tolerence used to judge
  multiple smoothing parameter convergence (default 1e-6). \code{mgcv.half} is the
  maximum number of times to halve the step length in an iterative update of the
  smoothing parameters if the step is failing to decrease the GCV/UBRE
  score (default 15). Setting \code{trace} to \code{TRUE} will cause various 
  diagnostics to be produced as fitting progresses, including plots of the GCV function 
  and current minimum, against model EDF. \code{maxit} and \code{globit} control the maximum
  iterations of the IRLS algorithm, as follows: the algorithm will first execute up to
\code{globit} steps in which the GCV/UBRE algorithm performs a global search for the best overall 
smoothing parameter at every iteration, if convergence has not occured by then, then a further 
\code{maxit} steps are taken, in which the overall smoothing parameter estimate is taken as the 
one locally minimising the GCV/UBRE score and resulting in the lowest EDF change. The difference 
between the two phases is only significant if the GCV/UBRE function develops more than one minima. 
The reason for this approach is that the GCV/UBRE score for the IRLS problem can develop "phantom"
minimima for some models: these are minima which are not present in the GCV/UBRE score of the IRLS
problem resulting from moving the parameters to the minimum! Such minima can lead to convergence 
failures, which are usually fixed by the second phase. 

See \code{\link{glm.control}} for more information. 
}
}



\references{


Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428


\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}


\seealso{   \code{\link{gam}} \code{\link{gam.fit}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





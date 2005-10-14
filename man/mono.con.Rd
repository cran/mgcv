\name{mono.con}
\alias{mono.con}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Monotonicity constraints for a cubic regression spline}
\description{ Finds linear constraints sufficient for monotonicity (and
  optionally upper and/or lower boundedness) of a cubic regression
  spline. The basis representation assumed is that given by the
  \code{gam}, \code{"cr"} basis: that is the spline has a set of knots,
  which have fixed x values, but the y values of which constitute the
  parameters of the spline. }
\usage{
mono.con(x,up=TRUE,lower=NA,upper=NA)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{The array of knot locations.}
 \item{up}{If \code{TRUE} then the constraints imply increase, if
   \code{FALSE} then decrease. }
 \item{lower}{This specifies the lower bound on the spline unless it is
   \code{NA} in which case no lower bound is imposed.}
 \item{upper}{This specifies the upper bound on the spline unless it is
   \code{NA} in which case no upper bound is imposed.}
}
\details{ Consider the natural cubic spline passing through the points:
  \eqn{ \{x_i,p_i:i=1 \ldots n \} }{ (x_i,p_i), i=1..n}. Then it is possible
  to find a relatively small set of linear constraints on \eqn{ \bf p}{p}
  sufficient to ensure monotonicity (and bounds if required):
  \eqn{ {\bf Ap}\ge{\bf b} }{Ap>=b}. Details are given in Wood (1994).
  This function returns a list containing \code{A} and \code{b}. 
}
\value{ The function returns a list containing constraint matrix
  \code{A} and constraint vector \code{b}.   
}
\references{

Gill, P.E., Murray, W. and Wright, M.H. (1981) Practical Optimization. Academic
Press, London. 

Wood, S.N. (1994) Monotonic smoothing splines fitted by cross validation SIAM
Journal on Scientific Computing 15(5):1126-1133


\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org} } 

\seealso{  \code{\link{mgcv} } \code{\link{pcls}}  }

\examples{
## see ?pcls
}
\keyword{models} \keyword{smooth} \keyword{regression} %-- one or more ..







\name{null.space.dimension}
\alias{null.space.dimension}
\alias{null.space.basis.powers}
\alias{null.space.basis.labels}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{The basis of the space of un-penalized functions for a t.p.r.s.}
\description{ The thin plate spline penalties give zero penalty to some
  functions. The space of these functions is spanned by a set of
  polynomial terms. \code{null.space.dimension} finds the dimension of this space, \eqn{M}{M}, given
  the number of covariates that the smoother is a function of, \eqn{d}{d},
  and the order of the smoothing penalty, \eqn{m}{m}. If \eqn{m}{m} does not
  satisfy \eqn{2m>d}{2m>d} then the smallest possible dimension
  for the null space is found given \eqn{d}{d} and the requirement that
  the smooth should be visually smooth. \code{null.space.basis.powers()}
  produces an \code{M} by \code{d} array, the ith row of which gives the
  powers of each of the d covariates featuring in the ith basis
  term. \code{null.space.basis.names()} gives the equivalent information
  as a list of basis labels - the key feature of which is that the same
  term will always have the same label, irrespective of the order in
  which variables are supplied to the routine. For some models a smooth
  may be multiplied by a \code{by} variable: in such cases the basis
  labels can include this name as well, but note that the constant label
  will be returned as if the \code{by} variable didn't exist (to meet
  requirements of function \code{gam.side.conditions}).
 
}
\usage{
null.space.dimension(d,m)
null.space.basis.powers(m,d)
null.space.basis.labels(names,m,by="NA")
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{d}{ is a positive integer - the number of variables of which the
   t.p.s. is a function. }
 \item{m}{ a non-negative integer giving the order of the penalty
   functional, or signalling that the default order should be used.}
 \item{names}{is an array of \code{d} variable names}
 \item{by}{the name of any \code{by} variable, which mutliplies all the
   terms in the penalty null space basis.}
}
\details{ Thin plate splines are only visually smooth if the order of the
  wiggliness penalty, \eqn{m}{m}, satisfies \eqn{2m > d+1}{2m >
    d+1}. If \eqn{2m<d+1}{2m<d+1} then this routine finds the smallest
  \eqn{m}{m} giving visual smoothness
  for the given \eqn{d}{d}, otherwise the supplied \eqn{m}{m} is used. The null space dimension is given by:

  \eqn{M=(m+d+1)!/(d!(m-d)!)}{M=(m+d+1)!/(d!(m-d)!}

  which is the value returned.

  
}
\value{
\code{null.space.dimension} returns an integer (array), the null space dimension
\eqn{M}{M}: this is the only one of these functions that accepts array arguments. \cr  \code{null.space.basis.powers} produces a 2-d array,
\code{p}, say: each row corresponds to one of the \code{M} basis vectors, while each column
corresponds to one covariate. So if \code{x1}, \code{x2}, \code{x3}, etc
are the covariates, the ith basis vector is given by
\code{x1^p[i,1]*x2^p[i,2]*x3^p[i,3]...}. \cr
\code{null.space.basis.labels()} just produces uniques labels for the
basis terms.    
 
}

\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\references{
Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

\url{http://www.stats.gla.ac.uk/~simon/}
}

\section{WARNING}
{ This routine operates symbolically, not numerically: your model may
  still be numerically un-identifiable, and the routine can't tell if
  two variables with different names are really the same.
}  

\seealso{\code{\link{gam.side.conditions}}, \code{\link{gam}}}

\examples{
null.space.dimension(2,0)
null.space.basis.powers(2,2)
null.space.basis.labels(c("x","z"),m=2)
null.space.basis.labels(c("z","x"),m=2) # labels are invariant
}
\keyword{models} \keyword{regression}%-- one or more ..



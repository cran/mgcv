\name{print.gam}
\alias{print.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Print a Generalized Additive Model object.}
\description{ The default print method for a \code{gam} object. 
}
\usage{
\method{print}{gam}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x, ...}{ fitted model objects of class \code{gam} as produced by \code{gam()}.}
}
\details{ Prints out the family, model formula, effective degrees of freedom for each smooth term, and optimized 
value of the smoothness selection criterion used. See \code{\link{gamObject}} (or \code{names(x)}) for a listing 
of what the object contains. \code{\link{summary.gam}} provides more detail. 

Note that the optimized smoothing parameter selection criterion reported is one of GCV, UBRE(AIC), GACV, 
negative log marginal likelihood (ML), or negative log restricted likelihood (REML).

If rank deficiency of the model was detected then the apparent rank is reported, along with the length of the cofficient vector 
(rank in absense of rank deficieny). Rank deficiency occurs when not all coefficients are identifiable given the data. Although 
the fitting routines (except \code{gamm}) deal gracefully with rank deficiency, interpretation of rank deficient models may be difficult. 

}



\references{

Wood, S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). CRC/
Chapmand and Hall, Boca Raton, Florida.

\url{http://www.maths.bris.ac.uk/~sw15190/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org} }



\seealso{   \code{\link{gam}}, \code{\link{summary.gam}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










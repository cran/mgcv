\name{gamm.setup}
\alias{gamm.setup}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Mixed Model set up.}
\description{ This is an internal function of package \code{mgcv}. It is called 
by \code{gamm} to set up a generalized additive mixed model in a form suitable for fitting by calls to 
\code{glmmPQL} from the \code{MASS} library or \code{lme} from the \code{nlme} library. The main task is the representation of the smooth terms as random effects. 
}

\references{

Wood, S.N. (2004a) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. Journal of the American Statistical Association.

Wood, S.N. (2004b) Low rank scale invariant tensor product smooths for
Generalized Additive Mixed Models. Technical Report of the Department of
Statistics, University of Glasgow, UK.


\url{http://www.stats.gla.ac.uk/~simon/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}



\seealso{   \code{\link{gamm}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









\name{gamm.setup}
\alias{gamm.setup}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized additive mixed model set up}
\description{ This is an internal function of package \code{mgcv}. It is called 
by \code{gamm} to set up a generalized additive mixed model in a form suitable for fitting by calls to 
\code{gammPQL} or \code{lme} from the \code{nlme} library. The main task is the representation 
of the smooth terms as random effects. The routine first calls \code{\link{gam.setup}}, then checks that
a model estimable by \code{\link{gamm}} has been specified in the model formula, and then transforms the 
problem to a form suitable for mixed model estimation.
}

\references{

Wood, S.N. (2004) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. Journal of the American Statistical Association. 99:673-686


Wood, S.N. (2006a) Low rank scale invariant tensor product smooths for
generalized additive mixed models. Biometrics 62(4):1025-1036

Wood S.N. (2006b) Generalized Additive Models: An Introduction with R. Chapman
and Hall/CRC Press.


\url{http://www.maths.bath.ac.uk/~sw283/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}



\seealso{   \code{\link{gamm}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









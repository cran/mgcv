\name{gam.nbut}
\alias{gam.nbut}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Models using Negative Binomial errors with unknown theta}
\description{ This is an internal function of package \code{mgcv}. It is called by \code{gam()} 
for data where the error distribution is known to follow the Negative Binomial, but the 
value of theta is unknown. Based on \code{glm.nb()} from the MASS library, it calls
\code{gam.fit()} repeatedly while iteratively improving its estimate of
theta, starting from a poisson error model and uses the \code{neg.binom}
error family (based again on code from the MASS library).
}

\value{ returns a \code{gam} object. See \code{\link{gam}} for details.}

\author{ Mike Lonergan \email{mel@st-and.ac.uk} }

\seealso{  \code{\link{gam}} \code{\link{neg.binom}} \code{\link{gam.nbut} }}

\keyword{models}
\keyword{smooth}
\keyword{regression}
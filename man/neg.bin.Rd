% based on MASS/negative.binomial
%
\name{neg.bin}
\alias{neg.bin}
\title{
Family function for Negative Binomial GAMs 
}
\description{
Specifies the information required to fit a Negative Binomial Generalized
Additive Model, with known or unknown \code{theta} parameter, using
\code{gam()} and \code{gam.nbut()}. The code is a minor modification of
 \code{negative.binomial} from \code{library(MASS)}. } 
\usage{
neg.bin(theta=NA, link="log")
}
\arguments{
\item{theta}{
The known value of the additional parameter, \code{theta}, or \code{NA}
if this is unknown.
}
\item{link}{
The link function.  Currently must be one of \code{log}, \code{sqrt} or \code{identity}.
}}
\value{
A list of functions and expressions needed by \code{gam()} to fit a
generalised additive mode with Negative Binomial errors.
}
\seealso{
\code{\link{gam}}, \code{\link{gam.nbut}}
}

\keyword{regression}
\keyword{models}
\keyword{smooth}

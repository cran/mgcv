\name{full.score}
\alias{full.score}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GCV/UBRE score for use withn nlm}
\description{ \code{gam} allows the option to `polish' smoothing parameter 
estimates by minimising the GCV/UBRE score calculated at convergence of the IRLS 
algorithm, given a set of smoothing parameters, rather than using the `performance 
iteration' method which estimates smoothing parameters within the IRLS loop. The 
estimates are often slightly different, since the power iteration effectively 
neglects the depedence of the iterative weights on the smoothing parameters.

The `polishing' optimization is fairly crude and numerically costly. 
\code{\link{nlm}} is used to minimize the scores with respect to the smoothin 
parameters, and this routine is designed to be supplied to it as an argument. 

This is basically a service routine for \code{\link{gam}}, and is not usually 
called directly by users.
}
\usage{
full.score(sp,G,family,control,gamma)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{sp}{The logs of the smoothing parameters}

\item{G}{a list reurned by \code{\link{GAMsetup}}}

\item{family}{The family object for the GAM.}

\item{control}{a list returned be \code{\link{gam.control}}}

\item{gamma}{the degrees of freedom inflation factor (usually 1).}
}
\value{ The value of the GCV/UBRE score, with attribute \code{"full.gam.object"}
which is the full object returned by \code{\link{gam.fit}}.
}
 
}

\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










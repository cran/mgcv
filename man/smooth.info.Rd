\name{smooth.info}
\alias{smooth.info}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generic function to provide extra information about smooth specification}
\description{Takes a smooth specification object and adds extra basis specific information to it before smooth constructor called. Default method returns supplied object unmodified.
}

\usage{
smooth.info(object)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{object}{ is a smooth specification object} 
}

\value{ A smooth specification object, which may be modified in some way.}


\details{Sometimes it is necessary to know something about a smoother before it is constructed, beyond what is in the initial smooth specification object.
For example, some smooth terms could be set up as tensor product smooths and it is useful for \code{\link{bam}} to take advantage of this when discrete covariate methods are used. However, \code{\link{bam}} needs to know whether a smoother falls into this category before it is constructed in order to discretize its covariates marginally instead of jointly. Rather than \code{\link{bam}} having a hard coded list of such smooth classes it is preferable for the smooth specification object to report this themselves. \code{smooth.info} method functions are the means for achieving this. When interpreting a gam formula the \code{smooth.info} function is applied to each smooth specification object as soon as it is produced (in \code{interpret.gam0}). 
}


\references{ 

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman
and Hall/CRC Press.

}

\author{Simon N. Wood \email{simon.wood@r-project.org}}

\seealso{ \code{\link{bam}},
\code{\link{smooth.construct}}, \code{\link{PredictMat}} }

\examples{# See smooth.construct examples
spec <- s(a,bs="re")
class(spec)
spec$tensor.possible
spec <- smooth.info(spec)
spec$tensor.possible
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...


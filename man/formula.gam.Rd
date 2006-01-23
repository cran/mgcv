\name{formula.gam}
\alias{formula.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract the formula from a gam object}
\description{ Extracts the formula from a fitted \code{gam} object. 
}
\usage{
formula.gam(x,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x}{ fitted model objects of class \code{gam} (see \code{\link{gamObject}}) as produced by \code{gam()}.}
\item{...}{un-used in this case} 
}

\details{ Returns \code{x$formula}. Provided so that \code{anova} methods
print an appropriate description of the model.
}

\value{ A model formula.}


\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  \code{\link{gam}}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







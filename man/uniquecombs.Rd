\name{uniquecombs}
\alias{uniquecombs}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{find the unique rows in a matrix }
\description{
This routine returns a matrix containing all the unique rows of the
matrix supplied as its argument. That is, all the duplicate rows are
stripped out. Note that the ordering of the rows on exit is not the same
as on entry.
}
\usage{
uniquecombs(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ is an \R matrix }
}
\details{ Models with more parameters than unique combinations of
  covariates are not identifiable. This routine provides a means of
  evaluating the number of unique combinations of coavariates in a
  model. The routine calls compiled C code.
   
}
\value{
A matrix consisting of the unique rows of \code{x} (in arbitrary order).
 
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\examples{
X<-matrix(c(1,2,3,1,2,3,4,5,6,1,3,2,4,5,6,1,1,1),6,3,byrow=TRUE)
print(X)
uniquecombs(X)
}
\keyword{models} \keyword{regression}%-- one or more ..



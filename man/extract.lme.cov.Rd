\name{extract.lme.cov}
\alias{extract.lme.cov}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Extract the data covariance matrix from an lme object}
\description{ This is a service routine for \code{\link{gamm}}. It extracts 
the estimated covariance matrix of the data from an \code{lme} object, allowing the 
user control about which levels of random effects to include in this 
calculation. 
}
\usage{
extract.lme.cov(b,data,start.level=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{b}{ A fitted model object returned by a call to \code{\link{lme}}}

\item{data}{ The data frame/ model frame that was supplied to
 \code{\link{lme}}.}

\item{start.level}{The level of nesting at which to start including random 
effects in the calculation. This is used to allow smooth terms to be estimated
as random effects, but treated like fixed effects for variance calculations.}

}

\details{ The random effects, correlation structure and variance structure used
for a linear mixed model combine to imply a covariance matrix for the 
response data being modelled. This routine extracts that covariance matrix.
The process is slightly complicated, because different components of the 
fitted model object are stored in different orders (see function code for 
details!).  

The calculation is not optimally efficient, since it forms the full matrix,
which may in fact be sparse. In applications in which the main objective is
to allow non-independent `errors' in GAMs this is unlikely to cause great 
computational losses.
}



\value{ An estimated covariance matrix.}

\references{

For \code{lme} see:

Pinheiro J.C. and Bates, D.M. (2000) Mixed effects Models in S and S-PLUS. Springer

For details of how GAMMs are set up here for estimation using \code{lme} see:
 
Wood, S.N. (manuscript) Tensor product smooth interaction terms in 
Generalized Additive Mixed Models.


\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}


\seealso{  
\code{\link{gamm}}
}

\examples{
library(nlme)
data(Rail)
b <- lme(travel~1,Rail,~1|Rail)
extract.lme.cov(b,Rail)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..






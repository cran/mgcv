\name{SANtest}
\alias{SANtest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Example of simple additive GAM using penalized regression splines.}
\description{ Provides an example of the use of \code{mgcv()} and \code{GAMsetup()} for 
simple GAM modelling with penalized regression splines. Truth is simulated as the sum of 
3 univariate functions of 3 covariates: normal errors are added to produce simulated 
data. A GAM is then fitted using \code{mgcv()} (including one spurious covariate). 
Scatter plots of data against the 4 covariates and plots of true and reconstructed 
functions are produced. 

Of little interest for practical analysis (see \code{gam} instead).  
}
\usage{
SANtest(n=100, sig2=-1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{n}{ \code{n} is the number of data to simulate. }
 \item{sig2}{ The magnitude of \code{sig2} gives the error variance to use for 
simulation of data. If \code{sig2} is greater than zero the UBRE is used with the error 
variance, \eqn{\sigma^2}{sigma^2}, given by  \code{sig2}, otherwise GCV is used.}
}
\value{
 An estimate of the error variance used for simulation (GCV), or the error variance 
(UBRE).
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk} }


\seealso{ \code{\link{mgcv}} \code{\link{gam}} }

\examples{
set.seed(0)
SANtest(300,-1) # use GCV with 300 data point example, variance 1
SANtest(300,0.5) # use UBRE 300 data, variance 0.5

}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









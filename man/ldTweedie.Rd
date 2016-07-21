\name{ldTweedie}
\alias{ldTweedie}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Log Tweedie density evaluation}
\description{A function to evaluate the log of the Tweedie density for variance powers between 1 and 2, inclusive.
Also evaluates first and second derivatives of log density w.r.t. its scale parameter, \code{phi}, and \code{p}, 
or w.r.t. \code{rho=log(phi)} and \code{theta} where \code{p = (a+b*exp(theta))/(1+exp(theta))}.
}

\usage{
ldTweedie(y,mu=y,p=1.5,phi=1,rho=NA,theta=NA,a=1.001,b=1.999,all.derivs=FALSE)
}
\arguments{
\item{y}{values at which to evaluate density.}

\item{mu}{corresponding means (either of same length as \code{y} or a single value).}

\item{p}{the variance of \code{y} is proportional to its mean to the power \code{p}. \code{p} must
be between 1 and 2. 1 is Poisson like (exactly Poisson if \code{phi=1}), 2 is gamma. }

\item{phi}{The scale parameter. Variance of \code{y} is \code{phi*mu^p}.}

\item{rho}{optional log scale parameter. Over-rides \code{phi} if \code{theta} also supplied.}

\item{theta}{parameter such that  \code{p = (a+b*exp(theta))/(1+exp(theta))}. Over-rides \code{p} if \code{rho} 
also supplied.}

\item{a}{lower limit parameter (>1) used in definition of \code{p} from \code{theta}.}

\item{b}{upper limit parameter (<2) used in definition of \code{p} from \code{theta}.}

\item{all.derivs}{if \code{TRUE} then derivatives w.r.t. \code{mu} are also returned. Only available with \code{rho} and
\code{phi} parameterization.}
}
\value{ A matrix with 6 columns, or 10 if \code{all.derivs=TRUE}. The first is the log density of \code{y} (log probability if \code{p=1}). 
The second and third are the first and second derivatives of the log density w.r.t. \code{phi}. 4th and 5th 
columns are first and second derivative w.r.t. \code{p}, final column is second derivative w.r.t. \code{phi} and \code{p}.

If \code{rho} and \code{theta} were supplied then derivatives are w.r.t. these. In this case, and if \code{all.derivs=TRUE} then the 7th colmn is the derivative w.r.t. \code{mu}, the 8th is the 2nd derivative w.r.t. \code{mu}, the 9th is the mixed derivative w.r.t. \code{theta} and\code{mu} and the 10th is the mixed derivative w.r.t. \code{rho} and \code{mu}.
}

\details{ A Tweedie random variable with 1<p<2 is a sum of \code{N} gamma random variables 
where \code{N} has a Poisson distribution. The p=1 case is a generalization of a Poisson distribution and is a discrete 
distribution supported on integer multiples of the scale parameter. For 1<p<2 the distribution is supported on the 
positive reals with a point mass at zero. p=2 is a gamma distribution. As p gets very close to 1 the continuous 
distribution begins to converge on the discretely supported limit at p=1. 

\code{ldTweedie} is based on the series evaluation method of Dunn and Smyth (2005). Without 
the restriction on \code{p} the calculation of Tweedie densities is less straightforward. If you really need this 
case then the \code{tweedie} package is the place to start. 

The \code{rho}, \code{theta} parameterization is useful for optimization of \code{p} and \code{phi}, in order to keep \code{p}
bounded well away from 1 and 2, and \code{phi} positive. The derivatives near \code{p=1} tend to infinity.

Note that if \code{p} and \code{phi} (or \code{theta} and \code{rho}) both contain only a single unique value, then the underlying
code is able to use buffering to avoid repeated calls to expensive log gamma, di-gamma and tri-gamma functions (\code{mu} can still be a vector of different values). This is much faster than is possible when these parameters are vectors with different values. 

}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\references{
Dunn, P.K. and G.K. Smith (2005) Series evaluation of Tweedie exponential dispersion model densities. 
Statistics and Computing 15:267-280

Tweedie, M. C. K. (1984). An index which distinguishes between
     some important exponential families. Statistics: Applications and
     New Directions. Proceedings of the Indian Statistical Institute
     Golden Jubilee International Conference (Eds. J. K. Ghosh and J.
     Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
}


\examples{
  library(mgcv)
  ## convergence to Poisson illustrated
  ## notice how p>1.1 is OK
  y <- seq(1e-10,10,length=1000)
  p <- c(1.0001,1.001,1.01,1.1,1.2,1.5,1.8,2)
  phi <- .5
  fy <- exp(ldTweedie(y,mu=2,p=p[1],phi=phi)[,1])
  plot(y,fy,type="l",ylim=c(0,3),main="Tweedie density as p changes")
  for (i in 2:length(p)) {
    fy <- exp(ldTweedie(y,mu=2,p=p[i],phi=phi)[,1])
    lines(y,fy,col=i)
  }


}
\keyword{models} \keyword{regression}%-- one or more ..



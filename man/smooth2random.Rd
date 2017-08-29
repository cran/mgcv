\name{smooth2random}
\alias{smooth2random}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convert a smooth to a form suitable for estimating as random effect}
\description{A generic function for converting \code{mgcv} smooth objects to forms suitable for estimation as random effects by e.g. \code{lme}. Exported mostly for use by other package developers. 
}
\usage{
smooth2random(object,vnames,type=1)

%\method{summary}{gam}(object, dispersion=NULL, freq=FALSE, p.type = 0, ...)

}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{an \code{mgcv} smooth object.}

\item{vnames}{a vector of names to avoid as dummy variable names in the random effects form.} 

\item{type}{\code{1} for \code{lme}, otherwise \code{lmer}.}

}

\details{There is a duality between smooths and random effects which means that smooths can be estimated using mixed modelling software. This function converts standard \code{mgcv} smooth objects to forms suitable for estimation by \code{lme}, for example. A service routine for \code{\link{gamm}} exported for use by package developers. 
}

\value{A list.

\item{rand}{ a list of random effects, including grouping factors, and 
       a fixed effects matrix. Grouping factors, model matrix and             model matrix name attached as attributes, to each element.}

\item{trans.D}{A vector, trans.D, that transforms coefs, in order [rand1, rand2,... fix] back to original parameterization. If null, then taken as vector of ones. \code{b.original = trans.U \%*\% (trans.D*b.fit)}.}

\item{trans.U}{A matrix, trans.U, that transforms coefs, in order [rand1, rand2,... fix] back to original parameterization. If null, then not needed. If null then taken as identity.}

\item{Xf}{A matrix for the fixed effects, if any.}

\item{fixed}{\code{TRUE/FALSE}, indicating if term was unpenalized or not. If unpenalized then other stuff may not be returned (it's not a random effect).}

\item{rind}{an index vector such that if br is the vector of 
random coefficients for the term, br[rind] is the coefs in 
order for this term. }

\item{pen.ind}{index of which penalty penalizes each coefficient: 0 for unpenalized.}
}


\references{

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman
and Hall/CRC Press.

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}.}


\seealso{  \code{\link{gamm}} }

\examples{
library(mgcv)
x <- runif(30)
sm <- smoothCon(s(x),data.frame(x=x))[[1]]
smooth2random(sm,"")

}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







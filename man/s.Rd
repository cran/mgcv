\name{s}
\alias{s}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Defining smooths in GAM formulae}
\description{ Function used in definition of smooth terms within
  \code{gam} model formulae. The function does not evaluate a (spline)
  smooth - it exists purely to help set up a model using spline based smooths.
}

\usage{s(..., k=-1,fx=FALSE,bs="tp",m=0,by=NA)}
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{...}{ a list of variables that are the covariates that this
    smooth is a function of.}
  \item{k}{ the dimension of the basis used to represent the smooth term.
    The default depends on the number of variables that the smooth is a
    function of. \code{k} should not be less than the dimension of the
    null space of the penalty for the term (see
    \code{\link{null.space.dimension}}), but will be reset if
  it is.}
  \item{fx}{indicates whether the term is a fixed d.f. regression
  spline (\code{TRUE}) or a penalized regression spline (\code{FALSE}).}
  \item{bs}{this can be \code{"cr"} for a cubic regression spline or
    \code{"tp"} for a thin plate regression spline. Only thin plate
    regression splines can be used for multidimensional smooths, so this
    argument only has an effect for univariate smooths. Note that the
    \code{"cr"} basis is faster than the \code{"tp"} basis, particularly
  on large data sets.}
  \item{m}{The order of the penalty for this t.p.r.s. term (e.g. 2 for
    normal cubic spline penalty with 2nd derivatives). O signals
    autoinitialization, which sets the order to the lowest value
    satisfying 2m>d+1, where d is the number of covariates: this choise
    ensures visual smoothness. In addition, m must satisfy the technical
    restriction 2m>d, otherwise it will be autoinitialized.}
  \item{by}{specifies a covariate by which the whole smooth term is to
    be multiplied. This is particularly useful for creating models in
    which a smooth interacts with a factor: in this case the \code{by}
    variable would usually be the dummy variable coding one level of the
  factor. See the examples below.}
  }
}
\details{The function does not evaluate the variable arguments. It will
  correctly interpret calls like \code{s(x,14|f)} (results in pure
  regression spline smooth of one variable with a 14 knot cubic
  regression spline basis), \code{s(x,z,20)} (a penalized regression
  spline of 2 covariates using a 20 dimensional t.p.r.s. basis),
  etc. but this feature is purely for back compatibility reasons, and
  may not be maintained indefinitely. 

}
\value{ A list with the following items:
  \item{term}{An array of text strings giving the names of the covariates that 
the term is a function of.}
  \item{bs.dim}{The dimension of the basis used to represent the smooth.}
  \item{bs.type}{The type of basis. 0 is cubic regression spline. 1 is
    thin plate regression spline. 0 can only be used for 1-d smooths.}
  \item{fixed}{TRUE if the term is to be treated as a pure regression
    spline (with fixed degrees of freedom); FALSE if it is to be treated
  as a penalized regression spline}
  \item{dim}{The dimension of the smoother - i.e. the number of
    covariates that it is a function of.}
  \item{p.order}{The order of the t.p.r.s. penalty, or 0 for
    auto-selection of the penalty order.}
  \item{by}{is the name of any \code{by} variable as text ("NA" for none).}
  \item{full.call}{Text for pasting into a string to be converted to a
    gam formula, which has the values of function options given explicitly -
  this is useful for constructing a fully expanded gam formula which can
be used without needing access to any variables that may have been used
to define k, fx, bs or m in the original call. i.e. this is text which
when parsed and evaluated generates a call to \code{s()} with all the
options spelled out explicitly.}
}


\author{ Simon N. Wood \email{snw@st-and.ac.uk}}


\seealso{   \code{\link{gam}}}

\examples{
# example utilising `by' variables
library(mgcv)
set.seed(0)
n<-200;sig2<-4
x1 <- runif(n, 0, 1);x2 <- runif(n, 0, 1);x3 <- runif(n, 0, 1)
fac<-c(rep(1,n/2),rep(2,n/2)) # create factor
fac.1<-rep(0,n)+(fac==1);fac.2<-1-fac.1 # and dummy variables
fac<-as.factor(fac)
f1 <-  exp(2 * x1) - 3.75887
f2 <-  0.2 * x1^11 * (10 * (1 - x1))^6 + 10 * (10 * x1)^3 * (1 - x1)^10 - 1.396
f<-f1*fac.1+f2*fac.2+x2
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- f + e
# NOTE: smooths will be centered, so need to include fac in model....
b<-gam(y~fac+s(x1,by=fac.1)+s(x1,by=fac.2)+x2) 
plot(b,pages=1)
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










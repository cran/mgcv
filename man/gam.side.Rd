\name{gam.side}
\alias{gam.side}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Identifiability side conditions for a GAM}
\description{ GAM formulae with repeated variables only correspond to
  identifiable models given some side conditions. This routine works 
  out appropriate side conditions, based on zeroing redundant parameters.
  It is called from \code{\link{gam.setup}} and is not intended to be called by users. 

  The method identifies nested and repeated variables by their names, but
  evaluates what constraints to impose, numerically. Constraints are always
  applied to smooths of more variables in preference to smooths of fewer
  variables. The numerical approach allows appropriate constraints to be
  applied to models constructed using any smooths, including user defined smooths.
}
\usage{
gam.side(sm,tol=.Machine$double.eps^.5)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{sm}{ A list of smooth objects as returned by
  \code{\link{smooth.construct}}.}
  \item{tol}{The tolernce to use when assessing linear dependence of smooths.}
}
\details{ Models such as  \code{y~s(x)+s(z)+s(x,z)} can be estimated by
  \code{\link{gam}}, but require identifiability constraints to be applied, to
  make them identifiable. This routine does this, effectively setting redundant parameters
  to zero. When the redundancy is between smooths of lower and higher numbers
  of variables, the constraint is always applied to the smooth of the higher
  number of variables. 

  Dependent smooths are identified symbolically, but which constraints are 
  needed to ensure identifiability of these smooths is determined numerically, using
  \code{\link{fixDependence}}. This makes the routine rather general, and not
  dependent on any particular basis.

}
\value{ A list of smooths, with model matrices and penalty matrices adjusted
to automatically impose the required constraints. Any smooth that has been
modified will have an attribute \code{"del.index"}, listing the columns of its
model matrix that were deleted. This index is used in the creation of
prediction matrices for the term.
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}



\examples{
set.seed(0)
n<-400
sig2<-4
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
pi <- asin(1) * 2
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f + 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10 - 1.396
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- f + e
b<-gam(y~s(x0)+s(x1)+s(x0,x1)+s(x2))
plot(b,pages=1)
test1<-function(x,z,sx=0.3,sz=0.4)  
{ (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
  0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
}
n<-500
old.par<-par(mfrow=c(2,2))
x<-runif(n);z<-runif(n);
y<-test1(x,z)+rnorm(n)*0.1
## a fully nested tensor product example
b<-gam(y~s(x,bs="cr",k=6)+s(z,bs="cr",k=6)+te(x,z,k=6))
plot(b)
par(old.par)
rm(list=c("f","x0","x1","x2","x","z","y","b","test1","n","sig2","pi","e"))
}
\keyword{models} \keyword{regression}%-- one or more ..



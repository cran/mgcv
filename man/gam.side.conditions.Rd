\name{gam.side.conditions}
\alias{gam.side.conditions}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Identifiability side conditions for a GAM.}
\description{ GAM formulae with repeated variables only correspond to
  identifiable models given some side conditions on the null space of
  the model penalties. This routine works out appropriate side conditions
  and returns the coefficient matrix that imposes them. It is called
  from \code{\link{gam}} and is not intended to be called by users.
}
\usage{
gam.side.conditions(G)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{G}{ the object obtained by calling \code{\link{gam.setup}} followed by
    \code{\link{GAMsetup}} }
}
\details{ Part of the basis of a thin plate regression spline (of the
  sort used to define GAMs in package \code{mgcv}) is a polynomial basis
  for the space of functions of zero wiggliness according to the spline
  wiggliness penalty. A GAM with repeated dependence on one or more
  variables will have multiple copies of these basis terms, and will
  hence not be identifiable without extra constraints. For example, the
  design matrix for the model 
  \code{y~s(x)+s(z)+s(x,z)} will feature two copies of the columns for
  \code{x} and \code{z}, but by imposing the constraint that the second
  parameter for each repeat variable is zero, the model becomes
  identifiable again.

  This routine automates the process of producing these constraints in
  general. The method used is to create labels for each null basis
  column that uniquely identify the term while allowing multiple copies
  of the column to be identified. The algorithm works upwards throught the
  dimensions - i.e. starts with 1-d terms, then does 2-d terms , then
  3-d and so on. A stack of null basis columns is maintained, and each
  new basis term tested against it - if the new term is a repeat then
  its parameter is constrained to be zero - otherwise it is added to
  the stack. The method is symbolic and not numerical - you
  can fool it by making two copies of the same covariate with different names.

  If the model can not be made identifiable then an error is generated.

  The code critically depends on equivalence between the R code
  null.space.basis.powers and the equivalent C code in tprs.c.

  Note that the routine does nothing to the constant terms in the basis
  - these are already dealt with by the constraints that centre all
  smooths. 

}
\value{
A matrix \eqn{ \bf C}{C} defining the constraints on the model. The constraints
are that \eqn{ { \bf Cp}={ \bf 0}}{Cp=0} where \eqn{ \bf p}{p} is the
parameter vector. \code{C} includes the original sum to zero constraints from
  \code{G$C} and is hence suitable for directly replacing that matrix. 
}

\author{ Simon N. Wood \email{snw@st-and.ac.uk}}

\references{
Wood (2003) Thin plate regression splines JRSSB in press. 
}


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
b<-gam(y~s(x)+s(z)+s(x,z))
plot(b)
par(old.par)
rm(list=c("f","x0","x1","x2","x","z","y","b","test1","n","sig2","pi","e"))
}
\keyword{models} \keyword{regression}%-- one or more ..



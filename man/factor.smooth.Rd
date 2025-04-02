\name{factor.smooth}
\alias{factor.smooth.interaction}
\alias{factor.smooth}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Factor smooth interactions in GAMs}

\description{The interaction of one or more factors with a smooth effect, produces a separate smooth for each factor level. These smooths can have different smoothing parameters, or all have the same smoothing parameter. There are several vays to set them up.
\describe{
\item{Factor \code{by} variables.}{If the \code{by} variables for a smooth (specified using \code{\link{s}}, \code{\link{te}}, \code{\link{ti}} or \code{\link{t2}}) is a factor, then a separate smooth is produced for each factor level. If the factor is ordered, then no smooth is produced for its first level: this is useful for setting up models which have a reference level smooth and then difference to reference smooths for each factor level except the first (which is the reference). Giving the smooth an \code{id} forces the same smoothing parameter to be used for all levels of the factor. For example \code{s(x,by=fac,id=1)} would produce a separate smooth of \code{x} for each level of \code{fac}, with each smooth having the same smoothing parameter. See \link{gam.models} for more.}

\item{Sum to zero smooth interactions}{\code{bs="sz"} These factor smooth interactions are specified using \code{s(...,bs="sz")}. There may be several factors supplied, and a smooth is produced for each combination of factor levels. The smooths are constructed to exclude the `main effect' smooth, or the effects of individual smooths produced for lower order combinations of factor levels. For example, with a single factor, the smooths for the different factor levels are so constrained that the sum over all factor levels of equivalent spline coefficients are all zero. This allows the meaningful and identifiable construction of models with a main effect smooth plus smooths for the difference between each factor level and the main effect. Such a construction is often more natural than the \code{by} variable with ordered factors construction. See \code{\link{smooth.construct.sz.smooth.spec}}.}

\item{Random wiggly curves}{\code{bs="fs"} This approach produces a smooth curve for each level of a single factor, treating the curves as entirely random. This means that in principle a model can be constructed with a main effect plus factor level smooth deviations from that effect. However the model is not forced to make the main effect do as much of the work as possible, in the way that the \code{"sz"} approach does. This approach can be very efficient with \code{\link{gamm}} as it exploits the nested estimation available in \code{lme}. See \code{\link{smooth.construct.fs.smooth.spec}}.
}
}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org} with input from Matteo Fasiolo.}

\seealso{\code{\link{smooth.construct.fs.smooth.spec}}, \code{\link{smooth.construct.sz.smooth.spec}}}

\examples{
library(mgcv)
set.seed(0)
## simulate data...
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x,a=2,b=-1) exp(a * x)+b
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
n <- 500;nf <- 25
fac <- sample(1:nf,n,replace=TRUE)
x0 <- runif(n);x1 <- runif(n);x2 <- runif(n)
a <- rnorm(nf)*.2 + 2;b <- rnorm(nf)*.5
f <- f0(x0) + f1(x1,a[fac],b[fac]) + f2(x2)
fac <- factor(fac)
y <- f + rnorm(n)*2
## so response depends on global smooths of x0 and 
## x2, and a smooth of x1 for each level of fac.

## fit model...
bm <- gamm(y~s(x0)+ s(x1,fac,bs="fs",k=5)+s(x2,k=20))
plot(bm$gam,pages=1)
summary(bm$gam)

bd <- bam(y~s(x0)+ s(x1) + s(x1,fac,bs="sz",k=5)+s(x2,k=20),discrete=TRUE)
plot(bd,pages=1)
summary(bd)



## Could also use...
## b <- gam(y~s(x0)+ s(x1,fac,bs="fs",k=5)+s(x2,k=20),method="ML")
## ... but its slower (increasingly so with increasing nf)
## b <- gam(y~s(x0)+ t2(x1,fac,bs=c("tp","re"),k=5,full=TRUE)+
##        s(x2,k=20),method="ML"))
## ... is exactly equivalent. 
}
\keyword{models} \keyword{regression}%-- one or more ..



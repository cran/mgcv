\name{choose.k}
\alias{choose.k}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Basis dimension choice for smooths}
\description{Choosing the basis dimension, and checking the choice, when using
penalized regression smoothers.

Penalized regression smoothers gain computational efficiency by virtue of
being defined using a basis of relatively modest size, \code{k}. When setting
up models in the \code{mgcv} package, using \code{\link{s}} or \code{\link{te}}
terms in a model formula, \code{k} must be chosen. 

In practice \code{k-1} sets the upper limit on the degrees of freedom
associated with an \code{\link{s}} smooth (1 degree of freedom is lost to the identifiability
constraint on the smooth). For \code{\link{te}} smooths the upper limit of the
degrees of freedom is given by the product of the \code{k} values provided for
each marginal smooth less one, for the constraint. However the actual
effective degrees of freedom are controlled by the degree of penalization
selected during fitting, by GCV, AIC or whatever is specified. The exception
to this is if a smooth is specified using the \code{fx=TRUE} option, in which
case it is unpenalized.  

So, exact choice of \code{k} is not generally critical: it should be chosen to
be large enough that you are reasonably sure of having enough degrees of
freedom to represent the underlying `truth' reasonably well, but small enough
to maintain reasonable computational efficiency. Clearly `large' and `small'
are dependent on the particular problem being addressed. 

As with all model assumptions, it is useful to be able to check the choice of
\code{k} informally. If the effective degrees of freedom for a model term are
estimated to be much less than \code{k-1} then this is unlikely to be very
worthwhile, but as the EDF approach \code{k-1}, checking can be important. A
useful general purpose approach goes as follows: (i) fit your model and
extract the deviance residuals; (ii) for each smooth term in your model, fit
an equivalent, single, smooth to the residuals, using a substantially
increased \code{k} to see if there is pattern in the residuals that could
potentially be explained by increasing \code{k}. Examples are provided below.

More sophisticated approaches based on partial residuals are also possible.

One scenario that can cause confusion is this: a model is fitted with
\code{k=10} for a smooth term, and the EDF for the term is estimated as 7.6,
some way below the maximum of 9. The model is then refitted with \code{k=20}
and the EDF increases to 8.7 - what is happening - how come the EDF was not
8.7 the first time around? The explanation is that the function space with
\code{k=20} contains a larger subspace of functions with EDF 8.7 than did the
function space with \code{k=10}: one of the functions in this larger subspace
fits the data a little better than did any function in the smaller
subspace. These subtleties seldom have much impact on the statistical
conclusions to be drawn from a model fit, however. 

}


\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\references{

Wood, S.N. (2006) Generalized Additive Models: An Introduction with R. CRC.

\url{http://www.maths.bath.ac.uk/~sw283/}
}

\examples{
## Simulate some data ....
library(mgcv)
set.seed(0) 
n<-400;sig<-2
x0 <- runif(n, 0, 1);x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1);x3 <- runif(n, 0, 1)
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f+0.2*x2^11*(10*(1-x2))^6+10*(10*x2)^3*(1-x2)^10-1.396
e <- rnorm(n, 0, sig)
y <- f + e
## fit a GAM with quite low `k'
b<-gam(y~s(x0,k=6)+s(x1,k=6)+s(x2,k=6)+s(x3,k=6))
plot(b,pages=1)

## check for residual pattern, removeable by increasing `k'
## typically `k', below, chould be substantially larger than 
## the original, `k' but certainly less than n/2.
## Note use of cheap "cs" shrinkage smoothers, and gamma=1.4
## to reduce chance of overfitting...
rsd <- residuals(b)
gam(rsd~s(x0,k=40,bs="cs"),gamma=1.4) ## fine
gam(rsd~s(x1,k=40,bs="cs"),gamma=1.4) ## fine
gam(rsd~s(x2,k=40,bs="cs"),gamma=1.4) ## original `k' too low
gam(rsd~s(x3,k=40,bs="cs"),gamma=1.4) ## fine

## similar example with multi-dimensional smooth
b1 <- gam(y~s(x0)+s(x1,x2,k=15)+s(x3))
rsd <- residuals(b1)
gam(rsd~s(x0,k=40,bs="cs"),gamma=1.4) ## fine
gam(rsd~s(x1,x2,k=100,bs="ts"),gamma=1.4) ## original `k' too low
gam(rsd~s(x3,k=40,bs="cs"),gamma=1.4) ## fine
 
## and a `te' example
b2 <- gam(y~s(x0)+te(x1,x2,k=4)+s(x3))
rsd <- residuals(b2)
gam(rsd~s(x0,k=40,bs="cs"),gamma=1.4) ## fine
gam(rsd~te(x1,x2,k=10,bs="cs"),gamma=1.4) ## original `k' too low
gam(rsd~s(x3,k=40,bs="cs"),gamma=1.4) ## fine

## same approach works with other families in the original model
g<-exp(f/4)
y<-rpois(rep(1,n),g)
bp<-gam(y~s(x0,k=6)+s(x1,k=6)+s(x2,k=6)+s(x3,k=6),family=poisson)
rsd <- residuals(bp)
gam(rsd~s(x0,k=40,bs="cs"),gamma=1.4) ## fine
gam(rsd~s(x1,k=40,bs="cs"),gamma=1.4) ## fine
gam(rsd~s(x2,k=40,bs="cs"),gamma=1.4) ## original `k' too low
gam(rsd~s(x3,k=40,bs="cs"),gamma=1.4) ## fine
 
}

\keyword{models} \keyword{regression}%-- one or more ..


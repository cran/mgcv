\name{identifiability}
\alias{identifiability}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Identifiability constraints}
\description{Smooth terms are generally only identifiable up to an additive constant. In consequence sum-to-zero identifiability constraints are imposed on most smooth terms. The exceptions are terms with \code{by} variables which cause the smooth to be identifiable without constraint (that doesn't include factor \code{by} variables), and random effect terms. Alternatively smooths can be set up to pass through zero at a user specified point.    
}


\details{ By default each smooth term is subject to the sum-to-zero constraint \deqn{\sum_i f(x_i) = 0.}{sum_i f(x_i) = 0.}
The constraint is imposed by reparameterization. The sum-to-zero constraint causes the term to be orthogonal to the intercept: alternative constraints lead to wider confidence bands for the constrained smooth terms.

No constraint is used for random effect terms, since the penalty (random effect covariance matrix) anyway ensures identifiability in this case. Also if a \code{by} variable means that the smooth is anyway identifiable, then no extra constraint is imposed. Constraints are imposed for factor \code{by} variables, so that the main effect of the factor must usually be explicitly added to the model (the example below is an exception). 

Occasionally it is desirable to substitute the constraint that a particular smooth curve should pass through zero at a particular point: the \code{pc} argument to \code{\link{s}}, \code{\link{te}}, \code{\link{ti}} and  \code{\link{t2}} allows this: if specified then such constraints are always applied.

}


%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood (s.wood@r-project.org)
}
\examples{

## Example of three groups, each with a different smooth dependence on x
## but each starting at the same value...
require(mgcv)
set.seed(53)
n <- 100;x <- runif(3*n);z <- runif(3*n)
fac <- factor(rep(c("a","b","c"),each=100))
y <- c(sin(x[1:100]*4),exp(3*x[101:200])/10-.1,exp(-10*(x[201:300]-.5))/
       (1+exp(-10*(x[201:300]-.5)))-0.9933071) + z*(1-z)*5 + rnorm(100)*.4

## 'pc' used to constrain smooths to 0 at x=0...
b <- gam(y~s(x,by=fac,pc=0)+s(z)) 
plot(b,pages=1)
}


\keyword{models} \keyword{regression}%-- one or more ..



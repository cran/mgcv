\name{gam.selection}
\alias{gam.selection}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Model Selection}
\description{ This page is intended to provide some more information on
  how to select GAMs. Given a model structure specified by a gam model formula,
\code{gam()} attempts to find the appropriate smoothness for each applicable model 
term using Generalized Cross Validation (GCV) or an Un-Biased Risk Estimator (UBRE),
the latter being used in cases in which the scale parameter is assumed known. GCV and 
UBRE are covered in Craven and Wahba (1979) and Wahba (1990), see
\code{\link{gam.method}} for more detail about the numerical optimization
approaches available. 


Automatic smoothness selection is unlikely to be successful with few data, particularly with
multiple terms to be selected. In addition GCV and UBRE score can occasionally
display local minima that can trap the minimisation algorithms. GCV/UBRE
scores become constant with changing smoothing parameters at very low or very
high smoothng parameters, and on occasion these `flat' regions can be
separated from regions of lower score by a small `lip'. This seems to be the
most common form of local minimum, but is usually avoidable by avoiding
extreme smoothing parameters as starting values in optimization, and by
avoiding big jumps in smoothing parameters while optimizing. Never the less, if you are
suspicious of smoothing parameter estimates, try changing fit method (see
\code{\link{gam.method}}) and see if the estimates change, or try changing
some or all of the smoothing parameters `manually' (argument \code{sp} of \code{\link{gam}}). 


In general the most logically consistent method to use for deciding which terms to include in the model is to compare GCV/UBRE scores
for models with and without the term. More generally the score for the model with a smooth term can be compared to the score for the model with 
the smooth term replaced by appropriate parametric terms. Candidates for removal can be identified by reference to the approximate p-values provided by \code{summary.gam}. Candidates for replacement by parametric terms are smooth terms with estimated degrees of freedom close to their 
minimum possible.

One appealing approach to model selection is via shrinkage. Smooth classes
\code{cs.smooth} and \code{tprs.smooth} (specified by \code{"cs"} and
\code{"ts"} respectively) have smoothness penalties which include a small
shrinkage component, so that for large enough smoothing parameters the smooth 
becomes identically zero. This allows automatic smoothing parameter selection
methods to effectively remove the term from the model altogether. The
shrinkage component of the penalty is set at a level that usually makes negligable
contribution to the penalization of the model, only becoming effective when
the term is effectively `completely smooth' according to the conventional penalty.


}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\references{
Craven and Wahba (1979) Smoothing Noisy Data with Spline Functions. Numer. Math. 31:377-403

Venables and Ripley (1999) Modern Applied Statistics with S-PLUS

Wahba (1990) Spline Models of Observational Data. SIAM.

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

\url{http://www.maths.bath.ac.uk/~sw283/}
}

\examples{
## an example of GCV based model selection
library(mgcv)
set.seed(0) 
n<-400;sig<-2
x0 <- runif(n, 0, 1);x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1);x3 <- runif(n, 0, 1)
x4 <- runif(n, 0, 1);x5 <- runif(n, 0, 1)
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f+0.2*x2^11*(10*(1-x2))^6+10*(10*x2)^3*(1-x2)^10-1.396
e <- rnorm(n, 0, sig)
y <- f + e
## Note the increased gamma parameter below to favour
## slightly smoother models (See Kim & Gu, 2004, JRSSB)...
b<-gam(y~s(x0,bs="ts")+s(x1,bs="ts")+s(x2,bs="ts")+
   s(x3,bs="ts")+s(x4,bs="ts")+s(x5,bs="ts"),gamma=1.4)
summary(b)
plot(b,pages=1)

}

\keyword{models} \keyword{regression}%-- one or more ..



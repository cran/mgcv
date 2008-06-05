\name{step.gam}
\alias{step.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Alternatives to step.gam}
\description{ There is no \code{step.gam} in package \code{mgcv}. The
\code{mgcv} default for model selection is to use MSE/KL-distance criteria
such as GCV or UBRE/AIC. Since the smoothness estimation part of model
selection is done in this way it is logically most consistent to perform model
selection on the basis of such criteria: i.e. to decide which terms to include
or omit by looking at changes in GCV/UBRE/AIC score. 

To facilitate fully automatic model selection the package includes 2 classes
of smoothers (\code{"cs"} and \code{"ts"}: see \code{\link{s}}) which can be
penalized to zero for sufficiently high smoothing parameter estimates: use of
such smooths provides an effective alternative to step-wise model
selection. The example below shows an example of the application of this
approach, where selection is a fully integrated part of model estimation.

}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}



\examples{
## an example of GCV based model selection as
## an alternative to stepwise selection
library(mgcv)
set.seed(0);n <- 400
dat <- gamSim(1,n=n,scale=2)
dat$x4 <- runif(n, 0, 1)
dat$x5 <- runif(n, 0, 1)
attach(dat)
## Note the increased gamma parameter below to favour
## slightly smoother models...
b<-gam(y~s(x0,bs="ts")+s(x1,bs="ts")+s(x2,bs="ts")+
   s(x3,bs="ts")+s(x4,bs="ts")+s(x5,bs="ts"),gamma=1.4)
summary(b)
plot(b,pages=1)
detach(dat);rm(dat)
}

\keyword{models} \keyword{regression}%-- one or more ..



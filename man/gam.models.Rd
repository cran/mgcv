\name{gam.models}
\alias{gam.models}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Specifying generalized additive models}
\description{ This page is intended to provide some more information on
  how to specify GAMs. Assume that in general we are interested in
  modelling some response variable \eqn{y_i}{y} with an exponential
  family distribution (e.g. Normal, Poisson, binomial, etc.) using predictor variables
  \eqn{x_{1i}, x_{2i}, \ldots}{x1,x2,...}, and let \eqn{\mu_i \equiv
    E(y_i)}{m=E(y)}. A typical GAM might be:
  \deqn{g(\mu_i) = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} +
    f_1(x_{3i})+f_2(x_{4i},x_{4i})}{g(m)=b0 + b1x1 + b2x2 + f1(x3) +
    f2(x4,x5) }
  where \eqn{g} is a smooth monotonic `link' function, and \eqn{f_1}{f1}
  and \eqn{f_2}{f2} are smooth functions. The key idea here is that the
  dependence of the response on the predictors can be represented as a
  parametric sub-model plus the sum of some smooth fucntions of one or
  more of the predictor variables. Thus the model is quite flexible
  relative to strictly parametric linear or generalized linear models,
  but still has much more structure than the completely general model
  that says that the response is just some smooth function of all the
  covariates.

  Note one important point. In order for the model to be identifiable
  the smooth functions have to be constrained to have zero mean (usually
  taken over the set of covariate values). Such constraints are always
  applied by \code{gam}.

  Specification of the distribution and link function is done using the \code{\link{family}}
  argument to \code{\link{gam}} and works in the same way as for
  \code{\link{glm}}. This page hence assumes the default identity link
  normal error family, since the generalizations are easy.

  Starting with the model given above then, the \code{gam} formula would
  be \cr
  \code{y~x1+x2+s(x3)+s(x4,x5)}. \cr
  This would use the default basis for the smooths (a thin plate
  regression spline basis for each), with automatic selection of the
  effective degrees of freedom for both smooths. The dimension of the
  smoothing basis is given a default value as well (the dimension of the
  basis sets an upper limit on the maximum possible degrees of
  freedom for the basis - the limit is typically one less than basis
  dimension). Full details of how to control smooths are given in
  \code{\link{s}} and \code{\link{te}}, and further discussion of basis
  dimension choice can be found in \code{\link{choose.k}}. For the moment suppose that we would like to change
  the basis of the first smooth to a cubic regression spline basis with
  a dimension of 20, while fixing the second term at 25 degrees of
  freedom. The appropriate formula would be:\cr
  \code{y~x1+x2+s(x3,bs="cr",k=20)+s(x4,x5,k=26,fx=T)}.


  Now consider some more unusual models. Consider a model in which
  \eqn{y}{y} is a smooth function of \eqn{x}{x} except at a point
  \eqn{x^*}{x*} where the function jumps discontinuously. This model can
  be written as:
  \deqn{E(y_i)= \beta_0+\beta_1h(x^*,x_i)+f_1(x_i)}{E(y) = b0 + b1
    h(x*,x) + f1(x)}
  where \eqn{h} is a step function jumping from 0 to 1 at
  \eqn{x^*}{x*}. The way to fit this model is to create a variable
  \code{h} which is zero for all \code{x} less than \eqn{x^*}{x*} and
  one otherwise. Then the model formula is:\cr
  \code{y~h+s(x)}. \cr

  Another situation that occurs quite often is the one in which we
  would like to find out if the model:
  \deqn{E(y_i) = f(x_i,z_i)}{E(y)=f(x,z)}
  is really necessary or whether:
  \deqn{E(y_i)=f_1(x_i)+f_2(z_i)}{E(y)=f1(x)+f2(z)}
  wouldn't do just as well. One way to do this is to look at the results
  of fitting:\cr
  \code{y~s(x)+s(z)+s(x,z)}.\cr
  \code{gam} automatically generates side conditions to make this model
  identifiable. You can also estimate `overlapping' models like:\cr
  \code{y~s(x,z)+s(z,v)}.
  
  Sometimes models of the form:
\deqn{E(y_i) = \beta_0+ f(x_i)z_i}{E(y)=b0+f(x)z}
  need to be estimated (where \eqn{f} is a smooth function, as usual.)
  The appropriate formula is:\cr
  \code{y~z+s(x,by=z)}\cr
  - the \code{by} argument ensures that the smooth function gets multiplied by
  covariate \code{z}, but GAM smooths are centred (average value zero),
  so the \code{z+} term is needed as well (\eqn{f} is being
  represented by a constant plus a centred smooth). If we'd wanted:
  \deqn{E(y_i) = f(x_i)z_i}{E(y)=f(x)z}
  then the appropriate formula would be:\cr
   \code{y~z+s(x,by=z)-1}.\cr
  
  The \code{by} mechanism also allows models to be estimated in which
  the form of a smooth depends on the level of a factor, but to do this
  the user must generate the dummy variables for each level of the
  factor. Suppose for example that \code{fac} is a factor with 3 levels
  \code{1}, \code{2}, \code{3},
  and at each level of this factor ther response depends smoothly on a
  variable \code{x} in a manner that is level dependent. Three dummy
  variables \code{fac.1}, \code{fac.2}, \code{fac.3}, can be generated for the factor   
  (e.g. \code{fac.1<-as.numeric(fac==1)}). Then the model formula would
  be:\cr
  \code{y~fac+s(x,by=fac.1)+s(x,by=fac.2)+s(x,by=fac.3)}.\cr

In the above examples the smooths of more than one covariate have all employed 
single penalty thin plate regression splines. These isotropic smooths are not 
alway appropriate: if variables are not naturally `well scaled' relative to each 
other then it is often preferable to use tensor product smooths, with a wiggliness 
penalty for each covariate of the term. See \code{\link{te}} for examples.

}

%- maybe also `usage' for other objects documented here.
\section{WARNING}{
There are no identifiability checks made between the smooth and
parametric parts of a gam formula, and any such lack of identifiability will
cause problems if the underlying fitting routine is \code{\link{mgcv}} (not
the default). 
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
set.seed(10)
n<-400
sig2<-4
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
pi <- asin(1) * 2
f1 <- 2 * sin(pi * x2)
f2 <-  exp(2 * x2) - 3.75887
f3 <-  0.2 * x2^11 * (10 * (1 - x2))^6 + 
       10 * (10 * x2)^3 * (1 - x2)^10 - 1.396
e <- rnorm(n, 0, sqrt(abs(sig2)))
# A continuous `by' variable example.... 
y <- f3*x1 + e
b<-gam(y~x1-1+s(x2,by=x1))
plot(b,pages=1)
summary(b)
# A dummy `by' variable example (with a spurious covariate x0)
fac<-as.factor(c(rep(1,100),rep(2,100),rep(3,200)))
fac.1<-as.numeric(fac==1);fac.2<-as.numeric(fac==2);
fac.3<-as.numeric(fac==3)
y<-f1*fac.1+f2*fac.2+f3*fac.3+ e
b<-gam(y~fac+s(x2,by=fac.1)+s(x2,by=fac.2)+s(x2,by=fac.3)+s(x0))
plot(b,pages=1)
summary(b)
}
\keyword{models} \keyword{regression}%-- one or more ..



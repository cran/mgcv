\name{gam}
\alias{gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Models using penalized regression splines and 
GCV}
\description{ Fits the specified  generalized additive model (GAM) to
  data. \code{gam()} is not a clone of what Splus provides.
Smooth terms are represented using penalized regression splines
with smoothing parameters selected by GCV or by regression splines with
fixed degrees of freedom (mixtures of the two are
permitted). Multi-dimensional smooths are available using penalized thin plate
regression splines, but the user must make sure that covariates are sensibly scaled
relative to each other when using such terms. For a general overview see
Wood (2001).
}
}
\usage{
 gam(formula,family=gaussian(),data=list(),weights=NULL,control=gam.control,scale=0)
}
%- maybe also `usage' for other objects documented here.
\details{ Two alternative bases are available for representing model
  terms. Univariate smooth terms can be represented using conventional
  cubic regression splines - which are very efficient
  computationally - or thin plate regression splines. Multivariate terms
  must be represented using thin plate regression splines. For either
  basis the user specifies the dimension of the basis for each smooth
  term. The dimension of the basis is one more than the maximum degrees of freedom that the 
  term can have, but usually the term will be fitted by penalized
  maximum likelihood estimation and the actual degrees of freedom will be
  chosen by GCV. However, the user can choose to fix the degrees of
  freedom of a term, in which case the actual degrees of freedom will be
  one less than the basis dimension.

  Thin plate regression splines are constructed by starting with the
  basis for a full thin plate spline and then truncating this basis in
  an optimal manner, to obtain a low rank smoother. Details are given in
  Wood (MS submitted). One key advantage of the approach is that it avoids
  the knot placement problems of conventional regression spline
  modelling, but it also has the advantage that smooths of lower rank
  are nested within smooths of higher rank, so that it is legitimate to
  use conventional hypothesis testing methods to compare models based on
  pure regression splines.  
  
  In the case of the cubic regression spline basis, knots  of the spline are placed evenly
  throughout the covariate values to which the term refers:  For
  example, if fitting 101 data with an 11 knot spline of \code{x} then
  there would be a knot at every 10th (ordered)  \code{x} value. The
  parameterization used represents the spline in terms of its
   values at the knots. Connection of these values at neighbouring knots
     by sections of  cubic polynomial constrainted to join at the knots so as to be 
     continuous up to and including second derivative yields a natural cubic 
spline through the values at the knots (given two extra conditions specifying 
     that the second derivative of the curve should be zero at the two end 
     knots). This parameterization gives the parameters a nice interpretability. 

  Given a basis for each smooth term, it easy to obtain a wiggliness
  penalty for each, and to construct a penalized likelihood, which
  balances the fit of the overall model against it's
  complexity. This consists of the log likelihood for the model minus a
  sum of wiggliness penalties (one for each smooth) each multiplied by a
  smoothing parameter. The smoothing parameters control the trade-off
  between fit and smoothness. 
     
  So, the gam fitting problem has become  a penalized glm fitting problem, which can be fitted using a
 slight modification of  \code{glm.fit} : \code{gam.fit}.  The penalized
glm approach also allows smoothing parameters for all smooth terms to
     be  selected  simultaneously by GCV or UBRE. This is achieved as
     part of fitting by calling \code{mgcv}  within \code{gam.fit}.        
     
     Details of the GCV/UBRE minimization method are given in Wood (2000).
 }

\arguments{ 
\item{formula}{ A GAM formula. This is exactly like the formula for a
glm except that smooth terms can be added to the right hand side of the
formula (and a formula of the form \code{y ~ .} is not allowed).
Smooth terms are specified by expressions of the form:
\code{s(var1,var2,...,k=12,fx=FALSE,bs="tp")} where \code{var1},
\code{var2}, etc. are the covariates which the smooth
is a function of and \code{k} is the dimension of the basis used to
represent the smooth term. If \code{k} is not
 specified then \code{k=10*3^(d-1)} is used where \code{d} is the number
 of covariates for this term. \code{fx} is used to indicate whether or
 not this term has a fixed muber of degrees of freedom (\code{fx=FALSE}
 to select d.f. by GCV/UBRE). \code{bs} indicates the basis to use, with
 \code{"cr"} indicating cubic regression spline, and \code{"tp"}
 indicating thin plate regression spline: \code{"cr"} can only be used
 with 1-d smooths.

For backwards compatibility the formula may also include terms like 
\code{s(x,12|f)}, which specifies a regression spline which is not to be penalized
and has 12 knots, or \code{s(x,z,25)} indicating a rank 25 penalized
t.p.r.s. In such cases arguements \code{k}, \code{fx} and \code{bs} are
ignored if supplied and a one dimensional term will always use a cubic
regression spline basis. Note that a term of the form \code{s(x)} will
result in a term with a \code{"tp"} basis.} 
 \item{family}{
This is a family object specifying the distribution and link to use on
fitting etc. See \code{\link{glm}} and \code{\link{family}} for more
details. Where the family is \code{neg.bin} then a negative binomial
family is used based on the implementation in the MASS library.
In this case, if the value of \code{theta} is not given, a version of
\code{glm.nb} :\code{\link{gam.nbut}} is used to estimate theta iteratively, starting from a 
Poisson distribution. This extra layer of iteration  slows  down fitting.
} 
\item{data}{ A data frame containing the model response variable and covariates required by the
formula. If this is missing then the frame from which \code{gam} was called is
searched for the variables specified in the formula.} 
\item{weights}{ 
prior weights on the data.
}
\item{control}{A list as returned by \code{gam.control}, with five user controllable elements:
              \code{maxit} controls maximum iterations in
	      \code{gam.fit}, convergence tolerance in \code{gam.fit} is controlled by \code{epsilon}  
              and the third item is \code{trace}. The smoothing
	      parameter selection method is controlled by two further
	      items: \code{mgcv.tol} controls the convergence tolerance
	      to use in smoothing parameter estimation, while
	      \code{mgcv.half.max} controls the maximum number of step
	      halvings to try in each Newton step if the step fails to
	      reduce the GCV score.}
\item{scale}{ If this is zero then GCV is used for all distributions
except Poisson binomial and negative binomial where UBRE is used with scale parameter
assumed to be 1. If this is greater than 1 it is assumed to be the scale
parameter/variance and UBRE is used. If \code{scale} is negative  GCV 
is always used (for binomial models in particular, it is probably worth 
comparing UBRE and GCV results; for ``over-dispersed Poisson'' GCV is
probably more appropriate than UBRE.)} 
}

\value{ 
 The function returns an object of class \code{"gam"} which has the following elements: 

\item{coefficients}{the coefficients of the fitted model. Parametric
          coefficients are  first, followed  by coefficients for each
          spline term in turn.}

\item{residuals}{the deviance residuals for the fitted model.}

\item{fitted.values}{fitted model predictions of expected value for each
          datum.}
\item{family}{family object specifying distribution and link used.
}
\item{linear.predictor}{fitted model prediction of link function of
expected value for  each datum.}
\item{deviance}{(unpenalized)}

\item{null.deviance}{deviance for single parameter model.} 

\item{df.null}{null degrees of freedom} 

\item{iter}{number of iterations of IRLS taken to get convergence.}

\item{weights}{final weights used in IRLS iteration.}

\item{prior.weights}{prior weights on observations.} 

\item{df.null}{number of data}

\item{y}{response data.}

\item{converged}{indicates whether or not the iterative fitting method converged.} 

\item{sig2}{estimated or supplied variance/scale parameter.}

\item{edf}{estimated degrees of freedom for each smooth.}

\item{boundary}{did parameters end up at boundary of parameter space?} 
\item{sp}{smoothing parameter for each smooth.}

\item{df}{number of knots for each smooth (one more than maximum
          degrees of freedom).}

\item{nsdf}{number of parametric, non-smooth, model terms including the
          intercept.}

\item{Vp}{estimated covariance matrix for parameters.}

\item{covariate.shift}{covariates get shifted so that they are centred
  around zero - this is by how much.}

\item{xp}{knot locations for each cubic regression spline based smooth.
  \code{xp[i,]} are the locations for the ith smooth.}

\item{UZ}{array storing the matrices for transforming from t.p.r.s. basis to
  equivalent t.p.s. basis - see \link{GAMsetup} for details of how the
  matrices are packed in this array.}
\item{Xu}{The set of unique covariate locations used to define t.p.s. from which
  t.p.r.s. basis was derived. Again see \link{GAMsetup} for details of
  the packing algorithm.}
\item{xu.length}{The number of unique covariate combinations in the
  data. }

\item{formula}{the model formula.}

\item{full.formula}{the model formula with each smooth term fully
  expanded and with option arguments given explicitly (i.e. not with reference
  to other variables) - useful for later prediction from the model.}

\item{x}{parametric design matrix columns (including intercept)
          followed by the data that form arguments of the smooths.}
\item{s.type}{type of spline basis used: 0 for cubic, 1 for t.p.r.s.}
\item{p.order}{the order of the penalty used for each term. 0 signals
  auto-selection.}
\item{dim}{number of covariates of which term is a function}
\item{call}{a mode \code{call} object containing the call to \code{gam()} that produced 
            this \code{gam} object (useful for constructing model frames).}

}

\references{

  Hastie and Tibshirani (1990) Generalized Additive Models. Chapman and
  Hall.

  Green and Silverman (1994) Nonparametric Regression and Generalized
  Linear Models. Chapman and Hall.
  
  
Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428

Wood (2001) mgcv:GAMs and Generalized Ridge Regression for R. R News 1(2):20-25
   
   Wood (MS submitted) Thin Plate Regression Splines

  Wahba (1990) Spline Models of Observational Data. SIAM 

  \url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}}

\section{WARNINGS }{The code does not check for rank defficiency of the
model matrix -it will likely just fail instead!

You must have more unique combinations of covariates than the model has total
parameters. (Total parameters is sum of knots plus sum of non-spline terms
less the number of spline terms). 

Automatic smoothing parameter selection is not likely to work well when 
fitting models to very few response data.

Relative scaling of covariates to a multi-dimensional smooth term has an
affect on the results: make sure that relative scalings are
sensible. For example, measuring one spatial co-ordinate in millimetres
and the other in lightyears will usually produce poor results.   

With large datasets (more than a few thousand data) the \code{"tp"}
basis gets very slow to use. In this case use \code{"cr"} for 1-d
smooths. If you need to use multi-dimensional terms with large datasets
and find \code{gam} too slow, please let me know - and I'll up the
priority for fixing this!
 } 

\seealso{  \code{\link{predict.gam}} \code{\link{plot.gam}}}

\examples{
library(mgcv)
set.seed(8)
n<-200
sig2<-4
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
pi <- asin(1) * 2
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f + 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10 - 1.396
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- f + e
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3))
plot(b,pages=1)
# now a GAM with 3df regression spline term & 2 penalized terms
b1<-gam(y~s(x0,k=4,fx=T,bs="tp")+s(x1,k=12)+s(x2,15))
plot(b1,pages=1)
# now fit a 2-d term to x0,x1
b3<-gam(y~s(x0,x1)+s(x2)+s(x3))
par(mfrow=c(2,2))
plot(b3)
par(mfrow=c(1,1))
# now simulate poisson data
g<-exp(f/5)
y<-rpois(rep(1,n),g)
b2<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=poisson)
plot(b2,pages=1)
# negative binomial data
y<-rnbinom(g,size=2,mu=g)
b3<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=neg.bin)
plot(b3,pages=1)
# and a pretty 2-d smoothing example....
test1<-function(x,z,sx=0.3,sz=0.4)  
{ (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
  0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
}
n<-500
old.par<-par(mfrow=c(2,2))
x<-runif(n);z<-runif(n);
xs<-seq(0,1,length=30);zs<-seq(0,1,length=30)
pr<-data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth<-matrix(test1(pr$x,pr$z),30,30)
contour(xs,zs,truth)
y<-test1(x,z)+rnorm(n)*0.1
b4<-gam(y~s(x,z))
fit1<-matrix(predict.gam(b4,pr,se=F),30,30)
contour(xs,zs,fit1)
persp(xs,zs,truth)
persp(b4)
par(old.par)
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
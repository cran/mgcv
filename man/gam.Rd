\name{gam}
\alias{gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Models using penalized regression splines and 
GCV}
\description{ Fits the specified  generalized additive model (GAM) to
  data. \code{gam()} is not a clone of what Splus provides.
Smooth terms are represented using penalized regression splines
with smoothing parameters selected by GCV/UBRE or by regression splines with
fixed degrees of freedom (mixtures of the two are
permitted). Multi-dimensional smooths are available using penalized thin plate
regression splines, but the user must make sure that covariates are sensibly scaled
relative to each other when using such terms. For a general overview see
Wood (2001). For more on specifying models see \code{\link{gam.models}}. For more on model 
selection see \code{\link{gam.selection}}.
}
\usage{

gam(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,
    na.action,control=gam.control(),scale=0,knots=NULL,sp=NULL,
    min.sp=NULL,H=NULL,gamma=1,...)
}
%- maybe also `usage' for other objects documented here.
\details{ 
A generalized additive model (GAM) is a generalized linear model (GLM) in which the linear 
predictor is given by a user specified sum of smooth functions of the covariates plus a 
conventional parametric component of the linear predictor. A simple example is:
\deqn{\log(E(y_i)) = f_1(x_{1i})+f_2(x_{2i})}{log(E(y_i))=f_1(x_1i)+f_2(x_2i)}
where the (idependent) response variables \eqn{y_i \sim {\rm Poi }}{y_i~Poi}, and
\eqn{f_1}{f_1} and \eqn{f_2}{f_2} are smooth functions of covariates \eqn{x_1}{x_1} and 
\eqn{x_2}{x_2}. The log is an example of a link function. 

If absolutely any smooth functions were allowed in model fitting then maximum likelihood 
estimation of such models would invariably result in complex overfitting estimates of 
\eqn{f_1}{f_1}  and \eqn{f_2}{f_2}. For this reason the models are usually fit by 
penalized likelihood 
maximization, in which the model (negative log) likelihood is modified by the addition of 
a penalty for each smooth function, penalizing its `wiggliness'. To control the tradeoff 
between penalizing wiggliness and penalizing badness of fit each penalty is multiplied by 
an associated smoothing parameter: how to estimate these parameters, and 
how to practically represent the smooth functions are the main statistical questions 
introduced by moving from GLMs to GAMs. 

The \code{mgcv} implementation of \code{gam} represents the smooth functions using 
penalized regression splines, and by default uses basis functions for these splines that 
are designed to be optimal, given the number basis functions used. The smooth terms can be 
functions of any number of covariates and the user has some control over how smoothness of 
the functions is measured. 

\code{gam} in \code{mgcv} solves the smoothing parameter estimation problem by using the 
Generalized Cross Validation (GCV) criterion or an Un-Biased Risk Estimator criterion 
(UBRE) which is in practice an approximation to AIC. Smoothing parameters are chosen to 
minimize the GCV or UBRE score for the model, and the main computational challenge solved 
by the \code{mgcv} package is to do this efficiently and reliably. Two alternative 
numerical methods are provided, see \code{\link{mgcv}}, \code{\link{magic}} and 
\code{\link{gam.control}}. 

Broadly \code{gam} works by first constructing basis functions and a quadratic penalty 
coefficient matrix for each smooth term in the model formula, obtaining a model matrix for 
the strictly parametric part of the model formula, and combining these to obtain a 
complete model matrix (/design matrix) and a set of penalty matrices for the smooth terms. 
Some linear identifiability constraints are also obtained at this point. The model is 
fit using \code{\link{gam.fit}}, a modification of \code{\link{glm.fit}}. The GAM 
penalized likelihood maximization problem is solved by penalized Iteratively 
Reweighted  Least Squares (IRLS) (see e.g. Wood 2000). At each iteration a penalized 
weighted least squares problem is solved, and the smoothing parameters of that problem are 
estimated by GCV or UBRE. Eventually both model parameter estimates and smoothing 
parameter estimates converge. 

The fitting approach just described, in which the smoothing parameters are estimated for 
each approximating linear model of the IRLS process was suggested by Chong Gu (see, e.g. 
Gu 2002), and is very computationally efficient. However, because the approach neglects 
the dependence of the iterative weights on the smoothing parameters, it is usually 
possible to find smoothing parameters which actually yield slightly lower GCV/UBRE score 
estimates than those resulting from this `performance iteration'. \code{gam} therefore also 
allows the user to `improve' the smoothing parameter estimates, by using O'Sullivan's 
(1986) suggested method, in which for each trial set of smoothing parameters the IRLS is 
iterated to convergance before the UBRE/GCV score is evaluated. This requires much less 
efficient minimisation of the power iteration based on \code{\link{nlm}}, and is 
therefore quite slow. 





Two alternative bases are available for representing model
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
  Wood (2003). One key advantage of the approach is that it avoids
  the knot placement problems of conventional regression spline
  modelling, but it also has the advantage that smooths of lower rank
  are nested within smooths of higher rank, so that it is legitimate to
  use conventional hypothesis testing methods to compare models based on
  pure regression splines. The t.p.r.s. basis can become expensive to
  calculate for large datasets. In this case the user can supply a reduced 
  set of knots to use in basis construction (see knots, in the argument list).
  
  In the case of the cubic regression spline basis, knots  of the spline are placed evenly
  throughout the covariate values to which the term refers:  For
  example, if fitting 101 data with an 11 knot spline of \code{x} then
  there would be a knot at every 10th (ordered)  \code{x} value. The
  parameterization used represents the spline in terms of its
   values at the knots. The values at neighbouring knots
     are connected by sections of  cubic polynomial constrainted to be 
     continuous up to and including second derivative at the knots. The resulting curve
     is a natural cubic  spline through the values at the knots (given two extra conditions specifying 
     that the second derivative of the curve should be zero at the two end 
     knots). This parameterization gives the parameters a nice interpretability. 

     
     Details of \code{"mgcv"} GCV/UBRE minimization method are given in Wood (2000): 
the basis of the approach 
     is to alternate efficient global optimization with respect to one overall smoothing 
     parameter with Newton updates of a set of relative smoothing parameters for each smooth term. The Newton updates are
backed up by steepest descent, since the GCV/UBRE score functions are not positive definite everywhere.

 }

\arguments{ 
\item{formula}{ A GAM formula (see also \code{\link{gam.models}}). This is exactly like the formula for a
glm except that smooth terms can be added to the right hand side of the
formula (and a formula of the form \code{y ~ .} is not allowed).
Smooth terms are specified by expressions of the form: \cr
\code{s(var1,var2,...,k=12,fx=FALSE,bs="tp",by=a.var)} where \code{var1},
\code{var2}, etc. are the covariates which the smooth
is a function of and \code{k} is the dimension of the basis used to
represent the smooth term. If \code{k} is not
 specified then \code{k=10*3^(d-1)} is used where \code{d} is the number
 of covariates for this term. \code{fx} is used to indicate whether or
 not this term has a fixed muber of degrees of freedom (\code{fx=FALSE}
 to select d.f. by GCV/UBRE). \code{bs} indicates the basis to use, with
 \code{"cr"} indicating cubic regression spline, and \code{"tp"}
 indicating thin plate regression spline: \code{"cr"} can only be used
 with 1-d smooths. \code{by} can be used to specify a variable by which
 the smooth should be multiplied. For example \code{gam(y~z+s(x,by=z))}
   would specify a model \eqn{ E(y) = f(x)z}{E(y)=f(x)z} where
   \eqn{f(\cdot)}{f(.)} is a smooth function (the formula is
   \code{y~x+s(x,by=z)} rather than \code{y~s(x,by=z)} because in GAMs
   the smooths are always set up to sum to zero over the covariate
   values). The \code{by} option is particularly useful for models in
   which different functions of the same variable are required for
   each level of a factor: see \code{\link{s}}.

For backwards compatibility the formula may also include terms like 
\code{s(x,12|f)}, which specifies a regression spline which is not to be penalized
and has 12 knots, or \code{s(x,z,25)} indicating a rank 25 penalized
t.p.r.s. In such cases arguements \code{k}, \code{fx} and \code{bs} are
ignored if supplied and a one dimensional term will always use a cubic
regression spline basis. Note that a term of the form \code{s(x)} will
result in a term with a \code{"tp"} basis.

Formulae can involve nested or "overlapping" terms such as \cr
\code{y~s(x)+s(z)+s(x,z)} or \code{y~s(x,z)+s(z,v)}: see
\code{\link{gam.side.conditions}} for further details and examples. } 
 \item{family}{
This is a family object specifying the distribution and link to use in
fitting etc. See \code{\link{glm}} and \code{\link{family}} for more
details. The negative binomial families provided by the MASS library 
can be used, with or without known \eqn{\theta}{theta} parameter: see
\code{\link{gam.neg.bin}} for details.
} 

\item{data}{ A data frame containing the model response variable and covariates required by the
formula. If this is missing then the frame from which \code{gam} was called is
searched for the variables specified in the formula.} 

\item{weights}{ prior weights on the data.}

\item{subset}{ an optional vector specifying a subset of observations to be
          used in the fitting process.}

\item{na.action}{ a function which indicates what should happen when the data
          contain `NA's.  The default is set by the `na.action' setting
          of `options', and is `na.fail' if that is unset.  The
          ``factory-fresh'' default is `na.omit'.}

\item{control}{A list of fit control parameters returned by 
\code{\link{gam.control}}.}

\item{scale}{ If this is zero then GCV is used for all distributions
except Poisson and binomial where UBRE is used with scale parameter
assumed to be 1. If this is greater than 1 it is assumed to be the scale
parameter/variance and UBRE is used: to use the negative binomial in this case
\eqn{\theta}{theta} must be known. If \code{scale} is negative  GCV 
is always used, which means that the scale parameter will be estimated by GCV and the Pearson 
estimator, or in the case of the negative binomial \eqn{\theta}{theta} will be estimated 
in order to force the GCV/Pearson scale estimate to unity (if this is possible). For binomial models in 
particular, it is probably worth  comparing UBRE and GCV results; for ``over-dispersed Poisson'' GCV is
probably more appropriate than UBRE.} 

\item{knots}{this is an optional list containing user specified knot values to be used for basis construction. 
For the \code{cr} basis the user simply supplies the knots to be used, and there must be the same number as the basis
dimension, \code{k}, for the smooth concerned. For the \code{tp} basis \code{knots} has two uses. Firstly, for large datasets 
the calculation of the \code{tp} basis can be time-consuming. The user can retain most of the advantages of the t.p.r.s. 
approach by supplying  a reduced set of covariate values from which to obtain the basis - 
typically the number of covariate values used will be substantially 
smaller than the number of data, and substantially larger than the basis dimension, \code{k}. The second possibility 
is to avoid the eigen-decomposition used to find the t.p.r.s. basis altogether and simply use 
the basis implied by the chosen knots: this will happen if the number of knots supplied matches the 
basis dimension, \code{k}. For a given basis dimension the second option is 
faster, but gives poorer results (and the user must be quite careful in choosing knot locations). 
Different terms can use different 
numbers of knots, unless they share a covariate.
}

\item{sp}{A vector of smoothing parameters for each term can be provided here.
 Smoothing parameters must 
be supplied in the order that the smooth terms appear in the model 
formula. With fit method \code{"magic"} (see \code{\link{gam.control}} 
and \code{\link{magic}}) then negative elements indicate that the 
parameter should be estimated, and hence a mixture of fixed and estimated 
parameters is possible. With fit method \code{"mgcv"}, if \code{sp} is 
supplied then all its elements must be positive.  }

\item{min.sp}{for fit method \code{"magic"} only, lower bounds can be 
supplied for the smoothing parameters. Note that if this option is used then
the smoothing parameters, \code{sp}, in the returned object will need to be added to
what is supplied here to get the actual smoothing parameters. Lower bounds on the smoothing 
parameters can sometimes help stabilize otherwise divergent P-IRLS iterations.}

\item{H}{a user supplied fixed quadratic penalty on the parameters of the 
GAM can be supplied, with this as its coefficient matrix. A common use of this term is 
to add a ridge penalty to the parameters of the GAM in circumstances in which the model
is close to un-identifiable on the scale of the linear predictor, but perfectly well
defined on the response scale.}

\item{gamma}{It is sometimes useful to inflate the model degrees of 
freedom in the GCV or UBRE score by a constant multiplier. This allows 
such a multiplier to be supplied if fit method is \code{"magic"}.} 
\item{...}{further arguments for 
passing on e.g. to \code{gam.fit}} }

\value{ 
 The function returns an object of class \code{"gam"} which has the following elements: 

\item{boundary}{did parameters end up at boundary of parameter space?} 

\item{by}{a 2-d array of \code{by} variables (i.e. covariates that multiply a
smooth term) \code{by[i,j]} is the jth value for the ith \code{by}
variable. There are only as many rows of this array as there are
\code{by} variables in the model (often 0). The rownames of \code{by}
give the \code{by} variable names.}

\item{by.exists}{an array of logicals: \code{by.exists[i]} is
  \code{TRUE} if the ith smooth has a \code{by} variable associated with
it, \code{FALSE} otherwise.}


\item{coefficients}{the coefficients of the fitted model. Parametric
          coefficients are  first, followed  by coefficients for each
          spline term in turn.}

\item{converged}{indicates whether or not the iterative fitting method converged.} 

\item{covariate.shift}{covariates get shifted so that they are centred
  around zero - this is by how much.}

\item{deviance}{(unpenalized)}

\item{df}{The maximum degrees of freedom for each smooth term plus one. In the \code{"cr"} basis case 
this is the number of knots used. For the \code{"tp"} basis it is the rank of the smooth (one greater than the 
maximum degrees of freedom because of the centering constraint).}

\item{df.null}{null degrees of freedom.} 

\item{dim}{number of covariates of which term is a function.}

\item{edf}{estimated degrees of freedom for each smooth.}

\item{family}{family object specifying distribution and link used.}

\item{fit.method}{ The underlying multiple GCV/UBRE method used: \code{"magic"} 
for the new more stable method, \code{"mgcv"} for the Wood (2000) method.}

\item{fitted.values}{fitted model predictions of expected value for each
          datum.}

\item{formula}{the model formula.}

\item{full.formula}{the model formula with each smooth term fully
  expanded and with option arguments given explicitly (i.e. not with reference
  to other variables) - useful for later prediction from the model.}

\item{gcv.ubre}{The minimized GCV or UBRE score.}

\item{gcv.used}{TRUE if GCV used for smoothing parameter selection, FALSE if UBRE used.}

\item{hat}{array of elements from the leading diagonal of the `hat' (or `influence') matrix. Same length as response data vector.}

\item{iter}{number of iterations of P-IRLS taken to get convergence.}

\item{linear.predictor}{fitted model prediction of link function of
expected value for  each datum.}

\item{mgcv.conv}{ A list of covergence diagnostics relating to smoothing
parameter estimation. Differs for method \code{"magic"} and \code{"mgcv"}. Here is 
the \code{"mgcv"} version:

\item{score}{corresponding to edf, an array of GCV or UBRE scores for the model given the final 
     estimated relative smoothing parameters.}

\item{g}{the gradient of the GCV/UBRE score w.r.t. the relative smoothing parameters at termination.}

\item{h}{the second derivatives corresponding to \code{g} above - i.e. the leading diagonal of the Hessian.}

\item{e}{the eigen-values of the Hessian. All non-negative indicates a positive definite Hessian.}

\item{iter}{the number of iterations taken.}

\item{in.ok}{\code{TRUE} if the second smoothing parameter guess improved the GCV/UBRE score.}

\item{step.fail}{\code{TRUE} if the algorithm terminated by failing to improve the GCV/UBRE score rather than by "converging". 
Not necessarily a problem, but check the above derivative information quite carefully.}

In the case of \code{"magic"} the items are:

\item{full.rank}{The apparent rank of the problem given the model matrix and 
constraints.}

\item{rank}{The numerical rank of the problem.}

\item{fully.converged}{\code{TRUE} is multiple GCV/UBRE converged by meeting 
convergence criteria. \code{FALSE} if method stopped with a steepest descent step 
failure.}

\item{hess.pos.def}{Was the hessian of the gCV/UBRE score positive definite at 
smoothing parameter estimation convergence?}

\item{iter}{How many iterations were required to find the smoothing parameters?}

\item{score.calls}{and how many times did the GCV/UBRE score have to be evaluated?}

\item{rms.grad}{root mean square of the gradient of the GCV/UBRE score at 
convergence.}

} % end of mgcv.conv listing

\item{min.edf}{Minimum possible degrees of freedom for whole model.}

\item{model}{model frame containing all variables needed in original model fit.}

\item{nsdf}{number of parametric, non-smooth, model terms including the
          intercept.}

\item{null.deviance}{deviance for single parameter model.} 

\item{p.order}{the order of the penalty used for each term. 0 signals
  auto-selection.}

\item{prior.weights}{prior weights on observations.} 

\item{residuals}{the deviance residuals for the fitted model.}

\item{sig2}{estimated or supplied variance/scale parameter.}

\item{sp}{smoothing parameter for each smooth.}

\item{s.type}{type of spline basis used: 0 for conventional cubic regression spline, 1 for t.p.r.s.}

\item{UZ}{array storing the matrices for transforming from t.p.r.s. basis to
  equivalent t.p.s. basis - see \link{GAMsetup} for details of how the
  matrices are packed in this array.}

\item{Vp}{estimated covariance matrix for parameter estimators.}

\item{weights}{final weights used in IRLS iteration.}

\item{xp}{knot locations for each cubic regression spline based smooth.
  \code{xp[i,]+covariate.shift[i]} are the locations for the ith smooth.}

\item{Xu}{The set of unique covariate locations used to define t.p.s. from which
  t.p.r.s. basis was derived. Again see \link{GAMsetup} for details of
  the packing algorithm.}

\item{xu.length}{The number of unique covariate combinations in the
  data. }

\item{y}{response data.}

}

\references{

Key References:

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

Background References:

Green and Silverman (1994) Nonparametric Regression and Generalized  Linear Models. Chapman and Hall.
  
Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Gu (2002) Smoothing Spline ANOVA Models, Springer.

Hastie and Tibshirani (1990) Generalized Additive Models. Chapman and Hall.

O'Sullivan, Yandall and Raynor (1986) Automatic smoothing of regression functions in generalized linear models.
J. Am. Statist.Ass. 81:96-103 

Wahba (1990) Spline Models of Observational Data. SIAM 

Wood (2001) mgcv:GAMs and Generalized Ridge Regression for R. R News 1(2):20-25
   
Wood and Augustin (2002) GAMs with integrated model selection using penalized regression splines and applications 
to environmental modelling. Ecological Modelling 157:157-177

\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\section{WARNINGS }{

The \code{"mgcv"} code does not check for rank defficiency of the
model matrix that may result from lack of identifiability between the
parametric and smooth components of the model. 

You must have more unique combinations of covariates than the model has total
parameters. (Total parameters is sum of basis dimensions plus sum of non-spline 
terms less the number of spline terms). 

Automatic smoothing parameter selection is not likely to work well when 
fitting models to very few response data.

Relative scaling of covariates to a multi-dimensional smooth term has an
affect on the results: make sure that relative scalings are
sensible. For example, measuring one spatial co-ordinate in millimetres
and the other in lightyears will usually produce poor results.   

With large datasets (more than a few thousand data) the \code{"tp"}
basis gets very slow to use: use the \code{knots} argument as discussed above and 
shown in the examples. Alternatively, for 1-d smooths  you can use the \code{"cr"} basis.

For data with many  zeroes clustered together in the covariate space it is quite easy to set up 
GAMs which suffer from identifiability problems, particularly when using poisson or binomial
families. The problem is that with e.g. log or logit links, mean value zero corresponds to
an infinite range on the linear predictor scale. Some regularization is possible in such cases: see 
\code{\link{gam.control}} for details.
} 

\seealso{\code{\link{gam.models}}, \code{\link{s}}, \code{\link{predict.gam}},
  \code{\link{plot.gam}}, \code{\link{summary.gam}}, \code{\link{gam.side.conditions}},
\code{\link{gam.selection}},\code{\link{mgcv}}, \code{\link{gam.control}}
\code{\link{gam.check}}, \code{\link{gam.neg.bin}}, \code{\link{magic}},\code{\link{vis.gam}}
 }

\examples{
library(mgcv)
set.seed(0) 
n<-400
sig2<-4
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
pi <- asin(1) * 2
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f+0.2*x2^11*(10*(1-x2))^6+10*(10*x2)^3*(1-x2)^10-1.396
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- f + e
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3))
summary(b)
plot(b,pages=1)
# an extra ridge penalty (useful with convergence problems) ....
bp<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),H=diag(0.5,41)) 
print(b);print(bp);rm(bp)
# set the smoothing parameter for the first term, estimate rest ...
bp<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),sp=c(0.01,-1,-1,-1))
plot(bp,pages=1);rm(bp)
# set lower bounds on smoothing parameters ....
bp<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),min.sp=c(0.001,0.01,0,10)) 
print(b);print(bp);rm(bp)

# now a GAM with 3df regression spline term & 2 penalized terms
b0<-gam(y~s(x0,k=4,fx=TRUE,bs="tp")+s(x1,k=12)+s(x2,15))
plot(b0,pages=1)
# now fit a 2-d term to x0,x1
b1<-gam(y~s(x0,x1)+s(x2)+s(x3))
par(mfrow=c(2,2))
plot(b1)
par(mfrow=c(1,1))
# now simulate poisson data
g<-exp(f/5)
y<-rpois(rep(1,n),g)
b2<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=poisson)
plot(b2,pages=1)
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
fit1<-matrix(predict.gam(b4,pr,se=FALSE),30,30)
contour(xs,zs,fit1)
persp(xs,zs,truth)
vis.gam(b4)
par(old.par)
# very large dataset example using knots
n<-10000
x<-runif(n);z<-runif(n);
y<-test1(x,z)+rnorm(n)
ind<-sample(1:n,1000,replace=FALSE)
b5<-gam(y~s(x,z,k=50),knots=list(x=x[ind],z=z[ind]))
vis.gam(b5)
# and a pure "knot based" spline of the same data
b6<-gam(y~s(x,z,k=100),knots=list(x= rep((1:10-0.5)/10,10),
        z=rep((1:10-0.5)/10,rep(10,10))))
vis.gam(b6,color="heat")
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..








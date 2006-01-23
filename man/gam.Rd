\name{gam}
\alias{gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized additive models with integrated smoothness estimation}
\description{ Fits a generalized additive model (GAM) to
  data. The degree of smoothness of model terms is estimated as part of
  fitting; isotropic or scale invariant smooths of any number of variables are
  available as model terms; confidence/credible intervals are readily
  available for any quantity predicted using a fitted model; \code{gam} is
  extendable: i.e. users can add smooths. 

 Smooth terms are represented using penalized regression splines (or similar smoothers)
with smoothing parameters selected by GCV/UBRE or by regression splines with
fixed degrees of freedom (mixtures of the two are
permitted). Multi-dimensional smooths are available using penalized thin plate
regression splines (isotropic) or tensor product splines (when an isotropic smooth is inappropriate). 
For more on specifying models see \code{\link{gam.models}}. For more on model 
selection see \code{\link{gam.selection}}. For faster fits use the \code{"cr"}
bases for smooth terms, \code{\link{te}} smooths for smooths of several
variables, and performance iteration for smoothing parameter estimation (see
\code{\link{gam.method}}).  For large datasets see warnings.

 \code{gam()} is not a clone of what S-PLUS provides: the major
  differences are (i) that by default estimation of the
  degree of smoothness of model terms is part of model fitting, (ii) a
  Bayesian approach to variance estimation is employed that makes for easier
  confidence interval calculation (with good coverage probabilites) and (iii) the
  facilities for incorporating smooths of more than one variable are
  different: specifically there are no \code{lo} smooths, but instead (a) \code{s}
  terms can have more than one argument, implying an isotropic smooth and (b) \code{te} smooths are
  provided as an effective means for modelling smooth interactions of any
  number of variables via scale invariant tensor product smooths. If you want
  a clone of what S-PLUS provides use \link[gam]{gam} from package \code{gam}.
}
\usage{

gam(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,
    na.action,offset=NULL,control=gam.control(),method=gam.method(),
    scale=0,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,
    fit=TRUE,G=NULL,...)
}
%- maybe also `usage' for other objects documented here.
\details{ 
A generalized additive model (GAM) is a generalized linear model (GLM) in which the linear 
predictor is given by a user specified sum of smooth functions of the covariates plus a 
conventional parametric component of the linear predictor. A simple example is:
\deqn{\log(E(y_i)) = f_1(x_{1i})+f_2(x_{2i})}{log(E(y_i))=f_1(x_1i)+f_2(x_2i)}
where the (independent) response variables \eqn{y_i \sim {\rm Poi }}{y_i~Poi}, and
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
Generalized Cross Validation (GCV) criterion
\deqn{n D / (n - DoF)^2}{n D/(n - DoF)^2}
 or an Un-Biased Risk Estimator (UBRE )criterion
\deqn{D/n + 2 s DoF / n - s }{D/n + 2 s DoF / n -s} 
where \eqn{D}{D} is the deviance, \eqn{n}{n} the number of data, \eqn{s}{s}
the scale parameter and 
\eqn{DoF}{DoF} the effective degrees of freedom of the model. Notice that UBRE is effectively
just AIC rescaled, but is only used when \eqn{s}{s} is known. It is also
possible to replace D by the Pearson statistic (see \code{\link{gam.method}}),
but this can lead to over smoothing. Smoothing parameters are chosen to 
minimize the GCV or UBRE score for the model, and the main computational challenge solved 
by the \code{mgcv} package is to do this efficiently and reliably. Various
alternative numerical methods are provided: see \code{\link{gam.method}}.


Broadly \code{gam} works by first constructing basis functions and one or more quadratic penalty 
coefficient matrices for each smooth term in the model formula, obtaining a model matrix for 
the strictly parametric part of the model formula, and combining these to obtain a 
complete model matrix (/design matrix) and a set of penalty matrices for the smooth terms. 
Some linear identifiability constraints are also obtained at this point. The model is 
fit using \code{\link{gam.fit}}, a modification of \code{\link{glm.fit}}. The GAM 
penalized likelihood maximization problem is solved by penalized Iteratively 
Reweighted  Least Squares (IRLS) (see e.g. Wood 2000). At each iteration a penalized 
weighted least squares problem is solved, and the smoothing parameters of that problem are 
estimated by GCV or UBRE. Eventually both model parameter estimates and smoothing 
parameter estimates converge. Alternatively the P-IRLS scheme is iterated to
convergence for each trial set of smoothing parameters, and GCV or UBRE scores
are only evaluated on convergence - optimization is then `outer' to the P-IRLS
loop: in this case extra derivatives have to be carried along with the P-IRLS
iteration, to facilitate optimization, and \code{\link{gam.fit2}} is used.


Five alternative basis-penalty types  are built in for representing model
smooths, but alternatives can easily be added (see \code{\link{smooth.construct}} which
uses p-splines to illustrate how to add new smooths). 
The built in alternatives for univariate smooths terms are: a conventional penalized
cubic regression spline basis, parameterized in terms of the function values at the knots; 
a cyclic cubic spline with a similar parameterization and thin plate regression splines. 
The cubic spline bases are computationally very efficient, but require `knot' locations to be 
chosen (automatically by default). The thin plate regression splines are optimal low rank 
smooths which do not have knots, but are more computationally costly to set
up. Smooths of several variables can be represented using thin plate
regression splines, or tensor products of any available basis 
including user defined bases (tensor product penalties are obtained automatically form 
the marginal basis penalties). The t.p.r.s. basis is isotropic, so if this is not appropriate tensor 
product terms should be used. Tensor product smooths have one penalty and smoothing parameter per marginal 
basis, which means that the relative scaling of covariates is essentially determined automatically by GCV/UBRE. 
The t.p.r.s. basis and cubic regression spline bases are both available with
either conventional `wiggliness penalties' or penalties augmented with a
shrinkage component: the conventional penalties treat some space of functions
as `completely smooth' and do not penalize such functions at all; the
penalties with extra shrinkage will zero a term altogether for high enough
smoothing parameters: \code{\link{gam.selection}} has an example of the use of
such terms.


For any  basis the user specifies the dimension of the basis for each smooth
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
  set of knots to use in basis construction (see knots, in the argument list), or 
use tensor products of cheaper bases.
  
  In the case of the cubic regression spline basis, knots  of the spline are placed evenly
  throughout the covariate values to which the term refers:  For
  example, if fitting 101 data with an 11 knot spline of \code{x} then
  there would be a knot at every 10th (ordered)  \code{x} value. The
  parameterization used represents the spline in terms of its
   values at the knots. The values at neighbouring knots
     are connected by sections of  cubic polynomial constrained to be 
     continuous up to and including second derivative at the knots. The resulting curve
     is a natural cubic  spline through the values at the knots (given two extra conditions specifying 
     that the second derivative of the curve should be zero at the two end 
     knots). This parameterization gives the parameters a nice interpretability. 

     
     Details of the underlying fitting methods are given in Wood (2000, 2004a). 

 }

\arguments{ 
\item{formula}{ A GAM formula (see also \code{\link{gam.models}}). This is exactly like the formula for a
GLM except that smooth terms can be added to the right hand side of the
formula (and a formula of the form \code{y ~ .} is not allowed).
Smooth terms are specified by expressions of the form: \cr
\code{s(var1,var2,...,k=12,fx=FALSE,bs="tp",by=a.var)} where \code{var1},
\code{var2}, etc. are the covariates which the smooth
is a function of and \code{k} is the dimension of the basis used to
represent the smooth term. If \code{k} is not
 specified then \code{k=10*3^(d-1)} is used where \code{d} is the number
 of covariates for this term. \code{fx} is used to indicate whether or
 not this term has a fixed number of degrees of freedom (\code{fx=FALSE}
 to select d.f. by GCV/UBRE). \code{bs} indicates the basis to use for the smooth: for a full
 list see \code{\link{s}}, but note that the default \code{"tp"}, while
 it possesses nice optimality properties is slow and memory hungry for very
 large datasets (but see examples for how to get around this). \code{by} can be used to specify a variable by which
 the smooth should be multiplied. For example \code{gam(y~z+s(x,by=z))}
   would specify a model \eqn{ E(y) = f(x)z}{E(y)=f(x)z} where
   \eqn{f(\cdot)}{f(.)} is a smooth function (the formula is
   \code{y~x+s(x,by=z)} rather than \code{y~s(x,by=z)} because
   the smooths are always set up to sum to zero over the covariate
   values). The \code{by} option is particularly useful for models in
   which different functions of the same variable are required for
   each level of a factor and for `variable parameter models': see \code{\link{s}}. 

An alternative for specifying smooths of more than one covariate is e.g.: \cr
\code{te(x,z,bs=c("tp","tp"),m=c(2,3),k=c(5,10))} which would specify a tensor product 
smooth of the two covariates \code{x} and \code{z} constructed from marginal t.p.r.s. bases 
of dimension 5 and 10 with marginal penalties of order 2 and 3. Any combination of basis types is 
possible, as is any number of covariates.


Formulae can involve nested or ``overlapping'' terms such as \cr
\code{y~s(x)+s(z)+s(x,z)} or \code{y~s(x,z)+s(z,v)}: see
\code{\link{gam.side}} for further details and examples.} 

\item{family}{
This is a family object specifying the distribution and link to use in
fitting etc. See \code{\link{glm}} and \code{\link{family}} for more
details. The negative binomial families provided by the MASS library 
can be used, with or without known \eqn{\theta}{theta} parameter: see
\code{\link{gam.neg.bin}} for details.
} 

\item{data}{ A data frame containing the model response variable and 
covariates required by the formula. By default the variables are taken 
from \code{environment(formula)}: typically the environment from 
which \code{gam} is called.} 

\item{weights}{ prior weights on the data.}

\item{subset}{ an optional vector specifying a subset of observations to be
          used in the fitting process.}

\item{na.action}{ a function which indicates what should happen when the data
          contain `NA's.  The default is set by the `na.action' setting
          of `options', and is `na.fail' if that is unset.  The
          ``factory-fresh'' default is `na.omit'.}

\item{offset}{Can be used to supply a model offset for use in fitting. Note
that this offset will always be completely ignored when predicting, unlike an offset 
included in \code{formula}: this conforms to the behaviour of
\code{lm} and \code{glm}.}

\item{control}{A list of fit control parameters returned by 
\code{\link{gam.control}}.}
\item{method}{A list controlling the fitting methods used. This can make a big
difference to computational speed, and, in some cases, reliability of
convergence: see \code{\link{gam.method}} for details.}
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
For the \code{cr} and \code{cc} bases the user simply supplies the knots to be used, and there must be the same number as the basis
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
supplied then all its elements must be positive. Note that \code{fx=TRUE} 
in a smooth term over-rides what is supplied here effectively setting the 
smoothing parameter to zero.}

\item{min.sp}{for fit method \code{"magic"} only, lower bounds can be 
supplied for the smoothing parameters. Note that if this option is used then
the smoothing parameters \code{sp}, in the returned object, will need to be added to
what is supplied here to get the actual smoothing parameters. Lower bounds on the smoothing 
parameters can sometimes help stabilize otherwise divergent P-IRLS iterations.}

\item{H}{With fit method \code{"magic"} a user supplied fixed quadratic 
penalty on the parameters of the 
GAM can be supplied, with this as its coefficient matrix. A common use of this term is 
to add a ridge penalty to the parameters of the GAM in circumstances in which the model
is close to un-identifiable on the scale of the linear predictor, but perfectly well
defined on the response scale.}

\item{gamma}{It is sometimes useful to inflate the model degrees of 
freedom in the GCV or UBRE score by a constant multiplier. This allows 
such a multiplier to be supplied if fit method is \code{"magic"}.} 

\item{fit}{If this argument is \code{TRUE} then \code{gam} sets up the model and fits it, but if it is
\code{FALSE} then the model is set up and an object \code{G} is returned which is the output from 
\code{\link{gam.setup}} plus some extra items required to complete the GAM fitting process.}
\item{G}{Usually \code{NULL}, but may contain the object returned by a previous call to \code{gam} with 
\code{fit=FALSE}, in which case all other arguments are ignored except for \code{gamma}, \code{family}, \code{control} and \code{fit}.}
\item{...}{further arguments for 
passing on e.g. to \code{gam.fit}} }

\value{ 
If \code{fit == FALSE} the function returns a list \code{G} of items needed to
fit a GAM, but doesn't actually fit it. 

Otherwise the function returns an object of class \code{"gam"} as described in \code{\link{gamObject}}.
}

\references{

Key References on this implementation:

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

Wood, S.N. (2004a) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass. 99:673-686

Wood, S.N. (2004b) On confidence intervals for GAMs based on penalized
regression splines. Technical Report 04-12 Department of Statistics,
University of Glasgow.

Wood, S.N. (2006) Low rank scale invariant tensor product smooths for
generalized additive mixed models. Biometrics

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. CRC Press.

Key Reference on GAMs and related models:

Hastie (1993) in Chambers and Hastie (1993) Statistical Models in S. Chapman
and Hall.

Hastie and Tibshirani (1990) Generalized Additive Models. Chapman and Hall.

Wahba (1990) Spline Models of Observational Data. SIAM 

Background References:

Green and Silverman (1994) Nonparametric Regression and Generalized  Linear Models. Chapman and Hall.
  
Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Gu (2002) Smoothing Spline ANOVA Models, Springer.

O'Sullivan, Yandall and Raynor (1986) Automatic smoothing of regression
functions in generalized linear models.
J. Am. Statist.Ass. 81:96-103 

Wood (2001) mgcv:GAMs and Generalized Ridge Regression for R. R News 1(2):20-25
   
Wood and Augustin (2002) GAMs with integrated model selection using penalized regression splines and applications 
to environmental modelling. Ecological Modelling 157:157-177

\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}

Front end design inspired by the S function of the same name based on the work
of Hastie and Tibshirani (1990). Underlying methods owe much to the work of
Wahba (e.g. 1990) and Gu (e.g. 2002).

}

\section{WARNINGS }{

If fit method \code{"mgcv"} is selected, the code does not check for rank deficiency of the
model matrix that may result from lack of identifiability between the
parametric and smooth components of the model. 

You must have more unique combinations of covariates than the model has total
parameters. (Total parameters is sum of basis dimensions plus sum of non-spline 
terms less the number of spline terms). 

Automatic smoothing parameter selection is not likely to work well when 
fitting models to very few response data.

With large datasets (more than a few thousand data) the \code{"tp"}
basis gets very slow to use: use the \code{knots} argument as discussed above and 
shown in the examples. Alternatively, for 1-d smooths  you can use the \code{"cr"} basis and 
for multi-dimensional smooths use \code{te} smooths.

For data with many  zeroes clustered together in the covariate space it is quite easy to set up 
GAMs which suffer from identifiability problems, particularly when using Poisson or binomial
families. The problem is that with e.g. log or logit links, mean value zero corresponds to
an infinite range on the linear predictor scale. Some regularization is possible in such cases: see 
\code{\link{gam.control}} for details.
} 

\seealso{\code{\link{gamObject}}, \code{\link{gam.models}}, \code{\link{s}}, \code{\link{predict.gam}},
  \code{\link{plot.gam}}, \code{\link{summary.gam}}, \code{\link{gam.side}},
\code{\link{gam.selection}},\code{\link{mgcv}}, \code{\link{gam.control}}
\code{\link{gam.check}}, \code{\link{gam.neg.bin}}, \code{\link{magic}},\code{\link{vis.gam}}
 }

\examples{
library(mgcv)
set.seed(0) 
n<-400
sig<-2
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f3 <- function(x) 0*x
f <- f0(x0) + f1(x1) + f2(x2)
e <- rnorm(n, 0, sig)
y <- f + e
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3))
summary(b)
plot(b,pages=1,residuals=TRUE)
# same fit in two parts .....
G<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),fit=FALSE)
b<-gam(G=G)
# an extra ridge penalty (useful with convergence problems) ....
bp<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),H=diag(0.5,37)) 
print(b);print(bp);rm(bp)
# set the smoothing parameter for the first term, estimate rest ...
bp<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),sp=c(0.01,-1,-1,-1))
plot(bp,pages=1);rm(bp)
# set lower bounds on smoothing parameters ....
bp<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),min.sp=c(0.001,0.01,0,10)) 
print(b);print(bp);rm(bp)

# now a GAM with 3df regression spline term & 2 penalized terms
b0<-gam(y~s(x0,k=4,fx=TRUE,bs="tp")+s(x1,k=12)+s(x2,k=15))
plot(b0,pages=1)
# now fit a 2-d term to x0,x1
b1<-gam(y~s(x0,x1)+s(x2)+s(x3))
par(mfrow=c(2,2))
plot(b1)
par(mfrow=c(1,1))

# now simulate poisson data
g<-exp(f/4)
y<-rpois(rep(1,n),g)
b2<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=poisson)
plot(b2,pages=1)
# repeat fit using performance iteration
gm <- gam.method(gam="perf.magic")
b3<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=poisson,method=gm)
plot(b3,pages=1)

# a binary example 
g <- (f-5)/3
g <- binomial()$linkinv(g)
y <- rbinom(g,1,g)
lr.fit <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=binomial)
## plot model components with truth overlaid in red
op <- par(mfrow=c(2,2))
for (k in 1:4) {
  plot(lr.fit,residuals=TRUE,select=k)
  xx <- sort(eval(parse(text=paste("x",k-1,sep=""))))
  ff <- eval(parse(text=paste("f",k-1,"(xx)",sep="")))
  lines(xx,(ff-mean(ff))/3,col=2)
}
par(op)
anova(lr.fit)
lr.fit1 <- gam(y~s(x0)+s(x1)+s(x2),family=binomial)
lr.fit2 <- gam(y~s(x1)+s(x2),family=binomial)
AIC(lr.fit,lr.fit1,lr.fit2)

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








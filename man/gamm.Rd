\name{gamm}
\alias{gamm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Mixed Models}
\description{ Fits the specified  generalized additive mixed model (GAMM) to
data, by a call to \code{lme} in the normal errors identity link case, or by 
a call to \code{glmmPQL} from the \code{MASS} library otherwise. 
In the latter case estimates are only approximately MLEs. The routine is typically 
much slower than \code{gam}, and not quite as numerically robust.

Smooths are specified as in a call to \code{\link{gam}} as part of the fixed 
effects model formula, but the wiggly components of the smooth are treated as 
random effects. The random effects structures and correlation structures 
availabale for \code{lme} are used to specify other random effects and 
correlations. 

It is assumed that the random effects and correlation structures are employed 
primarily to model residual correlation in the data and that the prime interest
is in inference about the terms in the fixed effects model formula including 
the smooths. For this reason the routine calculates a posterior covariance 
matrix for the coefficients of all the terms in the fixed effects formula, 
including the smooths.

}
\usage{
gamm(formula,random=NULL,correlation=NULL,family=gaussian(),
data=list(),weights=NULL,subset=NULL,na.action,knots=NULL,
control=lmeControl(niterEM=3),...)
}

\arguments{ 
\item{formula}{ A GAM formula (see also \code{\link{gam.models}}). 
This is exactly like the formula for a
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
 \code{"cr"} indicating cubic regression spline, \code{"cc"} a periodic 
cubic regression spline, and \code{"tp"}
 indicating thin plate regression spline: \code{"cr"} and \code{"cc"} can only be used
 with 1-d smooths. Tensor product smooths are specified using \code{\link{te}} terms.
} 
\item{random}{The (optional) random effects structure as specified in a call to 
\code{\link{lme}}: only the \code{list} form is allowed, to facilitate 
manipulation of the random effects structure within \code{gamm} in order to deal
with smooth terms.}
\item{correlation}{An optional \code{corStruct} object (see \code{\link{corClasses}}) as used to define correlation 
structures in \code{\link{lme}}.}
\item{family}{A \code{family} as used in a call to \code{\link{glm}} or \code{\link{gam}}.
The default \code{gaussian} with identity link causes \code{gamm} to fit by a direct call to
\code{\link{lme}}, otherwise \code{glmmPQL} from the \code{MASS} library is used.}
\item{data}{ A data frame containing the model response variable and 
covariates required by the formula.  By default the variables are taken 
from \code{environment(formula)}, typically the environment from 
which \code{gamm} is called.} 
\item{weights}{ prior weights on the data. Read the documentation for \code{lme} and \code{glmmPQL} very 
carefully before even thinking about using this argument.}
\item{subset}{ an optional vector specifying a subset of observations to be
          used in the fitting process.}
\item{na.action}{ a function which indicates what should happen when the data
          contain `NA's.  The default is set by the `na.action' setting
          of `options', and is `na.fail' if that is unset.  The
          ``factory-fresh'' default is `na.omit'.}



\item{knots}{this is an optional list containing user specified knot values to be used for basis construction. 
For the \code{cr} basis the user simply supplies the knots to be used, and there must be the same number as the basis
dimension, \code{k}, for the smooth concerned. For the \code{tp} basis \code{knots} has two uses. 
Firstly, for large datasets 
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
\item{control}{A list of fit control parameters for \code{\link{lme}} returned by 
\code{\link{lmeControl}}. Note the default setting for the number of EM iterations 
used by \code{lme}: high settings tend to lead to numerical problems because variance components
for smooth terms can legitimately be non-finite.}

\item{...}{further arguments for passing on e.g. to \code{lme}} 
}
%- maybe also `usage' for other objects documented here.
\details{ The Bayesian model of spline smoothing introduced by Wahba (1983) and Silverman (1985) opens 
up the possibility of estimating the degree of smoothness of terms in a generalized additive model
as variances of the wiggly components of the smooth terms treated as random effects. Several authors 
have recognised this (see Wang 1998; Ruppert, Wand and Carroll, 2003) and in the normal errors, identity link case estimation can 
be performed using general linear mixed effects modelling software such as \code{lme}. In the generalized case only 
approximate inference is so far available, for example using the Penalized Quasi-Likelihood approach of Breslow 
and Clayton (1993) as implemented in \code{glmmPQL} by Venables and Ripley (2002). 
One advantage of this approach is that it allows correlated errors to be dealt with via random effects 
or the correlation structures available in the \code{nlme} library. 


Some brief details of how GAMs are represented as mixed models and estimated using \code{lme} or \code{glmmPQL} in \code{gamm} can be found in Wood (manuscript). In addition \code{gamm} obtains a posterior covariance matrix for the parameters of all the fixed effects and the smooth terms. The approach is similar to that described in (Lin & Zhang, 1999) - the covariance matrix of the data (or pseudodata in the generalized case) implied by the weights, correlation and random effects structure is obtained, based on the estimates of the parameters of these terms and this is used to obtain the posterior covariance matrix of the fixed and smooth effects. 

The bases used to represent smooth terms are the same as those used in \code{\link{gam}}.

 }


\value{ Returns a list with two items:
\item{gam}{an object of class \code{gam}, less information relating to GCV/UBRE model selection.}
\item{lme}{the fitted model object returned by \code{lme} or \code{glmmPQL}. Note that the model formulae and grouping 
structures may appear to be rather bizarre, because of the manner in which the GAMM is split up and the calls to 
\code{lme} and \code{glmmPQL} are constructed.}

}

\references{

Breslow, N. E. and Clayton, D. G. (1993) Approximate inference in generalized linear 
mixed models. Journal of the American Statistical Association 88, 9-25.

Lin, X and Zhang, D. (1999) Inference in generalized additive mixed models by using smoothing 
splines. JRSSB. 55(2):381-400

Pinheiro J.C. and Bates, D.M. (2000) Mixed effects Models in S and S-PLUS. Springer

Ruppert, D., Wand, M.P. and Carroll, R.J. (2003) Semiparametric Regression. 
Cambridge

Silverman, B.W. (1985) Some aspects of the spline smoothing approach to nonparametric regression.
JRSSB 47:1-52 

Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics
with S. Fourth edition.  Springer.

Wahba, G. (1983) Bayesian confidence intervals for the cross validated smoothing spline. 
JRSSB 45:133-150

Wood, S.N. (in press) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. Journal of the American Statistical Association.

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

Wood, S.N. (manuscript) Tensor product smooth interactions in Generalized Additive Mixed Models.

Wang, Y. (1998) Mixed effects smoothing spline analysis of variance. J.R. Statist. Soc. B 60, 159-174

\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\section{WARNINGS }{
\code{lme} and \code{glmmPQL} will not deal with offsets, so neither can \code{gamm}.

Models like \code{s(z)+s(x)+s(x,z)} are not currently supported.

\code{gamm} is not as numerically stable as \code{gam}: an \code{lme} call will occasionally fail. Experimenting with 
\code{niterEM} in the \code{control} argument can sometimes help.

\code{gamm} is usually much slower than \code{gam}.
 
} 

\seealso{\code{\link{te}}, \code{\link{s}}, \code{\link{predict.gam}},
  \code{\link{plot.gam}}, \code{\link{summary.gam}}, \code{\link{gam.neg.bin}}, 
\code{\link{vis.gam}},\code{\link{pdTens}},\code{\link{gamm.setup}}
 }

\examples{
library(mgcv)
set.seed(0) 
n <- 400
sig <- 2
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f+0.2*x2^11*(10*(1-x2))^6+10*(10*x2)^3*(1-x2)^10-1.396
e <- rnorm(n, 0, sig)
y <- f + e
b <- gamm(y~s(x0)+s(x1)+s(x2)+s(x3))
plot(b$gam,pages=1)

b <- gamm(y~te(x0,x1)+s(x2)+s(x3)) 
op <- par(mfrow=c(2,2))
plot(b$gam)
par(op)

g<-exp(f/5)
y<-rpois(rep(1,n),g)
b2<-gamm(y~s(x0)+s(x1)+s(x2)+s(x3),family=poisson)
plot(b2$gam,pages=1)

# now an example with autocorrelated errors....
x <- 0:(n-1)/(n-1)
f <- 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10-1.396
e <- rnorm(n,0,sig)
for (i in 2:n) e[i] <- 0.6*e[i-1] + e[i]
y <- f + e
op <- par(mfrow=c(2,2))
b <- gamm(y~s(x,k=20),correlation=corAR1())
plot(b$gam);lines(x,f-mean(f),col=2)
b <- gamm(y~s(x,k=20))
plot(b$gam);lines(x,f-mean(f),col=2)
b <- gam(y~s(x,k=20))
plot(b);lines(x,f-mean(f),col=2)
par(op)

# and a "spatial" example
library(nlme);set.seed(1)
test1<-function(x,z,sx=0.3,sz=0.4)
{ (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
  0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
}
n<-200
old.par<-par(mfrow=c(2,2))
x<-runif(n);z<-runif(n);
xs<-seq(0,1,length=30);zs<-seq(0,1,length=30)
pr<-data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth <- matrix(test1(pr$x,pr$z),30,30)
contour(xs,zs,truth)  # true function
f <- test1(x,z)  # true expectation of response

cstr <- corGaus(.1,form = ~x+z)  
cstr <- Initialize(cstr,data.frame(x=x,z=z))
V <- corMatrix(cstr) # correlation matrix for data
Cv <- chol(V)
e <- t(Cv) \%*\% rnorm(n)*0.05 # correlated errors
y <- f + e 
b<- gamm(y~s(x,z,k=50),correlation=corGaus(.1,form=~x+z))
plot(b$gam) # gamm fit accounting for correlation
# overfits when correlation ignored.....  
b1 <- gamm(y~s(x,z,k=50));plot(b1$gam) 
b2 <- gam(y~s(x,z,k=50));plot(b2)
par(old.par)

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..








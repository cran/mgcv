\name{jagam}
\alias{jagam}
\alias{sim2jam}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Just Another Gibbs Additive Modeller: JAGS support for mgcv.}

\description{Facilities to auto-generate model specification code and associated data to simulate with GAMs in JAGS (or BUGS). This is useful for inference about models with complex random effects structure best coded in JAGS. It is a very innefficient approach to making inferences about standard GAMs. The idea is that \code{jagam} generates template JAGS code, and associated data, for the smooth part of the model. This template is then directly edited to include other stochastic components. After simulation with the resulting model, facilities are provided for plotting and prediction with the model smooth components. 
}
\usage{
jagam(formula,family=gaussian,data=list(),file,weights=NULL,na.action,
offset=NULL,knots=NULL,sp=NULL,drop.unused.levels=TRUE,
control=gam.control(),centred=TRUE,sp.prior = "gamma",diagonalize=FALSE)

sim2jam(sam,pregam,edf.type=2,burnin=0)
}
%- maybe also `usage' for other objects documented here.

\arguments{ 
\item{formula}{ A GAM formula (see \code{\link{formula.gam}} and also \code{\link{gam.models}}). 
This is exactly like the formula for a GLM except that smooth terms,  \code{\link{s}}, \code{\link{te}}, \code{\link{ti}} 
and \code{\link{t2}} can be added to the right hand side to specify that the linear predictor depends on smooth functions of predictors (or linear functionals of these).
} 

\item{family}{
This is a family object specifying the distribution and link function to use. 
See \code{\link{glm}} and \code{\link{family}} for more
details. Currently only gaussian, poisson, binomial and Gamma families are supported, 
but the user can easily modify the assumed distribution in the JAGS code.
} 

\item{data}{ A data frame or list containing the model response variable and 
covariates required by the formula. By default the variables are taken 
from \code{environment(formula)}: typically the environment from 
which \code{jagam} is called.} 

\item{file}{Name of the file to which JAGS model specification code should be written. See \code{\link{setwd}} for setting and querying the current working directory.}

\item{weights}{ prior weights on the data.}

\item{na.action}{ a function which indicates what should happen when the data
          contain `NA's.  The default is set by the `na.action' setting
          of `options', and is `na.fail' if that is unset.  The
          ``factory-fresh'' default is `na.omit'.}

\item{offset}{Can be used to supply a model offset for use in fitting. Note
that this offset will always be completely ignored when predicting, unlike an offset 
included in \code{formula}: this conforms to the behaviour of
\code{lm} and \code{glm}.}


\item{control}{A list of fit control parameters to replace defaults returned by 
\code{\link{gam.control}}. Any control parameters not supplied stay at their default values. 
little effect on \code{jagam}.}


\item{knots}{this is an optional list containing user specified knot values to be used for basis construction. 
For most bases the user simply supplies the knots to be used, which must match up with the \code{k} value
supplied (note that the number of knots is not always just \code{k}). 
See \code{\link{tprs}} for what happens in the \code{"tp"/"ts"} case. 
Different terms can use different numbers of knots, unless they share a covariate.
}

\item{sp}{A vector of smoothing parameters can be provided here.
 Smoothing parameters must be supplied in the order that the smooth terms appear in the model 
formula (without forgetting null space penalties). Negative elements indicate that the parameter should be estimated, and hence a mixture 
of fixed and estimated parameters is possible. If smooths share smoothing parameters then \code{length(sp)} 
must correspond to the number of underlying smoothing parameters.}

\item{drop.unused.levels}{by default unused levels are dropped from factors before fitting. For some smooths 
involving factor variables you might want to turn this off. Only do so if you know what you are doing.}

\item{centred}{Should centring constraints be applied to the smooths, as is usual with GAMS? Only set 
this to \code{FALSE} if you know exactly what you are doing. If \code{FALSE} there is a (usually global) 
intercept for each smooth.}

\item{sp.prior}{\code{"gamma"} or \code{"log.uniform"} prior for the smoothing parameters? Do check that the 
default parameters are appropriate for your model in the JAGS code.}

\item{diagonalize}{Should smooths be re-parameterized to have i.i.d. Gaussian priors (where possible)? For Gaussian data this 
allows efficient conjugate samplers to be used, and it can also work well with GLMs if the JAGS \code{"glm"} module is loaded, but otherwise it is often better to update smoothers blockwise, and not do this.}

\item{sam}{jags sample object, containing at least fields \code{b} (coefficients) and \code{rho} (log 
smoothing parameters). May also contain field \code{mu} containing monitored expected response.}

\item{pregam}{standard \code{mgcv} GAM setup data, as returned in \code{jagam} return list.}

\item{edf.type}{Since EDF is not uniquely defined and may be affected by the stochastic structure added to the JAGS model file, 3 options are offered. See details.}

\item{burnin}{the amount of burn in to discard from the simulation chains. Limited to .9 of the chain length.}

}

\value{ 
For \code{jagam} a three item list containing 
\item{pregam}{standard \code{mgcv} GAM setup data.}
\item{jags.data}{list of arguments to be supplied to JAGS containing information referenced in model specification.}
\item{jags.ini}{initialization data for smooth coefficients and smoothing parameters.}

For \code{sim2jam} an object of class \code{"jam"}: a partial version of an \code{mgcv} \code{\link{gamObject}}, suitable for 
plotting and predicting.
}


\details{Smooths are easily incorportated into JAGS models using multivariate normal priors on the smooth coefficients. The smoothing parameters and smoothing penalty matrices directly specifiy the prior multivariate normal precision matrix. Normally a smoothing penalty does not correspond to a full rank precision matrix, implying an improper prior inappropriate for Gibbs sampling. To rectify this problem the null space penalties suggested in Marra and Wood (2011) are added to the usual penalties. 

In an additive modelling context it is usual to centre the smooths, to avoid the identifiability issues associated with having an intercept for each smooth term (in addition to a global intercept). Under Gibbs sampling with JAGS it is technically possible to omit this centring, since we anyway force propriety on the priors, and this propiety implies formal model identifiability. However, in most situations this formal identifiability is rather artificial and does not imply statistically meaningfull identifiability. Rather it serves only to massively inflate confidence intervals, since the multiple intercept terms are not identifiable from the data, but only from the prior. By default then, \code{jagam} imposes standard GAM identifiability constraints on all smooths. The \code{centred} argument does allow you to turn this off, but it is not recommended. If you do set \code{centred=FALSE} then chain convergence and mixing checks should be particularly stringent. 

The final technical issue for model setup is the setting of initial conditions for the coefficients and smoothing parameters. The approach taken is to take the default initial smoothing parameter values used elsewhere by \code{mgcv}, and to take a single PIRLS fitting step with these smoothing parameters in order to obtain starting values for the smooth coefficients. In the setting of fully conjugate updating the initial values of the coefficients are not critical, and good results are obtained without supplying them. But in the usual setting in which slice sampling is required for at least some of the updates then very poor results can sometimes be obtained without initial values, as the sampler simply fails to find the region of the posterior mode. 

The \code{sim2jam} function takes the partial \code{gam} object (\code{pregam}) from \code{jagam} along with simulation output in standard \code{rjags} form and creates a reduced version of a \code{gam} object, suitable for plotting and prediction of the model's smooth components. \code{sim2gam} computes effective degrees of freedom for each smooth, but it should be noted that there are several possibilites for doing this in the context of a model with a complex random effects structure. The simplest approach (\code{edf.type=0}) is to compute the degrees of freedom that the smooth would have had if it had been part of an unweighted Gaussian additive model. One might choose to use this option if the model has been modified so that the response distribution and/or link are not those that were specified to \code{jagam}. The second option is (\code{edf.type=1}) uses the edf that would have been computed by \code{\link{gam}} had it produced these estimates - in the context in which the JAGS model modifications have all been about modifying the random effects structure, this is equivalent to simply setting all the random effects to zero for the effective degrees of freedom calculation. The default option (\code{edf.type=2}) is to base the EDF on the sample covariance matrix, \code{Vp}, of the model coefficients. If the simulation output (\code{sim}) includes a \code{mu} field, then this will be used to form the weight matrix \code{W} in \code{XWX = t(X)\%*\%W\%*\%X}, where the EDF is computed from \code{rowSums(Vp*XWX)*scale}. If \code{mu} is not supplied then it is estimated from the the model matrix \code{X} and the mean of the simulated coefficients, but the resulting \code{W} may not be strictly comaptible with the \code{Vp} matrix in this case. In the situation in which the fitted model is very different in structure from the regression model of the template produced by \code{jagam} then the default option may make no sense, and indeed it may be best to use option 0.

}


\references{
Wood, S.N. (2016) Just Another Gibbs Additive Modeller: Interfacing JAGS and mgcv. Journal of Statistical Software 75(7):1-15 doi:10.18637/jss.v075.i07)

Marra, G. and S.N. Wood (2011) Practical variable selection for generalized additive models.
Computational Statistics & Data Analysis 55(7): 2372-2387

Here is a key early reference to smoothing using BUGS (although the approach and smooths used are a bit different to jagam)

Crainiceanu, C. M. D Ruppert, & M.P. Wand (2005) Bayesian Analysis for Penalized Spline Regression Using WinBUGS Journal of Statistical Software 14.
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}

\section{WARNINGS }{
Gibb's sampling is a very slow inferential method for standard GAMs. It is only likely to be worthwhile when complex random effects structures are required above what is possible with direct GAMM methods.

Check that the parameters of the priors on the parameters are fit for your purpose.
} 

\seealso{\code{\link{gam}}, \code{\link{gamm}}, \code{\link{bam}}
}

\examples{
## the following illustrates a typical workflow. To run the 
## 'Not run' code you need rjags (and JAGS) to be installed.
require(mgcv)
  
set.seed(2) ## simulate some data... 
n <- 400
dat <- gamSim(1,n=n,dist="normal",scale=2)
## regular gam fit for comparison...
b0 <- gam(y~s(x0)+s(x1) + s(x2)+s(x3),data=dat,method="REML")

## Set directory and file name for file containing jags code.
## In real use you would *never* use tempdir() for this. It is
## only done here to keep CRAN happy, and avoid any chance of
## an accidental overwrite. Instead you would use
## setwd() to set an appropriate working directory in which
## to write the file, and just set the file name to what you
## want to call it (e.g. "test.jags" here). 

jags.file <- paste(tempdir(),"/test.jags",sep="") 

## Set up JAGS code and data. In this one might want to diagonalize
## to use conjugate samplers. Usually call 'setwd' first, to set
## directory in which model file ("test.jags") will be written.

jd <- jagam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat,file=jags.file,
            sp.prior="gamma",diagonalize=TRUE)

## In normal use the model in "test.jags" would now be edited to add 
## the non-standard stochastic elements that require use of JAGS....

\dontrun{
require(rjags)
load.module("glm") ## improved samplers for GLMs often worth loading
jm <-jags.model(jags.file,data=jd$jags.data,inits=jd$jags.ini,n.chains=1)
list.samplers(jm)
sam <- jags.samples(jm,c("b","rho","scale"),n.iter=10000,thin=10)
jam <- sim2jam(sam,jd$pregam)
plot(jam,pages=1)
jam
pd <- data.frame(x0=c(.5,.6),x1=c(.4,.2),x2=c(.8,.4),x3=c(.1,.1))
fv <- predict(jam,newdata=pd)
## and some minimal checking...
require(coda)
effectiveSize(as.mcmc.list(sam$b))
}

## a gamma example...
set.seed(1); n <- 400
dat <- gamSim(1,n=n,dist="normal",scale=2)
scale <- .5; Ey <- exp(dat$f/2)
dat$y <- rgamma(n,shape=1/scale,scale=Ey*scale)
jd <- jagam(y~s(x0)+te(x1,x2)+s(x3),data=dat,family=Gamma(link=log),
            file=jags.file,sp.prior="log.uniform")

## In normal use the model in "test.jags" would now be edited to add 
## the non-standard stochastic elements that require use of JAGS....

\dontrun{
require(rjags)
## following sets random seed, but note that under JAGS 3.4 many
## models are still not fully repeatable (JAGS 4 should fix this)
jd$jags.ini$.RNG.name <- "base::Mersenne-Twister" ## setting RNG
jd$jags.ini$.RNG.seed <- 6 ## how to set RNG seed
jm <-jags.model(jags.file,data=jd$jags.data,inits=jd$jags.ini,n.chains=1)
list.samplers(jm)
sam <- jags.samples(jm,c("b","rho","scale","mu"),n.iter=10000,thin=10)
jam <- sim2jam(sam,jd$pregam)
plot(jam,pages=1)
jam
pd <- data.frame(x0=c(.5,.6),x1=c(.4,.2),x2=c(.8,.4),x3=c(.1,.1))
fv <- predict(jam,newdata=pd)
}

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
\concept{Varying coefficient model}
\concept{Functional linear model}
\concept{Penalized GLM}
\concept{Generalized Additive Model}
\concept{Penalized regression}
\concept{Spline smoothing}
\concept{Penalized regression spline}
\concept{Generalized Cross Validation}
\concept{Smoothing parameter selection}
\concept{tensor product smoothing}
\concept{thin plate spline}
\concept{P-spline}
\concept{Generalized ridge regression}







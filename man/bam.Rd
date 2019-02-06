\name{bam}
\alias{bam}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized additive models for very large datasets}

\description{ Fits a generalized additive model (GAM) to a very large
data set, the term `GAM' being taken to include any quadratically penalized GLM (the extended families
listed in \code{\link{family.mgcv}} can also be used).  
The degree of smoothness of model terms is estimated as part of
fitting. In use the function is much like \code{\link{gam}}, except that the numerical methods
are designed for datasets containing upwards of several tens of thousands of data (see Wood, Goude and Shaw, 2015). The advantage 
of \code{bam} is much lower memory footprint than \code{\link{gam}}, but it can also be much faster, 
for large datasets. \code{bam} can also compute on a cluster set up by the \link[parallel]{parallel} package.

An alternative fitting approach (Wood et al. 2017) is provided by the \code{discrete==TRUE} method. In this case a method based on discretization of covariate values and C code level parallelization (controlled by the \code{nthreads} argument instead of the \code{cluster} argument) is used. This extends both the data set and model size that are practical.

}
\usage{
bam(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,
    na.action=na.omit, offset=NULL,method="fREML",control=list(),
    select=FALSE,scale=0,gamma=1,knots=NULL,sp=NULL,min.sp=NULL,
    paraPen=NULL,chunk.size=10000,rho=0,AR.start=NULL,discrete=FALSE,
    cluster=NULL,nthreads=1,gc.level=1,use.chol=FALSE,samfrac=1,
    coef=NULL,drop.unused.levels=TRUE,G=NULL,fit=TRUE,drop.intercept=NULL,...)
}
%- maybe also `usage' for other objects documented here.

\arguments{ 
\item{formula}{ A GAM formula (see \code{\link{formula.gam}} and also \code{\link{gam.models}}). 
This is exactly like the formula for a GLM except that smooth terms, \code{s} and \code{te} can be added 
to the right hand side to specify that the linear predictor depends on smooth functions of predictors 
(or linear functionals of these).
} 

\item{family}{
This is a family object specifying the distribution and link to use in
fitting etc. See \code{\link{glm}} and \code{\link{family}} for more
details. The extended families listed in \code{\link{family.mgcv}} can also be used.
} 

\item{data}{ A data frame or list containing the model response variable and 
covariates required by the formula. By default the variables are taken 
from \code{environment(formula)}: typically the environment from 
which \code{gam} is called.} 

\item{weights}{  prior weights on the contribution of the data to the log likelihood. Note that a weight of 2, for example, 
                is equivalent to having made exactly the same observation twice. If you want to reweight the contributions 
                of each datum without changing the overall magnitude of the log likelihood, then you should normalize the weights
                (e.g. \code{weights <- weights/mean(weights)}).}

\item{subset}{ an optional vector specifying a subset of observations to be
          used in the fitting process.}

\item{na.action}{ a function which indicates what should happen when the data
          contain `NA's.  The default is set by the `na.action' setting
          of `options', and is `na.fail' if that is unset.  The
          ``factory-fresh'' default is `na.omit'.}

\item{offset}{Can be used to supply a model offset for use in fitting. Note
that this offset will always be completely ignored when predicting, unlike an offset 
included in \code{formula} (this used to conform to the behaviour of
\code{lm} and \code{glm}).}



\item{method}{The smoothing parameter estimation method. \code{"GCV.Cp"} to use GCV for unknown scale parameter and
Mallows' Cp/UBRE/AIC for known scale. \code{"GACV.Cp"} is equivalent, but using GACV in place of GCV. \code{"REML"} 
for REML estimation, including of unknown scale, \code{"P-REML"} for REML estimation, but using a Pearson estimate 
of the scale. \code{"ML"} and \code{"P-ML"} are similar, but using maximum likelihood in place of REML. Default 
\code{"fREML"} uses fast REML computation.}

\item{control}{A list of fit control parameters to replace defaults returned by 
\code{\link{gam.control}}. Any control parameters not supplied stay at their default values.}

\item{select}{Should selection penalties be added to the smooth effects, so that they can in principle be 
penalized out of the model? See \code{gamma} to increase penalization.  Has the side effect that smooths no longer have a fixed effect component (improper prior from a Bayesian perspective) allowing REML comparison of models with the same fixed effect structure. 
}

\item{scale}{ If this is positive then it is taken as the known scale parameter. Negative signals that the 
scale paraemter is unknown. 0 signals that the scale parameter is 1  for Poisson and binomial and unknown otherwise. 
Note that (RE)ML methods can only work with scale parameter 1 for the Poisson and binomial cases.    
} 

\item{gamma}{Increase above 1 to force smoother fits. \code{gamma} is used to multiply the effective degrees of freedom in the GCV/UBRE/AIC score (so \code{log(n)/2} is BIC like). \code{n/gamma} can be viewed as an effective sample size, which allows it to play a similar role for RE/ML smoothing parameter estimation.} 


\item{knots}{this is an optional list containing user specified knot values to be used for basis construction. 
For most bases the user simply supplies the knots to be used, which must match up with the \code{k} value
supplied (note that the number of knots is not always just \code{k}). 
See \code{\link{tprs}} for what happens in the \code{"tp"/"ts"} case. 
Different terms can use different numbers of knots, unless they share a covariate.
}

\item{sp}{A vector of smoothing parameters can be provided here.
 Smoothing parameters must be supplied in the order that the smooth terms appear in the model 
formula. Negative elements indicate that the parameter should be estimated, and hence a mixture 
of fixed and estimated parameters is possible. If smooths share smoothing parameters then \code{length(sp)} 
must correspond to the number of underlying smoothing parameters.}

\item{min.sp}{Lower bounds can be supplied for the smoothing parameters. Note
that if this option is used then the smoothing parameters \code{full.sp}, in the 
returned object, will need to be added to what is supplied here to get the 
 smoothing parameters actually multiplying the penalties. \code{length(min.sp)} should 
always be the same as the total number of penalties (so it may be longer than \code{sp},
if smooths share smoothing parameters).}


\item{paraPen}{optional list specifying any penalties to be applied to parametric model terms. 
\code{\link{gam.models}} explains more.}


\item{chunk.size}{The model matrix is created in chunks of this size, rather than ever being formed whole. 
Reset to \code{4*p} if \code{chunk.size < 4*p} where \code{p} is the number of coefficients.}

\item{rho}{An AR1 error model can be used for the residuals (based on dataframe order), of Gaussian-identity 
           link models. This is the AR1 correlation parameter. Standardized residuals (approximately 
           uncorrelated under correct model) returned in 
           \code{std.rsd} if non zero. Also usable with other models when \code{discrete=TRUE}, in which case the AR model
           is applied to the working residuals and corresponds to a GEE approximation.}

\item{AR.start}{logical variable of same length as data, \code{TRUE} at first observation of an independent
section of AR1 correlation. Very first observation in data frame does not need this. If \code{NULL} then 
there are no breaks in AR1 correlaion.}

\item{discrete}{with \code{method="fREML"} it is possible to discretize covariates for storage and efficiency reasons.
If \code{discrete} is \code{TRUE}, a number or a vector of numbers for each smoother term, then discretization happens. If numbers are supplied they give the number of discretization bins.}

\item{cluster}{\code{bam} can compute the computationally dominant QR decomposition in parallel using \link[parallel]{parLapply}
from the \code{parallel} package, if it is supplied with a cluster on which to do this (a cluster here can be some cores of a 
single machine). See details and example code. 
}

\item{nthreads}{Number of threads to use for non-cluster computation (e.g. combining results from cluster nodes).
If \code{NA} set to \code{max(1,length(cluster))}. See details.}

\item{gc.level}{to keep the memory footprint down, it helps to call the garbage collector often, but this takes 
a substatial amount of time. Setting this to zero means that garbage collection only happens when R decides it should. Setting to 2 gives frequent garbage collection. 1 is in between.}

\item{use.chol}{By default \code{bam} uses a very stable QR update approach to obtaining the QR decomposition
of the model matrix. For well conditioned models an alternative accumulates the crossproduct of the model matrix
and then finds its Choleski decomposition, at the end. This is somewhat more efficient, computationally.}


\item{samfrac}{For very large sample size Generalized additive models the number of iterations needed for the model fit can 
be reduced by first fitting a model to a random sample of the data, and using the results to supply starting values. This initial fit is run with sloppy convergence tolerances, so is typically very low cost. \code{samfrac} is the sampling fraction to use. 0.1 is often reasonable. }

\item{coef}{initial values for model coefficients}

\item{drop.unused.levels}{by default unused levels are dropped from factors before fitting. For some smooths 
involving factor variables you might want to turn this off. Only do so if you know what you are doing.}

\item{G}{if not \code{NULL} then this should be the object returned by a previous call to \code{bam} with 
\code{fit=FALSE}. Causes all other arguments to be ignored except \code{sp}, \code{chunk.size}, \code{gamma},\code{nthreads}, \code{cluster}, \code{rho}, \code{gc.level}, \code{samfrac}, \code{use.chol}, \code{method} and \code{scale} (if >0).}

\item{fit}{if \code{FALSE} then the model is set up for fitting but not estimated, and an object is returned, suitable for passing as the \code{G} argument to \code{bam}.}

\item{drop.intercept}{Set to \code{TRUE} to force the model to really not have the a constant in the parametric model part,
even with factor variables present.}

\item{...}{further arguments for 
passing on e.g. to \code{gam.fit} (such as \code{mustart}). }

}

\value{ 
An object of class \code{"gam"} as described in \code{\link{gamObject}}.
}


\details{ When \code{discrete=FALSE}, \code{bam} operates by first setting up the basis characteristics for the smooths, using a representative subsample of the data. Then the model matrix is constructed in blocks using \code{\link{predict.gam}}. For each block the factor R,
from the QR decomposition of the whole model matrix is updated, along with Q'y. and the sum of squares of y. At the end of 
block processing, fitting takes place, without the need to ever form the whole model matrix. 

In the generalized case, the same trick is used with the weighted model matrix and weighted pseudodata, at each step of the PIRLS.
Smoothness selection is performed on the working model at each stage (performance oriented iteration), to maintain the 
small memory footprint. This is trivial to justify in the case of GCV or Cp/UBRE/AIC based model selection, and 
for REML/ML is justified via the asymptotic multivariate normality of Q'z where z is the IRLS pseudodata. 

For full method details see Wood, Goude and Shaw (2015).

Note that POI is not as stable as the default nested iteration used with \code{\link{gam}}, but that for very large, information rich,
datasets, this is unlikely to matter much. 

Note also that it is possible to spend most of the computational time on basis evaluation, if an expensive basis is used. In practice this means that the default \code{"tp"} basis should be avoided: almost any other basis (e.g. \code{"cr"} or \code{"ps"}) 
can be used in the 1D case, and tensor product smooths (\code{te}) are typically much less costly in the multi-dimensional case. 

If \code{cluster} is provided as a cluster set up using \code{\link[parallel]{makeCluster}} (or \code{\link[parallel]{makeForkCluster}}) from the \code{parallel} package, then the rate limiting QR decomposition of the model matrix is performed in parallel using this cluster. Note that the speed ups are often not that great. On a multi-core machine it is usually best to set the cluster size to the number of physical cores, which is often less than what is reported by \code{\link[parallel]{detectCores}}. Using more than the number of physical cores can result in no speed up at all (or even a slow down). Note that a highly parallel BLAS may negate all advantage from using a cluster of cores. Computing in parallel of course requires more memory than computing in series. See examples.

When \code{discrete=TRUE} the covariate data are first discretized. Discretization takes place on a smooth by smooth basis, or in the case of tensor product smooths (or any smooth that can be represented as such, such as random effects), separately for each marginal smooth. The required spline bases are then evaluated at the discrete values, and stored, along with index vectors indicating which original observation they relate to. Fitting is by a version of performance oriented iteration/PQL using REML smoothing parameter selection on each iterative working model (as for the default method). The iteration is based on the derivatives of the REML score, without computing the score itself, allowing the expensive computations to be reduced to one parallel block Cholesky decomposition per iteration (plus two basic operations of equal cost, but easily parallelized). Unlike standard POI/PQL, only one step of the smoothing parameter update for the working model is taken at each step (rather than iterating to the optimal set of smoothing parameters for each working model). At each step a weighted model matrix crossproduct of the model matrix is required - this is efficiently computed from the pre-computed basis functions evaluated at the discretized covariate values. Efficient computation with tensor product terms means that some terms within a tensor product may be re-ordered for maximum efficiency.

When \code{discrete=TRUE} parallel computation is controlled using the \code{nthreads} argument. For this method no cluster computation is used, and the \code{parallel} package is not required. Note that actual speed up from parallelization depends on the BLAS installed and your hardware. With the (R default) reference BLAS using several threads can make a substantial difference, but with a single threaded tuned BLAS, such as openblas, the effect is less marked (since cache use is typically optimized for one thread, and is then sub optimal for several). However the tuned BLAS is usually much faster than using the reference BLAS, however many threads you use. If you have a multi-threaded BLAS installed then you should leave \code{nthreads} at 1, since calling a multi-threaded BLAS from multiple threads usually slows things down: the only exception to this is that you might choose to form discrete matrix cross products (the main cost in the fitting routine) in a multi-threaded way, but use single threaded code for other computations: this can be achieved by e.g. \code{nthreads=c(2,1)}, which would use 2 threads for discrete inner products, and 1 for most code calling BLAS. Not that the basic reason that multi-threaded performance is often disappointing is that most computers are heavily memory bandwidth limited, not flop rate limited. It is hard to get data to one core fast enough, let alone trying to get data simultaneously to several cores.   

\code{discrete=TRUE} will often produce identical results to the methods without discretization, since covariates often only take a modest number of discrete values anyway, so no approximation at all is involved in the discretization process. Even when some approximation is involved, the differences are often very small as the algorithms discretize marginally whenever possible. For example each margin of a tensor product smooth is discretized separately, rather than discretizing onto a grid of covariate values (for an equivalent isotropic smooth we would have to discretize onto a grid). The marginal approach allows quite fine scale discretization and hence very low approximation error. Note that when using the smooth \code{id} mechanism to link smoothing parameters, the discrete method cannot force the linked bases to be identical, so some differences to the none discrete methods will be noticable.   

The extended families given in \code{\link{family.mgcv}} can also be used. The extra parameters of these are estimated by maximizing the penalized likelihood, rather than the restricted marginal likelihood as in \code{\link{gam}}. So estimates may differ slightly from those returned by \code{\link{gam}}. Estimation is accomplished by a Newton iteration to find the extra parameters (e.g. the theta parameter of the negative binomial or the degrees of freedom and scale of the scaled t) maximizing the log likelihood given the model coefficients at each iteration of the fitting procedure. 

}


\references{
Wood, S.N., Goude, Y. & Shaw S. (2015) Generalized additive models for large datasets. Journal of the Royal Statistical Society, Series C 64(1): 139-155.
\url{http://dx.doi.org/10.1111/rssc.12068}

Wood, S.N., Li, Z., Shaddick, G. & Augustin N.H. (2017) Generalized additive models for gigadata: modelling the UK black smoke network daily data. Journal of the American Statistical Association. 112(519):1199-1210
\url{http://dx.doi.org/10.1080/01621459.2016.1195744}
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}


}

\section{WARNINGS }{

The routine will be slow if the default \code{"tp"} basis is used. 

Unless discrete=TRUE, you must have more unique combinations of covariates than the model has total
parameters. (Total parameters is sum of basis dimensions plus sum of non-spline 
terms less the number of spline terms). 

This routine is less stable than `gam' for the same dataset.

With \code{discrete=TRUE}, \code{te} terms are efficiently computed, but \code{t2} are not.

The negbin family is only supported for the *known theta* case.
} 

\seealso{\code{\link{mgcv.parallel}}, 
\code{\link{mgcv-package}}, \code{\link{gamObject}}, \code{\link{gam.models}}, \code{\link{smooth.terms}},
\code{\link{linear.functional.terms}}, \code{\link{s}},
\code{\link{te}} \code{\link{predict.gam}},
\code{\link{plot.gam}}, \code{\link{summary.gam}}, \code{\link{gam.side}},
\code{\link{gam.selection}}, \code{\link{gam.control}}
\code{\link{gam.check}}, \code{\link{linear.functional.terms}} \code{\link{negbin}}, \code{\link{magic}},\code{\link{vis.gam}}
 }

\examples{
library(mgcv)
## See help("mgcv-parallel") for using bam in parallel

## Sample sizes are small for fast run times.

set.seed(3)
dat <- gamSim(1,n=25000,dist="normal",scale=20)
bs <- "cr";k <- 12
b <- bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+
           s(x3,bs=bs),data=dat)
summary(b)
plot(b,pages=1,rug=FALSE)  ## plot smooths, but not rug
plot(b,pages=1,rug=FALSE,seWithMean=TRUE) ## `with intercept' CIs

\donttest{ 
ba <- bam(y ~ s(x0,bs=bs,k=k)+s(x1,bs=bs,k=k)+s(x2,bs=bs,k=k)+
            s(x3,bs=bs,k=k),data=dat,method="GCV.Cp") ## use GCV
summary(ba)}

## A Poisson example...

k <- 15
dat <- gamSim(1,n=21000,dist="poisson",scale=.1)

system.time(b1 <- bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
            data=dat,family=poisson()))
b1

## Similar using faster discrete method...

\donttest{ 
system.time(b2 <- bam(y ~ s(x0,bs=bs,k=k)+s(x1,bs=bs,k=k)+s(x2,bs=bs,k=k)+
            s(x3,bs=bs,k=k),data=dat,family=poisson(),discrete=TRUE))
b2
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







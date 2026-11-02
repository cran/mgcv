\name{scasm}
\alias{scasm}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Shape constrained additive smooth models}

\description{ Fits a generalized additive model (GAM, GAMM, GAMLSS) subject to shape constraints on the components, using quadratic and linear programming and the extended Fellner-Schall method for smoothness estimation. Usable with all the families listed in \code{\link{family.mgcv}}.
}
\usage{

scasm(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,
      na.action, offset=NULL,control=list(),scale=0,select=FALSE,
      knots=NULL,sp=NULL,gamma=1,fit=TRUE,G=NULL,drop.unused.levels=TRUE,
      drop.intercept=NULL,bs=0,...)

}
%- maybe also `usage' for other objects documented here.

\arguments{
\item{formula}{ A GAM formula, or a list of formulae (see \code{\link{formula.gam}} and also \code{\link{gam.models}}). 
These are exactly like the formula for a GLM except that smooth terms, \code{\link{s}}, \code{\link{te}}, \code{\link{ti}} 
and \code{\link{t2}}, can be added to the right hand side to specify that the linear predictor depends on smooth functions of predictors (or linear functionals of these).
} 

\item{family}{
This is a family object specifying the distribution and link to use in
fitting etc (see \code{\link{glm}} and \code{\link{family}}). See 
\code{\link{family.mgcv}} for a full list of what is available, which goes well beyond exponential family.
Note that \code{quasi} families actually result in the use of extended quasi-likelihood (McCullagh and Nelder, 1989, 9.6).
} 


\item{data}{ A data frame or list containing the model response variable and 
covariates required by the formula. By default the variables are taken 
from \code{environment(formula)}: typically the environment from 
which \code{scasm} is called.} 

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


\item{control}{A list of fit control parameters to replace defaults returned by 
\code{\link{gam.control}}. Any control parameters not supplied stay at their default values.}

\item{select}{Should selection penalties be added to the smooth effects, so that they can in principle be 
penalized out of the model? See \code{gamma} to increase penalization.  Has the side effect that smooths no longer have a fixed effect component (improper prior from a Bayesian perspective) allowing REML comparison of models with the same fixed effect structure. 
}

\item{scale}{ If this is positive then it is taken as the known scale parameter. Negative signals that the 
scale parameter is unknown. 0 signals that the scale parameter is 1  for Poisson and binomial and unknown otherwise. 
Note that (RE)ML methods can only work with scale parameter 1 for the Poisson and binomial cases.    
} 

\item{gamma}{Increase above 1 to force smoother fits. Effectively optimizes estimate of what REML would have been for sample size \code{n/gamma}. This is equivalent to simply multiplying the EFS ratio by \code{gamma}.} 


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
must correspond to the number of underlying smoothing parameters. Note that \code{discrete=TRUE}may result in
re-ordering of variables in tensor product smooths for improved efficiency, and \code{sp} must be supplied in re-ordered order.}

\item{drop.unused.levels}{by default unused levels are dropped from factors before fitting. For some smooths 
involving factor variables you might want to turn this off. Only do so if you know what you are doing.}

\item{G}{if not \code{NULL} then this should be the object returned by a previous call to \code{bam} with 
\code{fit=FALSE}. Causes all other arguments to be ignored except \code{sp}, \code{chunk.size}, \code{gamma},\code{nthreads}, \code{cluster}, \code{rho}, \code{gc.level}, \code{samfrac}, \code{use.chol}, \code{method} and \code{scale} (if >0).}

\item{fit}{if \code{FALSE} then the model is set up for fitting but not estimated, and an object is returned, suitable for passing as the \code{G} argument to \code{bam}.}

\item{drop.intercept}{Set to \code{TRUE} to force the model to really not have the a constant in the parametric model part,
even with factor variables present.}

\item{bs}{number of bootstrap replicates to use to approximate sampling distribution of model coefficients. If not zero or \code{FALSE} then results in a matrix \code{bs} being returned in the fit object, whose columns are the bootstrap replicate coefficients.}

\item{...}{further arguments for 
passing on e.g. to \code{gam.fit} (such as \code{mustart}). }



}

\value{ 
An object of class \code{"gam"} as described in \code{\link{gamObject}}.
}


\details{Operates much as \code{\link{gam}}, but allows linear constraints on smooths to facilitate shape constraints. Smoothness selection uses the extended Fellner-Schall method of Wood and Fasiolo (2017). Model coefficient estimation uses linear programming to find initial feasible coefficients, and then sequential quadratic programming to maximize the penalized likelihood given smoothing parameters.

For simple shape constraints on univariate smooths (and hence tensor product smooths) see \code{\link{smooth.construct.sc.smooth.spec}}. General linear constraints can be imposed on any smooth term via the \code{pc} argument to \code{\link{s}} etc. To do this \code{pc} is provided as a list of two lists. The first list component specifies the linear equality constraints, and the second, if present specifies linear inequality constraints.

The first list has a matrix element for each covariate of the smooth, named as the covariates. Two further elements are un-named: a matrix of weights and a vector specifying the right hand side of the condition to be met. For example, suppose that a constraint is to be applied to a smooth term term \eqn{f(X,Z)}{f(X,Z)}. The list would contain named matrices \eqn{\bf X}{X} and \eqn{\bf Z}{Z} and a third unnamed matrix, \eqn{\bf W}{W}, say. All are of the same dimension. The list also contains an unnamed vector with the same number of rows as the matrices, \eqn{\bf b}{b}, say. The constraints encoded are
\deqn{\sum_j f(X_{ij},Z_{ij}) W_{ij} > b_i}{sum_j f(X_ij,Z_ij) W_ij = b_i}
where the summation is over all columns of the matrices. There is one constraint for each row of \eqn{\bf b}{b}. The second list, if present, has the same structure, but encodes inequality constraints (\eqn{>b_i}{>b_i}). If the first list is empty then there are no equality constraints. 

This structure is quite general - derivative and integral constraints are easily implemented numerically. Also note that the simple constraints of \code{\link{smooth.construct.sc.smooth.spec}} can be combined with general \code{pc} constraints.

Note that \code{pc} constraints replace the usual sum to zero identifiability constraints, so if these are required, they need to be explicitly included as equality constraints.

}


\references{
Wood, S.N. (2026) Shape constrained additive smooth models (with discussion). Statistical Modelling. 26(4):382-409
\doi{10.1177/1471082X261436505}

Wood, S.N. and M. Fasiolo (2017) A generalized Fellner-Schall method for smoothing
parameter optimization with application to Tweedie location, scale and shape models.
Biometrics 73 (4), 1071-1081 \doi{10.1111/biom.12666}
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}


}

\section{WARNINGS }{
This is generally slower than \code{gam}, and can only use EFS smoothness estimation
} 

\seealso{
\code{\link{smooth.construct.sc.smooth.spec}},
\code{\link{mgcv-package}}, \code{\link{gamObject}}, \code{\link{gam.models}}, \code{\link{smooth.terms}},
 \code{\link{predict.gam}},
\code{\link{plot.gam}}, \code{\link{summary.gam}},
\code{\link{gam.check}}
 }

\examples{
library(mgcv)
## a density estimation example, illustrating general constraint use... 
n <- 400; set.seed(4)
x <- c(rnorm(n*.8,.4,.13),runif(n*.2))
y <- rep(1e-10,n)
xp=(1:n-.5)/n
options(warn=-1)
## Use an exponential distribution, inverse link and zero pseudoresponse
## data to get correct likelihood. Use a positive smooth and impose
## integration to 1 via a general linear 'pc' constraint...
b <- scasm(y ~ s(x,bs="sc",xt="+",k=15,pc=list(list(x=matrix(xp,1,n),
      matrix(1/n,1,n),1)))-1, family=Gamma,scale=1,knots=list(x=c(0,1)))
options(warn=0)
plot(b,scheme=2)
lines(xp,.2+dnorm(xp,.4,.13)*.8,col=2)

## some more standard examples 

fs <- function(x,k) { ## some simulation functions
  if (k==2) 0.2*x^11*(10*(1-x))^6 + 10*(10*x)^3* (1-x)^10 else
  if (k==0) 2 * sin(pi * x) else
  if (k==1) exp(2 * x) else
  if (k==3) exp((x-.4)*20)/(1+exp((x-.4)*20)) else
  if (k==4) (x-.55)^4 else x
} ## fs

## Simple Poisson example... 
n <- 400; set.seed(4)
x0 <- runif(n);x1 <- runif(n);x2 <- runif(n);x3 <- runif(n);
fac <- factor(sample(1:40,n,replace=TRUE))
f <- fs(x0,k=0) + 2*fs(x1,k=3) + fs(x2,k=2) + 10*fs(x3,k=4)
mu <- exp((f-1)/4)
y <- rpois(n,mu)
k <- 10
b0 <- gam(y~s(x0,bs="bs",k=k)+s(x1,bs="bs",k=k)+s(x2,bs="bs",k=k)+
          s(x3,bs="bs",k=k),method="REML",family=poisson)
	  
b1 <- scasm(y~s(x0,bs="bs",k=k)+s(x1,bs="bs",k=k)+s(x2,bs="bs",k=k)+
            s(x3,bs="bs",k=k),family=poisson)
b0;b1 ## note similarity when no constraints	    

plot(b1,pages=1,scheme=2) ## second term not monotonic, so impose this... 

b2 <- scasm(y~s(x0,bs="bs",k=k)+s(x1,bs="sc",xt="m+",k=k)+
            s(x2,bs="bs",k=k)+s(x3,bs="bs",k=k),family=poisson)
plot(b2,pages=1,scheme=2)

## constraint respecting bootstrap CIs...

b3 <- scasm(y~s(x0,bs="bs",k=k)+s(x1,bs="sc",xt="m+",k=k)+
            s(x2,bs="bs",k=k)+s(x3,bs="bs",k=k),family=poisson,bs=200)
plot(b3,pages=1,scheme=2)
fv <- predict(b3,se=TRUE)
fv1 <- predict(b3,se=.95)
plot(fv$se*3.92,fv1$ul-fv1$ll) ## compare CI widths
abline(0,1,col=2)

## gamma distribution location-scale example

## simulate data...
x0 <- runif(n);x1 <- runif(n);x2 <- runif(n);x3 <- runif(n);
fac <- factor(sample(1:40,n,replace=TRUE))
f <- fs(x0,k=0) + 2*fs(x1,k=3) + fs(x2,k=2) #+ 10*fs(x3,k=4)
mu <- exp((fs(x0,k=0)+ fs(x2,k=2))/5)
th <- exp((fs(x1,k=1)+ 10*fs(x3,k=4))/2-2)
y <- rgamma(n,shape=1/th,scale=mu*th)

## fit model
b <- scasm(list(y~s(x0,bs="sc",xt="c-")+s(x2),~s(x1,bs="sc",xt="m+")
       +s(x3)),family=gammals,bs=200)
plot(b,scheme=2,pages=1,scale=0) ## bootstrap CIs
b$bs <- NULL; plot(b,scheme=2,pages=1,scale=0) ## constraint free CI

## A survival example...

library(survival) ## for data
col1 <- colon[colon$etype==1,] ## focus on single event
col1$differ <- as.factor(col1$differ)
col1$sex <- as.factor(col1$sex)

## set up the AFT response... 
logt <- cbind(log(col1$time),log(col1$time))
logt[col1$status==0,2] <- Inf ## right censoring
col1$logt <- -logt ## -ve conventional for AFT versus Cox PH comparison

## fit the model...
b0 <- gam(logt~s(age,by=sex)+sex+s(nodes)+perfor+rx+obstruct+adhere,
         family=cnorm(),data=col1)
plot(b0,pages=1)

## nodes effect should be increasing...	  
b <- scasm(logt~s(age,by=sex)+sex+s(nodes,bs="sc",xt="m+")+perfor+rx+
           obstruct+adhere,family=cnorm(),data=col1)
plot(b,pages=1)
b 

## same with logistic distribution...
b <- scasm(logt~s(age,by=sex)+sex+s(nodes,bs="sc",xt="m+")+perfor+rx+
           obstruct+adhere,family=clog(),data=col1)
plot(b,pages=1)

## and with Cox PH model... 
b <- scasm(time~s(age,by=sex)+sex+s(nodes,bs="sc",xt="m+")+perfor+rx+
           obstruct+adhere,family=cox.ph(),data=col1,weights=status)
summary(b) 
plot(b,pages=1) ## plot effects

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
\concept{Varying coefficient model}
\concept{Functional linear model}
\concept{Penalized GLM}
\concept{Generalized Additive Model}
\concept{Penalized regression}
\concept{Spline smoothing}
\concept{Penalized regression spline}
\concept{Shape constraint}
\concept{Smoothing parameter selection}








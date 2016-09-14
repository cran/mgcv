\name{cox.ph}
\alias{cox.ph}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Additive Cox Proportional Hazard Model}
\description{The \code{cox.ph} family implements the Cox Proportional Hazards model with Peto's 
correction for ties, and estimation by penalized partial likelihood maximization, for use with 
\code{\link{gam}}. In the model formula, event time is the response. The \code{weights} vector provides 
the censoring information (0 for censoring, 1 for event). Deals with the case in which each subject has
one event/censoring time and one row of covariate values. When each subject has several time dependent
covariates see \code{\link{cox.pht}}.

}

\usage{
cox.ph(link="identity")
}
\arguments{
\item{link}{currently (and possibly for ever) only \code{"identity"} supported.}

}
\value{
 An object inheriting from class \code{general.family}.
}

\details{Used with \code{\link{gam}} to fit Cox Proportional Hazards models to survival data. The model formula will 
have event/censoring times on the left hand side and the linear predictor specification on the right hand side. Censoring
information is provided by the \code{weights} argument to \code{gam}, with 1 indicating an event and 0 indicating 
censoring. 

Prediction from the fitted model object (using the \code{predict} method) with \code{type="response"} will predict on the 
survivor function scale. See example code below for extracting the cumulative baseline hazard/survival directly.
Martingale or deviance 
residuals can be extracted. The \code{fitted.values} stored in the model object are survival function estimates for each 
subject at their event/censoring time. 

Estimation of model coefficients is by maximising the log-partial likelihood penalized by the smoothing penalties. See e.g. 
Hastie and Tibshirani, 1990, section 8.3. for the partial likelihood used (with Peto's approximation for ties), but note that 
optimization of the partial likelihood does not follow Hastie and Tibshirani. See Klein amd Moeschberger (2003) for 
estimation of residuals, the cumulative baseline hazard, survival function and associated standard errors (the survival standard error expression has a typo).  

The percentage deviance explained reported for Cox PH models is based on the sum of squares of the deviance residuals, as the model deviance, and the sum of squares of the deviance residuals when the covariate effects are set to zero, as the null deviance. The same baseline hazard estimate is used for both.

This family deals efficiently with the case in which each subject has one event/censoring time and one row of covariate values. For studies in which there are multiple time varying covariate measures for each subject then the equivalent Poisson model should be fitted to suitable pseudodata using \code{bam(...,discrete=TRUE)}. 
}

\references{
Hastie and Tibshirani (1990) Generalized Additive Models, Chapman and Hall.

Klein, J.P and Moeschberger, M.L. (2003) Survival Analysis: Techniques for
   Censored and Truncated Data (2nd ed.) Springer.

Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association.
\url{http://arxiv.org/abs/1511.03864}
}


\examples{
library(mgcv)
library(survival) ## for data
col1 <- colon[colon$etype==1,] ## concentrate on single event
col1$differ <- as.factor(col1$differ)
col1$sex <- as.factor(col1$sex)

b <- gam(time~s(age,by=sex)+sex+s(nodes)+perfor+rx+obstruct+adhere,
         family=cox.ph(),data=col1,weights=status)

summary(b) 

plot(b,pages=1,all.terms=TRUE) ## plot effects

plot(b$linear.predictors,residuals(b))

## plot survival function for patient j...

np <- 300;j <- 6
newd <- data.frame(time=seq(0,3000,length=np))
dname <- names(col1)
for (n in dname) newd[[n]] <- rep(col1[[n]][j],np)
newd$time <- seq(0,3000,length=np)
fv <- predict(b,newdata=newd,type="response",se=TRUE)
plot(newd$time,fv$fit,type="l",ylim=c(0,1),xlab="time",ylab="survival")
lines(newd$time,fv$fit+2*fv$se.fit,col=2)
lines(newd$time,fv$fit-2*fv$se.fit,col=2)

## crude plot of baseline survival...

plot(b$family$data$tr,exp(-b$family$data$h),type="l",ylim=c(0,1),
     xlab="time",ylab="survival")
lines(b$family$data$tr,exp(-b$family$data$h + 2*b$family$data$q^.5),col=2)
lines(b$family$data$tr,exp(-b$family$data$h - 2*b$family$data$q^.5),col=2)
lines(b$family$data$tr,exp(-b$family$data$km),lty=2) ## Kaplan Meier

## Simple simulated known truth example...
ph.weibull.sim <- function(eta,gamma=1,h0=.01,t1=100) { 
  lambda <- h0*exp(eta)
  n <- length(eta)
  U <- runif(n)
  t <- (-log(U)/lambda)^(1/gamma)
  d <- as.numeric(t <= t1)
  t[!d] <- t1
  list(t=t,d=d)
}
n <- 500;set.seed(2)
x0 <- runif(n, 0, 1);x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1);x3 <- runif(n, 0, 1)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f3 <- function(x) 0*x
f <- f0(x0) + f1(x1) + f2(x2)
g <- (f-mean(f))/5
surv <- ph.weibull.sim(g)
surv$x0 <- x0;surv$x1 <- x1;surv$x2 <- x2;surv$x3 <- x3

b <- gam(t~s(x0)+s(x1)+s(x2,k=15)+s(x3),family=cox.ph,weights=d,data=surv)

plot(b,pages=1)

}
\keyword{models} \keyword{regression}%-- one or more ..



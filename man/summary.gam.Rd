\name{summary.gam}
\alias{summary.gam}
\alias{print.summary.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summary for a GAM fit}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and produces various useful
summaries from it. (See \code{\link{sink}} to divert output to a file.)
}
\usage{
\method{summary}{gam}(object, dispersion=NULL, freq=FALSE, re.test=TRUE, ...)

\method{print}{summary.gam}(x,digits = max(3, getOption("digits") - 3), 
                  signif.stars = getOption("show.signif.stars"),...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{ a fitted \code{gam} object as produced by \code{gam()}.}

\item{x}{a \code{summary.gam} object produced by \code{summary.gam()}.} 

\item{dispersion}{A known dispersion parameter. \code{NULL} to use estimate or
                  default (e.g. 1 for Poisson).}

\item{freq}{By default p-values for parametric terms are calculated using the Bayesian estimated
covariance matrix of the parameter estimators. If this is set to \code{TRUE} then
the frequentist covariance matrix of the parameters is used instead. }

\item{re.test}{Should tests be performed for random effect terms (including any term with a zero dimensional null space)?
For large models these tests can be computationally expensive. }

\item{digits}{controls number of digits printed in output.}

\item{signif.stars}{Should significance stars be printed alongside output.}

\item{...}{ other arguments.}
}

\details{ Model degrees of freedom are taken as the trace of the influence (or
hat) matrix \eqn{ {\bf A}}{A} for the model fit.
Residual degrees of freedom are taken as number of data minus model degrees of
freedom. 
Let \eqn{ {\bf P}_i}{P_i} be the matrix 
giving the parameters of the ith smooth when applied to the data (or pseudodata in the generalized case) and let \eqn{ {\bf X}}{X} 
be the design matrix of the model. Then \eqn{ tr({\bf XP}_i )}{tr(XP_i)} is the edf for the ith term. Clearly this 
definition causes the edf's to add up properly! An alternative version of EDF is more appropriate for p-value computation, and is based on the trace of  \eqn{ 2{\bf A} - {\bf AA}}{2A - AA}.

\code{print.summary.gam} tries to print various bits of summary information useful for term selection in a pretty way.

P-values for smooth terms are usually based on a 
test statistic motivated by an extension of Nychka's (1988) analysis of the frequentist properties 
of Bayesian confidence intervals for smooths (Marra and Wood, 2012). 
These have better frequentist performance (in terms of power and distribution under the null) 
than the alternative strictly frequentist approximation. When the Bayesian intervals have good 
across the function properties then the p-values have close to the correct null distribution 
and reasonable power (but there are no optimality results for the power). Full details are in Wood (2013b), 
although what is computed is actually a slight variant in which the components of the test statistic are weighted by the iterative fitting weights. 

Note that for terms with no unpenalized terms (such as Gaussian random effects) the Nychka (1988) requirement for smoothing bias to be substantially less than variance breaks down (see e.g. appendix of Marra and Wood, 2012), and this results in incorrect null distribution for p-values computed using the above approach. In this case it is necessary to use an alternative 
approach designed for random effects variance components, and this is done. See Wood (2013a) for details: the test is based on a likelihood ratio statistic (with the reference distribution appropriate for the null hypothesis on the boundary of the parameter space).

All p-values are computed without considering uncertainty in the smoothing parameter estimates.

In simulations the p-values have best behaviour under ML smoothness selection, with REML coming second. In general the p-values behave well, but neglecting smoothing parameter uncertainty means that they may be somewhat too low when smoothing parameters are highly uncertain. High uncertainty happens in particular when smoothing parameters are poorly identified, which can occur with nested smooths or highly correlated covariates (high concurvity).

By default the p-values for parametric model terms are also based on Wald tests using the Bayesian 
covariance matrix for the coefficients. This is appropriate when there are "re" terms present, and is 
otherwise rather similar to the results using the frequentist covariance matrix (\code{freq=TRUE}), since 
the parametric terms themselves are usually unpenalized. Default P-values for parameteric terms that are 
penalized using the \code{paraPen} argument will not be good. However if such terms represent conventional 
random effects with full rank penalties, then setting \code{freq=TRUE} is appropriate. 

}

\value{\code{summary.gam} produces a list of summary information for a fitted \code{gam} object. 

\item{p.coeff}{is an array of estimates of the strictly parametric model coefficients.}

\item{p.t}{is an array of the \code{p.coeff}'s divided by their standard errors.}

\item{p.pv}{is an array of p-values for the null hypothesis that the corresponding parameter is zero. 
Calculated with reference to the t distribution with the estimated residual
degrees of freedom for the model fit if the dispersion parameter has been
estimated, and the standard normal if not.}

\item{m}{The number of smooth terms in the model.}

\item{chi.sq}{An array of test statistics for assessing the significance of
model smooth terms. See details.}

\item{s.pv}{An array of approximate p-values for the null hypotheses that each
smooth term is zero. Be warned, these are only approximate.}

\item{se}{array of standard error estimates for all parameter estimates.}

\item{r.sq}{The adjusted r-squared for the model. Defined as the proportion of variance explained, where original variance and 
residual variance are both estimated using unbiased estimators. This quantity can be negative if your model is worse than a one 
parameter constant model, and can be higher for the smaller of two nested models! The proportion null deviance 
explained is probably more appropriate for non-normal errors. Note that \code{r.sq} does not include any offset in the one parameter model.}

\item{dev.expl}{The proportion of the null deviance explained by the model. The null deviance is computed taking account of any offset, so 
\code{dev.expl} can be substantially lower than \code{r.sq} when an offset is present.}

\item{edf}{array of estimated degrees of freedom for the model terms.}

\item{residual.df}{estimated residual degrees of freedom.}

\item{n}{number of data.}

\item{np}{number of model coefficients (regression coefficients, not smoothing parameters or other parameters of likelihood).}

\item{rank}{apparent model rank.}

\item{method}{The smoothing selection criterion used.}

\item{sp.criterion}{The minimized value of the smoothness selection criterion. Note that for ML and REML methods, 
what is reported is the negative log marginal likelihood or negative log restricted likelihood. }

\item{scale}{estimated (or given) scale parameter.}

\item{family}{the family used.}

\item{formula}{the original GAM formula.}

\item{dispersion}{the scale parameter.}

\item{pTerms.df}{the degrees of freedom associated with each parametric term
(excluding the constant).}

\item{pTerms.chi.sq}{a Wald statistic for testing the null hypothesis that the
each parametric term is zero.}

\item{pTerms.pv}{p-values associated with the tests that each term is
zero. For penalized fits these are approximate. The reference distribution 
is an appropriate chi-squared when the
scale parameter is known, and is based on an F when it is not.}

\item{cov.unscaled}{The estimated covariance matrix of the parameters (or
estimators if \code{freq=TRUE}), divided
by scale parameter.}

\item{cov.scaled}{The estimated covariance matrix of the parameters
(estimators if \code{freq=TRUE}).}

\item{p.table}{significance table for parameters}

\item{s.table}{significance table for smooths}

\item{p.Terms}{significance table for parametric model terms}
}

\references{


Marra, G and S.N. Wood (2012) Coverage Properties of Confidence Intervals for Generalized Additive
Model Components. Scandinavian Journal of Statistics, 39(1), 53-74.

Nychka (1988) Bayesian Confidence Intervals for Smoothing Splines. 
Journal of the American Statistical Association 83:1134-1143.

Wood, S.N. (2013a) A simple test for random effects in regression models. Biometrika 100:1005-1010

Wood, S.N. (2013b) On p-values for smooth components of an extended generalized additive model. Biometrika 100:221-228

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman
and Hall/CRC Press.
}
\author{ Simon N. Wood \email{simon.wood@r-project.org} with substantial
improvements by Henric Nilsson.}

\section{WARNING }{ 
The p-values are approximate and neglect smoothing parameter uncertainty. They are likely to be somewhat too low 
when smoothing parameter estimates are highly uncertain: do read the details section. If the exact values matter,
read Wood (2013a or b). 

P-values for terms penalized via `paraPen' are unlikely to be correct.
} 

\seealso{  \code{\link{gam}}, \code{\link{predict.gam}},
\code{\link{gam.check}}, \code{\link{anova.gam}}, \code{\link{gam.vcomp}}, \code{\link{sp.vcov}} }

\examples{
library(mgcv)
set.seed(0)

dat <- gamSim(1,n=200,scale=2) ## simulate data

b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
plot(b,pages=1)
summary(b)

## now check the p-values by using a pure regression spline.....
b.d <- round(summary(b)$edf)+1 ## get edf per smooth
b.d <- pmax(b.d,3) # can't have basis dimension less than 3!
bc<-gam(y~s(x0,k=b.d[1],fx=TRUE)+s(x1,k=b.d[2],fx=TRUE)+
        s(x2,k=b.d[3],fx=TRUE)+s(x3,k=b.d[4],fx=TRUE),data=dat)
plot(bc,pages=1)
summary(bc)

## Example where some p-values are less reliable...
dat <- gamSim(6,n=200,scale=2)
b <- gam(y~s(x0,m=1)+s(x1)+s(x2)+s(x3)+s(fac,bs="re"),data=dat)
## Here s(x0,m=1) can be penalized to zero, so p-value approximation
## cruder than usual...
summary(b) 

## p-value check - increase k to make this useful!
k<-20;n <- 200;p <- rep(NA,k)
for (i in 1:k)
{ b<-gam(y~te(x,z),data=data.frame(y=rnorm(n),x=runif(n),z=runif(n)),
         method="ML")
  p[i]<-summary(b)$s.p[1]
}
plot(((1:k)-0.5)/k,sort(p))
abline(0,1,col=2)
ks.test(p,"punif") ## how close to uniform are the p-values?

## A Gamma example, by modify `gamSim' output...
 
dat <- gamSim(1,n=400,dist="normal",scale=1)
dat$f <- dat$f/4 ## true linear predictor 
Ey <- exp(dat$f);scale <- .5 ## mean and GLM scale parameter
## Note that `shape' and `scale' in `rgamma' are almost
## opposite terminology to that used with GLM/GAM...
dat$y <- rgamma(Ey*0,shape=1/scale,scale=Ey*scale)
bg <- gam(y~ s(x0)+ s(x1)+s(x2)+s(x3),family=Gamma(link=log),
          data=dat,method="REML")
summary(bg)

}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







\name{summary.gam}
\alias{summary.gam}
\alias{print.summary.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summary for a GAM fit}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and produces various useful
summaries from it. 
}
\usage{
summary.gam(object,...)
print.summary.gam(x,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{ a fitted \code{gam} object as produced by \code{gam()}.}
\item{x}{a \code{summary.gam} object produced by \code{summary.gam()}.} 
\item{...}{ other arguments.}
}
\details{ Model degrees of freedom are taken as the trace of the influence (or hat) matrix \eqn{ {\bf A}}{A} for the model fit.
Residual degrees of freedom are taken as number of data minus model degrees of freedom. Let \eqn{ {\bf P}_i}{P_i} be the matrix 
giving the parameters of the ith smooth when applied to the data (or pseudodata in the generalized case) and let \eqn{ {\bf X}}{X} 
be the design matrix of the model. Then \eqn{ tr({\bf XP}_i )}{tr(XP_i)} is the edf for the ith term. Clearly this definition causes 
the edf's to add up properly! 

\code{print.summary.gam} tries to print various bits of summary information useful for term selection in a pretty way.
}

\value{\code{summary.gam} produces a list of summary information for a fitted \code{gam} object. 
\item{p.coeff}{is an array of estimates of the strictly parametric model coefficients.}
\item{p.t}{is an array of the \code{p.coeff}'s divided by their standard errors.}
\item{p.pv}{is an array of p-values for the null hypothesis that the corresponding parameter is zero. 
Calculated with reference to the t distribution with the estimated residual degrees of freedom for the model fit.}
\item{m}{The number of smooth terms in the model.}
\item{chi.sq}{An array of test statistics for assessing the significance of model smooth terms. If \eqn{ {\bf p}_i}{p_i} 
is the parameter vector for the ith smooth term, and this term has estimated covariance matrix \eqn{ {\bf V}_i}{V_i} then the 
statistic is \eqn{ {\bf p}_i^\prime {\bf V}_i^{k-} {\bf p}_i}{p_i'V_i^{k-}p_i}, where \eqn{ {\bf V}^{k-}_i}{V_i^{k-}} is the rank k-1 
pseudo-inverse of \eqn{ {\bf V_i}}{V_i}, and k is the basis dimension.}
\item{s.pv}{An array of approximate p-values for the null hypotheses that each smooth term is zero. Be warned, these are only 
approximate. In the case in which UBRE has been used, they are obtained by comparing the chi.sq statistic given above to the 
chi-squared distribution with degrees of  freedom given by the estimated degrees of freedom for the term. In the GCV case (in 
which the scale parameter will have been estimated) the statistic is compared to an F distribution with upper d.f. given by
the estimate degrees of freedom for the term, and lower d.f. given by the residual degrees of freedom for the model . 
Use at your own risk! Typically the p-values will be somewhat inaccurate
because they are conditional on the smoothing parameters and the distributional assumption doesn't have a firm theoretical basis. 
A pragmatic approach to the latter issue is to check p-values by refitting the model using regression splines
with each basis dimension set to one more than the rounded edf for
the term (see example, below). In this latter case the distributional
assumption would be fine if the smoothing parameters were
known, but of course they are not, and conditioning on the smoothing
parameters will always be problematic. The difficulty is that GCV has
used the data to select the most plausible smoothing parameters for each
term - not surprisingly this tends to mean that the terms appear
a little more `significant' than they should. Of course this problem
is no different to the standard difficulty in interpreting p-values for
model terms when the model has been selected by hypothesis testing
methods such as backwards elimination.   
}
\item{se}{array of standard error estimates for all parameter estimates.}
\item{r.sq}{The adjusted r-squared for the model. Defined as the proportion of variance explained, where original variance and 
residual variance are both estimated using unbiased estimators. This quantity can be negative if your model is worse than a one 
parameter constant model, and can be higher for the smaller of two nested models! Note that proportion null deviance 
explained is probably more appropriate for non-normal errors.}
\item{dev.expl}{The proportion of the null deviance explained by the model.}
\item{edf}{array of estimated degrees of freedom for the model terms.}
\item{residual.df}{estimated residual degrees of freedom.}
\item{n}{number of data.}
\item{gcv}{minimized GCV score for the model, if GCV used.}
\item{ubre}{minimized UBRE score for the model, if UBRE used.}
\item{scale}{estimated (or given) scale parameter.}
\item{family}{the family used.}
\item{formula}{the original GAM formula.}
\item{dispersion}{the scale parameter.}
\item{pTerms.df}{the degrees of freedom associated with each parameteric term
(excluding the constant).}
\item{pTerms.chi.sq}{a Wald statistic for testing the null hypothesis that the
each parametric term is zero.}
\item{pTerms.pv}{p-values associated with the tests that each term is
zero. For penalized fits these are approximate, being conditional on the smoothing
parameters. The reference distribution is an appropriate chi-squared when the
scale parameter is known, and is based on an F when it is not.}
}

\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

Wood, S.N. (2004a) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass.

Wood, S.N. (2004b) On confidence intervals for GAMs based on penalized
regression splines. Technical Report 04-12 Department of Statistics,
University of Glasgow.


\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\section{WARNING }{ The supplied p-values are only approximate and should be
treated with scepticism, particularly for the smooth terms.
} 

\seealso{  \code{\link{gam}}, \code{\link{predict.gam}},
\code{\link{gam.check}}, \code{\link{anova.gam}} }

\examples{
library(mgcv)
set.seed(0)
n<-200
sig2<-4
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
pi <- asin(1) * 2
y <- 2 * sin(pi * x0)
y <- y + exp(2 * x1) - 3.75887
y <- y + 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10 - 1.396
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- y + e
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3))
plot(b,pages=1)
summary(b)
# now check the p-values by using a pure regression spline.....
b.d<-round(b$edf)+1 
b.d<-pmax(b.d,3) # can't have basis dimension less than this!
bc<-gam(y~s(x0,k=b.d[1],fx=TRUE)+s(x1,k=b.d[2],fx=TRUE)+
        s(x2,k=b.d[3],fx=TRUE)+s(x3,k=b.d[4],fx=TRUE))
plot(bc,pages=1)
summary(bc)
# p-value check - increase k to make this useful!
n<-200;p<-0;k<-20
for (i in 1:k)
{ b<-gam(y~s(x,z),data=data.frame(y=rnorm(n),x=runif(n),z=runif(n)))
  p[i]<-summary(b)$s.p[1]
}
plot(((1:k)-0.5)/k,sort(p))
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







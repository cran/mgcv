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
\item{x}{a \code{summary.gam} object prooduced by \code{summary.gam()}.} 
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
Calculated with reference to the t distribution with the estimated residual degrees of freedom for the model fit (rounded, of course).}
\item{m}{The number of smooth terms in the model.}
\item{chi.sq}{An array of test statistics for assessing the significance of model smooth terms. If \eqn{ {\bf p}_i}{p_i} 
is the parameter vector for the ith smooth term, and this term has estimated covariance matrix \eqn{ {\bf V}_i}{V_i} then the 
statistic is \eqn{ {\bf p}_i^\prime {\bf V}_i^{k-} {\bf p}_i}{p_i'V_i^{k-}p_i}, where \eqn{ {\bf V}^{k-}_i}{V_i^{k-}} is the rank k-1 
pseudo-inverse of \eqn{ {\bf V_i}}{V_i}, and k is the basis dimension for the smooth. The pseudoinverse is truncated to the numerical rank of 
\eqn{ {\bf V_i}}{V_i} if this is less than k-1.}
\item{s.pv}{An array of approximate p-values for the null hypotheses that each smooth term is zero. Be warned, these are only 
approximate: they are obtained by comparing the chi.sq statistic given above to the chi-squared distribution with degrees of 
freedom closest to the estimated degrees of freedom for the term. Use at your own risk (although limited simulation studies 
seem to indicate that these are more reasonable than you might reasonably expect).}
\item{se}{array of standard error estimates for all parameter estimates.}
\item{r.sq}{The adjusted r-squared for the model. Defined as the proportion of variance explained, where original variance and 
residual variance are both estimated using unbiased estimators. This quantity can be negative if your model is worse than a one 
parameter constant model, and can be higher for the smaller of two nested models!}
\item{edf}{array of estimated degrees of freedom for the model terms.}
\item{residual.df}{estimated residual degrees of freedom.}
\item{n}{number of data.}
\item{gcv}{GCV score for the model.}
\item{scale}{estimated (or given) scale parameter.}
\item{family}{the family used.}
\item{formula}{the original GAM formula.}
}

\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428


\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}}

\section{WARNING }{ The supplied p-values are only approximate.
} 

\seealso{  \code{\link{gam}} \code{\link{predict.gam}}}

\examples{

library(mgcv)
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
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







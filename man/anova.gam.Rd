\name{anova.gam}
\alias{anova.gam}
\alias{print.anova.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Hypothesis tests related to  GAM fits}
\description{ Performs hypothesis tests relating to one or more fitted
\code{gam} objects. For a single fitted \code{gam} object, Wald tests of
the significance of each parametric and smooth term are performed. Otherwise
the fitted models are compared using an analysis of deviance table. The tests
are usually approximate, unless the models are un-penalized.

}
\usage{
anova.gam(object, ..., dispersion = NULL, test = NULL)
print.anova.gam(x, digits = max(3, getOption("digits") - 3),...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object,...}{ fitted model objects of class \code{gam} as produced by \code{gam()}.}
\item{x}{an \code{anova.gam} object produced by a single model call to \code{anova.gam()}.} 
\item{dispersion}{ a value for the dispersion parameter: not normally used.}
\item{test}{what sort of test to perform for a multi-model call. One of
\code{"Chisq"}, \code{"F"} or \code{"Cp"}. }
\item{digits}{number of digits to use when printing output.}
}
\details{ If more than one fitted model is provided than \code{anova.glm} is
used. If only one model is provided then the significance of each model term
is assessed using Wald tests: see \code{\link{summary.gam}} for details of the
actual computations. In the latter case \code{print.anova.gam} is used as the
printing method.

P-values are usually reliable if the smoothing parameters are known, or
the model is unpenalized. If smoothing parameters have been estimated then the
p-values are typically somewhat too low under the null. This occurs because the uncertainty associated
with the smoothing parameters is neglected in the calculations of the
distributions under the null, which tends to lead to underdispersion in these
distributions, and in turn to p-value estimates that are too low. (In
simulations where the null is correct, I have seen p-values that are as low as half of what they should
be.) Note however that tests can have low power if the estimated rank of the test
statistic is much higher than the EDF, so that p-values can also be too low in
some cases. 

If it is important to have  p-values that are as accurate as possible, then,
at least in the single model case, it is probably advisable to perform tests using unpenalized smooths
(i.e. \code{s(...,fx=TRUE)}) with the basis dimension, \code{k}, left at what would
have been used with penalization. Such tests are not as powerful, of
course, but the p-values are more accurate under the null. Whether or not extra accuracy is
required will usually depend on whether or not hypothesis testing is a key
objective of the analysis.
}

\value{In the multi-model case \code{anova.gam} produces output identical to
\code{\link{anova.glm}}, which it in fact uses.

In the single model case an object of class \code{anova.gam} is produced,
which is in fact an object returned from \code{\link{summary.gam}}.

\code{print.anova.gam} simply produces tabulated output.

}


\author{ Simon N. Wood \email{simon.wood@r-project.org} with substantial
improvements by Henric Nilsson.}

\section{WARNING}{ P-values are only approximate, particularly as a result of ignoring
smoothing parameter uncertainty.
} 

\seealso{  \code{\link{gam}}, \code{\link{predict.gam}},
\code{\link{gam.check}}, \code{\link{summary.gam}} }

\examples{
library(mgcv)
set.seed(0)
n<-200
sig<-2
x0 <- rep(1:4,50)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
y <- 2 * x0
y <- y + exp(2 * x1)
y <- y + 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10
e <- rnorm(n, 0, sig)
y <- y + e
x0<-as.factor(x0)
b<-gam(y~x0+s(x1)+s(x2)+s(x3))
anova(b)
b1<-gam(y~x0+s(x1)+s(x2))
anova(b,b1,test="F")
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







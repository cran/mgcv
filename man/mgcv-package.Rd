\name{mgcv.package}
\alias{mgcv.package}
\alias{mgcv-package}
\alias{mgcv}
\docType{package}
\title{Mixed GAM Computation Vehicle with GCV/AIC/REML smoothness estimation and GAMMs by REML/PQL
}
\description{
\code{mgcv} provides functions for generalized additive modelling (\code{\link{gam}} and \code{\link{bam}})  and
generalized additive mixed modelling (\code{\link{gamm}}, and \code{\link{random.effects}}). The term GAM is taken to include 
any model dependent on unknown smooth functions of predictors and estimated by quadratically penalized (possibly quasi-) likelihood maximization. Available distributions are covered in \code{\link{family.mgcv}} and available smooths in \code{\link{smooth.terms}}.

Particular features of the package are facilities for automatic smoothness selection (Wood, 2004, 2011), 
and the provision of a variety of smooths of more than one variable. User defined 
smooths can be added. A Bayesian approach to confidence/credible interval calculation is
provided. Linear functionals of smooths, penalization of parametric model terms and linkage 
of smoothing parameters are all supported. Lower level routines for generalized ridge 
regression and penalized linearly constrained least squares are also available.
}

\details{ \code{mgcv} provides generalized additive modelling functions \code{\link{gam}},
\code{\link{predict.gam}} and \code{\link{plot.gam}}, which are very similar
in use to the S functions of the same name designed by Trevor Hastie (with some extensions). 
However the underlying representation and estimation of the models is based on a
penalized regression spline approach, with automatic smoothness selection. A
number of other functions  such as \code{\link{summary.gam}} and \code{\link{anova.gam}} 
are also provided, for extracting information from a fitted \code{\link{gamObject}}.

Use of \code{\link{gam}} is much like use of \code{\link{glm}}, except that
within a \code{gam} model formula, isotropic smooths of any number of predictors can be specified using
\code{\link{s}} terms, while scale invariant smooths of any number of
predictors can be specified using \code{\link{te}},  \code{\link{ti}} or \code{\link{t2}}  terms. 
\code{\link{smooth.terms}} provides an 
overview of the built in smooth classes, and \code{\link{random.effects}} should be refered to for an overview 
of random effects terms (see also \code{\link{mrf}} for Markov random fields). Estimation is by
penalized likelihood or quasi-likelihood maximization, with smoothness
selection by GCV, GACV, gAIC/UBRE or (RE)ML. See \code{\link{gam}}, \code{\link{gam.models}}, 
\code{\link{linear.functional.terms}} and \code{\link{gam.selection}} for some discussion of model specification and
selection. For detailed control of fitting see \code{\link{gam.convergence}},
\code{\link{gam}} arguments \code{method} and \code{optimizer} and \code{\link{gam.control}}. For checking and
visualization see \code{\link{gam.check}}, \code{\link{choose.k}}, \code{\link{vis.gam}} and \code{\link{plot.gam}}.
While a number of types of smoother are built into the package, it is also
extendable with user defined smooths, see \code{\link{smooth.construct}}, for example.

A Bayesian approach to smooth modelling is used to derive standard errors on
predictions, and hence credible intervals (see Marra and Wood, 2012). The Bayesian covariance matrix for
the model coefficients is returned in \code{Vp} of the
\code{\link{gamObject}}. See \code{\link{predict.gam}} for examples of how
this can be used to obtain credible regions for any quantity derived from the
fitted model, either directly, or by direct simulation from the posterior
distribution of the model coefficients. Approximate p-values can also be obtained for testing 
individual smooth terms for equality to the zero function, using similar ideas (see Wood, 2013a,b). Frequentist
approximations can be used for hypothesis testing based model comparison. See \code{\link{anova.gam}} and
\code{\link{summary.gam}} for more on hypothesis testing.

For large datasets (that is large n) see \code{\link{bam}} which is a version of \code{\link{gam}} with 
a much reduced memory footprint.

The package also provides a generalized additive mixed modelling function,
\code{\link{gamm}}, based on a PQL approach and  
\code{lme} from the \code{nlme} library (for an \code{lme4} based version, see package \code{gamm4}). 
\code{gamm} is particularly useful
for modelling correlated data (i.e. where a simple independence model for the
residual variation is inappropriate). In addition, low level routine \code{\link{magic}}
can fit models to data with a known correlation structure.


Some underlying GAM fitting methods are available as low level fitting
functions: see \code{\link{magic}}. But there is little functionality 
that can not be more conventiently accessed via \code{\link{gam}} . 
Penalized weighted least squares with linear equality and inequality constraints is provided by 
\code{\link{pcls}}.

For a complete list of functions type \code{library(help=mgcv)}. See also \code{\link{mgcv.FAQ}}.

}
\author{
Simon Wood <simon.wood@r-project.org> 

with contributions and/or help from Natalya Pya, Thomas Kneib, Kurt Hornik, Mike Lonergan, Henric Nilsson,
Fabian Scheipl and Brian Ripley. 

Polish translation - Lukasz Daniel; German translation - Chris Leick, Detlef Steuer; 
French Translation - Philippe Grosjean

Maintainer: Simon Wood <simon.wood@r-project.org>

Part funded by EPSRC: EP/K005251/1
}
\references{
These provide details for the underlying mgcv methods, and fuller 
references to the large literature on which the methods are based. 

Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models (with discussion).
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}

Wood, S.N. (2011) Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

Wood, S.N. (2004) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass. 99:673-686. 

Marra, G and S.N. Wood (2012) Coverage Properties of Confidence Intervals for Generalized Additive
Model Components. Scandinavian Journal of Statistics, 39(1), 53-74.

Wood, S.N. (2013a) A simple test for random effects in regression models. Biometrika 100:1005-1010

Wood, S.N. (2013b) On p-values for smooth components of an extended generalized additive model. Biometrika 100:221-228

Wood, S.N. (2017) \emph{Generalized Additive Models: an introduction with R (2nd edition)},
CRC

Development of mgcv version 1.8 was part funded by EPSRC grants EP/K005251/1 and EP/I000917/1.

}
\keyword{ package }
\keyword{models} 
\keyword{smooth} 
\keyword{regression}


\examples{
## see examples for gam and gamm
}

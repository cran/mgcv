\name{gamObject}
\alias{gamObject}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fitted gam object}
\description{A fitted GAM object returned by function \code{gam} and of class
\code{"gam"} inheriting from classes \code{"glm"} and \code{"lm"}. Method
functions \code{anova}, \code{logLik}, \code{influence}, \code{plot},
\code{predict}, \code{print}, \code{residuals} and \code{summary} exist for
this class.

All compulsory elements of \code{"glm"} and \code{"lm"} objects are present,
but the fitting method for a GAM is different to a linear model or GLM, so
that the elements relating to the QR decomposition of the model matrix are
absent.
}

\value{ A \code{gam} object has the following elements:

\item{aic}{AIC of the fitted model: bear in mind that the degrees of freedom
used to calculate this are the effective degrees of freedom of the model, and
the likelihood is evaluated at the maximum of the penalized likelihood in most
cases, not at the MLE.}

\item{assign}{Array whose elements indicate which model term (listed in
\code{pterms}) each parameter relates to: applies only to non-smooth terms.}

\item{boundary}{did parameters end up at boundary of parameter space?} 

\item{call}{the matched call (allows \code{update} to be used with \code{gam} objects, for example). }

\item{coefficients}{the coefficients of the fitted model. Parametric
          coefficients are  first, followed  by coefficients for each
          spline term in turn.}

\item{control}{the \code{gam} control list used in the fit.}

\item{converged}{indicates whether or not the iterative fitting method converged.} 

\item{data}{the original supplied data argument (for class \code{"glm"} compatibility).}

\item{deviance}{model deviance (not penalized deviance).}

\item{df.null}{null degrees of freedom.} 

\item{df.residual}{effective residual degrees of freedom of the model.}

\item{edf}{estimated degrees of freedom for each model parameter. Penalization
means that many of these are less than 1.}

\item{family}{family object specifying distribution and link used.}

\item{fit.method}{ Character string describing the multiple GCV/UBRE smoothing
parameter estimation method used.}

\item{fitted.values}{fitted model predictions of expected value for each
          datum.}

\item{formula}{the model formula.}

\item{full.formula}{the model formula with each smooth term fully
  expanded and with option arguments given explicitly (i.e. not with reference
  to other variables) - useful for later prediction from the model.}

\item{gcv.ubre}{The minimized GCV or UBRE score.}


\item{hat}{array of elements from the leading diagonal of the `hat' (or `influence') matrix. 
Same length as response data vector.}

\item{iter}{number of iterations of P-IRLS taken to get convergence.}

\item{linear.predictors}{fitted model prediction of link function of
expected value for  each datum.}

\item{method}{One of \code{"GCV"} or \code{"UBRE"}, depending on the fitting
criterion used.}

\item{mgcv.conv}{ A list of convergence diagnostics relating to the
\code{"mgcv"} or \code{"magic"} parts of smoothing
parameter estimation - this will not be very meaningful for pure \code{"outer"}
estimation of smoothing parameters. \code{mgcv.conv} differs for method \code{"magic"} and \code{"mgcv"}. Here is 
the \code{"mgcv"} version:

\item{score}{corresponding to edf, an array of GCV or UBRE scores for the model given the final 
     estimated relative smoothing parameters.}

\item{g}{the gradient of the GCV/UBRE score w.r.t. the relative smoothing parameters at termination.}

\item{h}{the second derivatives corresponding to \code{g} above - i.e. the leading diagonal of the Hessian.}

\item{e}{the eigen-values of the Hessian. All non-negative indicates a positive definite Hessian.}

\item{iter}{the number of iterations taken.}

\item{in.ok}{\code{TRUE} if the second smoothing parameter guess improved the GCV/UBRE score.}

\item{step.fail}{\code{TRUE} if the algorithm terminated by failing to improve the GCV/UBRE score rather than by `converging'. 
Not necessarily a problem, but check the above derivative information quite carefully.}

In the case of \code{"magic"} the items are:

\item{full.rank}{The apparent rank of the problem given the model matrix and 
constraints.}

\item{rank}{The numerical rank of the problem.}

\item{fully.converged}{\code{TRUE} is multiple GCV/UBRE converged by meeting 
convergence criteria. \code{FALSE} if method stopped with a steepest descent step 
failure.}

\item{hess.pos.def}{Was the hessian of the GCV/UBRE score positive definite at 
smoothing parameter estimation convergence?}

\item{iter}{How many iterations were required to find the smoothing parameters?}

\item{score.calls}{and how many times did the GCV/UBRE score have to be evaluated?}

\item{rms.grad}{root mean square of the gradient of the GCV/UBRE score at 
convergence.}

} % end of mgcv.conv listing

\item{min.edf}{Minimum possible degrees of freedom for whole model.}

\item{model}{model frame containing all variables needed in original model fit.}

\item{na.action}{The \code{\link{na.action}} used in fitting.}

\item{nsdf}{number of parametric, non-smooth, model terms including the
          intercept.}

\item{null.deviance}{deviance for single parameter model.} 

\item{offset}{model offset.}

\item{outer.info}{If `outer' iteration has been used to fit the model (see
\code{\link{gam.method}}) then this is present and contains whatever was
returned by the optimization routine used (currently \code{\link{nlm}} or \code{\link{optim}}). }

\item{prior.weights}{prior weights on observations.} 

\item{pterms}{\code{terms} object for strictly parametric part of model.}

\item{rank}{apparent rank of fitted model.}

\item{residuals}{the working residuals for the fitted model.}

\item{sig2}{estimated or supplied variance/scale parameter.}

\item{smooth}{list of smooth objects, containing the basis information for each term in the 
model formula in the order in which they appear. These smooth objects are what gets returned by
the \code{\link{smooth.construct}} objects.}

\item{sp}{smoothing parameter for each smooth.}

\item{terms}{\code{terms} object of \code{model} model frame.}

\item{Vp}{estimated covariance matrix for the parameters. This is a Bayesian
posterior covariance matrix that results from adopting a particular Bayesian
model of the smoothing process. Paricularly useful for creating
credible/confidence intervals.}

\item{Ve}{frequentist estimated covariance matrix for the parameter
estimators. Particularly useful for testing whether terms are zero. Not so
useful for CI's as smooths are usually biased.}

\item{weights}{final weights used in IRLS iteration.}

\item{y}{response data.}

}

\references{

Key References on this implementation:

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

Wood, S.N. (in press) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass.

Wood, S.N. (2004) On confidence intervals for GAMs based on penalized
regression splines. Technical Report 04-12 Department of Statistics,
University of Glasgow.

Wood, S.N. (2004) Low rank scale invariant tensor product smooths for
generalized additive mixed models. Technical Report 04-13 Department of Statistics,
University of Glasgow.

Key Reference on GAMs and related models:

Hastie (1993) in Chambers and Hastie (1993) Statistical Models in S. Chapman
and Hall.

Hastie and Tibshirani (1990) Generalized Additive Models. Chapman and Hall.

Wahba (1990) Spline Models of Observational Data. SIAM 

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\section{WARNINGS }{ This model object is different to that described in
Chambers and Hastie (1993) in order to allow smoothing parameter estimation etc.

} 

\seealso{\code{\link{gam}}}




\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..








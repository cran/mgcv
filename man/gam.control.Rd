\name{gam.control}
\alias{gam.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Setting GAM fitting defaults}
\description{ This is an internal function of package \code{mgcv} which allows 
control of the numerical options for fitting a GAM. 
Typically users will want to modify the defaults if model fitting fails to
converge, or if the warnings are generated which suggest a 
loss of numerical stability during fitting.  To change the default
choise of fitting method, see \code{\link{gam.method}}.
}
\usage{
gam.control(irls.reg=0.0,epsilon = 1e-06, maxit = 100,globit = 20,
            mgcv.tol=1e-7,mgcv.half=15,nb.theta.mult=10000, trace = FALSE,
            rank.tol=.Machine$double.eps^0.5,absorb.cons=TRUE,
            max.tprs.knots=5000,nlm=list(),optim=list(),newton=list(),
            outerPIsteps=1) 
}
\arguments{ 
\item{irls.reg}{For most models this should be 0. The iteratively re-weighted least squares method
by which GAMs are fitted  can fail to converge in some circumstances. For example, data with many zeroes can cause 
problems in a model with a log link, because a mean of zero corresponds to an infinite range of linear predictor 
values. Such convergence problems are caused by a fundamental lack of identifiability, but do not show up as 
lack of identifiability in the penalized linear model problems that have to be solved at each stage of iteration.
In such circumstances it is possible to apply a ridge regression penalty to the model to impose identifiability, and 
\code{irls.reg} is the size of the penalty. The penalty can not be used if the
underlying fitting method is \code{\link{mgcv}} (not the default - see
\code{\link{gam.method}} for details).
}

\item{epsilon}{This is used for judging conversion of the GLM IRLS loop in
\code{\link{gam.fit}} or \code{\link{gam.fit2}} .}

\item{maxit}{Maximum number of IRLS iterations to perform using cautious
GCV/UBRE optimization, after \code{globit} IRLS iterations with normal GCV 
optimization have been performed. Note that only fitting based on
\code{\link{mgcv}} (not default) makes 
any distinction between cautious and global optimization.}

\item{globit}{Maximum number of IRLS iterations to perform with normal 
GCV/UBRE optimization. If convergence is not achieved after these 
iterations then a further \code{maxit} iterations will be performed 
using cautious GCV/UBRE optimization. } 

\item{mgcv.tol}{The convergence tolerance parameter to use in GCV/UBRE optimization.}

\item{mgcv.half}{If a step of  the GCV/UBRE optimization method leads 
to a worse GCV/UBRE score, then the step length is halved. This is
the number of halvings to try before giving up.}

\item{nb.theta.mult}{Controls the limits on theta when negative binomial 
parameter is to be estimated. Maximum theta is set to the initial value 
multiplied by \code{nb.theta.mult}, while the minimum value is set to
the initial value divided by \code{nb.theta.mult}. }

\item{trace}{Set this to \code{TRUE} to turn on diagnostic output.}


\item{rank.tol}{The tolerance used to estimate the rank of the fitting
problem, for methods which deal with rank deficient cases (basically all
except those based on \code{\link{mgcv}}).}

\item{absorb.cons}{If \code{TRUE} then the GAM is set up using a
parameterization which requires no further constraint. Usually this means that
all the smooths are automatically centered (i.e. they sum to zero over the
covariate values). If \code{FALSE} then the ordinary parameterizations of the
smooths are used, which require constraints to be imposed during fitting.}

\item{max.tprs.knots}{This is the default initial maximum number of knots to allow
when constructing a t.p.r.s bases (\code{bs="tp"}). The set up cost (and
storage) for these smooths scales as the square of the number of initial knots, so if it's too
high you can appear to freeze R. Usually one would want to use an alternative
smoothing basis (or \code{\link{te}} terms), or the approach illustrated in the
examples in \code{\link{gam}}, rather than simply increasing this default.}

\item{nlm}{list of control parameters to pass to \code{\link{nlm}} if this is
used for outer estimation of smoothing parameters. See details.}

\item{optim}{list of control parameters to pass to \code{\link{optim}} if this
is used for outer estimation of smoothing parameters. See details.}

\item{newton}{list of control parameters to pass to default Newton optimizer
used for outer estimation of log smoothing parameters. See details.}

\item{outerPIsteps}{The number of performance interation steps used to
initialize outer iteration. Less than 1 means
that only one performance iteration step is taken to get the function scale,
but the corresponding smoothing parameter estimates are discarded. }
}

\details{ When outer iteration is used for fitting then the control list
\code{nlm} stores control arguments for calls to routine
\code{\link{nlm}}. The list has the following named elements: (i) \code{ndigit} is
the number of significant digits in the GCV/UBRE score - by default this is
worked out from \code{epsilon}; (ii) \code{gradtol} is the tolerance used to
judge convergence of the gradient of the GCV/UBRE score to zero - by default
set to \code{100*epsilon}; (iii) \code{stepmax} is the maximum allowable log
smoothing parameter step - defaults to 2; (iv) \code{steptol} is the minimum
allowable step length - defaults to 1e-4; (v) \code{iterlim} is the maximum
number of optimization steps allowed - defaults to 200; (vi)
\code{check.analyticals} indicates whether the built in exact derivative
calculations should be checked numerically - defaults to \code{FALSE}. Any of
these which are not supplied and named in the list are set to their default
values.

Outer iteration using \code{newton} is controlled by the list \code{newton}
with the following elements : \code{conv.tol} (default
1e-6) is the relative convergence tolerance; \code{maxNstep} is the maximum
length allowed for an element of the Newton search direction (default 5);
\code{maxSstep} is the maximum length allowed for an element of the steepest
descent direction (only used if Newton fails - default 2); \code{maxHalf} is
the maximum number of step halvings to permit before giving up (default 30); 
\code{use.svd} determines whether to use an SVD step or QR step for the part
of fitting that determines problem rank, etc. (defaults to \code{FALSE}).

Outer iteration using \code{\link{optim}} is controlled using list
\code{optim}, which currently has one element: \code{factr} which takes
default value 1e7.


When fitting is been done by calls to routine \code{\link{mgcv}}, 
\code{maxit} and \code{globit} control the maximum iterations of the IRLS algorithm, as follows: 
the algorithm will first execute up to
\code{globit} steps in which the GCV/UBRE algorithm performs a global search for the best overall 
smoothing parameter at every iteration. If convergence is not achieved within \code{globit} iterations, then a further 
\code{maxit} steps are taken, in which the overall smoothing parameter estimate is taken as the 
one locally minimising the GCV/UBRE score and resulting in the lowest EDF change. The difference 
between the two phases is only significant if the GCV/UBRE function develops more than one minima. 
The reason for this approach is that the GCV/UBRE score for the IRLS problem can develop `phantom'
minimima for some models: these are minima which are not present in the GCV/UBRE score of the IRLS
problem resulting from moving the parameters to the minimum! Such minima can lead to convergence 
failures, which are usually fixed by the second phase. 
}


\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 
65(1):95-114

Wood, S.N. (2004) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass.


\url{http://www.maths.bath.ac.uk/~sw283/}


}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{\code{\link{gam.method}}   \code{\link{gam}}, \code{\link{gam.fit}}, \code{\link{glm.control}} }

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





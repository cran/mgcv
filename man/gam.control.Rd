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
            max.tprs.knots=5000) 
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
}

\details{ 
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


\url{http://www.stats.gla.ac.uk/~simon/}


}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{\code{\link{gam.method}}   \code{\link{gam}}, \code{\link{gam.fit}}, \code{\link{glm.control}} }

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





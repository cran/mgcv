\name{gam.control}
\alias{gam.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Setting GAM fitting defaults}
\description{ This is an internal function of package \code{mgcv} which allows control of the numerical options for fitting a GAM. 
Typically users will want to modify the defaults if model fitting fails to converge, or if the warnings are generated which suggest a 
loss of numerical stability during fitting.  
}
\usage{
gam.control(irls.reg=0.0,epsilon = 1e-04, maxit = 20,globit = 20,
            mgcv.tol=1e-6,mgcv.half=15,nb.theta.mult=10000, trace = FALSE,
            fit.method="magic",perf.iter=TRUE,rank.tol=.Machine$double.eps^0.5) 
}
\arguments{ 
\item{irls.reg}{For most models this should be 0. The iteratively re-weighted least squares method
by which GAMs are fitted  can fail to converge in some circumstances. For example data with many zeroes can cause 
problems in a model with a log link, because a mean of zero corresponds to an infinite range of linear predictor 
values. Such convergence problems are caused by a fundamental lack of identifiability, but do not show up as 
lack of identifiability in the penalized linear model problems that have to be solved at each stage of iteration.
In such circumstances it is possible to apply a ridge regression penalty to the model to impose identifiability, and 
\code{irls.reg} is the size of the penalty. The penalty can only be used if \code{fit.method=="magic"}.
}

\item{epsilon}{This is used for judging conversion of the GLM IRLS loop in \code{gam.fit}.}

\item{maxit}{Maximum number of IRLS iterations to perform using cautious GCV/UBRE opttimization, after \code{globit} 
IRLS iterations with normal GCV optimization have been performed. Note that fit method \code{"magic"} makes no distinction 
between cautious and global optimization.}

\item{globit}{Maximum number of IRLS iterations to perform with normal GCV/UBRE optimization. If convergence is not achieved after these 
iterations then a further \code{maxit} iterations will be performed using cautious GCV/UBRE optimization. } 

\item{mgcv.tol}{The convergence tolerance parameter to use in GCV/UBRE optimization.}

\item{mgcv.half}{If a step of  the GCV/UBRE optimization method leads to a worse GCV/UBRE score, then the step length is halved. This is
the number of halvings to try before giving up.}

\item{nb.theta.mult}{Controls the limits on theta when negative binomial parameter is to be estimated. Maximum theta is set to 
the initial value multiplied by \code{nb.theta.mult}, while the minimum value is set to the initial value divided by 
\code{nb.theta.mult}. }

\item{trace}{Set this to \code{TRUE} to turn on diagnostic output.}

\item{fit.method}{set to \code{"mgcv"} to use the method described in Wood (2000). Set to \code{"magic"} to use a 
newer numerically more stable method, which allows regularization and mixtures of fixed and estimated smoothing parameters. 
Set to \code{"fastest"} to use \code{"mgcv"} for one smooth models and \code{"magic"} otherwise.}

\item{perf.iter}{set to \code{TRUE} to use Gu's approach to finding smoothing 
parameters in which GCV or UBRE is
applied to the penalized linear modelling problem produced at each IRLS iteration. This method is very fast, but means that it 
is often possible to find smoothing parameters yielding a slightly lower GCV/UBRE score by perturbing some of the estimated 
smoothing parameters a little. Technically this occurs because the performance iteration effectively neglects the dependence of the 
iterative weights on the smoothing parameters. Setting \code{perf.iter} to \code{FALSE} uses O'Sullivan's approach in which
the IRLS is run to convergence for each trial set of smoothing parameters. This is much slower, since it uses \code{nlm} with 
finite differences and requires many more IRLS steps. } 

\item{rank.tol}{The tolerance used to estimate rank when using \code{fit.method="magic"}.}
}

\details{ 
With fit method \code{"mgcv"}, 
\code{maxit} and \code{globit} control the maximum iterations of the IRLS algorithm, as follows: 
the algorithm will first execute up to
\code{globit} steps in which the GCV/UBRE algorithm performs a global search for the best overall 
smoothing parameter at every iteration. If convergence is not achieved within \code{globit} iterations, then a further 
\code{maxit} steps are taken, in which the overall smoothing parameter estimate is taken as the 
one locally minimising the GCV/UBRE score and resulting in the lowest EDF change. The difference 
between the two phases is only significant if the GCV/UBRE function develops more than one minima. 
The reason for this approach is that the GCV/UBRE score for the IRLS problem can develop "phantom"
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

\url{http://www.stats.gla.ac.uk/~simon/}


}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}


\seealso{   \code{\link{gam}}, \code{\link{gam.fit}}, \code{\link{glm.control}} }

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





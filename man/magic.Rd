\name{magic}
\alias{magic}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Stable Multiple Smoothing Parameter Estimation by GCV or UBRE, 
with optional fixed penalty}
\description{
Function to efficiently estimate smoothing parameters in generalized
ridge regression problems with multiple (quadratic) penalties, by GCV 
or UBRE. The function uses Newton's method in multi-dimensions, backed up by 
steepest descent to iteratively adjust the smoothing parameters for each penalty 
(one penalty may have a smoothing parameter fixed at unity ). 

For maximal numerical stability the method is based on orthogonal decomposition methods, 
and attempts to deal with numerical rank deficiency gracefully using a truncated singular 
value decomposition approach.
} %- end description

\usage{
magic(y,X,sp,S,off,rank=NULL,H=NULL,C=NULL,w=NULL,
      gamma=1,scale=1,gcv=TRUE,ridge.parameter=NULL,
      control=list(maxit=50,tol=1e-6,step.half=25,
              rank.tol=.Machine$double.eps^0.5))
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{y}{is the response data vector.}

\item{X}{is the model matrix.}

\item{sp}{is the array of smoothing parameters multiplying the penalty matrices stored in 
\code{S}. Any that are negative are autoinitialized, otherwise they are taken as supplying 
starting values. A supplied starting value will be reset to a default starting value if the
gradient of the GCV/UBRE score is too small at the supplied value.  }
 
\item{S}{ is a list of of penalty matrices. \code{S[[i]]} is the ith penalty matrix, but note
that it is not stored as a full matrix, but rather as the smallest square matrix including all 
the non-zero elements of the penalty matrix. Element 1,1 of \code{S[[i]]}  occupies 
element \code{off[i]}, \code{off[i]} of the ith penalty matrix. Each \code{S[[i]]} must be 
positive semi-definite.  }

\item{off}{is an array indicating the first parameter in the parameter vector that is 
penalized by the penalty involving \code{S[[i]]}.}

\item{rank}{ is an array specifying the ranks of the penalties. This is useful, but not 
essential, for forming square roots of the penalty matrices.}

\item{H}{ is the optional offset penalty - i.e. a penalty with a smoothing parameter fixed at 
1. This is useful for allowing regularization of the estimation process, fixed smoothing 
penalties etc.}

 \item{C}{ is the optional matrix specifying any linear equality constraints on the fitting 
problem. If \eqn{\bf b}{b} is the parameter vector then the parameters are forced to satisfy 
\eqn{ {\bf Cb} = {\bf 0} }{Cb=0}. }

\item{w}{ the regression weights. If this is a matrix then it is taken as being the 
square root of the inverse of the covariance matrix of \code{y}. If \code{w} is an array then 
it is taken as the diagonal of this matrix, or simply the weight for each element of 
\code{y}.}

\item{gamma}{is an inflation factor for the model degrees of freedom in the GCV or UBRE 
score.}

\item{scale}{ is the scale parameter for use with UBRE.}

\item{gcv}{ should be set to \code{TRUE} if GCV is to be used, \code{FALSE} for UBRE.}

\item{ridge.parameter}{It is sometimes useful to apply a ridge penalty to the fitting problem, 
penalizing the parameters in the constrained space directly. Setting this parameter to a value 
greater than zero will cause such a penalty to be used, with the magnitude given by the 
parameter value.}

\item{control}{ is a list of iteration control constants with the following elements:

\item{maxit}{The maximum number of iterations of the magic algorithm to allow.}

\item{tol}{The tolerance to use in judging convergence.}

\item{step.half}{If a trial step fails then the method tries halving it up to a maximum of 
\code{step.half} times.}

\item{rank.tol}{is a constant used to test for numerical rank deficiency of the problem. 
Basically any singular value less than \code{rank_tol} multiplied by the largest singular value of 
the  problem is set to zero.}

} %- end of control
}

\details{ 

The method is a computationally efficient means of applying GCV or UBRE (often approximately 
AIC) to the 
problem of smoothing parameter selection in generalized ridge regression problems 
of the form:
\deqn{ minimise~ \| { \bf W} ({ \bf Xb - y} ) \|^2 + {\bf b}^\prime {\bf Hb}  \sum_{i=1}^m
\theta_i {\bf b^\prime S}_i{\bf b} }{ minimise || W (Xb-y) ||^2 + b'Hb+
theta_1 b'S_1 b + theta_2 b'S_2 b + . . .}
possibly subject to constraints \eqn{ {\bf Cb}={\bf 0}}{Cb=0}. 
\eqn{ {\bf X}}{X} is a design matrix, \eqn{\bf b}{b} a parameter vector, 
\eqn{\bf y}{y} a data vector, \eqn{\bf W}{W} a weight matrix,
\eqn{ {\bf S}_i}{S_i} a positive semi-definite matrix  of coefficients
defining the ith penalty with associated smoothing parameter \eqn{\theta_i}{theta_i}, 
\eqn{\bf H}{H} is the positive semi-definite offset penalty matrix  and \eqn{\bf C}{C} a 
matrix of coefficients defining any linear equality constraints on the problem. 
\eqn{ {\bf X}}{X} need not be of full column rank.

The \eqn{\theta_i}{theta_i} are chosen to minize either the GCV score:

\deqn{V_g = \frac{n\|{\bf W}({\bf y} - {\bf Ay})\|^2}{[tr({\bf I} - \gamma {\bf A})]^2}}{V_g = n ||W(y-Ay)||^2/[tr(I - g A)]^2}

or the UBRE score:

\deqn{V_u=\|{\bf W}({\bf y}-{\bf Ay})\|^2/n-2 \phi tr({\bf I}-\gamma {\bf A})/n + \phi}{
V_u =||W(y-Ay||^2/n - 2 s tr(I - g A)/n + s }

where \eqn{\gamma}{g} is \code{gamma} the inflation factor for degrees of freedom (usually set to 1) and \eqn{\phi}{s} 
is \code{scale}, the scale parameter. \eqn{\bf A}{A} is the hat matrix (influence matrix) for the fitting problem (i.e
the matrix mapping data to fitted values). Dependence of the scores on the smoothing parameters is through \eqn{\bf A}{A}. 

The method operates by  Newton or steepest descent updates of the logs of the 
\eqn{\theta_i}{theta_i}. A key aspect of the method is stable and economical calculation of the 
first and second derivatives of the scores w.r.t. the log smoothing parameters. 
Because the GCV/UBRE scores are flat w.r.t. very large or very small \eqn{\theta_i}{theta_i}, 
it's important to get good starting parameters, and to be careful not to step into a flat region
of the smoothing parameter space. For this reason the algorithm rescales any Newton step that 
would result in a \eqn{log(\theta_i)}{log(theta_i)} change of more than 5. Newton steps are 
only used if the Hessian of the GCV/UBRE is postive definite, otherwise steepest descent is 
used. Similarly steepest descent is used if the Newton step has to be contracted too far 
(indicating that the quadratic model underlying Newton is poor). All initial steepest descent 
steps are scaled so that their largest component is 1. However a step is calculated, 
it is never expanded if it is successful (to avoid flat portions of the objective), 
but steps are successively halved if they do not decrease the GCV/UBRE score, until 
they do, or the direction is deemed to have failed. (Given the smoothing parameters the optimal 
\eqn{\bf b}{b} parameters are easily found.)


The method is coded in \code{C} with matrix factorizations performed using LINPACK and LAPACK routines.
}
\value{The function returns a list with the following items:

\item{b}{The best fit parameters given the estimated smoothing parameters.}

\item{scale}{the estimated (GCV) or supplied (UBRE) scale parameter.}

\item{score}{the minimized GCV or UBRE score.}

\item{sp}{an array of the estimated smoothing parameters.}

\item{rV}{a factored form of the parameter covariance matrix. The (Bayesian)  covariance
matrix of the parametes \code{b} is given by \code{rV\%*\%t(rV)*scale}. }

\item{gcv.info}{is a list of information about the performance of the method with the following elements:
\item{full.rank}{The apparent rank of the problem: number of parameters less number of equality constraints.}
\item{rank}{The estimated actual rank of the problem (at the final iteration of the method).}
\item{fully.converged}{is \code{TRUE} if the method converged by satisfying the convergence criteria, and \code{FALSE} if it coverged 
by failing to decrease the score along the search direction.}
\item{hess.pos.def}{is \code{TRUE} if the hessian of the UBRE or GCV score was positive definite at convergence.}
\item{iter}{is the number of Newton/Steepest descent iterations taken.}
\item{score.calls}{is the number of times that the GCV/UBRE score had to be evaluated.}
\item{rms.grad}{is the root mean square of the gradient of the UBRE/GCV score w.r.t. the smoothing parameters.}
}
Note that some further useful quantities can be obtained using \code{\link{magic.post.proc}}.

}
\references{
Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}


\seealso{  
\code{\link{magic.post.proc}},
\code{\link{mgcv}},
\code{\link{gam}},
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..






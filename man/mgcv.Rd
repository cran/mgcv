\name{mgcv}
\alias{mgcv}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Multiple Smoothing Parameter Estimation by GCV or UBRE}
\description{
Function to efficiently estimate smoothing parameters in Generalized
Ridge Regression Problem with multiple (quadratic) penalties, by GCV 
or UBRE. The function uses Newton's method in multi-dimensions, backed up by steepest descent to iteratively 
adjust a set of relative smoothing parameters for each penalty. To ensure that the overall level of smoothing
is optimal, and to guard against trapping by local minima, a highly efficient global minimisation with respect to 
one overall smoothing parameter is also made at each iteration.

The other functions in the \code{mgcv} package are listed under `See Also'.
}
\usage{
mgcv(M)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{M}{ is the single argument to \code{mgcv} and should be a list
    with the elements listed below:}
 \item{M$y}{The response data vector.}
 \item{M$w}{A vector of weights for the data (often proportional to the 
           reciprocal of the variance). }
 \item{M$X}{The design matrix for the problem, note that \code{ncol(M$X)}
            must give the number of model parameters, while \code{nrow(M$X)} 
            should give the number of data.}
\item{M$C}{Matrix containing any linear equality constraints 
            on the problem (i.e. \eqn{\bf C}{C} in \eqn{ {\bf Cp}={\bf 0} }{Cp=0}).}
 \item{M$S}{ A one dimensional array containing the non-zero elements of the
   penalty matrices. Let \code{start<-sum(M$df[1:(k-1)]^2)} or
   \code{start<-0} if \code{k==1}. Then penalty matrix \code{k} has
   \code{M$S[start+i+M$df[i]*(j-1)} on its ith row and jth column.
   The kth full penalty matrix is constructed by placing the matrix just
   derived into an appropriate matrix of zeroes with it's 1,1 element at
   \code{M$off[k],M$off[k]}
 }
\item{M$m}{The number of penalty terms that are to be applied.}
 \item{M$off}{ Offset values locating the elements of \code{M$S[]} in the correct location within each penalty coefficient matrix. indexing starts at zero.}
 \item{M$fix}{ \code{M$fix[i]} is FALSE if the corresponding smoothing parameter 
is to be estimated, or TRUE if it is to be fixed at zero.}
 \item{M$df}{ \code{M$df[i]} gives the number of rows and columns of
               \code{M$S[i]} that are non-zero.}
 \item{M$sig2}{ This serves 2 purposes. If it is strictly positive then UBRE is 
used with sig2 taken as the known error variance. If it is 
non-positive then GCV is used (and on return it will contain an estimate of the error variance). 
Note that it is assumed that the variance of \eqn{y_i}{y_i} is 
given
                by \eqn{\sigma^2/w_i}{\code{sig2}/w_i}.}   
\item{M$sp}{ An array of initial smoothing parameter estimates if any of 
these is not strictly positive then automatic generation of initial values
is used. If \code{M$fixed.sp} is non-zero then these parameters will used as known smoothing parameters.}
\item{M$fixed.sp}{Non-zero to indicate that all smoothing parameters are supplied as known, and are not 
      to be estimated.}
\item{M$conv.tol}{The convergence tolerence to use for the multiple
  smoothing parameter search. 1e-6 is usually ok.}
\item{M$max.half}{Each step of the GCV minimisation is repeatedly halved if it fails to
  improve the GCV score. This is the maximum number of halvings to try
  at each step. 15 is usually enough.}
\item{M$min.edf}{The minimum estimated degrees of freedom for the model: set to less than or equal to 0 to ignore.
                 Supplying this quantity helps to maintain numerical stability in cases in which the best model is 
                 actually very smooth. }
\item{M$target.edf}{Set this to negative for it to be ignored. If it is positive then this is taken as the target edf
to use for the model - if the GCV/UBRE function has multiple minima w.r.t. the overall smoothing parameter then the 
minimum closest to the target edf will be selected. This feature is primarily useful when selecting smoothing parameters as part 
of iterative least squares methods - the GCV/UBRE score can develop spurious minima in these cases that are 
not representative of the problem being approximated.
}
}
\details{ 

This is documentation for the code implementing the method described in section 
4 of 
Wood (2000) . The method is a computationally efficient means of applying GCV to 
the 
problem of smoothing parameter selection in generalized ridge regression problems 
of 
the form:
\deqn{ minimise~ \| { \bf W}^{1/2} ({ \bf Xp - y} ) \|^2 \rho +  \sum_{i=1}^m
\lambda_i {\bf p^\prime S}_i{\bf p} }{ minimise || W^0.5 (Xp-y) ||^2 rho + 
lambda_1 p'S_1 p + lambda_1 p'S_2 p + . . .}
possibly subject to constraints \eqn{ {\bf Cp}={\bf 0}}{Cp=0}. 
\eqn{ {\bf X}}{X} is a design matrix, \eqn{\bf p}{p} a parameter vector, 
\eqn{\bf y}{y} a data vector, \eqn{\bf W}{W} a diagonal weight matrix,
\eqn{ {\bf S}_i}{S_i} a positive semi-definite matrix  of coefficients
defining the ith penalty and \eqn{\bf C}{C} a matrix of coefficients 
defining any linear equality constraints on the problem. The smoothing
parameters are the \eqn{\lambda_i}{lambda_i} but there is an overall
smoothing parameter \eqn{\rho}{rho} as well. Note that \eqn{ {\bf X}}{X}
must be of full column rank, at least when projected  into the null space
of any equality constraints.  

The method operates by alternating very efficient direct searches for 
\eqn{\rho}{rho}
with Newton or steepest descent updates of the logs of the \eqn{\lambda_i}{lambda_i}. 
Because the GCV/UBRE scores are flat w.r.t. very large or very small \eqn{\lambda_i}{lambda_i}, 
it's important to get good starting parameters, and to be careful not to step into a flat region
of the smoothing parameter space. For this reason the algorithm rescales any Newton step that 
would result in a \eqn{log(\lambda_i)}{log(lambda_i)} change of more than 5. Newton steps are only used
if the Hessian of the GCV/UBRE is postive definite, otherwise steepest descent is used. Similarly steepest 
descent is used if the Newton step has to be contracted too far (indicating that the quadratic model 
underlying Newton is poor). All initial steepest descent steps are scaled so that their largest component is
1. However a step is calculated, it is never expanded if it is successful (to avoid flat portions of the objective), 
but steps are successively halved if they do not decrease the GCV/UBRE score, until they do, or the direction is deemed to have 
failed. \code{M$conv} provides some convergence diagnostics.

The method is coded in \code{C} and is intended to be portable. It should be 
noted that seriously ill conditioned problems (i.e. with close to column rank 
deficiency in the design matrix) may cause problems, especially if weights vary 
wildly between observations.  
}
\value{
  
  The function returns a list of the same form as the input list with the 
following
  elements added/modified:
  \item{M$p}{The best fit parameters given the estimated smoothing parameters.}
  \item{M$sp}{The estimated smoothing parameters 
(\eqn{\lambda_i/\rho}{lambda_i/rho})}
  \item{M$sig2}{The estimated (GCV) or supplied (UBRE) error variance.}
  \item{M$edf}{The estimated degrees of freedom associated with each penalized term. If \eqn{ {\bf p}_i}{p_i} is the vector of parameters penalized by the ith penalty then we can write \eqn{ {\bf p}_i = {\bf P}_i {\bf y}}{p_i = P_i y}, and the edf for the ith term is \eqn{ tr({\bf XP}_i)}{tr(XP_i)}... while this is sensible in the non-overlapping penalty case (e.g. GAMs) it makes alot less sense if penalties overlap!}
  \item{M$Vp}{Estimated covariance matrix of model parameters.}
\item{M$gcv.ubre}{The minimized GCV/UBRE score.}
  \item{M$conv}{A list of convergence diagnostics, with the following elements:
\item{edf}{Array of whole model estimated degrees of freedom.}
\item{hat}{Array of elements from leading diagonal of `hat' matrix (`influence' matrix). Same length as \code{M$y}.}
\item{score}{Array of ubre/gcv scores at the edfs for the final set of relatinve smoothing parameters.}
\item{g}{the gradient of the GCV/UBRE score w.r.t. the smoothing parameters at termination.}
\item{h}{the second derivatives corresponding to \code{g} above - i.e. the leading diagonal of the Hessian.}
\item{e}{the eigen-values of the Hessian. These should all be non-negative!}
\item{iter}{the number of iterations taken.}
\item{in.ok}{\code{TRUE} if the second smoothing parameter guess improved the GCV/UBRE score. (Please report examples 
where this is \code{FALSE})}
\item{step.fail}{\code{TRUE} if the algorithm terminated by failing to improve the GCV/UBRE score rather than by "converging". 
Not necessarily a problem, but check the above derivative information quite carefully.}

}
}
\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428


\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\section{WARNING }{ The method may not behave well with near column rank deficient \eqn{ {\bf
X}}{X}
especially in contexts where the weights vary wildly. } 

\seealso{  
\code{\link{gam}},
\code{\link{s}},
\code{\link{plot.gam}},
\code{\link{print.gam}},
\code{\link{summary.gam}},
\code{\link{predict.gam}},
\code{\link{persp.gam}},
\code{\link{residuals.gam}},
\code{\link{gam.models}},
\code{\link{gam.selection}},
\code{\link{gam.control}}
\code{\link{gam.check}}
\code{\link{pcls}},
\code{\link{mono.con}}
\code{\link{GAMsetup}},
\code{\link{gam.parser}},
\code{\link{gam.fit}},
\code{\link{gam.nbut}},
\code{\link{gam.setup}},
\code{\link{QT}},
\code{\link{get.family}},
\code{\link{neg.binom}},
\code{\link{null.space.dimension}},
\code{\link{theta.maxl}},
\code{\link{uniquecombs}}
}

\examples{
set.seed(0)    
n<-100 # number of observations to simulate
x <- runif(4 * n, 0, 1) # simulate covariates
x <- array(x, dim = c(4, n)) # put into array for passing to GAMsetup
pi <- asin(1) * 2  # begin simulating some data
y <- 2 * sin(pi * x[1, ])
y <- y + exp(2 * x[2, ]) - 3.75887
y <- y + 0.2 * x[3, ]^11 * (10 * (1 - x[3, ]))^6 + 10 * (10 *
        x[3, ])^3 * (1 - x[3, ])^10 - 1.396
sig2<- -1    # set magnitude of variance
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- y + e          # simulated data
w <- matrix(1, n, 1) # weight matrix
par(mfrow = c(2, 2)) # scatter plots of simulated data   
plot(x[1, ], y);plot(x[2, ], y)
plot(x[3, ], y);plot(x[4, ], y)
# create list for passing to GAMsetup
G <- list(m = 4, n = n, nsdf = 0, df = c(15, 15, 15, 15),dim=c(1,1,1,1),
          s.type=c(0,0,0,0),by=0,by.exists=c(FALSE,FALSE,FALSE,FALSE),
          p.order=c(0,0,0,0), x = x,n.knots=rep(0,4)) 
H <- GAMsetup(G)
H$y <- y    # add data to H
H$sig2 <- sig2  # add variance (signalling GCV use in this case) to H
H$w <- w       # add weights to H
H$sp<-array(-1,H$m)
H$fix<-array(FALSE,H$m)
H$conv.tol<-1e-6;H$max.half<-15
H$min.edf<-5;H$fixed.sp<-0
H <- mgcv(H)  # select smoothing parameters and fit model 
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..






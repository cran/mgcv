\name{mgcv}
\alias{mgcv}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Multiple Smoothing Parameter Estimation by GCV or UBRE}
\description{
Function to efficiently estimate smoothing parameters in Generalized
Ridge Regression Problem with multiple (quadratic) penalties, by GCV 
or UBRE. The function uses Newton's method in multi-dimensions.
}
\usage{
mgcv(M)
}
%- maybe also `usage' for other objects documented here.
\arguments{ The function takes list argument \code{M} with the 
  following elements:
 \item{M$y}{The response data vector}
 \item{M$X}{The design matrix for the problem, note that \code{ncol(M$X)}
            must give the number of model parameters, while \code{nrow(M$X)} 
            should give the number of data.}
\item{M$Z}{Matrix containing the null space of the linear equality constraints 
            on the problem, stored as Householder rotations as produced by 
routine
            \code{QT()}. Number of rows of \code{M$Z} must be number of 
constraints.}
 \item{M$S}{ A 3 dimensional array. \code{M$S[i]} is a 2 dimensional array 
             containing the elements of the ith penalty matrix. Note that to
             economize on storage space, rows and columns containing only zeroes 
             need not be stored. It is assumed that \code{M$S[i][j][k]} contains
             element \code{j+M$off[i],k+M$off[i]} of the ith penalty matrix, 
             \eqn{ {\bf S}_i}{S_i}. It is also assumed that only the first 
            \code{M$df[i]} rows and columns of \code{M$S[i]} contain non-zero
            elements. \code{dim(M$S)[1]} must evaluate to the number of 
penalties.}
 \item{M$off}{ Offset values locating the elements of \code{M$S[]} in the correct
               location within each penalty coefficient matrix.}
 \item{M$df}{ \code{M$df[i]} gives the number of rows and columns of
               \code{M$S[i]} that are non-zero.}
 \item{M$sig2}{ This serves 2 purposes. If it is strictly positive then UBRE is 
used
                with sig2 taken as the known error variance. If it is 
non-positive 
                then GCV is used (and an estimate of the error variance 
returned).
                Note that it is assumed that the variance of \eqn{y_i}{y_i} is 
given
                by \eqn{\sigma^2/w_i}{\code{sig2}/w_i}.}            
}
\details{ 

This is documentation for the code implementing the method described in section 
4. of 
Wood (2000) . The method is a computationally efficient means of applying GCV to 
the 
problem of smoothing parameter selection in generalized ridge regression problems 
of 
the form:
\deqn{ minimise~ \| { \bf W}^{1/2} ({ \bf Xp - y} ) \|^2 \rho +  \sum_{i=1}^m
\theta_i {\bf p^\prime S}_i{\bf p} }{ minimise || W^0.5 (Xp-y) ||^2 rho + 
lambda_1 p'S_1 p + lambda_1 p'S_2 p + . . .}
possibly subject to constraints \eqn{ {\bf Cp}={\bf 0}}{Cp=0}. 
\eqn{ {\bf X}}{X} is a design matrix, \eqn{\bf p}{p} a parameter vector, \eqn{\bf 
y}{y}
 a data vector, \eqn{\bf W}{W} a diagonal weight matrix, \eqn{ {\bf S}_i}{S_i} a 
matrix 
 of coefficients defining the ith penalty and \eqn{\bf C}{C} a matrix of 
coefficients
defining any linear equality constraints on the problem. The smoothing parameters 
are
the \eqn{\lambda_i}{lambda_i} but there is an overall smoothing parameter 
\eqn{\rho}{rho}
as well. 

The method operates by alternating very efficient direct searches for 
\eqn{\rho}{rho}
with Newton updates of the logs of the \eqn{\lambda_i}{lambda_i}

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

}
\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2)


\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}}

\section{WARNING }{ The method may not behave well with near column rank defficient \eqn{ {\bf
X}}{X}
especially in contexts where the weights vary wildly. } 

\seealso{  \code{\link{GAMsetup}} }

\examples{

 # This example modified from routine SANtest()
    
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
    plot(x[1, ], y)
    plot(x[2, ], y)
    plot(x[3, ], y)
    plot(x[4, ], y)
    G <- list(m = 4, n = n, nsdf = 0, df = c(15, 15, 15, 15),
        x = x) # creat list for passing to GAMsetup
    H <- GAMsetup(G)
    H$y <- y    # add data to H
    H$sig2 <- sig2  # add variance (signalling GCV use in this case) to H
    H$w <- w       # add weights to H
    H <- mgcv(H)  # select smoothing parameters and fit model
 
}
\keyword{GCV, AIC, UBRE, GAM }%-- one or more ...

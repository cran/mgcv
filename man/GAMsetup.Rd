\name{GAMsetup}
\alias{GAMsetup}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Set up GAM using penalized regression splines} 
\description{ Sets up design matrix \eqn{\bf X}{X}, penalty matrices \eqn{ {\bf
S}_i }{S_i} and linear equality constraint matrix \eqn{\bf C}{C} for a GAM defined in terms of 
penalized regression splines. Various other information characterising the bases used is also returned.
 The output is such that the model can be fitted and 
smoothing parameters estimated by the method of Wood (2000) as implemented in routine
\code{mgcv()}. This is usually called by \code{gam}.

}
\usage{
GAMsetup(G)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{G}{ is the single argument to this function: it is a list
  containing at least the elements listed below:}
\item{m}{the number of smooth terms in the model}
\item{df}{an array of \code{G$m} integers specifying the maximum d.f. for each spline 
  term.}
\item{n}{the number of data to be modelled}
\item{nsdf}{the number of user supplied columns of the design matrix for any parametric 
model parts}
\item{dim}{An array of dimensions for the smooths. \code{dim[i]} is the
  number of covariates that smooth \code{i} is a function of.}
\item{fix}{An array of logicals indicating whether each smooth term has fixed degrees of freedom or not.}
\item{s.type}{An array giving the type of basis used for each term. 0 for
  cubic regression spline, 1 for t.p.r.s }
\item{p.order}{An array giving the order of the penalty for each term. 0
for auto selection.}
\item{x}{an array of \code{G$n} element arrays of data and (optionally) design matrix 
columns. The first \code{G$nsdf} elements of \code{G$x} should contain the elements of 
the columns of the design matrix corresponding to the parametric part of the model. The 
remaining \code{G$m} elements of \code{G$x} are the values of the covariates that are 
arguments of the spline terms. Note that the smooths will be centred and no intercept term 
will be added unless an array of 1's is supplied as part of in
\code{G$x}}
\item{vnames}{Array of variable names, including the constant, if present.}
\item{w}{prior weights on response data.}
\item{by}{a 2-d array of \code{by} variables (i.e. covariates that multiply a
smooth term) \code{by[i,j]} is the jth value for the ith \code{by}
variable. There are only as many rows of this array as there are
\code{by} variables in the model (often 0). The rownames of \code{by}
give the \code{by} variable names.}
\item{by.exists}{an array of logicals: \code{by.exists[i]} is
\code{TRUE} if the ith smooth has a \code{by} variable associated with
it, \code{FALSE} otherwise.}
\item{knots}{ a compact array of user supplied knot locations for each smooth, in the order corresponding 
to the  row order in \code{G$x}. There are \code{G$dim[i]} arrays of length \code{G$n.knots[i]} for the ith
         smooth - all these arrays are packed end to end in 1-d array \code{G$knots} - 
        zero length 1 for no knots.}
\item{n.knots}{ array giving number of user supplied knots of basis for each smooth term 0's for none supplied.}
}
\value{
  A list \code{H}, containing the elements of \code{G} (the input list) plus the 
following:  
\item{X}{the full design matrix.}
\item{S}{ A one dimensional array containing the non-zero elements of the
   penalty matrices. Let \code{start[k+1]<-start[k]+H$df[1:(k-1)]^2} and
   \code{start[1]<-0}. Then penalty matrix \code{k} has
   \code{H$S[start[k]+i+H$df[i]*(j-1)} on its ith row and jth column.
   To get the kth full penalty matrix the matrix so obtained would be
   inserted into a full matrix of zeroes with it's 1,1 element at \code{H$off[k],H$off[k]}. 
 }
\item{off}{is an array of offsets, used to facilitate efficient storage of the penalty 
matrices and to indicate where in the overall parameter vector the parameters of the ith 
spline reside (e.g. first parameter of ith spline is at \code{p[off[i]+1]}).}
\item{C}{ a matrix defining the linear equality constraints on the parameters used to define the
  the model (i.e. \eqn{\bf C}{C} in \eqn{ {\bf Cp } ={\bf 0}}{Cp=0}). }
\item{UZ}{Array containing matrices, which transform from a t.p.r.s. basis to the
  equivalent t.p.s. basis (for t.p.r.s. terms only). The packing method
  is as follows: \cr
set \code{start[1]<-0} and
  \code{start[k+1]<-start[k]+(M[k]+n)*tp.bs[k]} where \code{n} is number
  of data, \code{M[k]} is penalty null space dimension and
  \code{tp.bs[k]} is zero for a cubic regression spline and the basis
  dimension for a t.p.r.s. Then element \code{i,j} of the UZ matrix for
  model term \code{k} is:\cr
 \code{UZ[start[k]+i+(j=1)*(M[k]+n)]}.}
\item{Xu}{Set of unique covariate combinations for each term.  The packing method
  is as follows: \cr
set \code{start[1]<-0} and
  \code{start[k+1]<-start[k]+(xu.length[k])*tp.dim[k]} where \code{xu.length[k]} is number
  of unique covariate combinations and \code{tp.dim[k]} is zero for a
  cubic regression spline
  and the dimension of the smooth (i.e. number of covariates it is a
  function of) for a t.p.r.s. Then element \code{i,j} of the Xu matrix for
  model term \code{k} is:\cr
\code{Xu[start[k]+i+(j=1)*(xu.length[k])]}.}
\item{xu.length}{Number of unique covariate combinations for each
  t.p.r.s. term.}
\item{covariate.shift}{All covariates are centred around zero before
  bases are constructed - this is an array of the applied shifts.}
\item{xp}{matrix whose rows contain the covariate values corresponding to the 
parameters  of each cubic regression spline - the cubic regression splines are parameterized using their \eqn{y}{y}- 
values at a series of \eqn{x}{x} values - these vectors contain those \eqn{x}{x} 
values! Note that these will be covariate shifted.}
}
\references{ Wood, S.N. (2000) "Modelling and smoothing parameter estimation with multiple
quadratic penalties" JRSSB 62(2):413-428 }
\author{Simon N. Wood \email{snw@st-and.ac.uk}}

\seealso{ \code{\link{mgcv}} \code{\link{gam}} }

\examples{
# This example modified from routine SANtest()
set.seed(0)
n<-100 # number of observations to simulate
x <- runif(5 * n, 0, 1) # simulate covariates
x <- array(x, dim = c(5, n)) # put into array for passing to GAMsetup
pi <- asin(1) * 2  # begin simulating some data
y <- 2 * sin(pi * x[2, ])
y <- y + exp(2 * x[3, ]) - 3.75887
y <- y + 0.2 * x[4, ]^11 * (10 * (1 - x[4, ]))^6 + 10 * (10 * 
     x[4, ])^3 * (1 - x[4, ])^10 - 1.396
sig2<- -1    # set magnitude of variance 
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- y + e          # simulated data
w <- matrix(1, n, 1) # weight matrix
par(mfrow = c(2, 2)) # scatter plots of simulated data
plot(x[2, ], y);plot(x[3, ], y);plot(x[4, ], y);plot(x[5, ], y)
x[1,]<-1
# create list for passing to GAMsetup....
G <- list(m = 4, n = n, nsdf = 0, df = c(15, 15, 15, 15),dim=c(1,1,1,1),
     s.type=c(0,0,0,0),by=0,by.exists=c(FALSE,FALSE,FALSE,FALSE),
     p.order=c(0,0,0,0),x = x,n.knots=rep(0,4))
H <- GAMsetup(G)
H$y <- y    # add data to H
H$sig2 <- sig2  # add variance (signalling GCV use in this case) to H
H$w <- w       # add weights to H
H$sp<-array(-1,H$m)
H$fix<-array(FALSE,H$m)
H$conv.tol<-1e-6;H$max.half<-15
H$min.edf<-5
H <- mgcv(H)  # select smoothing parameters and fit model    
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...


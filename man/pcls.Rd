\name{pcls}
\alias{pcls}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Penalized Constrained Least Squares Fitting}
\description{
Solves least squares problems with quadratic penalties subject to linear
equality and inequality constraints using quadratic programming.
}
\usage{
pcls(M)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{M}{is the single list argument to \code{pcls}. It should have  the 
  following elements (last 3 are not in argument for \link{mgcv}) :}
 \item{M$y}{The response data vector.}
 \item{M$w}{A vector of weights for the data (often proportional to the 
           reciprocal of the variance). }
 \item{M$X}{The design matrix for the problem, note that \code{ncol(M$X)}
            must give the number of model parameters, while \code{nrow(M$X)} 
            should give the number of data.}
\item{M$C}{Matrix containing any linear equality constraints 
            on the problem (i.e. \eqn{ \bf C}{C} in \eqn{ {\bf Cp}={\bf
		0} }{Cp=0}). If you have no equality constraints
	    initialize this to a zero by zero matrix.}
 \item{M$S}{ A one dimensional array containing the non-zero elements of the
   penalty matrices. Let \code{start<-sum(M$df[1:(k-1)]^2)} or
   \code{start<-0} if \code{k==1}. Then penalty matrix \code{k} has
   \code{M$S[start+i+M$df[i]*(j-1)} on its ith row and jth column.
   The kth full penalty matrix is constructed by placing the matrix just
   derived into an appropriate matrix of zeroes with it's 1,1 element at
   \code{M$off[k],M$off[k]}
 }
 \item{M$off}{ Offset values locating the elements of \code{M$S[]} in
   the correct location within each penalty coefficient matrix. (Zero
   offset implies starting in first location)}

 \item{M$df}{ \code{M$df[i]} gives the number of rows and columns of
               \code{M$S[i]} that are non-zero.}
  
\item{M$sp}{ An array of  smoothing parameter estimates.}
\item{M$p}{An array of feasible initial parameter estimates - these must
satisfy the constraints, but should avoid satisfying the inequality
constraints as equality constraints.}
\item{M$Ain}{Matrix for the inequality constraints \eqn{ {\bf A}_{in}
    {\bf p} > {\bf b}_{in}}{A_in p > b}. }
\item{M$bin}{vector in the inequality constraints. }
}
\details{ 

  This solves the problem:
 
\deqn{ minimise~ \| { \bf W}^{1/2} ({ \bf Xp - y} ) \|^2  +  \sum_{i=1}^m
\lambda_i {\bf p^\prime S}_i{\bf p} }{ minimise || W^0.5 (Xp-y) ||^2 + 
lambda_1 p'S_1 p + lambda_1 p'S_2 p + . . .}
subject to constraints \eqn{ {\bf Cp}={\bf 0}}{Cp=0} and \eqn{ {\bf
    A}_{in}{\bf p}>{\bf b}_{in}}{A_in p > b_in}, w.r.t. \eqn{\bf p}{p} given the
smoothing parameters \eqn{\lambda_i}{lambda_i}.
\eqn{ {\bf X}}{X} is a design matrix, \eqn{\bf p}{p} a parameter vector, 
\eqn{\bf y}{y} a data vector, \eqn{\bf W}{W} a diagonal weight matrix,
\eqn{ {\bf S}_i}{S_i} a positive semi-definite matrix  of coefficients
defining the ith penalty and \eqn{\bf C}{C} a matrix of coefficients 
defining the linear equality constraints on the problem. The smoothing
parameters are the \eqn{\lambda_i}{lambda_i}. Note that \eqn{ {\bf X}}{X}
must be of full column rank, at least when projected  into the null space
of any equality constraints. \eqn{ {\bf A}_{in}}{A_in} is a matrix of
coefficients defining the inequality constraints, while \eqn{ {\bf
    b}_{in}}{b_in} is a vector involved in defining the inequality constraints.  

Quadratic programming is used to perform the solution. The method used
is designed for maximum stability with least squares problems:
i.e. \eqn{ {\bf X}^\prime {\bf X}}{X'X} is not formed explicitly. See
Gill et al. 1981.

}
\value{ The function returns an array containing the estimated parameter
  vector. 
   
}
\references{

Gill, P.E., Murray, W. and Wright, M.H. (1981) Practical Optimization. Academic
Press, London. 

Wood, S.N. (1994) Monotonic smoothing splines fitted by cross validation SIAM
Journal on Scientific Computing 15(5):1126-1133

\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}} 

\seealso{  \code{\link{mgcv}} \code{\link{mono.con}}  }

\examples{
# first an un-penalized example - fit E(y)=a+bx subject to a>0
n<-100
x<-runif(n);y<-x-0.2+rnorm(n)*0.1
M<-list(X=matrix(0,n,2),p=c(0.1,0.5),off=array(0,0),df=array(0,0),S=0,
Ain=matrix(0,1,2),bin=0,C=matrix(0,0,0),sp<-0,y=y,w=y*0+1)
M$X[,1]<-1;M$X[,2]<-x;M$Ain[1,]<-c(1,0)
pcls(M)->M$p
plot(x,y);abline(M$p,col=2);abline(coef(lm(y~x)),col=3)

# and now a penalized example: a monotonic penalized regression spline .....

# Generate data from a monotonic truth.
set.seed(10);x<-runif(100)*4-1;x<-sort(x);
f<-exp(4*x)/(1+exp(4*x));y<-f+rnorm(100)*0.1;plot(x,y)
# Show regular spline fit (and save fitted object)
f.ug<-gam(y~s(x,k=10,bs="cr"));lines(x,fitted(f.ug))
# Create Design matrix, constriants etc. for monotonic spline....
gam.setup(y~s(x,k=10,bs="cr")-1)->G;GAMsetup(G)->G;F<-mono.con(G$xp);
G$Ain<-F$A;G$bin<-F$b;G$C<-matrix(0,0,0);G$sp<-f.ug$sp;
G$p<-G$xp;G$y<-y;G$w<-y*0+1

p<-pcls(G);  # fit spline (using s.p. from unconstrained fit)

# now modify the gam object from unconstrained fit a little, to use it
# for predicting and plotting constrained fit. 
p<-c(0,p);f.ug$coefficients<-p; 
lines(x,predict.gam(f.ug,newdata=data.frame(x=x)),col=2)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





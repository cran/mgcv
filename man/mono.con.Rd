\name{mono.con}
\alias{mono.con}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Monotonicity constraints for a cubic regression spline.}
\description{ Finds linear constraints sufficient for monotonicity (and
  optionally upper and/or lower boundedness) of a cubic regression
  spline. The basis representation assumed is that given by the
  \code{gam}, \code{"cr"} basis: that is the spline has a set of knots,
  which have fixed x values, but the y values of which constitute the
  parameters of the spline. }
\usage{
mono.con(x,up=TRUE,lower=NA,upper=NA)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{The array of knot locations.}
 \item{up}{If \code{TRUE} then the constraints imply increase, if
   \code{FALSE} then decrease. }
 \item{lower}{This specifies the lower bound on the spline unless it is
   \code{NA} in which case no lower bound is imposed.}
 \item{upper}{This specifies the upper bound on the spline unless it is
   \code{NA} in which case no upper bound is imposed.}
}
\details{ Consider the natural cubic spline passing through the points:
  \eqn{ \{x_i,p_i:i=1 \ldots n \} }{ (x_i,p_i), i=1..n}. Then it is possible
  to find a relatively small set of linear constraints on \eqn{ \bf p}{p}
  sufficient to ensure monotonicity (and bounds if required):
  \eqn{ {\bf Ap}\ge{\bf b} }{Ap>=b}. Details are given in Wood (1994).
  This function returns a list containing \code{A} and \code{b}. 
}
\value{ The function returns a list containing constraint matrix
  \code{A} and constraint vector \code{b}.   
}
\references{

Gill, P.E., Murray, W. and Wright, M.H. (1981) Practical Optimization. Academic
Press, London. 

Wood, S.N. (1994) Monotonic smoothing splines fitted by cross validation SIAM
Journal on Scientific Computing 15(5):1126-1133

\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk} } 

\seealso{  \code{\link{mgcv} } \code{\link{pcls}}  }

\examples{
# Fit a monotonic penalized regression spline .....

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
\keyword{models} \keyword{smooth} \keyword{regression} %-- one or more ..







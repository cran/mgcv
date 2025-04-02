\name{XWXd}
\alias{XWXd}
\alias{XWyd}
\alias{Xbd}
\alias{diagXVXd}
\alias{ijXVXd}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Internal functions for discretized model matrix handling}
\description{Routines for computing with discretized model matrices as described in Wood et al. (2017) and Li and Wood (2019).
}
\usage{
XWXd(X,w,k,ks,ts,dt,v,qc,nthreads=1,drop=NULL,ar.stop=-1,ar.row=-1,ar.w=-1,
     lt=NULL,rt=NULL)
XWyd(X,w,y,k,ks,ts,dt,v,qc,drop=NULL,ar.stop=-1,ar.row=-1,ar.w=-1,lt=NULL)
Xbd(X,beta,k,ks,ts,dt,v,qc,drop=NULL,lt=NULL)
diagXVXd(X,V,k,ks,ts,dt,v,qc,drop=NULL,nthreads=1,lt=NULL,rt=NULL)
ijXVXd(i,j,X,V,k,ks,ts,dt,v,qc,drop=NULL,nthreads=1,lt=NULL,rt=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{X}{A list of the matrices containing the unique rows of model matrices for terms of a full model matrix, or the model matrices of the terms margins.
 if term subsetting arguments \code{lt} and \code{rt} are non-NULL then this requires an \code{"lpip"} attribute: see details. The elements of \code{X} may
 be sparse matrices of class \code{"dgCMatrix"}, in which case the list requires attributes \code{"r"} and \code{"off"} defining reverse indices (see details).}
 \item{w}{An n-vector of weights}
 \item{y}{n-vector of data.}
 \item{beta}{coefficient vector.}
 \item{k}{A matrix whose columns are index n-vectors each selecting the rows of an X[[i]] required to create the full matrix.}
\item{ks}{The ith term has index vectors \code{ks[i,1]:(ks[i,2]-1)}. The corresponing full model matrices are summed over.} 
\item{ts}{The element of \code{X} at which each model term starts.}
\item{dt}{How many elements of \code{X} contribute to each term.}
\item{v}{\code{v[[i]]} is Householder vector for ith term, if \code{qc[i]>0}.}
\item{qc}{if \code{qc[i]>0} then term has a constraint.}
\item{nthreads}{number of threads to use}
\item{drop}{list of columns of model matrix/parameters to drop}
\item{ar.stop}{Negative to ignore. Otherwise sum rows \code{(ar.stop[i-1]+1):ar.stop[i]} of the rows selected by \code{ar.row} and weighted by \code{ar.w} to get ith row of model matrix to use.}
\item{ar.row}{extract these rows...}
\item{ar.w}{weight by these weights, and sum up according to \code{ar.stop}. Used to implement AR models.}
\item{lt}{use only columns of X corresponding to these model matrix terms (for left hand \code{X} in \code{XWXd}). If \code{NULL} set to \code{rt}.}
\item{rt}{as \code{lt} for right hand \code{X}. If \code{NULL} set to \code{lt}. If \code{lt} and \code{rt} are \code{NULL} use all columns.}
\item{V}{Coefficient covariance matrix.}
\item{i}{vector of rows of \eqn{XVX^T}{XVX'} at which to evaluate.}
\item{j}{vector of corresponding columns.}
}

\details{These functions are really intended to be internal, but are exported so that they can be used in the initialization code of families without problem. They are primarily used by \code{\link{bam}} to implement the methods given in the references. \code{XWXd} produces \eqn{X^TWX}{X'WX}, \code{XWy} produces \eqn{X^TWy}{X'Wy}, \code{Xbd} produces \eqn{X\beta}{Xb} and \code{diagXVXd} produces the diagonal of \eqn{XVX^T}{XVX'}, while \code{ijXVXd} evaluates the scattered \eqn{i,j}{i,j} elements indexed in \code{i} and \code{j}. 

The \code{"lpip"} attribute of \code{X} is a list of the coefficient indices for each term. Required if subsetting via \code{lt} and \code{rt}.

\code{X} can be a list of sparse matrices of class \code{"dgCMatrix"}, in which case reverse indices are needed, mapping stored matrix rows to rows in the full matrix (that is the reverse of \code{k} which maps full matrix rows to the stored unique matrix rows). \code{r} is the same dimension as \code{k} while \code{off} is a list with as many elements as \code{k} has columns. \code{r} and  \code{off} are supplied as attributes to \code{X} . For simplicity let \code{r} and \code{off} denote a single column and element corresponding to each other: then \code{r[off[j]:(off[j+1]-1)]} contains the rows of the full matrix corresponding to row \code{j} of the stored matrix. The reverse indices are essential for efficient computation with sparse matrices. See the example code for how to create them efficiently from the forward index matrix, \code{k}.    
}
 

\references{
Wood, S.N., Li, Z., Shaddick, G. & Augustin N.H. (2017) Generalized additive models for gigadata: modelling the UK black smoke network daily data. Journal of the American Statistical Association. 112(519):1199-1210
\doi{10.1080/01621459.2016.1195744}

Li, Z & S.N. Wood (2019) Faster model matrix crossproducts for large generalized linear models with discretized covariates. Statistics and Computing.
\doi{10.1007/s11222-019-09864-2}

}

\examples{
  library(mgcv);library(Matrix)
  ## simulate some data creating a marginal matrix sequence...
  set.seed(0);n <- 4000
  dat <- gamSim(1,n=n,dist="normal",scale=2)
  dat$x4 <- runif(n)
  dat$y <- dat$y + 3*exp(dat$x4*15-5)/(1+exp(dat$x4*15-5))
  dat$fac <- factor(sample(1:20,n,replace=TRUE))
  G <- gam(y ~ te(x0,x2,k=5,bs="bs",m=1)+s(x1)+s(x4)+s(x3,fac,bs="fs"),
           fit=FALSE,data=dat,discrete=TRUE)
  p <- ncol(G$X)
  ## create a sparse version...
  Xs <- list(); r <- G$kd*0; off <- list()
  for (i in 1:length(G$Xd)) Xs[[i]] <- as(G$Xd[[i]],"dgCMatrix")
  for (j in 1:nrow(G$ks)) { ## create the reverse indices...
    nr <- nrow(Xs[[j]]) ## make sure we always tab to final stored row 
    for (i in G$ks[j,1]:(G$ks[j,2]-1)) {
      r[,i] <- (1:length(G$kd[,i]))[order(G$kd[,i])]
      off[[i]] <- cumsum(c(1,tabulate(G$kd[,i],nbins=nr)))-1
    }
  }
  attr(Xs,"off") <- off;attr(Xs,"r") <- r 

  par(mfrow=c(2,3))

  beta <- runif(p)
  Xb0 <- Xbd(G$Xd,beta,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  Xb1 <- Xbd(Xs,beta,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  range(Xb0-Xb1);plot(Xb0,Xb1,pch=".")

  bb <- cbind(beta,beta+runif(p)*.3)
  Xb0 <- Xbd(G$Xd,bb,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  Xb1 <- Xbd(Xs,bb,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  range(Xb0-Xb1);plot(Xb0,Xb1,pch=".")

  p <- length(beta) ## extract full model matrix...
  X <- matrix(Xbd(G$Xd,diag(p),G$kd,G$ks,G$ts,G$dt,G$v,G$qc),ncol=p)

  w <- runif(n)
  XWy0 <- XWyd(G$Xd,w,y=dat$y,G$kd,G$ks,G$ts,G$dt,G$v,G$qc) 
  XWy1 <- XWyd(Xs,w,y=dat$y,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  range(XWy1-XWy0);plot(XWy1,XWy0,pch=".")

  yy <- cbind(dat$y,dat$y+runif(n)-.5)
  XWy0 <- XWyd(G$Xd,w,y=yy,G$kd,G$ks,G$ts,G$dt,G$v,G$qc) 
  XWy1 <- XWyd(Xs,w,y=yy,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  range(XWy1-XWy0);plot(XWy1,XWy0,pch=".")

  A <- XWXd(G$Xd,w,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  B <- XWXd(Xs,w,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  D <- crossprod(X,w*X) ## direct computation
  range(A-D)
  range(A-B);plot(A,B,pch=".")
  ## compute some cross product terms only...
  A <- XWXd(G$Xd,w,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,lt=1:3,rt=4:5)
  range(A-D[1:nrow(A),(nrow(A)+1):ncol(D)])

  V <- crossprod(matrix(runif(p*p),p,p))
  ii <- c(20:30,100:200)
  jj <- c(50:90,150:160)
  V[ii,jj] <- 0;V[jj,ii] <- 0
  d1 <- diagXVXd(G$Xd,V,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  Vs <- as(V,"dgCMatrix")
  d2 <- diagXVXd(Xs,Vs,G$kd,G$ks,G$ts,G$dt,G$v,G$qc)
  range(d1-d2);plot(d1,d2,pch=".")
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..










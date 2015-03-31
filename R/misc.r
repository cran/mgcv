## (c) Simon N. Wood 2011-2014
## Many of the following are simple wrappers for C functions, used largely 
## for testing purposes

rmvn <- function(n,mu,V) {
## generate multivariate normal deviates. e.g.
## V <- matrix(c(2,1,1,2),2,2); mu <- c(1,1);n <- 1000;z <- rmvn(n,mu,V);crossprod(sweep(z,2,colMeans(z)))/n
  p <- ncol(V)
  R <- mroot(V,rank=ncol(V)) ## RR' = V
  if (is.matrix(mu)) {
    if (ncol(mu)!=p||nrow(mu)!=n) stop("mu dimensions wrong")
    z <- matrix(rnorm(p*n),n,p)%*%t(R) + mu
  } else { 
    if (length(mu)!=p) stop("mu dimensions wrong")
    z <- t(R%*% matrix(rnorm(p*n),p,n) + mu)
    if (n==1) z <- as.numeric(z)
  }
  z
} ## rmvn

mgcv.omp <- function() {
## does open MP appear to be available?
  oo <- .C(C_mgcv_omp,a=as.integer(-1))
  if (oo$a==1) TRUE else FALSE
}

mvn.ll <- function(y,X,beta,dbeta=NULL) {
## to facilitate testing of MVN routine mvn_ll.
## X is a sequence of m model matrices bound columnwise, with m dim attribute lpi
##   indicating where the next starts in all cases.
## beta is parameter vector - last m*(m+1)/2 elements are chol factor of precision params.
## y is m by n data matrix.
  lpi <- attr(X,"lpi")-1;m <- length(lpi)
  nb <- length(beta)
  if (is.null(dbeta)) {
    nsp = 0;dbeta <- dH <- 0
  } else {
    nsp = ncol(dbeta)
    dH = rep(0,nsp*nb*nb)
  }
  oo <- .C(C_mvn_ll,y=as.double(y),X=as.double(X),XX=as.double(crossprod(X)),beta=as.double(beta),n=as.integer(nrow(X)),
                  lpi=as.integer(lpi),m=as.integer(m),ll=as.double(0),lb=as.double(beta*0),
                  lbb=as.double(rep(0,nb*nb)), dbeta = as.double(dbeta), dH = as.double(dH), 
                  deriv = as.integer(nsp>0),nsp = as.integer(nsp),nt=as.integer(1))
  if (nsp==0) dH <- NULL else {
    dH <- list();ind <- 1:(nb*nb)
    for (i in 1:nsp) { 
      dH[[i]] <- matrix(oo$dH[ind],nb,nb)
      ind <- ind + nb*nb
    }
  }
  list(l=oo$ll,lb=oo$lb,lbb=matrix(oo$lbb,nb,nb),dH=dH)
}

pinv <- function(X,svd=FALSE) {
## a pseudoinverse for n by p, n>p matrices
  qrx <- qr(X,tol=0,LAPACK=TRUE)
  R <- qr.R(qrx);Q <- qr.Q(qrx) 
  rr <- Rrank(R) 
  if (svd&&rr<ncol(R)) {
    piv <- 1:ncol(X); piv[qrx$pivot] <- 1:ncol(X)
    er <- svd(R[,piv])
    d <- er$d*0;d[1:rr] <- 1/er$d[1:rr]
    X <- Q%*%er$u%*%(d*t(er$v))
  } else {
    Ri <- R*0 
    Ri[1:rr,1:rr] <- backsolve(R[1:rr,1:rr],diag(rr))
    X[,qrx$pivot] <- Q%*%t(Ri)
  }
  X
} ## end pinv

pqr2 <- function(x,nt=1,nb=30) {
## Function for parallel pivoted qr decomposition of a matrix using LAPACK
## householder routines. Currently uses a block algorithm.
## library(mgcv); n <- 4000;p<-3000;x <- matrix(runif(n*p),n,p)
## system.time(qrx <- qr(x,LAPACK=TRUE))
## system.time(qrx2 <- mgcv:::pqr2(x,2)) 
## system.time(qrx3 <- mgcv:::pqr(x,2)) 
## range(qrx2$qr-qrx$qr)
  p <- ncol(x)
  beta <- rep(0.0,p)
  piv <- as.integer(rep(0,p))
  ## need to force a copy of x, otherwise x will be over-written 
  ## by .Call *in environment from which function is called*
  x <- x*1  
  rank <- .Call(C_mgcv_Rpiqr,x,beta,piv,nt,nb)
  ret <- list(qr=x,rank=rank,qraux=beta,pivot=piv+1)
  attr(ret,"useLAPACK") <- TRUE
  class(ret) <- "qr"
  ret
} ## pqr2

pbsi <- function(R,nt=1,copy=TRUE) {
## parallel back substitution inversion of upper triangular R
## library(mgcv); n <- 5000;p<-4000;x <- matrix(runif(n*p),n,p)
## qrx <- mgcv:::pqr2(x,2);R <- qr.R(qrx)
## system.time(Ri <- mgcv:::pbsi(R,2))
## system.time(Ri2 <- backsolve(R,diag(p)));range(Ri-Ri2)
  if (copy) R <- R * 1 ## ensure that R modified only within pbsi
 .Call(C_mgcv_Rpbsi,R,nt)
 R
} ## pbsi

pchol <- function(A,nt=1,nb=30) {
## parallel Choleski factorization.
## library(mgcv);
## set.seed(2);n <- 200;r <- 190;A <- tcrossprod(matrix(runif(n*r),n,r))
## system.time(R <- chol(A,pivot=TRUE));system.time(L <- mgcv:::pchol(A));range(R[1:r,]-L[1:r,])
## system.time(L <- mgcv:::pchol(A,nt=2,nb=30))
## piv <- attr(L,"pivot");attr(L,"rank");range(crossprod(L)-A[piv,piv])
## should nb be obtained from 'ILAENV' as page 23 of Lucas 2004??
  piv <- as.integer(rep(0,ncol(A)))
  A <- A*1 ## otherwise over-write in calling env!
  rank <- .Call(C_mgcv_Rpchol,A,piv,nt,nb)
  attr(A,"pivot") <- piv+1;attr(A,"rank") <- rank
  A
}

pRRt <- function(R,nt=1) {
## parallel RR' for upper triangular R
## following creates index of lower triangular elements...
## n <- 4000;a <- rep(1:n,n);b <- rep(1:n,each=n)-1;which(a-b>0) -> ii;a[ii]+b[ii]*n->ii
## library(mgcv);R <- matrix(0,n,n);R[ii] <- runif(n*(n+1)/2)
## Note: A[a-b<=0] <- 0 zeroes upper triangle 
## system.time(A <- mgcv:::pRRt(R,2))
## system.time(A2 <- tcrossprod(R));range(A-A2)
  n <- nrow(R)
  A <- matrix(0,n,n)
  .Call(C_mgcv_RPPt,A,R,nt)
  A
}

block.reorder <- function(x,n.blocks=1,reverse=FALSE) {
## takes a matrix x divides it into n.blocks row-wise blocks, and re-orders 
## so that the blocks are stored one after the other. 
## e.g. library(mgcv); x <- matrix(1:18,6,3);xb <- mgcv:::block.reorder(x,2)
## x;xb;mgcv:::block.reorder(xb,2,TRUE)

 r = nrow(x);cols = ncol(x);
 if (n.blocks <= 1) return(x);
 if (r%%n.blocks) { 
   nb = ceiling(r/n.blocks)
 } else nb = r/n.blocks;
 oo <- .C(C_row_block_reorder,x=as.double(x),as.integer(r),as.integer(cols),
          as.integer(nb),as.integer(reverse));
 matrix(oo$x,r,cols)
} ## block.reorder


pqr <- function(x,nt=1) {
## parallel QR decomposition, using openMP in C, and up to nt threads (only if worthwhile)
## library(mgcv);n <- 20;p<-4;X <- matrix(runif(n*p),n,p);er <- mgcv:::pqr(X,nt=2)
  x.c <- ncol(x);r <- nrow(x)
  oo <- .C(C_mgcv_pqr,x=as.double(c(x,rep(0,nt*x.c^2))),as.integer(r),as.integer(x.c),
           pivot=as.integer(rep(0,x.c)), tau=as.double(rep(0,(nt+1)*x.c)),as.integer(nt)) 
  list(x=oo$x,r=r,c=x.c,tau=oo$tau,pivot=oo$pivot+1,nt=nt)
}

pqr.R <- function(x) {
## x is an object returned by pqr. This extracts the R factor...
## e.g. as pqr then...
## R <- mgcv:::pqr.R(er); R0 <- qr.R(qr(X,tol=0))
## svd(R)$d;svd(R0)$d
  oo <- .C(C_getRpqr,R=as.double(rep(0,x$c^2)),as.double(x$x),as.integer(x$r),as.integer(x$c),
           as.integer(x$c),as.integer(x$nt))
  matrix(oo$R,x$c,x$c)
}

pqr.qy <- function(x,a,tr=FALSE) {
## x contains a parallel QR decomp as computed by pqr. a is a matrix. computes
## Qa or Q'a depending on tr.
## e.g. as above, then...
## a <- diag(p);Q <- mgcv:::pqr.qy(er,a);crossprod(Q)
## X[,er$pivot+1];Q%*%R
## Qt <- mgcv:::pqr.qy(er,diag(n),TRUE);Qt%*%t(Qt);range(Q-t(Qt))
## Q <- qr.Q(qr(X,tol=0));z <- runif(n);y0<-t(Q)%*%z
## mgcv:::pqr.qy(er,z,TRUE)->y
## z <- runif(p);y0<-Q%*%z;mgcv:::pqr.qy(er,z)->y
  if (is.matrix(a)) a.c <- ncol(a) else a.c <- 1
  if (tr) {
    if (is.matrix(a)) { if (nrow(a) != x$r) stop("a has wrong number of rows") }
    else if (length(a) != x$r) stop("a has wrong number of rows")
  } else {
    if (is.matrix(a)) { if (nrow(a) != x$c) stop("a has wrong number of rows") }
    else if (length(a) != x$c)  stop("a has wrong number of rows")
    a <- c(a,rep(0,a.c*(x$r-x$c)))
  }
  oo <- .C(C_mgcv_pqrqy,a=as.double(a),as.double(x$x),as.double(x$tau),as.integer(x$r),
                         as.integer(x$c),as.integer(a.c),as.integer(tr),as.integer(x$nt))
  if (tr) return(matrix(oo$a[1:(a.c*x$c)],x$c,a.c)) else
  return(matrix(oo$a,x$r,a.c))
}

pmmult <- function(A,B,tA=FALSE,tB=FALSE,nt=1) {
## parallel matrix multiplication (not for use on vectors or thin matrices)
## library(mgcv);r <- 10;c <- 5;n <- 8
## A <- matrix(runif(r*n),r,n);B <- matrix(runif(n*c),n,c);range(A%*%B-mgcv:::pmmult(A,B,nt=1))
## A <- matrix(runif(r*n),n,r);B <- matrix(runif(n*c),n,c);range(t(A)%*%B-mgcv:::pmmult(A,B,TRUE,FALSE,nt=1))
## A <- matrix(runif(r*n),n,r);B <- matrix(runif(n*c),c,n);range(t(A)%*%t(B)-mgcv:::pmmult(A,B,TRUE,TRUE,nt=1))
## A <- matrix(runif(r*n),r,n);B <- matrix(runif(n*c),c,n);range(A%*%t(B)-mgcv:::pmmult(A,B,FALSE,TRUE,nt=1))

 if (tA) { n = nrow(A);r = ncol(A)} else {n = ncol(A);r = nrow(A)}
 if (tB) { c = nrow(B)} else {c = ncol(B)}
 C <- rep(0,r * c) 
 oo <- .C(C_mgcv_pmmult,C=as.double(C),as.double(A),as.double(B),as.integer(tA),as.integer(tB),as.integer(r),
          as.integer(c),as.integer(n),as.integer(nt));
 matrix(oo$C,r,c)
}
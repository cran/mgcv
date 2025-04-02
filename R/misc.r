## (c) Simon N. Wood 2011-2023
## Many of the following are simple wrappers for C functions

mchol <- function(A) {
## Simple wrapper for Matrix sparse Cholesky routine. Basically restores the
## functionality of Matrix::chol that vanished when the maintainers decided
## not to return the pivot sequence from Matrix::chol(foo, pivot=TRUE) (?!)
  if (inherits(A,"matrix")) suppressWarnings(return(chol(A,pivot=TRUE)))
  cha <- suppressWarnings(try(Matrix::Cholesky(A,perm=TRUE,super=NA),silent=TRUE))
  if (inherits(cha,"try-error")) {
    R <- -1;attr(R,"rank") <- -1 ## signal rank deficient
  } else { 
    R <- Matrix::triu(Matrix::expand1(cha,"L."))
    p <- ncol(A);
    attr(R,"pivot") <- if (length(cha@perm)==0) 1:p else cha@perm+1
    attr(R,"rank") <- p
  }  
  R ## R'R = H[pivot,pivot]
} ## mchol

dpnorm <- function(x0,x1) {
  ## Cancellation avoiding evaluation of pnorm(x1)-pnorm(x0) 
  ## first avoid 1-1 problems by exchanging and changing sign of double +ve
  ii <- x1>0&x0>0
  d <- x0[ii];x0[ii] <- -x1[ii];x1[ii] <- -d
  ## now deal with points that are so close that cancellation error
  ## too large - might as well use density times interval width
  ii <- abs(x1-x0) < sqrt(.Machine$double.eps)*dnorm((x1+x0)/2)
  p <- x0; d <- x1[ii]-x0[ii]; m <- (x1[ii]+x0[ii])/2
  p[ii] <- dnorm(m)*d
  p[!ii] <- pnorm(x1[!ii]) - pnorm(x0[!ii])
  p
} ## dpnorm


"%.%" <- function(a,b) {
  if (inherits(a,"dgCMatrix")||inherits(b,"dgCMatrix"))
  tensor.prod.model.matrix(list( ## following is coercion to double, general (no special structure), compressed column 
  as(as(as(a, "dMatrix"), "generalMatrix"), "CsparseMatrix"), 
  as(as(as(b, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  )) else tensor.prod.model.matrix(list(as.matrix(a),as.matrix(b)))
#  tensor.prod.model.matrix(list(as(a,"dgCMatrix"),as(b,"dgCMatrix"))) else - deprecated
#  tensor.prod.model.matrix(list(as.matrix(a),as.matrix(b)))
}

blas.thread.test <- function(n=1000,nt=4) {
  inter <- interactive()
  if (inter) prg <- txtProgressBar(min = 0, max = n, initial = 0,
              char = "=",width = NA, title="Progress", style = 3)
  for (i in 1:n) {
    m <- sample(213:654,1)
    p <- sample(3:153,1)
    X <- matrix(runif(m*p),m,p)
    er <- pqr(X,nt=nt) 
    rr <- range(pqr.qy(er,pqr.R(er))-X[,er$pivot])
    if (rr[1] < -1e-4||rr[2] > 1e-4) {
      break;
    }
    if (inter) setTxtProgressBar(prg, i)
  }
  if (inter) close(prg)
  if (rr[1] < -1e-4||rr[2] > 1e-4) {
    cat("BLAS thread safety problem at iteration",i,"\n")
  } else cat("No problem encountered in",i,"iterations\n")
} ## blas.thread.test



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

sdiag <- function(A,k=0) {
## extract sub or super diagonal of matrix (k=0 is leading)  
 p <- ncol(A)
 n <- nrow(A)
 if (k>p-1||-k > n-1) return()
 if (k >= 0) {
   i <- 1:n
   j <- (k+1):p
 } else {
   i <- (-k+1):n
   j <- 1:p
 }
 if (length(i)>length(j)) i <- i[1:length(j)] else j <- j[1:length(i)]
 ii <- i + (j-1) * n 
 A[ii]
} ## sdiag

"sdiag<-" <- function(A,k=0,value) {
 p <- ncol(A)
 n <- nrow(A)
 if (k>p-1||-k > n-1) return()
 if (k >= 0) {
   i <- 1:n
   j <- (k+1):p
 } else {
   i <- (-k+1):n
   j <- 1:p
 }
 if (length(i)>length(j)) i <- i[1:length(j)] else j <- j[1:length(i)]
 ii <- i + (j-1) * n 
 A[ii] <- value
 A
} ## "sdiag<-"

bandchol <- function(B) {
## obtain R such that R'R = A. Where A is banded matrix contained in R.
  n <- ncol(B)
  k <- 0
  if (n==nrow(B)) { ## square matrix. Extract the diagonals
    A <- B*0
    for (i in 1:n) {
      b <- sdiag(B,i-1)
      if (sum(b!=0)!=0) {
        k <- i ## largest index of a non-zero band
        A[i,1:length(b)] <- b
      }
    } 
    B <- A[1:k,]
  }
  oo <- .C(C_band_chol,B=as.double(B),n=as.integer(n),k=as.integer(nrow(B)),info=as.integer(0))
  if (oo$info<0) stop("something wrong with inputs to LAPACK routine")
  if (oo$info>0) stop("not positive definite")
  B <- matrix(oo$B,nrow(B),n)
  if (k>0) { ## was square on entry, so also on exit...
    A <- A * 0
    for (i in 1:k) sdiag(A,i-1) <- B[i,1:(n-i+1)]
    B <- A
  }
  B
} ## bandchol

trichol <- function(ld,sd) {
## obtain chol factor R of symm tridiag matrix, A, with leading diag
## ld and sub/super diags sd. R'R = A. On exit ld is diag of R and
## sd its super diagonal.
  n <- length(ld)
  if (n<2) stop("don't be silly")
  if (n!=length(sd)+1) stop("sd should have exactly one less entry than ld")
  oo <- .C(C_tri_chol,ld=as.double(ld),sd=as.double(sd),n=as.integer(n),info=as.integer(0))
  if (oo$info<0) stop("something wrong with inputs to LAPACK routine")
  if (oo$info>0) stop("not positive definite")
  ld <- sqrt(oo$ld)
  sd <- oo$sd*ld[1:(n-1)]
  list(ld=ld,sd=sd)
}

mgcv.omp <- function() {
## does open MP appear to be available?
  oo <- .C(C_mgcv_omp,a=as.integer(-1))
  if (oo$a==1) TRUE else FALSE
}

## discretized covariate routines...

XWXd <- function(X,w,k,ks,ts,dt,v,qc,nthreads=1,drop=NULL,ar.stop=-1,ar.row=-1,ar.w=-1,lt=NULL,rt=NULL) {
## Form X'WX given weights in w and X in compressed form in list X.
## each element of X is a (marginal) model submatrix. Full version 
## is given by X[[i]][k[,i],] (see below for summation convention).
## list X relates to length(ts) separate
## terms. ith term starts at matrix ts[i] and has dt[i] marginal matrices.
## For summation convention, k[,ks[j,1]:ks[j,2]] gives index columns
## for matrix j, thereby allowing summation over matrix covariates....
## i.e. for q in ks[j,1]:ks[j,2] sum up X[[j]][k[,q],] 
## Terms with several marginals are tensor products and may have 
## constraints (if qc[i]>1), stored as a householder vector in v[[i]]. 
## check ts and k index start (assumed 1 here)
## if drop is non-NULL it contains index of rows/cols to drop from result
## * lt is array of terms to include in left matrix (assumed in ascending coef index order)
## * rt is array of terms to include in right matrix (assumed in ascending coef index order)
## * if both NULL all terms are included, if only one is NULL then used for left and right. 
  m <- unlist(lapply(X,nrow));p <- unlist(lapply(X,ncol))
  nx <- length(X);nt <- length(ts)
  n <- length(w);ptfull <- pt <- 0;
  for (i in 1:nt) {
    fullsize <- prod(p[ts[i]:(ts[i]+dt[i]-1)])
    ptfull <- ptfull + fullsize
    pt <- pt + fullsize - if (qc[i]>0) 1 else if (qc[i]<0) v[[i]][v[[i]][1]+2] else 0
  } 
  if (inherits(X[[1]],"dgCMatrix")) { ## the marginals are sparse
    if (length(ar.stop)>1||ar.stop!=-1) warning("AR not available with sparse marginals")
    ## create list for passing to C
    if (any(qc<0)) stop("sparse method for Kronecker product contrasts not implemented")
    m <- list(Xd=X,kd=k,ks=ks,v=v,ts=ts,dt=dt,qc=qc)
    m$off <- attr(X,"off"); m$r <- attr(X,"r")
    if (is.null(m$off)||is.null(m$r)) stop("reverse indices missing from sparse discrete marginals")
    ### code that could create the marginal reverse indices...
    #for (j in 1:nrow(m$ks)) {
    #nr <- nrow(m$Xd[[j]]) ## make sure we always tab to final stored row 
    #for (i in m$ks[j,1]:(m$ks[j,2]-1)) {
    #  m$r[,i] <- (1:length(m$kd[,i]))[order(m$kd[,i])]
    #  m$off[[i]] <- cumsum(c(1,tabulate(m$kd[,i],nbins=nr)))-1
    #}
    m$offstart <- cumsum(c(0,lapply(m$off,length)))
    m$off <- unlist(m$off)
    ## Now C base all indices...
    m$ks <- m$ks - 1; m$kd <- m$kd - 1; m$r <- m$r - 1; m$ts <- m$ts-1
    nthreads <- as.integer(nthreads)
    w <- as.double(w)
    if ((!is.null(lt)||!is.null(rt))&&!is.null(drop)) {
      lpip <- attr(X,"lpip") ## list of coefs for each term
      rpi <- unlist(lpip[rt])
      lpi <- unlist(lpip[lt])
      if (is.null(lt)) lpi <- rpi
      else if (is.null(rt)) rpi <- lpi
      ldrop <- which(lpi %in% drop)
      rdrop <- which(lpi %in% drop)
    } else rdrop <- ldrop <- drop 

    if (is.null(lt)&&is.null(rt)) {
      lt <- rt <- 1:nt
    } else if (is.null(lt)) {
      lt <- rt
    } else if (is.null(rt)) rt <- lt;
    lt <- as.integer(lt-1)
    rt <- as.integer(rt-1)
    XWX <- .Call(C_sXWXd,m,w,lt,rt,nthreads)
    if (!is.null(drop)) {
      Dl <- Diagonal(ncol(XWX),1)
      XWX <- Dl[-ldrop,] %*% XWX %*% t(Dl[-rdrop,])
    }
    return(XWX) ## note that this is sparse
  } ## sparse case 
  ## block oriented code...
  if (is.null(lt)&&is.null(lt)) {
    # old .C code - can't handle long vector k
    #oo <- .C(C_XWXd0,XWX =as.double(rep(0,ptfull^2)),X= as.double(unlist(X)),w=as.double(w),
    #       k=as.integer(k-1),ks=as.integer(ks-1),m=as.integer(m),p=as.integer(p), n=as.integer(n), 
    #       ns=as.integer(nx), ts=as.integer(ts-1), as.integer(dt), nt=as.integer(nt),
    #       v = as.double(unlist(v)),qc=as.integer(qc),nthreads=as.integer(nthreads),
    #       ar.stop=as.integer(ar.stop-1),ar.weights=as.double(ar.w))
    #XWX <- if (is.null(drop)) matrix(oo$XWX[1:pt^2],pt,pt) else matrix(oo$XWX[1:pt^2],pt,pt)[-drop,-drop]
    XWX <- numeric(ptfull^2)
    .Call(C_CXWXd0,XWX,as.double(unlist(X)),w,k-1L,as.integer(ks-1L),as.integer(m),as.integer(p),
          as.integer(ts-1L), as.integer(dt),as.double(unlist(v)),as.integer(qc),as.integer(nthreads),
	  as.integer(ar.stop-1L),as.double(ar.w))
    XWX <- if (is.null(drop)) matrix(XWX[1:pt^2],pt,pt) else matrix(XWX[1:pt^2],pt,pt)[-drop,-drop]
    
  } else {
    lpip <- attr(X,"lpip") ## list of coefs for each term
    rpi <- unlist(lpip[rt])
    lpi <- unlist(lpip[lt])
    if (is.null(lt)) { ## note nrs and ncs not used in current code 
      lpi <- rpi
      nrs <- lt <- 0 ## currently lt == 0 signals nrs = 0 !
      ncs <- length(rt)
    } else {
      nrs <- length(lt)
      if (is.null(rt)) {
        rt <- ncs <- 0
	rpi <- lpi
      } else ncs <- length(rt)
    }  
    #oo <- .C(C_XWXd1,XWX =as.double(rep(0,ptfull^2)),X= as.double(unlist(X)),w=as.double(w),
    #       k=as.integer(k-1),ks=as.integer(ks-1),m=as.integer(m),p=as.integer(p), n=as.integer(n), 
    #       ns=as.integer(nx), ts=as.integer(ts-1), dt=as.integer(dt), nt=as.integer(nt),
    #       v = as.double(unlist(v)),qc=as.integer(qc),nthreads=as.integer(nthreads),
    #       ar.stop=as.integer(ar.stop-1),ar.weights=as.double(ar.w),rs=as.integer(lt-1),
    # 	   cs=as.integer(rt-1),nrs=as.integer(nrs),ncs=as.integer(ncs))#)	   
    #XWX <- matrix(oo$XWX[1:(length(lpi)*length(rpi))],length(lpi),length(rpi))
    XWX <- numeric(ptfull^2)
    .Call(C_CXWXd1,XWX,as.double(unlist(X)),w,k-1L,as.integer(ks-1L),as.integer(m),as.integer(p),
          as.integer(ts-1L), as.integer(dt),as.double(unlist(v)),as.integer(qc),as.integer(nthreads),
	  as.integer(ar.stop-1L),as.double(ar.w),as.integer(lt-1L),as.integer(rt-1L))
    XWX <- matrix(XWX[1:(length(lpi)*length(rpi))],length(lpi),length(rpi))
    if (!is.null(drop)) {
      ldrop <- which(lpi %in% drop)
      rdrop <- which(lpi %in% drop)
      if (length(ldrop)>0||length(rdrop)>0) XWX <-
        if (length(ldrop==0)) XWX[,-rdrop] else if (length(rdrop)==0) XWX[-ldrop,] else XWX[-ldrop,-rdrop]
    }
  }
  XWX
} ## XWXd


XVXd <- function(X,e,k,ks,ts,dt,v,qc,nthreads=1,a,ma) {
## e is a residual vector. a[ma[i-1]+1:ma[i]] contains the indices of the neighbours of point i
## (ma[-1]=-1, by convention - C indexing assumed).
## This routine computes X'VX where V[i,j] = e[i]*e[j] if i is a neighbour of j and 0 otherwise.
  m <- unlist(lapply(X,nrow));p <- unlist(lapply(X,ncol))
  nx <- length(X);nt <- length(ts)
  n <- length(e);ptfull <- pt <- 0;
  for (i in 1:nt) {
    fullsize <- prod(p[ts[i]:(ts[i]+dt[i]-1)])
    ptfull <- ptfull + fullsize
    pt <- pt + fullsize - if (qc[i]>0) 1 else if (qc[i]<0) v[[i]][v[[i]][1]+2] else 0
  }
  XVX <- numeric(ptfull^2)
  .Call(C_CXVXd0,XVX,as.double(unlist(X)),e,k-1L,as.integer(ks-1L),as.integer(m),as.integer(p),
          as.integer(ts-1L), as.integer(dt),as.double(unlist(v)),as.integer(qc),as.integer(nthreads),
	  as.integer(a),ma)
  XVX <- matrix(XVX[1:pt^2],pt,pt)
  XVX
} ## XVXd


XWyd <- function(X,w,y,k,ks,ts,dt,v,qc,drop=NULL,ar.stop=-1,ar.row=-1,ar.w=-1,lt=NULL) {
## X'Wy...
## if lt if not NULL then it lists the discrete terms to include (from X)
## returned vector/matrix only includes rows for selected terms
  m <- unlist(lapply(X,nrow));p <- unlist(lapply(X,ncol))
  nx <- length(X);nt <- length(ts)
  n <- length(w);
  if (is.null(lt)) {
    pt <- 0
    for (i in 1:nt) pt <- pt + prod(p[ts[i]:(ts[i]+dt[i]-1)]) - if (qc[i]>0) 1 else if (qc[i]<0) v[[i]][v[[i]][1]+2] else 0
    lt <- 1:nt
  } else {
    lpip <- attr(X,"lpip") ## list of coefs for each term 
    lpi <- unlist(lpip[lt]) ## coefs corresponding to terms selected by lt
    if (!is.null(drop)) drop <- which(lpi %in% drop) ## rebase drop 
    pt <- length(lpi)
  }
  cy <- if (is.matrix(y)) ncol(y) else 1
  if (inherits(X[[1]],"dgCMatrix")) { ## the marginals are sparse
    ## create list for passing to C
    m <- list(Xd=X,kd=k,ks=ks,v=v,ts=ts,dt=dt,qc=qc)
    m$off <- attr(X,"off"); m$r <- attr(X,"r")
    if (is.null(m$off)||is.null(m$r)) stop("reverse indices missing from sparse discrete marginals")
    m$offstart <- cumsum(c(0,lapply(m$off,length)))
    m$off <- unlist(m$off)
    ## Now C base all indices...
    m$ks <- m$ks - 1; m$kd <- m$kd - 1; m$r <- m$r - 1; m$ts <- m$ts-1
    Wy <- as.double(w*y);lt <- as.integer(lt-1)
    XWy <- .Call(C_sXyd,m,Wy,lt)
    if (cy>1) XWy <- matrix(XWy,ncol=cy)
    if (!is.null(drop)) XWy <- if (cy>1) XWy[-drop,] else XWy[-drop]
  } else { ## dense marginals case  
    ## old .C code - can't handle long vector k
    #oo <- .C(C_XWyd,XWy=rep(0,pt*cy),y=as.double(y),X=as.double(unlist(X)),w=as.double(w),k=as.integer(k-1), 
    #       ks=as.integer(ks-1),
    #       m=as.integer(m),p=as.integer(p),n=as.integer(n),cy=as.integer(cy), nx=as.integer(nx), ts=as.integer(ts-1), 
    #       dt=as.integer(dt),nt=as.integer(nt),v=as.double(unlist(v)),qc=as.integer(qc),
    #       ar.stop=as.integer(ar.stop-1),ar.row=as.integer(ar.row-1),ar.weights=as.double(ar.w),
    # 	   cs=as.integer(lt-1),ncs=as.integer(length(lt)))
    #if (cy>1) XWy <- if (is.null(drop)) matrix(oo$XWy,pt,cy) else matrix(oo$XWy,pt,cy)[-drop,] else
    #XWy <- if (is.null(drop)) oo$XWy else oo$XWy[-drop]
    XWy <- numeric(pt*cy)
    .Call(C_CXWyd,XWy,as.double(y),as.double(unlist(X)),as.double(w),k-1L,as.integer(ks-1L),as.integer(m),as.integer(p),
          as.integer(cy),as.integer(ts-1L), as.integer(dt),as.double(unlist(v)),as.integer(qc),as.integer(ar.stop-1L),
	  as.integer(ar.row-1L),as.double(ar.w),as.integer(lt-1L))
    if (cy>1) { XWy <- if (is.null(drop)) matrix(XWy,pt,cy) else matrix(XWy,pt,cy)[-drop,] } else {
      XWy <- if (is.null(drop)) XWy else XWy[-drop]
    }
  }  
  XWy
} ## XWyd 

Xbd <- function(X,beta,k,ks,ts,dt,v,qc,drop=NULL,lt=NULL) {
## note that drop may contain the index of columns of X to drop before multiplying by beta.
## equivalently we can insert zero elements into beta in the appropriate places.
## if lt if not NULL then it lists the discrete terms to include (from X)
  n <- if (is.matrix(k)) nrow(k) else length(k) ## number of data
  m <- unlist(lapply(X,nrow)) ## number of rows in each discrete model matrix
  p <- unlist(lapply(X,ncol)) ## number of cols in each discrete model matrix
  nx <- length(X) ## number of model matrices
  if (length(p)!=nx) stop("something wrong with matrix list - not all matrices?")
  nt <- length(ts) ## number of terms
  if (!is.null(drop)) { 
    b <- if (is.matrix(beta)) matrix(0,nrow(beta)+length(drop),ncol(beta)) else rep(0,length(beta)+length(drop))
    if (is.matrix(beta)) b[-drop,] <- beta else b[-drop] <- beta
    beta <- b
  }
  if (is.null(lt)) lt <- 1:nt
 
  bc <- if (is.matrix(beta)) ncol(beta) else 1 ## number of columns in beta
  ## The C code mechanism for dealing with lt is very basic, and requires that beta is re-ordered and
  ## truncated to relate only to the selected terms, in the order they are selected.  
  lpip <- attr(X,"lpip")
  if (!is.null(lpip)) { ## then X list may not be in coef order...
    lpip <- unlist(lpip[lt])
    beta <- if (is.matrix(beta)) beta[lpip,] else beta[lpip] ## select params required in correct order
  }
  if (inherits(X[[1]],"dgCMatrix")) { ## the marginals are sparse
    ## create list for passing to C
    m <- list(Xd=X,kd=k,ks=ks,v=v,ts=ts,dt=dt,qc=qc)
    m$off <- attr(X,"off"); m$r <- attr(X,"r")
    if (is.null(m$off)||is.null(m$r)) stop("reverse indices missing from sparse discrete marginals")
    m$offstart <- cumsum(c(0,lapply(m$off,length)))
    m$off <- unlist(m$off)
    ## Now C base all indices...
    m$ks <- m$ks - 1; m$kd <- m$kd - 1; m$r <- m$r - 1; m$ts <- m$ts-1
    beta <- as.matrix(beta);storage.mode(beta) <- "double"
    lt <- as.integer(lt-1)
    Xb <- .Call(C_sXbd,m,beta,lt)
    if (bc>1) Xb <- matrix(Xb,ncol=bc)
  } else { ## dense marginals case
  
    #oo <- .C(C_Xbd,f=as.double(rep(0,n*bc)),beta=as.double(beta),X=as.double(unlist(X)),k=as.integer(k-1),
    #       ks = as.integer(ks-1), 
    #       m=as.integer(m),p=as.integer(p), n=as.integer(n), nx=as.integer(nx), ts=as.integer(ts-1), 
    #       as.integer(dt), as.integer(nt),as.double(unlist(v)),as.integer(qc),as.integer(bc),as.integer(lt-1),as.integer(length(lt)))
    #Xb <- if (is.matrix(beta)) matrix(oo$f,n,bc) else oo$f
    f <- numeric(n*bc)
    .Call(C_CXbd,f,beta,as.double(unlist(X)),k-1L, as.integer(ks-1L),as.integer(m),as.integer(p),
          as.integer(ts-1L),as.integer(dt),as.double(unlist(v)),as.integer(qc),as.integer(bc), as.integer(lt-1L))
    Xb <- if (is.matrix(beta)) matrix(f,n,bc) else f	  
  }
  return(Xb)
} ## Xbd

diagXVXd <- function(X,V,k,ks,ts,dt,v,qc,drop=NULL,nthreads=1,lt=NULL,rt=NULL) {
## discrete computation of diag(XVX')
  workXVXd(X,V,k,k1=NULL,ks,ts,dt,v,qc,drop,nthreads,lt,rt)
} ## diagXVXd

ijXVXd <- function(i,j,X,V,k,ks,ts,dt,v,qc,drop=NULL,nthreads=1,lt=NULL,rt=NULL) {
## discrete computation of scattered elements XVX'[i,j]. i and j are vectors of indices
  if (is.matrix(k)) {
    ki <- k[i,];kj <- k[j,]
  } else {
    ki <- k[i];kj <- k[j]
  }
  workXVXd(X,V,ki,kj,ks,ts,dt,v,qc,drop,nthreads,lt,rt)
} ## ijXVXd

workXVXd <- function(X,V,k,k1,ks,ts,dt,v,qc,drop=NULL,nthreads=1,lt=NULL,rt=NULL) {
## discrete computation of diag(XVX') or scattered elements of XWX'

  n <- if (is.matrix(k)) nrow(k) else length(k)
  m <- unlist(lapply(X,nrow));p <- unlist(lapply(X,ncol))
  nx <- length(X);nt <- length(ts)
  if (is.null(lt)&&is.null(rt)) rt <- lt <- 1:nt else {
    if (is.null(rt)) rt <- lt else if (is.null(lt)) lt <- rt
  }  
  if (is.null(rt)) rt <- 1:nt
  if (inherits(X[[1]],"dgCMatrix")) { ## the marginals are sparse
    if (!is.null(k1)) stop("scattered XVX computation not yet implemented") 
    ## create list for passing to C
    m <- list(Xd=X,kd=k,ks=ks,v=v,ts=ts,dt=dt,qc=qc)
    m$off <- attr(X,"off"); m$r <- attr(X,"r")
    if (is.null(m$off)||is.null(m$r)) stop("reverse indices missing from sparse discrete marginals")
    m$offstart <- cumsum(c(0,lapply(m$off,length)))
    m$off <- unlist(m$off)
    ## Now C base all indices...
    m$ks <- m$ks - 1; m$kd <- m$kd - 1; m$r <- m$r - 1; m$ts <- m$ts-1
    
    if (!is.null(drop)) {
      D <- Diagonal(ncol(V)+length(drop),1)[-drop,]
      V <- t(D) %*% V %*% D
    }
    ## The C code mechanism for dealing with rt and lt is very basic, and requires that V is
    ## re-ordered and truncated to relate only to the selected terms, in the order they are selected.  
    lpip <- attr(X,"lpip")
    if (!is.null(lpip)) { ## then X list may not be in coef order...
      lpi <- unlist(lpip[lt])
      rpi <- unlist(lpip[rt])
      V <- V[lpi,rpi,drop=FALSE] ## select part of V required in correct order
    }
    lt <- as.integer(lt-1);rt <- as.integer(rt-1)
    D <- .Call(C_sdiagXVXt,m , V, lt, rt)
  } else { ## dense marginals
    if (!is.null(drop)) { 
      pv <- ncol(V)+length(drop)
      V0 <- matrix(0,pv,pv)
      V0[-drop,-drop] <- V
      V <- V0;rm(V0)
    } else pv <- ncol(V)
    ## The C code mechanism for dealing with rt and lt is very basic, and requires that V is
    ## re-ordered and truncated to relate only to the selected terms, in the order they are selected.  
    lpip <- attr(X,"lpip")
    if (!is.null(lpip)) { ## then X list may not be in coef order...
      lpi <- unlist(lpip[lt])
      rpi <- unlist(lpip[rt])
      V <- V[lpi,rpi,drop=FALSE] ## select part of V required in correct order
    }
    D <- numeric(n)
    if (is.null(k1)) .Call(C_CdiagXVXt,D,V,as.double(unlist(X)),k-1L,as.integer(ks-1L),as.integer(m),as.integer(p),
          as.integer(ts-1L), as.integer(dt),as.double(unlist(v)),as.integer(qc),as.integer(nthreads),
	  as.integer(lt-1L),as.integer(rt-1L)) else
    .Call(C_CijXVXt,D,V,as.double(unlist(X)),k-1L,k1-1L,as.integer(ks-1L),as.integer(m),as.integer(p),
          as.integer(ts-1L), as.integer(dt),as.double(unlist(v)),as.integer(qc),as.integer(nthreads),
	  as.integer(lt-1L),as.integer(rt-1L))	  
  }
  D
} ## workXVXd

dchol <- function(dA,R) {
## if dA contains matrix dA/dx where R is chol factor s.t. R'R = A
## then this routine returns dR/dx...
  p <- ncol(R)
  oo <- .C(C_dchol,dA=as.double(dA),R=as.double(R),dR=as.double(R*0),p=as.integer(ncol(R)))
  return(matrix(oo$dR,p,p))
} ## dchol

choldrop <- function(R,k) {
## routine to update Cholesky factor R of A on dropping row/col k of A.
## R can be upper triangular, in which case (R'R=A) or lower triangular in
## which case RR'=A...
  n <- as.integer(ncol(R))
  k1 <- as.integer(k-1)
  ut <- as.integer(as.numeric(R[1,2]!=0))
  if (k<1||k>n) return(R)
  Rup <- matrix(0,n-1,n-1)
  #oo <- .C(C_chol_down,R=as.double(R),Rup=as.double(Rup),n=as.integer(n),k=as.integer(k-1),ut=as.integer(ut))
  .Call(C_mgcv_chol_down,R,Rup,n,k1,ut)
  #matrix(oo$Rup,n-1,n-1)
  Rup
} ## choldrop

cholup <- function(R,u,up=TRUE) {
## routine to update Cholesky factor R to the factor of R'R + uu' (up == TRUE)
## or R'R - uu' (up=FALSE). 
  n <- as.integer(ncol(R))
  up <- as.integer(up)
  eps <- as.double(.Machine$double.eps)
  R1 <- R * 1.0
  .Call(C_mgcv_chol_up,R1,u,n,up,eps)
  if (up==0) if ((n>1 && R1[2,1] < -1)||(n==1&&u[1]>R[1])) stop("update not positive definite")
  R1
} ## cholup


vcorr <- function(dR,Vr,trans=TRUE) {
## Suppose b = sum_k op(dR[[k]])%*%z*r_k, z ~ N(0,Ip), r ~ N(0,Vr). vcorr returns cov(b).
## dR is a list of p by p matrices. 'op' is 't' if trans=TRUE and I() otherwise.
  p <- ncol(dR[[1]])
  M <- if (trans) ncol(Vr) else -ncol(Vr) ## sign signals transpose or not to C code
  if (abs(M)!=length(dR)) stop("internal error in vcorr, please report to simon.wood@r-project.org")
  oo <- .C(C_vcorr,dR=as.double(unlist(dR)),Vr=as.double(Vr),Vb=as.double(rep(0,p*p)),
           p=as.integer(p),M=as.integer(M))
  return(matrix(oo$Vb,p,p))
} ## vcorr


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
  ret <- list(qr=x,rank=rank,qraux=beta[1:rank],pivot=piv+1)
  attr(ret,"useLAPACK") <- TRUE
  class(ret) <- "qr"
  ret
} ## pqr2

pbsi <- function(R,nt=1,copy=TRUE) {
## parallel back substitution inversion of upper triangular R
## library(mgcv); n <- 500;p<-400;x <- matrix(runif(n*p),n,p)
## qrx <- qr(x);R <- qr.R(qrx)
## system.time(Ri <- mgcv:::pbsi(R,2))
## system.time(Ri2 <- backsolve(R,diag(p)));range(Ri-Ri2)
  if (copy) R <- R * 1 ## ensure that R modified only within pbsi
 .Call(C_mgcv_Rpbsi,R,nt)
 R
} ## pbsi

pchol <- function(A,nt=1,nb=40) {
## parallel Choleski factorization.
## library(mgcv);
## set.seed(2);n <- 200;r <- 190;A <- tcrossprod(matrix(runif(n*r),n,r))
## system.time(R <- chol(A,pivot=TRUE));system.time(L <- mgcv:::pchol(A));range(R[1:r,]-L[1:r,])
## k <- 30;range(R[1:k,1:k]-L[1:k,1:k])
## system.time(L <- mgcv:::pchol(A,nt=2,nb=30))
## piv <- attr(L,"pivot");attr(L,"rank");range(crossprod(L)-A[piv,piv])
## should nb be obtained from 'ILAENV' as page 23 of Lucas 2004??
  piv <- as.integer(rep(0,ncol(A)))
  A <- A*1 ## otherwise over-write in calling env!
  rank <- .Call(C_mgcv_Rpchol,A,piv,nt,nb)
  attr(A,"pivot") <- piv+1;attr(A,"rank") <- rank
  A
}

pforwardsolve <- function(R,B,nt=1) {
## parallel forward solve via simple col splitting...
 if (!is.matrix(B)) B <- as.matrix(B)
 .Call(C_mgcv_Rpforwardsolve,R,B,nt)

}

pcrossprod <- function(A,trans=FALSE,nt=1,nb=30) {
## parallel cross prod A'A or AA' if trans==TRUE...
 if (!is.matrix(A)) A <- as.matrix(A)
 if (trans) A <- t(A)
 .Call(C_mgcv_Rpcross,A,nt,nb)
}

pRRt <- function(R,nt=1) {
## parallel RR' for upper triangular R
## following creates index of lower triangular elements...
## n <- 4000;a <- rep(1:n,n);b <- rep(1:n,each=n);which(a>=b) -> ii;a[ii]+(b[ii]-1)*n->ii ## lower
## n <- 4000;a <- rep(1:n,n);b <- rep(1:n,each=n);which(a<=b) -> ii;a[ii]+(b[ii]-1)*n->ii ## upper
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
## range(mgcv:::pqr.qy(er,mgcv:::pqr.R(er))-X[,er$pivot])
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

treig <- function(ld,sd,vec=FALSE,descend=FALSE) {
## eigen decomposition of tri-diagonal matrix with leading diagonal ld and
## sub-diagonal sd...
  n <- length(ld)
  v <- if (vec) rep(0,n*n) else 0
  oo <- .C(C_mgcv_trisymeig,d=as.double(ld),g=as.double(sd),v=as.double(v),n=as.integer(n),
                            get.vec=as.integer(vec),descending=as.integer(descend))
  v <- if (vec) matrix(oo$v,n,n) else NA
  list(values=oo$d,vectors=v)
} ## treig


lanczos <- function(A,v0,M,Av = function(A,v) A%*%v,n=ncol(A)) {
## Apply M steps of Lanczos starting at v0 for n by n +ve semi definite matrix A.
## Av is a function forming the product of matrix A with vector v.
## A can be a matrix, in which case the default Av applies, or
## it could be a list of arguments used to define the multiplication
## in some other way (e.g. as a sequence of matrix products) defined
## in a custom Av...
## The function is used to find an approximate eigen value CDF for A
## as described in Lin, Saad and Yang (2016) SIAM Review 58(1), 34-65
## section 3.2.1 in particular.
  v0.norm <- sqrt(sum(v0^2)) 
  gamma <- epsilon <- rep(0,M)
  q <- matrix(v0/v0.norm,n,M)
  for (j in 1:M) {
    c <- Av(A,q[,j])
    if (j==1) vAv <- sum(v0*c)*v0.norm ## v0'Av0 - useful for estimating tr(A)
    gamma[j] <- sum(c*q[,j])
    c <- c - gamma[j]*q[,j]  
    if (j>1) {
      c <- c - epsilon[j-1]*q[,j-1]
      cq <- drop(t(c) %*% q[,1:j,drop=FALSE]) 
      c <- c - colSums(cq*t(q[,1:j]))
      cq <- drop(t(c) %*% q[,1:j,drop=FALSE]) 
      c <- c - colSums(cq*t(q[,1:j]))
    }
    epsilon[j] <- sqrt(sum(c^2))
    if (j<M) q[,j+1] <- c/epsilon[j]
  }
  et <- treig(gamma,epsilon,descend=FALSE,vec=TRUE)
  ## compute error bounds on the eigenvalues of
  ## A using the method described in section 3.2 of
  ## Parlett, BN (1998) The Symmetric Eigenvalue Problem, SIAM
  err <- abs(et$vectors[M,])*epsilon[M]
  theta <- c(0,et$values)
  tau <- et$vectors[1,]
  eta <- c(0,cumsum(tau^2))
  ## theta is vector of eigenvalues at which CDF jumps
  ## tau^2 is jump size, eta is CDF at theta
  ## lam.ub is upper bound on larget eigenvalue
  list(theta=theta,tau=tau,eta=eta,err=err,vAv=vAv)
} ## lanczos

eigen.approx <- function(A,Av = function(A,v) A%*%v,M=20,n.rep=20,n=ncol(A),seed=1) {
## get the approximate eigenvalues of n by n +ve semi def matrix A. Av is
## the function for multiplying a vector by the matrix defined by A. If A
## is simply a matrix then the default Av is sufficient. M is the number of
## Lanczos steps to use, and n.rep the number of random replicates to average
## over. Routine restores RNG to pre-call state on exit.
## Based on Lin, Saad and Yang (2016) SIAM Review 58(1), 34-65
## section 3.2.1, but extended to only use this approximation for the
## eigen-values that have yet to converge.
## The CDF approximation is based on their Appendix C proposal, rather than
## using a Gausssian kernel approximation to the pdf and cdf and then
## inverting by tabulation. This is because the latter tends to oversmooth the
## CDF in a way that is unhelpful for rank deficient matrices.
## The Gaussian kernel approach would probably be prefereable for full rank matrices,
## since it then benefits from the extra stability of kernel smoothing.
   if (is.finite(seed)) a <- temp.seed(seed) ## seed RNG and store state
   eva <- rep(0,n)
   trA <- rep(0,n.rep)
   tol <- .Machine$double.eps^.5
   for (r in 1:n.rep) {
     v0 <- rnorm(n)
     lz <- lanczos(A,v0,M=M,Av=Av,n=n)
     trA[r] <- lz$vAv
     ## following is suggested in Appendix C of LSY, and is quite important
     ## to avoid slight downward bias...
     eta1 <- c(0,(lz$eta[1:M]*0.5+0.5*lz$eta[1:M+1])) 
     eta1[2] <-  lz$eta[2] ## correction to avoid over-estimation in lower tail if rank def
     conv <- lz$err<lz$theta[M+1]*tol ## these eigenvalues are converged
     upper.uconv <- if (any(!conv)) max(which(!conv)) else 0 ## last uncoverged
     n.conv <- M - upper.uconv ## number converged
     lz$theta[lz$theta<0] <- 0
     theta.conv <- if (n.conv) lz$theta[(upper.uconv+1):M+1] else rep(0,0)
     if (upper.uconv) {
       eta <- eta1[1:(upper.uconv+1)]
       theta <- lz$theta[1:(upper.uconv+1)]
       eta <- eta/max(eta)
       nri <- c(diff(eta)!=0,TRUE) ## strip out duplicates
       eva <- eva + c(approx(eta[nri],theta[nri],seq(0,1,length=n-n.conv),method="linear",rule=2)$y,theta.conv)
     } else {
       eva <- eva + c(rep(0,n-n.conv),theta.conv)
     }
   }
   if (is.finite(seed)) temp.seed(a) ## restore RNG state
   eva <- eva/n.rep
   trA.sd <- sd(trA)/sqrt(n.rep);trA <- mean(trA)
   if (abs(sum(eva)-trA)>2.5*trA.sd) { ## evidence for bias in eigen-spectrum
     eva <- eva*trA/sum(eva) ## correction
   }
   eva
} ## eigen.approx

temp.seed <- function(x) {
## when called with a numeric x stores the state of the RNG and sets its seed to
## x. Returns an object of class "rng.state". When called with an object of this
## class created on a previus call to this function, resets the RNG to the state
## on entry to that previous call. 
  if (inherits(x,"rng.state")) { 
    RNGkind(x$kind[1],x$kind[2])
    assign(".Random.seed",x$seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
  } else {
    seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
    if (inherits(seed,"try-error")) {
          runif(1)
          seed <- get(".Random.seed",envir=.GlobalEnv)
    }
    kind <- RNGkind(NULL)
    RNGkind("default","default")
    set.seed(x) ## ensure repeatability
    x <- list(seed=seed,kind=kind)
    class(x) <- "rng.state"
    return(x)
  }  
} ## temp.seed

mat.rowsum <- function(X,m,k) {
## Let X be n by p and m of length M. Produces an M by p matrix B
## where the ith row of B is the sum of the rows X[k[j],] where
## j = (m[i-1]+1):m[i]. m[0] is taken as 0.
## n <- 10;p <- 5;X <- matrix(runif(n*p),n,p)
## m <- c(3,5,8,11);k <- c(1,4,3,6,1,5,7,10,9,5,6)
## mgcv:::mat.rowsum(X,m,k)
  if (max(k)>nrow(X)||min(k)<1) stop("index vector has invalid entries")
  k <- k - 1 ## R to C index conversion
  .Call(C_mrow_sum,X,m,k)
} ## mat.rowsum

isa <- function(R,nt=1) {
## Finds the elements of (R'R)^{-1} on NZP(R+R').
  if (!inherits(R,c("dgCMatrix","dtCMatrix"))) stop("isa requires a dg/tCMatrix")
  nt <- round(nt)
  if (nt<1) nt = 1
  Hpi <- R + t(R)
  .Call(C_isa1p,t(R),Hpi,nt)
  Hpi
} ## isa

AddBVB <- function(A,Bt,VBt) {
## Add B %*% V %*% t(B) to calss 'dgCMatrix' A returning result on NZP(A) only
## (i.e. discarding elements of BVB' not in NZP(A)), Bt is the transpose
## of B. B and VBt are class 'matrix'
  A@x <- A@x * 1.0 ## force copy, otherwise A and return value modified
  .Call(C_AddBVB,A,Bt,VBt)
  A
} ## AddBVB

minres <- function(R,u,b) {
## routine to solve (R'R-uu')x = b using minres algorithm.
## set.seed(0);n <- 100;p <- 20;X <- matrix(runif(n*p)-.5,n,p);R <- chol(crossprod(X));b <- runif(p);k <- 1;
## solve(crossprod(X[-k,]),b);mgcv:::minres(R,t(X[k,]),b)
  x <- b; p <- length(b);
  m <- if (is.matrix(u)) ncol(u) else 1
  work <- rep(0,p*(m+7)+m)
  oo <- .C(C_minres,R=as.double(R), u=as.double(u),b=as.double(b), x=as.double(x), p=as.integer(p),m=as.integer(m),work=as.double(work))
  cat("\n niter : ",oo$m,"\n")
  oo$x
}

neicov <- function(Dd,D1=NULL,nei) {
## wrapper for nei_cov. Dd is n by p matrix of leave one out perturbations to
## coef vectors. nei is neighbourhood structure. If D1 is not NULL then it is
## perturbation matrix to be used on RHS of computation - used for correction
## terms.
  p <- ncol(Dd)
  V <- matrix(0,p,p)
  a <- nei$a-1
  if (is.null(D1)) .Call(C_nei_cov,V,Dd,Dd,nei$ma,a) else .Call(C_nei_cov,V,Dd,D1,nei$ma,a)
  (V+t(V))/2
} ## neicov


## Routines to force matrix to pdef

pdev <- function(A) {
## Force A to meet necessary conditions for +ve def from Thm 4.2.8. of Golub and van Load 4th ed.
## Non-positive diagonals are set to diagonal dominance. Off diagonals are then set to just meet
## necessary conditions.
## NOTE: A is modified directly in situ, so care needed to force a copy to be kept if it's needed
## A0 <- A will not work, as R only actually copies when R modifies!
## A0 <- A; A0[1,1] <- A0[1,1] + 1;A0[1,1] <- A0[1,1] - 1  works.  
  if (inherits(A,"Matrix")) {
    if (!inherits(A,"dgCMatrix")) A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    da <- diag(A); ii <- which(da==0)
    if (length(ii)) diag(A)[ii] <- -1e-30 ## Avoid C code having to insert extra non-zeroes
    mod <- .Call(C_spdev,A)
  } else { ## dense matrix
    mod <- .Call(C_dpdev,A)
  }
  if (mod>0) attr(A,"modified") <- mod
  return(A)
}

## following are wrappers for KP STZ constraints - intended for testing only

Zb <- function(b0,v,qc,p,w) {
  b1 <- rep(0,p)
  oo <- .C(C_Zb,b1=as.double(b1),as.double(b0),as.double(v),as.integer(qc),as.integer(p),as.double(w))
  oo$b1
}

Ztb <- function(b0,v,qc,di,p,w) {
  ## p is length(b0)/di
  w <- rep(0,2*p)
  M <- v[1]
  pp <- p
  for (i in 1:M) pp <- pp/v[i+1];
  p0 <- prod(v[1+1:M]-1)*pp
  b1 <- rep(0,p0*di)
  oo <- .C(C_Ztb,b1=as.double(b1),as.double(b0),as.double(v),as.integer(qc),as.integer(di),as.integer(p),as.double(w))
  oo$b1
}


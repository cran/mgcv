## routines for very large dataset generalized additive modelling.
## (c) Simon N. Wood 2009-2023


ls.size <- function(x) {
## If `x' is a list, return the size of its elements, in bytes, in a named array
## otherwise return the size of the object
 if (is.list(x)==FALSE) return(object.size(x))

 xn <- names(x)
 n <- length(x)
 sz <- rep(-1,n)
 for (i in 1:n) sz[i] <- object.size(x[[i]])
 names(sz) <- xn
 sz
} ## ls.size

rwMatrix <- function(stop,row,weight,X,trans=FALSE) {
## Routine to recombine the rows of a matrix X according to info in 
## stop, row and weight. Consider the ith row of the output matrix 
## ind <- 1:stop[i] if i==1 and ind <- (stop[i-1]+1):stop[i]
## otherwise. The ith output row is then X[row[ind],]*weight[ind]
  if (is.matrix(X)) { n <- nrow(X);p<-ncol(X);ok <- TRUE} else { n<- length(X);p<-1;ok<-FALSE}
  stop <- stop - 1;row <- row - 1 ## R indices -> C indices
  oo <-.C(C_rwMatrix,as.integer(stop),as.integer(row),as.double(weight),X=as.double(X),
          as.integer(n),as.integer(p),trans=as.integer(trans),work=as.double(rep(0,n*p)))
  if (ok) return(matrix(oo$X,n,p)) else
  return(oo$X) 
} ## rwMatrix

chol2qr <- function(XX,Xy,nt=1) {
## takes X'X and X'y and returns R and f
## equivalent to qr update.  
  op <- options(warn = -1) ## otherwise warns if +ve semidef
  R <- if (nt) pchol(XX,nt=nt) else chol(XX,pivot=TRUE)
  options(op)
  p <- length(Xy)
  ipiv <- piv <- attr(R,"pivot");ipiv[piv] <- 1:p
  rank <- attr(R,"rank");ind <- 1:rank
  if (rank<p) R[(rank+1):p,] <- 0 ## chol is buggy (R 3.1.0) with pivot=TRUE
  f <- c(forwardsolve(t(R[ind,ind]),Xy[piv[ind]]),rep(0,p-rank))[ipiv]
  R <- R[ipiv,ipiv]
  list(R=R,f=f)
} ## chol2qr

qr_update <- function(Xn,yn,R=NULL,f=rep(0,0),y.norm2=0,use.chol=FALSE,nt=1)
## Let X = QR and f = Q'y. This routine updates f and R
## when Xn is appended to X and yn appended to y. If R is NULL
## then initial QR of Xn is performed. ||y||^2 is accumulated as well.
## if use.chol==TRUE then quicker but less stable accumulation of X'X and
## X'y are used. Results then need post processing, to get R =chol(X'X)
## and f= R^{-1} X'y.
## if nt>1 and use.chol=FALSE then parallel QR is used 
{ p <- ncol(Xn)  
  y.norm2 <- y.norm2+sum(yn*yn)
  if (use.chol) { 
    if (is.null(R)) { 
      R <- crossprod(Xn)
      fn <- as.numeric(t(Xn)%*%yn) 
    } else {
      R <- R + crossprod(Xn)
      fn <- f + as.numeric(t(Xn)%*%yn)
    } 
    return(list(R=R,f=fn,y.norm2=y.norm2))
  } else { ## QR update
    if (!is.null(R)) {
      Xn <- rbind(R,Xn)
      yn <- c(f,yn)
    }
    qrx <- if (nt==1) qr(Xn,tol=0,LAPACK=TRUE) else pqr2(Xn,nt)
    fn <- qr.qty(qrx,yn)[1:min(p,nrow(Xn))]
    rp <- qrx$pivot;rp[rp] <- 1:p # reverse pivot
    return(list(R = qr.R(qrx)[,rp],f=fn,y.norm2=y.norm2))
  }
} ## qr_update


qr_up <- function(arg) {
## routine for parallel computation of the QR factorization of 
## a large gam model matrix, suitable for calling with parLapply.
  wt <- rep(0,0) 
  dev <- 0
  eta <- arg$eta
  efam <- !is.null(arg$family) ## extended family?
  for (b in 1:arg$n.block) {
    ind <- arg$start[b]:arg$stop[b]
    X <- predict(arg$G,newdata=arg$mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))
    rownames(X) <- NULL
    if (is.null(arg$coef)) eta1 <- arg$eta[ind] else eta[ind] <- eta1 <- drop(X%*%arg$coef) + arg$offset[ind]
    mu <- arg$linkinv(eta1) 
    y <- arg$G$y[ind] ## arg$G$model[[arg$response]] 
    weights <- arg$G$w[ind]
    if (efam) { ## extended family case
       dd <- dDeta(y,mu,weights,theta=arg$theta,arg$family,0)
       ## note: no handling of infinities and wz case yet
       w <- dd$EDeta2 * .5 
       #w <- w
       z <- (eta1-arg$offset[ind]) - dd$Deta.EDeta2
       good <- is.finite(z)&is.finite(w)     
    } else { ## regular exp fam case
      mu.eta.val <- arg$mu.eta(eta1)
      good <- (weights > 0) & (mu.eta.val != 0)
      z <- (eta1 - arg$offset[ind]) + (y - mu)/mu.eta.val
      w <- (weights * mu.eta.val^2)/arg$variance(mu)
    }
    w[!good] <- 0 ## drop if !good
    #z[!good] <- 0 ## irrelevant
    dev <- dev + if (efam) sum(arg$family$dev.resids(y,mu,weights,arg$theta)) else sum(arg$dev.resids(y,mu,weights))
    wt <- c(wt,w)
    z <- z[good];w <- w[good]
    w <- sqrt(w)
    ## note assumption that nt=1 in following qr_update - i.e. each cluster node is strictly serial
    if (b == 1) qrx <- qr_update(w*X[good,,drop=FALSE],w*z,use.chol=arg$use.chol) 
    else qrx <- qr_update(w*X[good,,drop=FALSE],w*z,qrx$R,qrx$f,qrx$y.norm2,use.chol=arg$use.chol)
    rm(X);if(arg$gc.level>1) gc() ## X can be large: remove and reclaim
  }
  qrx$dev <- dev;qrx$wt <- wt;qrx$eta <- eta
  if (arg$gc.level>1) { rm(arg,ind,mu,y,weights,mu.eta.val,good,z,w,wt,w);gc()}
  qrx
} ## qr_up

compress.df <- function(dat,m=NULL) {
## Takes dataframe in dat and compresses it by rounding and duplicate 
## removal. Returns index array as attribute. If covariates were matrices
## then the index is a matrix.
## For metric variables we first find the unique cases.
## If there are <= m of these then these are employed, otherwise 
## rounding is used. Factors are always reduced to the number of 
## levels present in the data. Idea is that this function is called 
## with columns of dataframes corresponding to single smooths or marginals.
## Note that this uses random sampling, so random seed manipulation
## is typically used before calling to force exact repeatability. 
  d <- ncol(dat) ## number of variables to deal with
  n <- nrow(dat) ## number of data/cases
  if (is.null(m)) {
    m <- if (d==1) 1000 else if (d==2) 100 else 25
  } else if (d>1) m <- round(m^{1/d}) + 1
  
  mf <- mm <- 1 ## total grid points for factor and metric
  for (i in 1:d) if (is.factor(dat[,i])) {  
    mf <- mf * length(unique(as.vector(dat[,i]))) 
  } else {
    mm <- mm * m 
  } 
  if (is.matrix(dat[[1]])) { ## must replace matrix terms with vec(dat[[i]])
    dat0 <- data.frame(as.vector(dat[[1]]))
    if (d>1) for (i in 2:d) dat0[[i]] <- as.vector(dat[[i]])
    names(dat0) <- names(dat)
    dat <- dat0;rm(dat0)
  }
  xu <- uniquecombs(dat,TRUE)
  if (nrow(xu)>mm*mf) { ## too many unique rows to use only unique
    for (i in 1:d) if (!is.factor(dat[,i])) { ## round the metric variables
      xl <- range(dat[,i])
      xu <- seq(xl[1],xl[2],length=m)
      dx <- xu[2]-xu[1]
      kx <- round((dat[,i]-xl[1])/dx)+1
      dat[,i] <- xu[kx] ## rounding the metric variables
    }
    xu <- uniquecombs(dat,TRUE)
  }  
  k <- attr(xu,"index")
  if (nrow(xu)==nrow(dat)) { ## might as well return original data
    k <- 1:nrow(dat)
    if (length(k)>n) k <- matrix(k,nrow=n) ## deal with matrix arguments
    attr(dat,"index") <- k 
    return(dat)
  }
  ## shuffle rows in order to avoid induced dependencies between discretized
  ## covariates (which can mess up gam.side)...
  ## Any RNG setting should be done in routine calling this one!!
  
  ii <- sample(1:nrow(xu),nrow(xu),replace=FALSE) ## shuffling index
  
  xu[ii,] <- xu  ## shuffle rows of xu
  k <- ii[k]     ## correct k index accordingly
  ## ... finished shuffle
  ## if arguments were matrices, then return matrix index
  if (length(k)>n) k <- matrix(k,nrow=n)
  storage.mode(k) <- "integer"
  k -> attr(xu,"index")
  xu
} ## compress.df

check.term <- function(term,rec) {
## utility function for discrete.mf. Checks whether variables in "term"
## have already been discretized, and if so whether this discretization 
## can be re-used for the current "term". Stops if term already discretized
## but we can't re-use discretization. Otherwise returns index of k index
## or 0 if the term is not in the existing list.
  ii <- which(rec$vnames%in%term)
  if (length(ii)) { ## at least one variable already discretized
    if (length(term)==rec$d[min(ii)]) { ## dimensions match previous discretization
      if (sum(!(term%in%rec$vnames[ii]))) stop("bam can not discretize with this nesting structure")
      else return(rec$ki[min(ii)]) ## all names match previous - return index of previous
    } else stop("bam can not discretize with this nesting structure")
  } else return(0) ## no match
} ## check.term


discrete.mf <- function(gp,mf,names.pmf,m=NULL,full=TRUE) {
## Attempt at a cleaner design, in which k.start and nr have an entry
## for each variable in mf. For jointly discretized covariates these
## will simply be identical for each covariate. 
## discretize the covariates for the terms specified in smooth.spec
## id not allowed. names.pmf gives the names of the parametric part
## of mf. If full is FALSE then the mf returned is a list where columns
## can be of different lengths.
## On exit... 
## * mf is a model frame containing the unique discretized covariate
##   values, in randomized order, padded to all be same length (if full=TRUE)
## * nr records the number of unique discretized covariate values
##   for each variable in mf i.e. the number of rows before the padding starts -
##   elements are labelled, corresponding to names in mf.
## * ks is a two column matrix. ks[i,1] is the first column in index
##   matrix k for the corresponding element of mf, ks[i,2]-1 is the 
##   final column in k. row names match names of mf.
## * k is the index matrix. The ith record of the 1st column of the 
##   jth variable is in row k[i,ks[j,1]] of the corresponding 
##   column of mf.
## ... there is an element of nr and ks for each variable of
## mf. Variables are only discretized and stored in mf once.
## Model terms can be matched to elements of nr and ks based on the
## 'xnames' returned from terms2tensor, or the 'terms' component of
## smooth objects.
## ks has an extra row, k an extra column and nr an extra entry for
## any "(Intercept)"
## Note that the margins component of a tensor product smooth spec
## is used to establish which covariates of the smooth must be
## jointly discretized.

  ## some sub sampling here... want to set and restore RNG state used for this
  ## to ensure strict repeatability.
  rngs <- temp.seed(8547) ## keep different to tps constructor!

  mf0 <- list()
  nk <- 0 ## count number of index vectors to avoid too much use of cbind
  nlp <- if (is.null(gp$nlp)) 1 else sum(unlist(lapply(gp,inherits,"split.gam.formula")))
  for (lp in 1:nlp) { ## loop over linear predictors
    smooth.spec <- if (is.null(gp$nlp)) gp$smooth.spec else gp[[lp]]$smooth.spec
    if (length(smooth.spec)>0) for (i in 1:length(smooth.spec)) nk <- nk + as.numeric(smooth.spec[[i]]$by!="NA") + length(smooth.spec[[i]]$term)
      ## if (inherits(smooth.spec[[i]],"tensor.smooth.spec")) length(smooth.spec[[i]]$margin) else 1
  }
  names.pmf <- names.pmf[names.pmf %in% names(mf)] ## drop names.pmf not in mf (usually response)
 
  nk <- nk + length(names.pmf)
  k <- matrix(0L,nrow(mf),nk) ## each column is an index vector
  ks <- matrix(NA,nk,2,dimnames=list(rep("",nk)))
  ik <- 0 ## index counter i.e. counts marginal smooths and their indices
  nr <- rep(0,nk) ## number of rows for marginal term
  ## structure to record terms already processed...
  rec <- list(vnames = rep("",0), ## variable names
              ki = rep(0,0),      ## index of original index vector var relates to  
              d = rep(0,0))       ## dimension of terms involving this var
  ## loop through the terms discretizing the covariates...
  ii <- 0
  for (lp in 1:nlp) { ## loop over linear predictors
    smooth.spec <- if (is.null(gp$nlp)) gp$smooth.spec else gp[[lp]]$smooth.spec
    if (length(smooth.spec)>0) for (i in 1:length(smooth.spec)) { ## loop over smooths in this lp
      nmarg <- if (inherits(smooth.spec[[i]],"tensor.smooth.spec")) length(smooth.spec[[i]]$margin) else 1
      maxj <- if (smooth.spec[[i]]$by=="NA") nmarg else nmarg + 1 
      ii <- ii + 1 ## all smooths counter (ii==i if nlp = 1)
      mi <- if (is.null(m)||length(m)==1) m else m[ii]
      j <- 0
      for (jj in 1:maxj) { ## loop through marginals
        if (jj==1&&maxj!=nmarg) termi <- smooth.spec[[i]]$by else {
          j <- j + 1
          termi <- if (inherits(smooth.spec[[i]],"tensor.smooth.spec")) smooth.spec[[i]]$margin[[j]]$term else 
                   smooth.spec[[i]]$term          
        } 
        ik.prev <- check.term(termi,rec) ## term already discretized?
        
        if (ik.prev==0) { ## new discretization required
	  ik <- ik + 1 ## increment index counter
          mfd <- compress.df(mf[termi],m=mi)
          ki <- attr(mfd,"index")
	  ks[ik,1] <- if (ik==1) 1 else ks[ik-1,2] 
          if (is.matrix(ki)) {
	     ks[ik,2] <- ks[ik,1] + ncol(ki)
             k <- cbind(k,matrix(0L,nrow(k),ncol(ki)-1)) ## extend index matrix
             ind <- ks[ik,1]:(ks[ik,2]-1) 
             k[,ind] <- ki 
           } else {
	     ks[ik,2] <- ks[ik,1] + 1
             k[,ks[ik,1]] <- ki
           }
	   names(nr)[ik] <- rownames(ks)[ik] <- names(mfd)[1]
           nr[ik] <- nrow(mfd)
	   if (length(termi)>1) for (ii in 2:length(termi)) { ## duplicate index info for each jointly discretized variable 
             ik <- ik+1
	     ks[ik,] <- ks[ik-1,]
	     nr[ik] <- nr[ik-1]
	     names(nr)[ik] <- rownames(ks)[ik] <- names(mfd)[ii]
           }
           mf0 <- c(mf0,mfd) 
           ## record variable discretization info, for later duplicate avoiding check...
           d <- length(termi)
           rec$vnames <- c(rec$vnames,termi)
           rec$ki <- c(rec$ki,rep(ik,d)) ### NOTE: OK, as ki contents essentially unused
           rec$d <- c(rec$d,rep(d,d))
         } else { ## re-use an earlier discretization...
           ## nothing to do in fact
         }
      } ## end marginal jj loop
    } ## term loop (i)
  } ## linear predictor, lp, loop

  ## obtain parametric terms and..
  ## pad mf0 so that all rows are the same length
  ## padding is necessary if gam.setup is to be used for setup

  if (full) {
    maxr <- max(nr)
    ## If NA's caused rows to be dropped in mf, then they should
    ## also be dropped in pmf, otherwise we can end up with factors
    ## with more levels than unique observations, for example.
    ## The next couple of lines achieve this.

    mfp <- mf[names.pmf] ## retain only variables needed in parametric part

    ## now discretize parametric covariates...
    mi <- if (is.null(m)) m else max(m)
    if (length(mfp)) for (i in 1:length(mfp)) {
      ii <- which(rec$vnames %in% names.pmf[i])
      ik.prev <- if (length(ii)>0) rec$ki[ii] else 0  ## term already discretized?
      if (ik.prev==0) { ## new discretization needed (no need to record - no repeat vars within para)
        ## note that compress.df(mfp[i]) is set up to deal with matrices in the manner appropriate
	## to the smooth summation convention, but not to matrices in the parametric part of the model
	## hence the following work around...
        if (is.matrix(mfp[[i]])) {
          mfd <- compress.df(mfp[[i]],m=mi)
	  mr <- nrow(mfd)
	  if (is.matrix(mfd)&&ncol(mfd)==1) mfd <- drop(mfd)
	  mf0[[names(mfp[i])]] <- mfd
	} else {
          mfd <- compress.df(mfp[i],m=mi);mf0 <- c(mf0,mfd)
	  mr <- length(mfd[[1]])
        }
        ki <- attr(mfd,"index")
        ik <- ik + 1
        ks[ik,1] <- if (ik==1) 1 else ks[ik-1,2]
	ks[ik,2] <- ks[ik,1] + 1
        k[,ks[ik,1]] <- ki
        if (maxr<mr) maxr <- mr
        nr[ik] <- mr
	names(nr)[ik] <- rownames(ks)[ik] <- names(mfp[i])
      } else { ## re-use previous discretization
        ## nothing needed
      }
    }

    for (i in 1:length(mf0)) { ## pad out columns to be same length
      if (is.matrix(mf0[[i]])) {
        me <- nrow(mf0[[i]])
	mf0[[i]] <- rbind(mf0[[i]],matrix(sample(mf0[[i]],(maxr-me)*ncol(mf0[[i]]),replace=TRUE),maxr-me,ncol(mf0[[i]])))
      } else {
        me <- length(mf0[[i]]) 
        if (me < maxr) mf0[[i]][(me+1):maxr] <- sample(mf0[[i]],maxr-me,replace=TRUE)
      }	
    }
    
    ## mf0 is the discretized model frame (actually a list), padded to have equal length rows
    ## now copy back into mf so terms unchanged
    mf <- mf[sample(1:nrow(mf),maxr,replace=TRUE),,drop=FALSE]
    for (na in names(mf0)) mf[[na]] <- mf0[[na]] 
  } else mf <- mf0
  nr <- nr[names(mf)] ## same order for both mf and nr

  ## reset RNG to old state...
  temp.seed(rngs)

  ## k can end up with too many columns if marginals themselves have
  ## multiple arguments, so drop these here...
  ks <- ks[!is.na(ks[,1]),,drop=FALSE]
  if (ncol(k)>max(ks)-1) k <- k[,1:(max(ks)-1),drop=FALSE] 
  k <- cbind(k,1L) ## add an intercept index column
  nk <- ncol(k)
  ks <- rbind(ks,c(nk,nk+1)) ## add intercept row to ks
  nr <- c(nr,1) ## and intercept entry to nr
  names(nr)[length(nr)] <- rownames(ks)[nrow(ks)] <- "(Intercept)"
  list(mf=mf,k=k,nr=nr,ks=ks)
} ## discrete.mf

mini.mf <-function(mf,chunk.size) {
## takes a model frame and produces a representative subset of it, suitable for 
## basis setup.
  ## first count the minimum number of rows required for representiveness
  mn <- 0
  for (j in 1:length(mf)) mn <- mn + if (is.factor(mf[[j]])) length(levels(mf[[j]])) else 2
  if (chunk.size < mn) chunk.size <- mn
  n <- nrow(mf)
  if (n <= chunk.size) return(mf)
  rngs <- temp.seed(66)
  ## randomly sample from original frame...
  ind <- sample(1:n,chunk.size)
  mf0 <- mf[ind,,drop=FALSE]
  ## ... now need to ensure certain sorts of representativeness

  ## work through elements collecting the rows containing 
  ## max and min for each variable, and a random row for each 
  ## factor level....

  ind <- sample(1:n,n,replace=FALSE) ## randomized index for stratified sampling w.r.t. factor levels
  fun <- function(X,fac,ind) ind[fac[ind]==X][1] ## stratified sampler
  k <- 0 
  for (j in 1:length(mf)) if (is.numeric(mf0[[j]])) {
    if (is.matrix(mf0[[j]])) { ## find row containing minimum
      j.min <- min((1:n)[as.logical(rowSums(mf[[j]]==min(mf[[j]])))])
      j.max <- min((1:n)[as.logical(rowSums(mf[[j]]==max(mf[[j]])))])
    } else { ## vector
      j.min <- min(which(mf[[j]]==min(mf[[j]])))
      j.max <- min(which(mf[[j]]==max(mf[[j]])))
    }
    k <- k + 1; mf0[k,] <- mf[j.min,]
    k <- k + 1; mf0[k,] <- mf[j.max,] 
  } else if (is.factor(mf[[j]])) { ## factor variable...
    ## randomly sample one row from each factor level...
    find <- apply(X=as.matrix(levels(mf[[j]])),MARGIN=1,FUN=fun,fac=mf[[j]],ind=ind)
    find <- find[is.finite(find)] ## in case drop.unused.levels==FALSE, so that there ar levels without rows
    nf <- length(find)
    mf0[(k+1):(k+nf),] <- mf[find,]
    k <- k + nf
  }
  temp.seed(rngs) ## reset RNG to initial state

  mf0
} ## mini.mf


bgam.fitd <- function (G, mf, gp ,scale , coef=NULL, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs),rho=0, control = gam.control(), intercept = TRUE, 
    gc.level=0,nobs.extra=0,npt=c(1,1),gamma=1,in.out=NULL,method="fREML",nei=NULL,...) {
## This is a version of bgam.fit designed for use with discretized covariates. 
## Difference to bgam.fit is that XWX, XWy and Xbeta are computed in C
## code using compressed versions of X. Parallelization of XWX formation
## is performed at the C level using openMP.
## Alternative fitting iteration using Cholesky only, including for REML.
## Basic idea is to take only one Newton step for parameters per iteration
## and to control the step length to ensure that at the end of the step we
## are not going uphill w.r.t. the REML criterion...
    
  y <- G$y
  weights <- G$w 
  conv <- FALSE
  nobs <- nrow(mf)
  offset <- G$offset 

  if (method=="NCV") {
    if (rho!=0) {
      warning("rho ignored with NCV"); rho <- 0
    }
    nei$af <- nei$a; nei$maf <- nei$ma ## for use in cov matrix computation
    if (!is.null(nei$sample)) { ## only a sub-sample of points used in optimization NCV
      if (length(nei$sample)==1) nei$sample <- sample(nei$ma,nei$sample)
      m1 <- nei$ma[nei$sample]; m0 <- c(0,nei$ma)[nei$sample]
      nei$a <- nei$a[sequence(m1-m0,m0+1L,1L)]
      nei$ma <- cumsum(m1-m0)-1L
      m1 <- nei$md[nei$sample]; m0 <- c(0,nei$md)[nei$sample]
      nei$d <- nei$d[sequence(m1-m0,m0+1L,1L)]
      nei$md <- cumsum(m1-m0)-1L
    }
  }
  
  if (inherits(G$family,"extended.family")) { ## preinitialize extended family
    efam <- TRUE
    if (!is.null(G$family$preinitialize) && !is.null(attr(G$family$preinitialize,"needG"))) attr(G$family,"G") <- G
    pini <- if (is.null(G$family$preinitialize)) NULL else G$family$preinitialize(y,G$family)
    if (is.null(G$family$preinitialize)) attr(G$family,"G") <- NULL
    if (!is.null(pini$family)) G$family <- pini$family
    if (!is.null(pini$Theta)) G$family$putTheta(pini$Theta)
    if (!is.null(pini$y)) y <- pini$y
    if (is.null(G$family$scale)) scale <- 1 else scale <- if (G$family$scale<0) scale else G$family$scale
    scale1 <- scale
    if (scale < 0) scale <- var(y) *.1 ## initial guess
  } else efam <- FALSE

  if (rho!=0) { ## AR1 error model
      
    ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
    sd <- -rho*ld         ## sub diagonal
    N <- nobs    
    ## see rwMatrix() for how following are used...
    ar.row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)]) ## index of rows to reweight
    ar.weight <- c(1,rep(c(sd,ld),N-1))     ## row weights
    ar.stop <- c(1,1:(N-1)*2+1)    ## (stop[i-1]+1):stop[i] are the rows to reweight to get ith row
    if (!is.null(mf$"(AR.start)")) { ## need to correct the start of new AR sections...
      ii <- which(mf$"(AR.start)"==TRUE)
      if (length(ii)>0) {
        if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
        ar.weight[ii*2-2] <- 0 ## zero sub diagonal
        ar.weight[ii*2-1] <- 1 ## set leading diagonal to 1
      }
    }
  } else { ## AR setup complete
    ar.row <- ar.weight <- ar.stop <- -1 ## signal no re-weighting
  }

  family <- G$family
  additive <- if (family$family=="gaussian"&&family$link=="identity") TRUE else FALSE
  linkinv <- family$linkinv;#dev.resids <- family$dev.resids
  if (!efam) {
    variance <- family$variance
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
  }
  valideta <- family$valideta
  if (is.null(valideta))
      valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu))
      validmu <- function(mu) TRUE
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }

  if (is.matrix(y)&&ncol(y)>1) stop("This family should not have a matrix response")
  eta <- if (!is.null(etastart)) etastart else family$linkfun(mustart)
    
  mu <- linkinv(eta)
  if (!(validmu(mu) && valideta(eta)))
    stop("cannot find valid starting values: please specify some")
  dev <- sum(family$dev.resids(y, mu, weights))*2 ## just to avoid converging at iter 1

  conv <- FALSE
  
  G$coefficients <- rep(0,ncol(G$X))
  class(G) <- "gam"  
    
  ## need to reset response and weights to post initialization values
  ## in particular to deal with binomial properly...
  G$y <- y
  G$w <- weights
  ## need to keep untransformed S matrices with NCV...
  Sl <- Sl.setup(G,keepS=(method=="NCV")) ## setup block diagonal penalty object
  rank <- 0
  if (length(Sl)>0) for (b in 1:length(Sl)) rank <- rank + Sl[[b]]$rank
  Mp <- ncol(G$X) - rank ## null space dimension
  Nstep <- 0
  if (efam) theta <- family$getTheta()

  if(!is.null(coef)) { ## use supplied coef if supplied, including for step halving (Matteo Fasiolo)
    # Define all quantities needed for step-halving
    coef0 <- coef
    b0 <- Sl.initial.repara(Sl, coef0, inverse = FALSE, both.sides = FALSE, cov = FALSE)
    eta0 <- Xbd(G$Xd,coef0,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) + offset
    mu0 <- linkinv(eta0)
    dev0 <- if (efam) sum(family$dev.resids(G$y,mu0,G$w,theta)) else sum(family$dev.resids(G$y,mu0,G$w))
    dev <- dev0 * 2  # just to avoid converging at iter 1 (see above)
    # Note: possibly over-writing arguments etastart and mustart, and corresponding dev, but we must be coherent with coef0
    etastart <- eta <- eta0
    mustart <- mu <- mu0
    c.iter <- 1 # Used later on to indicate whether to perform step-halving after 1st or...
  } else c.iter <- 2 # ... 2nd iteration 


  for (iter in 1L:control$maxit) { ## main fitting loop 
    devold <- dev
    dev <- 0
     
    if (iter==1||!additive) {
      qrx <- list()

      if (iter>1) {
        ## form eta = X%*%beta
        eta <- Xbd(G$Xd,coef,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) + offset
        lsp.full <- G$lsp0
	if (n.sp>0) lsp.full <- lsp.full + if (is.null(G$L)) lsp[1:n.sp] else G$L %*% lsp[1:n.sp]
	rSb <- Sl.rSb(Sl,lsp.full,prop$beta) ## store S beta to allow rapid step halving
	if (iter>c.iter) {
	  rSb0 <- Sl.rSb(Sl,lsp.full,b0)
	  bSb0 <- sum(rSb0^2) ## penalty at start of beta step
	  ## get deviance at step start, with current theta if efam
	  dev0 <- if (efam) sum(family$dev.resids(G$y,mu0,G$w,theta)) else
	                    sum(family$dev.resids(G$y,mu0,G$w))
	    
        }
      }
      kk <- 1
      repeat { ## step control
        mu <- linkinv(eta)
	dev <- if (efam) sum(family$dev.resids(G$y,mu,G$w,theta)) else
	                 sum(family$dev.resids(G$y,mu,G$w))
        if (iter>c.iter) { ## coef step length control
	  bSb <- sum(rSb^2) ## penalty at end of beta step
	  #stepr <- max(abs(coef-coef0))/(max(abs(coef0)))
	  #if (stepr<2)
	  stepr <- 2
	  istepr <- 1 - 1/stepr
          if ((!is.finite(dev) || dev0 + bSb0 < dev + bSb) && kk < 30) { ## beta step not improving current pen dev
            coef <- coef0*istepr + coef/stepr ## reduce the step
	    rSb <- rSb0*istepr + rSb/stepr
	    eta <- eta0*istepr + eta/stepr
	    prop$beta <- b0*istepr + prop$beta/stepr
	    kk <- kk + 1
          } else break
        } else break
      }	## step control	 

      if (iter>1) { ## save components of penalized deviance for step control
        coef0 <- coef ## original para
        eta0 <- eta
	mu0 <- mu
	b0 <- prop$beta ## beta repara
	dev <- dev + sum(rSb^2) ## add penalty to deviance
      } else crit <- dev ## for convergence checking
	
      if (efam) { ## extended family
	if (iter>1) { ## estimate theta
          if (family$n.theta>0||scale1<0) theta <- estimate.theta(theta,family,y,mu,scale=scale1,wt=G$w,tol=1e-7)
          if (!is.null(family$scale) && scale1<0) {
	    scale <- exp(theta[family$n.theta+1])
	    theta <- theta[1:family$n.theta]
	  }  
          family$putTheta(theta)
        } ## theta estimation
	  
        dd <- dDeta(y,mu,G$w,theta=theta,family,0)
	## note: no handling of infinities and wz case yet

        if (rho==0) {
          w <- dd$Deta2 * .5 
          z <- (eta-offset) - dd$Deta.Deta2
        } else { ## use fisher weights
	  w <- dd$EDeta2 * .5 
          z <- (eta-offset) - dd$Deta.EDeta2
	}
        good <- is.finite(z)&is.finite(w)
	w[!good] <- 0 ## drop if !good
	z[!good] <- 0 ## irrelevant
      } else { ## exponential family
        mu.eta.val <- mu.eta(eta)
        good <- mu.eta.val != 0
        mu.eta.val[!good] <- .1 ## irrelvant as weight is zero
        z <- (eta - offset) + (G$y - mu)/mu.eta.val
        w <- (G$w * mu.eta.val^2)/variance(mu)
      }
      
  
      qrx$y.norm2 <- if (rho==0) sum(w*z^2) else   ## AR mod needed
          sum(rwMatrix(ar.stop,ar.row,ar.weight,sqrt(w)*z,trans=FALSE)^2) 
       
      ## form X'WX efficiently...
      qrx$R <- XWXd(G$Xd,w,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,npt[1],G$drop,ar.stop,ar.row,ar.weight)
      ## form X'Wz efficiently...
      qrx$f <- XWyd(G$Xd,w,z,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop,ar.stop,ar.row,ar.weight)
      if (gc.level>1) gc()
     
      ## following reparameterizes X'X and f=X'y, according to initial reparameterizarion...
     
      qrx$XX <- Sl.initial.repara(Sl,qrx$R,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=npt[1])
      qrx$Xy <- Sl.initial.repara(Sl,qrx$f,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=npt[1])  
       
      G$n <- nobs
    } else {  ## end of if (iter==1||!additive)
      dev <- qrx$y.norm2 - sum(coef*qrx$f) ## actually penalized deviance
    }
  
    if (control$trace)
         message(gettextf("Deviance = %s Iterations - %d", dev, iter, domain = "R-mgcv"))

    if (!is.finite(dev)) stop("Non-finite deviance")

    ## preparation for working model fit is ready, but need to test for convergence first
    if (iter>2 && abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon && (scale>0 || method=="NCV" || abs(Nstep[n.sp+1])<control$epsilon*(abs(log.phi)+1))) {
      conv <- TRUE
      break
    }

    ## use fast REML code
    ## block diagonal penalty object, Sl, set up before loop

    if (iter==1) { ## need to get initial smoothing parameters 
      lambda.0 <- if (is.null(in.out)) initial.sp(qrx$R,G$S,G$off,XX=TRUE) else in.out$sp ## note that this uses the untransformed X'X in qrx$R
      ## convert intial s.p.s to account for L 
      lsp0 <- log(lambda.0) ## initial s.p.
      if (!is.null(G$L)) lsp0 <- 
        if (ncol(G$L)>0) as.numeric(coef(lm(lsp0 ~ G$L-1+offset(G$lsp0)))) else rep(0,0)
      n.sp <- length(lsp0) 
    }
     
    ## carry forward scale estimate if possible...
    if (scale>0) log.phi <- log(scale) else {
      if (iter==1&&method!="NCV") {
        if (is.null(in.out)) {
          if (is.null(coef)||qrx$y.norm2==0) lsp0[n.sp+1] <- log(var(as.numeric(G$y))*.05) else
              lsp0[n.sp+1] <- log(qrx$y.norm2/(nobs+nobs.extra))
        } else lsp0[n.sp+1] <- log(in.out$scale)		 
      }
    }

    ## get beta, grad and proposed Newton step... 
    repeat { ## Take a Newton step to update log sp and phi
      lsp <- lsp0 + Nstep ## Note that Nstep is 0 at iter==1
        
      if (method=="NCV") {
        ## NOTE: scale param, tol based on 'reml',  nthreads!!
	## basically only here for testing at the moment - not fully plumbed in
	#prop2 <- Sl.fitChol(Sl,qrx$XX,qrx$Xy,rho=lsp[1:n.sp],yy=qrx$y.norm2,L=G$L,rho0=G$lsp0,log.phi=log.phi,
        #        phi.fixed=scale>0,nobs=nobs,Mp=Mp,nt=npt,tol=abs(crit)*.Machine$double.eps^.5,gamma=gamma)
        prop <- Sl.ncv(z,G$Xd,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,nei=nei,Sl=Sl,XX=qrx$XX,w=w,f=qrx$Xy,rho=lsp[1:n.sp],
	               nt=npt,L=G$L,rho0=G$lsp0,drop=G$drop,tol=abs(crit)*.Machine$double.eps^.5,nthreads=npt[1])
	if (scale<=0) log.phi = log(sum((z-eta)^2*w)/nobs) ## not really needed		 
        deriv.check <- FALSE
        if (deriv.check) { ## for debug derivative testing
          Hd <- prop$hess; eps <- 1e-7 ;fd <- prop$grad
	  bfd <- prop$db
          for (i in 1:n.sp) {
            lspfd <- lsp; lspfd[i] <- lsp[i] + eps
	    prop1 <- Sl.ncv(z,G$Xd,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,nei=nei,Sl=Sl,XX=qrx$XX,w=w,f=qrx$Xy,rho=lspfd[1:n.sp],
	                    nt=npt,L=G$L,rho0=G$lsp0,drop=G$drop,tol=abs(crit)*.Machine$double.eps^.5,nthreads=1)
	    fd[i] <- (prop1$NCV - prop$NCV)/eps
	    Hd[,i] <- (prop1$grad - prop$grad)/eps
	    bfd[,i] <- (prop1$beta - prop$beta)/eps
          }
	  Hd <- 0.5 * (Hd + t(Hd))
        }
      } else { ## method is REML
        if (scale<=0) log.phi <- lsp[n.sp+1]
        prop <- Sl.fitChol(Sl,qrx$XX,qrx$Xy,rho=lsp[1:n.sp],yy=qrx$y.norm2,L=G$L,rho0=G$lsp0,log.phi=log.phi,
                phi.fixed=scale>0,nobs=nobs,Mp=Mp,nt=npt,tol=abs(crit)*.Machine$double.eps^.5,gamma=gamma)
        deriv.check <- FALSE
        if (deriv.check) { ## for debug derivative testing
          Hd <- prop$hess; eps <- 1e-7
          bfd <- prop$db
          for (i in 1:n.sp) {
            lspfd <- lsp; lspfd[i] <- lsp[i] + eps
            prop1 <- Sl.fitChol(Sl,qrx$XX,qrx$Xy,rho=lspfd[1:n.sp],yy=qrx$y.norm2,L=G$L,rho0=G$lsp0,log.phi=log.phi,
                  phi.fixed=scale>0,nobs=nobs,Mp=Mp,nt=npt,tol=abs(crit)*.Machine$double.eps^.5,gamma=gamma)
            Hd[,i] <- (prop1$grad - prop$grad)/eps
            bfd[,i] <- (prop1$beta - prop$beta)/eps
          }
          Hd <- 0.5 * (Hd + t(Hd))
        }	  
      }
      if (max(Nstep)==0) { 
        Nstep <- prop$step;lsp0 <- lsp;
        break 
      } else { ## step length control
        if (sum(prop$grad*Nstep)>dev*1e-7) Nstep <- Nstep/2 else {
          Nstep <- prop$step;lsp0 <- lsp;break;
        }
      }
    } ## end of sp update

    coef <- Sl.initial.repara(Sl,prop$beta,inverse=TRUE,both.sides=FALSE,cov=FALSE)

    if (any(!is.finite(coef))) {
      conv <- FALSE
      warning(gettextf("non-finite coefficients at iteration %d",
              iter))
      break
    }
    crit <- if (method=="NCV") prop$NCV else (dev/(exp(log.phi)*gamma) - prop$ldetS + prop$ldetXXS)/2
  } ## end fitting iteration

  if (!conv) warning("algorithm did not converge")
   
  eps <- 10 * .Machine$double.eps
  if (family$family == "binomial") {
    if (any(mu > 1 - eps) || any(mu < eps))
        warning("fitted probabilities numerically 0 or 1 occurred")
  }
  if (family$family == "poisson") {
    if (any(mu < eps))
      warning("fitted rates numerically 0 occurred")
  }

  object <- list(mgcv.conv=conv,rank=prop$r,
                 scale.estimated = scale<=0,
		 outer.info=NULL, optimizer=c("perf","chol")) 
  Mp <- G$nsdf
  if (length(G$smooth)>1) for (i in 1:length(G$smooth)) Mp <- Mp + G$smooth[[i]]$null.space.dim
  scale <- exp(log.phi)
  crit <- if (method=="NCV") prop$NCV else (dev/(scale*gamma) - prop$ldetS + prop$ldetXXS +
                                     (length(y)/gamma-Mp)*log(2*pi*scale)+Mp*log(gamma))/2
  if (rho!=0) { ## correct REML score for AR1 transform
    df <- if (is.null(mf$"(AR.start)")) 1 else sum(mf$"(AR.start)")
    crit <- crit - (nobs/gamma-df)*log(ld)
  }

  if (ncol(prop$db)>0) for (i in 1:ncol(prop$db)) prop$db[,i] <- ## d beta / d rho matrix
        Sl.initial.repara(Sl,as.numeric(prop$db[,i]),inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=npt[1]) 

  object$db.drho <- prop$db
  object$gcv.ubre <- crit
      
  object$coefficients <- coef
  object$family <- family
  ## form linear predictor efficiently...
  object$linear.predictors <- Xbd(G$Xd,coef,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) + G$offset
  object$fitted.values <- family$linkinv(object$linear.predictors)
  if (efam) { ## deal with any post processing
     if (!is.null(family$postproc)) {
      posr <- family$postproc(family=object$family,y=G$y,prior.weights=G$w,
              fitted=object$fitted.values,linear.predictors=object$linear.predictors,offset=G$offset,
	      intercept=G$intercept)
      if (!is.null(posr$family)) object$family$family <- posr$family
      if (!is.null(posr$deviance)) object$deviance <- posr$deviance
      if (!is.null(posr$null.deviance)) object$null.deviance <- posr$null.deviance
    }
    if (is.null(object$null.deviance)) object$null.deviance <- sum(family$dev.resids(G$y,weighted.mean(G$y,G$w),G$w,theta))   
  }
 
  PP <- Sl.initial.repara(Sl,prop$PP,inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=npt[1])
  F <- pmmult(PP,qrx$R,FALSE,FALSE,nt=npt[1])  ##crossprod(PP,qrx$R) - qrx$R contains X'WX in this case
 
  object$edf <- diag(F)
  object$edf1 <- 2*object$edf - rowSums(t(F)*F)
  lsp <- if (n.sp>0) lsp[1:n.sp] else rep(0,0)
  object$sp <- exp(lsp)
  object$full.sp <- if (is.null(G$L)) object$sp else exp(drop(G$L%*%lsp + G$lsp0))
  if (object$scale.estimated&&method=="NCV") {
    rss <- if (additive) sum(w*(y-object$fitted.values)^2) else  sum((z-eta)^2*w)
    scale <- rss/(nobs-sum(object$edf))
  }
  object$sig2 <- object$scale <- scale
  object$Vp <- PP * scale
  object$Ve <- pmmult(F,object$Vp,FALSE,FALSE,nt=npt[1]) ## F%*%object$Vp

  if (method=="NCV"&& max(diff(nei$ma))>1 && max(diff(nei$md))==1) { ## replace Vp with NCV estimate 
    e <- if (additive) w*(y-object$fitted.values) else (z-eta)*w ## weighted residuals
    XVX <- XVXd(G$Xd,e,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,nthreads=1,a=nei$af,ma=nei$maf)
    ## debug check of XVX
    debug <- FALSE
    if (debug) {
      X <- Xbd(G$Xd,diag(1,nrow=nrow(XVX)),G$kd,G$ks,G$ts,G$dt,G$v,G$qc,drop=NULL,lt=NULL)
      Ve <- matrix(0,nobs,nobs)
      j0 <- 1;af <- nei$af + 1; maf <- nei$maf + 1
      for (i in 1:nobs) {
        for (j in af[j0:maf[i]]) Ve[i,j] <- e[i]*e[j]
        j0 <- maf[i]+1
      }
      XVX1 <- t(X)%*%Ve%*%X
    }  
    ## debug end
    ## note: prop$rsd contains weighted cross validated residuals
    inflate <- max(1,1 + .6*(mad(prop$rsd)/mad(e)-1))
    #inflate <- max(1,(prop$NCV/sum(e[nei$d]^2/w[nei$d])-1)/2+1)
    V1 <- PP %*% XVX %*% PP *inflate^2
    object$Vp <- sum(diag(V1))/sum(diag(object$Ve))*(object$Vp-object$Ve) + V1
    ## NOTE: generally we will not have cross validated residuals for all data points.
    ##       This is because we usually base NCV on a sample for large data. But average
    ##       CV inflation factor is computable.
  }
  
  ## sp uncertainty correction... 
  if (!is.null(G$L)) prop$db <- prop$db%*%G$L
  M <- ncol(prop$db) 
  if (M>0&&method!="NCV") {
    ev <- eigen(prop$hess,symmetric=TRUE)
    ind <- ev$values <= 0
    ev$values[ind] <- 0;ev$values[!ind] <- 1/sqrt(ev$values[!ind])
    rV <- (ev$values*t(ev$vectors))[,1:M]
    V.sp <- crossprod(rV)
    attr(V.sp,"L") <- G$L
    Vc <- pcrossprod(rV%*%t(prop$db),nt=npt[1])
  } else {
    Vc <- 0
    V.sp <- NULL
  }  
  Vc <- object$Vp + Vc  ## Bayesian cov matrix with sp uncertainty
  object$edf2 <- rowSums(Vc*qrx$R)/scale
  object$Vc <- Vc
  object$V.sp <- V.sp
  object$outer.info <- list(grad = prop$grad,hess=prop$hess)  
  object$AR1.rho <- rho
  object$R <- if (npt[2]>1) pchol(qrx$R,npt) else suppressWarnings(chol(qrx$R,pivot=TRUE)) ## latter much faster under optimized BLAS
  piv <- attr(object$R,"pivot") 
  object$R[,piv] <- object$R   
  object$iter <- iter 
  object$wt <- w
  object$y <- G$y
  object$prior.weights <- G$w
  rm(G);if (gc.level>0) gc()
  object
} ## end bgam.fitd 


regular.Sb <- function(S,off,sp,beta) {
## form S %*% beta for a normal G list
  a <- beta*0
  if (length(S)>0) for (i in 1:length(S)) {
    ind <- off[i] - 1 + 1:ncol(S[[i]])
    a[ind] <- a[ind] + sp[i] * S[[i]] %*% beta[ind]
  }
  a
} ## regular.Sb


bgam.fit <- function (G, mf, chunk.size, gp ,scale ,gamma,method, coef=NULL,etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), control = gam.control(), intercept = TRUE, 
    cl = NULL,gc.level=0,use.chol=FALSE,nobs.extra=0,samfrac=1,npt=1,in.out=NULL,...) {
    y <- G$y
    weights <- G$w
    conv <- FALSE
    nobs <- nrow(mf)
    offset <- G$offset
    family <- G$family
    ## extended family may have non-standard y that requires careful subsetting (e.g. cnorm)
    subsety <- if (is.null(G$family$subsety)) function(y,ind) y[ind] else G$family$subsety
    if (inherits(G$family,"extended.family")) { ## preinitialize extended family
      efam <- TRUE
      pini <- if (is.null(G$family$preinitialize)) NULL else G$family$preinitialize(y,G$family)
      if (!is.null(pini$Theta)) G$family$putTheta(pini$Theta)
      if (!is.null(pini$y)) y <- pini$y
      if (is.null(G$family$scale)) scale <- 1 else scale <- if (G$family$scale<0) scale else G$family$scale
      scale1 <-scale
      if (scale < 0) scale <- var(y) *.1 ## initial guess
    } else efam <- FALSE

 
    G$family <- gaussian() ## needed if REML/ML used
    G$family$drop.intercept <- family$drop.intercept ## needed in predict.gam
    linkinv <- family$linkinv
    if (!efam) {
      variance <- family$variance
      mu.eta <- family$mu.eta
      if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
    }
    dev.resids <- family$dev.resids
    ## aic <- family$aic
   
    valideta <- family$valideta
    if (is.null(valideta))
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu))
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }

    if (is.matrix(y)&&ncol(y)>1) stop("This family should not have a matrix response")

    ##coefold <- NULL
    eta <- if (!is.null(etastart))
         etastart
    else family$linkfun(mustart)
    
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
       stop("cannot find valid starting values: please specify some")
    dev <- sum(dev.resids(y, mu, weights))*2 ## just to avoid converging at iter 1
    conv <- FALSE
   
    G$coefficients <- rep(0,ncol(G$X))
    class(G) <- "gam"  
    
    ## need to reset response and weights to post initialization values
    ## in particular to deal with binomial properly...
    G$y <- y
    G$w <- weights

  
    ## set up cluster for parallel computation...

    if (!is.null(cl)&&inherits(cl,"cluster")) {
      n.threads <- length(cl)
      while(nobs/n.threads < ncol(G$X)) n.threads <- n.threads - 1
      if (n.threads < length(cl)) { 
        warning("Too many cluster nodes to use all efficiently")
      }
    } else n.threads <- 1

    if (n.threads>1) { ## set up thread argument lists
      ## number of obs per thread
      nt <- rep(ceiling(nobs/n.threads),n.threads)
      nt[n.threads] <- nobs - sum(nt[-n.threads])
      arg <- list()
      n1 <- 0
      for (i in 1:n.threads) {
        n0 <- n1+1;n1 <- n1+nt[i]
        ind <- n0:n1 ## this thread's data block from mf
        n.block <- nt[i]%/%chunk.size ## number of full sized blocks
        stub <- nt[i]%%chunk.size ## size of end block
        if (n.block>0) {
          start <- (0:(n.block-1))*chunk.size+1
          stop <- (1:n.block)*chunk.size
          if (stub>0) {
            start[n.block+1] <- stop[n.block]+1
            stop[n.block+1] <- nt[i]
            n.block <- n.block+1
          } 
        } else {
          n.block <- 1
          start <- 1
          stop <- nt[i]
        }
        arg[[i]] <- list(nobs= nt[i],start=start,stop=stop,n.block=n.block,
                         linkinv=linkinv,dev.resids=dev.resids,gc.level=gc.level,
                         mf = mf[ind,],
                         eta = eta[ind],offset = offset[ind],G = G,use.chol=use.chol)
        if (efam) {
          arg[[i]]$family <- family
        } else {
          arg[[i]]$mu.eta <- mu.eta
	  arg[[i]]$variance <- variance
        }
        arg[[i]]$G$w <- G$w[ind];arg[[i]]$G$model <- NULL
        arg[[i]]$G$y <- subsety(G$y,ind) ##G$y[ind]
      }
    } else { ## single thread, requires single indices
      ## construct indices for splitting up model matrix construction... 
      n.block <- nobs%/%chunk.size ## number of full sized blocks
      stub <- nobs%%chunk.size ## size of end block
      if (n.block>0) {
        start <- (0:(n.block-1))*chunk.size+1
        stop <- (1:n.block)*chunk.size
        if (stub>0) {
          start[n.block+1] <- stop[n.block]+1
          stop[n.block+1] <- nobs
          n.block <- n.block+1
        } 
      } else {
        n.block <- 1
        start <- 1
        stop <- nobs
      }
   } ## single thread indices complete
 
    conv <- FALSE

    if (method=="fREML") Sl <- Sl.setup(G) ## setup block diagonal penalty object

    if (efam) theta <- family$getTheta()

    for (iter in 1L:control$maxit) { ## main fitting loop
       ## accumulate the QR decomposition of the weighted model matrix
       devold <- dev
       kk <- 0
       repeat { 
         dev <- 0;wt <- rep(0,0) 
         if (n.threads == 1) { ## use original serial update code
           wt <- G$y
           for (b in 1:n.block) {
             ind <- start[b]:stop[b]
	     if (!is.null(family$setInd)) family$setInd(ind) ## family may need to know ind - e.g. gfam
             X <- predict(G,newdata=mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))
             rownames(X) <- NULL
             if (is.null(coef)) eta1 <- eta[ind] else eta[ind] <- eta1 <- drop(X%*%coef) + offset[ind]
             mu <- linkinv(eta1) 
             y <- subsety(G$y,ind) #G$y[ind] ## G$model[[gp$response]] ## - G$offset[ind]
             weights <- G$w[ind]
             if (efam) { ## extended family case
                dd <- dDeta(y,mu,weights,theta=theta,family,0)
	        ## note: no handling of infinities and wz case yet
               
	        w <- dd$EDeta2 * .5 
                z <- (eta1-offset[ind]) - dd$Deta.EDeta2
	        good <- is.finite(z)&is.finite(w)
	       
             } else { ## regular exp fam case
               mu.eta.val <- mu.eta(eta1)
               good <- (weights > 0) & (mu.eta.val != 0)
               z <- (eta1 - offset[ind]) + (y - mu)/mu.eta.val
               w <- (weights * mu.eta.val^2)/variance(mu)
             }
             dev <- dev + if (efam) sum(dev.resids(y,mu,weights,theta)) else sum(dev.resids(y,mu,weights))
	     w[!good] <- 0 ## drop if !good
	     #z[!good] <- 0 ## irrelevant 
             wt[ind] <- w ## wt <- c(wt,w)
             w <- w[good];z <- z[good]
             w <- sqrt(w)
             ## note that QR may be parallel using npt>1, even under serial accumulation...
             if (b == 1) qrx <- qr_update(w*X[good,,drop=FALSE],w*z,use.chol=use.chol,nt=npt) 
             else qrx <- qr_update(w*X[good,,drop=FALSE],w*z,qrx$R,qrx$f,qrx$y.norm2,use.chol=use.chol,nt=npt)
             rm(X);if(gc.level>1) gc() ## X can be large: remove and reclaim
          }
          if (use.chol) { ## post proc to get R and f...
            y.norm2 <- qrx$y.norm2 
            qrx <- chol2qr(qrx$R,qrx$f,nt=npt)
            qrx$y.norm2 <- y.norm2
          }
        } else { ## use parallel accumulation
	  
          for (i in 1:length(arg)) {
	    arg[[i]]$coef <- coef
	    if (efam) arg[[i]]$theta <- theta
	  }
          res <- parallel::parLapply(cl,arg,qr_up) 
          ## now consolidate the results from the parallel threads...
          if (use.chol) {
            R <- res[[1]]$R;f <- res[[1]]$f;dev <- res[[1]]$dev
            wt <- res[[1]]$wt;y.norm2 <- res[[1]]$y.norm2
	    eta <- res[[1]]$eta
            for (i in 2:n.threads) {
              R <- R + res[[i]]$R; f <- f + res[[i]]$f
              wt <- c(wt,res[[i]]$wt);eta <- c(eta,res[[i]]$eta);
	      dev <- dev + res[[i]]$dev
              y.norm2 <- y.norm2 + res[[i]]$y.norm2
            }         
            qrx <- chol2qr(R,f,nt=npt)
            qrx$y.norm2 <- y.norm2
          } else { ## proper QR
            R <- res[[1]]$R;f <- res[[1]]$f;dev <- res[[1]]$dev
            wt <- res[[1]]$wt;y.norm2 <- res[[1]]$y.norm2; eta <- res[[1]]$eta
            for (i in 2:n.threads) {
              R <- rbind(R,res[[i]]$R); f <- c(f,res[[i]]$f)
              wt <- c(wt,res[[i]]$wt);eta <- c(eta,res[[i]]$eta)
	      dev <- dev + res[[i]]$dev
              y.norm2 <- y.norm2 + res[[i]]$y.norm2
            }         
            ## use parallel QR here if npt>1...
            qrx <- if (npt>1) pqr2(R,npt) else qr(R,tol=0,LAPACK=TRUE) 
            f <- qr.qty(qrx,f)[1:ncol(R)]
            rp <- qrx$pivot;rp[rp] <- 1:ncol(R) # reverse pivot
            qrx <- list(R=qr.R(qrx)[,rp],f=f,y.norm2=y.norm2)
          }
        } 

        ## if the routine has been called with only a random sample of the data, then 
        ## R, f and ||y||^2 can be corrected to estimate the full versions...
 
        qrx$R <- qrx$R/sqrt(samfrac)
        qrx$f <- qrx$f/sqrt(samfrac)
        qrx$y.norm2 <- qrx$y.norm2/samfrac

        G$n <- nobs
      
        rss.extra <- qrx$y.norm2 - sum(qrx$f^2)
      
        if (control$trace)
           message(gettextf("Deviance = %s Iterations - %d", dev, iter, domain = "R-mgcv"))

        if (!is.finite(dev)) stop("Non-finite deviance")

        ## preparation for working model fit is ready, but need to test for convergence first
        if (iter>2 && abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            coef <- start
            break
        }
        if (kk > 0) break ## already shrunk the step
        ## At this point it is worth checking that coef update actually improved the penalized
        ## deviance. If not try step halving, and redo the above once a suitable step has been
        ## found...
        if (iter>2) { ## can test divergence
          ## need to compute penalty at start and end of step
	  if (efam) {
	    dev0 <- sum(dev.resids(G$y,linkinv(eta0),G$w,theta0)) ## depends on theta, which will have changed
	    dev1 <- sum(dev.resids(G$y,linkinv(eta),G$w,theta0)) ## depends on theta, which will have changed
          } else { dev1 <- dev }
          if (method=="fREML") {
	    pcoef <- fit$beta
            Sb0 <- Sl.Sb(um$Sl,rho=log(object$full.sp),pcoef0)
	    Sb <- Sl.Sb(um$Sl,rho=log(object$full.sp),pcoef)
          } else {
	    pcoef <- coef
	    full.sp <- if (is.null(object$full.sp)) object$sp else object$full.sp
            Sb0 <- regular.Sb(G$S,G$off,full.sp,pcoef0)
	    Sb <- regular.Sb(G$S,G$off,full.sp,pcoef)
          }
	  while (dev0 + sum(pcoef0*Sb0) < dev1 + sum(pcoef * Sb) && kk < 6) {
            ## shrink step ...
            coef <- (coef0 + coef)/2
	    pcoef <- (pcoef0 + pcoef)/2
	    eta <- (eta0 + eta)/2
	    Sb <- (Sb0 + Sb)/2
	    ## recompute deviance ...
	    dev <- if (efam) sum(dev.resids(G$y,linkinv(eta),G$w,theta)) else sum(dev.resids(G$y,linkinv(eta),G$w)) 
            dev1 <- if (efam) sum(dev.resids(G$y,linkinv(eta),G$w,theta0)) else dev
            kk <- kk + 1
          }
        }
	if (kk == 0) break ## step was ok
      } ## repeat
      
      if (conv) break

      if (iter>1) { ## store coef and eta for divergence checking
        coef0 <- coef
	if (efam) theta0 <- theta ## theta used for determining step
	pcoef0 <- if (method=="fREML") fit$beta else coef
	eta0 <- eta
	dev0 <- dev
      }

      if (efam && iter>1) { ## estimate theta
        if (family$n.theta>0||scale1<0) theta <- estimate.theta(theta,family,G$y,linkinv(eta),scale=scale1,wt=G$w,tol=1e-7)
        if (!is.null(family$scale) && scale1<0) {
	   scale <- exp(theta[family$n.theta+1])
	   theta <- theta[1:family$n.theta]
	}  
        family$putTheta(theta)
      }
	  
      if (method=="GCV.Cp") {
         fit <- magic(qrx$f,qrx$R,G$sp,G$S,G$off,L=G$L,lsp0=G$lsp0,rank=G$rank,
                      H=G$H,C=matrix(0,0,ncol(qrx$R)),     ##C=G$C,
                      gamma=gamma,scale=scale,gcv=(scale<=0),
                      extra.rss=rss.extra,n.score=nobs+nobs.extra)
 
         post <- magic.post.proc(qrx$R,fit,qrx$f*0+1) 
      } else if (method=="fREML") { ## use fast REML code
        ## block diagonal penalty object, Sl, set up before loop
        um <- Sl.Xprep(Sl,qrx$R,nt=npt)
        lambda.0 <- if (is.null(in.out)) initial.sp(qrx$R,G$S,G$off) else in.out$sp
        lsp0 <- log(lambda.0) ## initial s.p.
        ## carry forward scale estimate if possible...
        if (scale>0) log.phi <- log(scale) else {
          if (iter>1) log.phi <- log(object$scale) else {
	    if (is.null(in.out)) {
              if (is.null(coef)||qrx$y.norm2==0) log.phi <- log(var(as.numeric(G$y))*.05) else
                 log.phi <- log(qrx$y.norm2/(nobs+nobs.extra))
	    } else log.phi <- log(in.out$scale)	 
          }
        }
        fit <- fast.REML.fit(um$Sl,um$X,qrx$f,rho=lsp0,L=G$L,rho.0=G$lsp0,
                             log.phi=log.phi,phi.fixed=scale>0,rss.extra=rss.extra,
                             nobs =nobs+nobs.extra,Mp=um$Mp,nt=npt,gamma=gamma)
        res <- Sl.postproc(Sl,fit,um$undrop,qrx$R,cov=FALSE,L=G$L,nt=npt)
        object <- list(coefficients=res$beta,db.drho=fit$d1b,
                       gcv.ubre=fit$reml,mgcv.conv=list(iter=fit$iter,
                       message=fit$conv),rank=ncol(um$X),
                       Ve=NULL,scale.estimated = scale<=0,outer.info=fit$outer.info,
                        optimizer=c("perf","newton"))
 
        if (scale<=0) { ## get sp's and scale estimate
          nsp <- length(fit$rho)
          object$sig2 <- object$scale <- exp(fit$rho[nsp])
          object$sp <- exp(fit$rho[-nsp]) 
          nsp <- length(fit$rho.full)
          object$full.sp <- exp(fit$rho.full[-nsp])
        } else { ## get sp's
          object$sig2 <- object$scale <- scale  
          object$sp <- exp(fit$rho)
          object$full.sp <- exp(fit$rho.full)
        }
        class(object)<-c("gam")               
      } else { ## method is one of "ML", "P-REML" etc...
        y <- G$y; w <- G$w; n <- G$n;offset <- G$offset
        G$y <- qrx$f
        G$w <- G$y*0+1
        G$X <- qrx$R
        G$n <- length(G$y)
        G$offset <- G$y*0
        G$dev.extra <- rss.extra
        G$pearson.extra <- rss.extra
        G$n.true <- nobs+nobs.extra
        object <- gam(G=G,method=method,gamma=gamma,scale=scale,control=gam.control(nthreads=npt))
        y -> G$y; w -> G$w; n -> G$n;offset -> G$offset
	object$deviance <- object$family <- object$null.deviance <- object$fitted.values <- NULL
      }
     
      if (method=="GCV.Cp") { 
        object <- list()
        object$coefficients <- fit$b
        object$edf <- post$edf 
        object$edf1 <- post$edf1
        object$full.sp <- fit$sp.full
        object$gcv.ubre <- fit$score
        object$hat <- post$hat
        object$mgcv.conv <- fit$gcv.info 
        object$optimizer="magic"
        object$rank <- fit$gcv.info$rank
        object$Ve <- post$Ve
        object$Vp <- post$Vb
        object$sig2 <- object$scale <- fit$scale
        object$sp <- fit$sp
        names(object$sp) <- names(G$sp)
        class(object)<-c("gam")
      }

      coef <- object$coefficients
        
      if (any(!is.finite(coef))) {
          conv <- FALSE
          warning(gettextf("non-finite coefficients at iteration %d",
                  iter))
          break
      }
    } ## end fitting iteration

    if (!is.null(family$setInd)) family$setInd(NULL) ## reset family to state before indexing

    if (method=="fREML") { ## do expensive cov matrix cal only at end
      res <- Sl.postproc(Sl,fit,um$undrop,qrx$R,cov=TRUE,scale=scale,L=G$L,nt=npt)
      object$edf <- res$edf
      object$edf1 <- res$edf1
      object$edf2 <- res$edf2
      object$hat <- res$hat
      object$Vp <- res$Vp
      object$Ve <- res$Ve
      object$Vc <- res$Vc
      object$V.sp <- res$V.sp
    }
    
    if (efam) { ## deal with any post processing
       if (!is.null(family$postproc)) {
         object$family <- family
         posr <- family$postproc(family=family,y=G$y,prior.weights=G$w,
              fitted=linkinv(eta),linear.predictors=eta,offset=G$offset,
	      intercept=G$intercept)
         if (!is.null(posr$family)) object$family$family <- posr$family
         if (!is.null(posr$deviance)) object$deviance <- posr$deviance
         if (!is.null(posr$null.deviance)) object$null.deviance <- posr$null.deviance
       }
      if (is.null(object$null.deviance)) object$null.deviance <- sum(family$dev.resids(G$y,weighted.mean(G$y,G$w),G$w,theta))   
    }

    if (!conv)
       warning("algorithm did not converge")
   
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
         if (any(mu > 1 - eps) || any(mu < eps))
                warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
            if (any(mu < eps))
                warning("fitted rates numerically 0 occurred")
    }
  object$R <- qrx$R    
  object$iter <- iter 
  object$wt <- wt
  object$y <- G$y
  object$prior.weights <- G$w
  rm(G);if (gc.level>0) gc()
  object
} ## end bgam.fit




ar.qr_up <- function(arg) {
## function to perform QR updating with AR residuals, on one execution thread
  if (arg$rho!=0) { ## AR1 error model
     ld <- 1/sqrt(1 - arg$rho^2) ## leading diagonal of root inverse correlation
     sd <- -arg$rho * ld         ## sub diagonal
  } 
  yX.last <- NULL
  qrx <- list(R=NULL,f=array(0,0),y.norm2=0) ## initial empty qr object
  for (i in 1:arg$n.block) {
    ind <- arg$start[i]:arg$end[i] 
    if (arg$rho!=0) { ## have to find AR1 transform...
       N <- arg$end[i]-arg$start[i]+1
       ## note first row implied by this transform
       ## is always dropped, unless really at beginning of data.
       row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)])
       weight <- c(1,rep(c(sd,ld),N-1))
       stop <- c(1,1:(N-1)*2+1)
       if (!is.null(arg$mf$"(AR.start)")) { ## need to correct the start of new AR sections...
           ii <- which(arg$mf$"(AR.start)"[ind]==TRUE)
           if (length(ii)>0) {
             if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
             weight[ii*2-2] <- 0 ## zero sub diagonal
             weight[ii*2-1] <- 1 ## set leading diagonal to 1
           }
       }
     } 
     w <- sqrt(arg$G$w[ind])
     X <- w*predict(arg$G,newdata=arg$mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))
     y <- w*(arg$mf[ind,arg$response] - arg$offset[ind]) ## w*(arg$G$model[[arg$response]] - arg$offset[ind])
     if (arg$rho!=0) {
       ## Apply transform...
       if (arg$last&&arg$end[i]==arg$nobs) yX.last <- 
           c(y[nrow(X)],X[nrow(X),]) ## store final row, in case of update
       if (arg$first&&i==1) {
          X <- rwMatrix(stop,row,weight,X)
          y <- rwMatrix(stop,row,weight,y)
       } else {
          X <- rwMatrix(stop,row,weight,X)[-1,]
          y <- rwMatrix(stop,row,weight,y)[-1]
       } 
     } ## dealt with AR1      
     qrx <- qr_update(X,y,qrx$R,qrx$f,qrx$y.norm2,use.chol=arg$use.chol)
     rm(X);if (arg$gc.level>1) {gc()} ## X can be large: remove and reclaim
  } ## all blocks dealt with
  qrx$yX.last <- yX.last
  if (arg$gc.level>1) {rm(arg,w,y,ind);gc()}
  qrx
} ## ar.qr_up

pabapr <- function(arg) {
## function for parallel calling of predict.gam
## QUERY: ... handling?
  predict.gam(arg$object,newdata=arg$newdata,type=arg$type,se.fit=arg$se.fit,terms=arg$terms,
                        block.size=arg$block.size,newdata.guaranteed=arg$newdata.guaranteed,
                        na.action=arg$na.action)
}

predict.bam <- function(object,newdata,type="link",se.fit=FALSE,terms=NULL,exclude=NULL,
                        block.size=50000,newdata.guaranteed=FALSE,na.action=na.pass,
                        cluster=NULL,discrete=TRUE,n.threads=1,gc.level=0,...) {
## function for prediction from a bam object, possibly in parallel
  
  #if (is.function(na.action)) na.action <- deparse(substitute(na.action)) ## otherwise predict.gam can't detect type
  if (discrete && !is.null(object$dinfo)) {
    return(predict.bamd(object,newdata,type,se.fit,terms,exclude,
                        block.size,newdata.guaranteed,na.action,n.threads,gc.level=gc.level,...))
  }
  ## remove some un-needed stuff from object
  object$Sl <- object$qrx <- object$R <- object$F <- object$Ve <-
  object$Vc <- object$G <- object$residuals <- object$fitted.values <-
  object$linear.predictors <- NULL
  if (gc.level>0) gc()
  if (!is.null(cluster)&&inherits(cluster,"cluster")) { 
     ## require(parallel)
     n.threads <- length(cluster)
  } else n.threads <- 1
  if (missing(newdata)) n <- nrow(object$model) else {
    n <- if (is.matrix(newdata[[1]])) nrow(newdata[[1]]) else length(newdata[[1]]) 
  }
  if (n < 100*n.threads) n.threads <- 1 ## not worth the overheads
  if (n.threads==1) { ## single threaded call
    if (missing(newdata)) return(
      predict.gam(object,newdata=object$model,type=type,se.fit=se.fit,terms=terms,exclude=exclude,
                        block.size=block.size,newdata.guaranteed=newdata.guaranteed,
                        na.action=na.action,...)
    ) else return(
      predict.gam(object,newdata=newdata,type=type,se.fit=se.fit,terms=terms,exclude=exclude,
                        block.size=block.size,newdata.guaranteed=newdata.guaranteed,
                        na.action=na.action,...))
  } else { ## parallel call...
    nt <- rep(floor(n/n.threads),n.threads)
    nt[1] <- n - sum(nt[-1])
    arg <- list()
    n1 <- 0
    for (i in 1:n.threads) { 
      n0 <- n1+1;n1 <- n1+nt[i]
      ind <- n0:n1 ## this thread's data block from mf
      arg[[i]] <- list(object=object,type=type,se.fit=se.fit,terms=terms,exclude=exclude,
                        block.size=block.size,newdata.guaranteed=newdata.guaranteed,
                        na.action=na.action)
      arg[[i]]$object$model <- object$model[1:2,] ## save space
      if (missing(newdata)) {
        arg[[i]]$newdata <- object$model[ind,]
      } else {
        arg[[i]]$newdata <- newdata[ind,]
      }
    } ## finished setting up arguments
    ## newdata and object no longer needed - all info in thread lists...
    if (!missing(newdata)) rm(newdata)
    rm(object)
    if (gc.level>0) gc()
    res <- parallel::parLapply(cluster,arg,pabapr) ## perform parallel prediction
    if (gc.level>0) gc()
    ## and splice results back together...
    if (type=="lpmatrix") {
      X <- res[[1]]
      for (i in 2:length(res)) X <- rbind(X,res[[i]])
      return(X)
    } else if (se.fit==TRUE) {
      rt <- list(fit = res[[1]]$fit,se.fit = res[[1]]$se.fit)
      if (type=="terms") {
        for (i in 2:length(res)) { 
          rt$fit <- rbind(rt$fit,res[[i]]$fit)
          rt$se.fit <- rbind(rt$se.fit,res[[i]]$se.fit)
        }
      } else {
        for (i in 2:length(res)) { 
          rt$fit <- c(rt$fit,res[[i]]$fit)
          rt$se.fit <- c(rt$se.fit,res[[i]]$se.fit)
        }
      }
      return(rt)
    } else { ## no se's returned
      rt <- res[[1]]
       if (type=="terms") {
        for (i in 2:length(res)) rt <- rbind(rt,res[[i]])
      } else {
        for (i in 2:length(res)) rt <- c(rt,res[[i]])
      }
      return(rt)
    } 
  }
} ## end predict.bam 


bam.fit <- function(G,mf,chunk.size,gp,scale,gamma,method,rho=0,
                    cl=NULL,gc.level=0,use.chol=FALSE,in.out=NULL,npt=1) {
## function that does big additive model fit in strictly additive case
   ## first perform the QR decomposition, blockwise....
   n <- nrow(mf)
   if (rho!=0) { ## AR1 error model
     ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
     sd <- -rho*ld         ## sub diagonal
   }

   if (n>chunk.size) { ## then use QR accumulation approach
     if (!is.null(cl)&&inherits(cl,"cluster")) { 
       n.threads <- length(cl)
       while(n/n.threads < ncol(G$X)) n.threads <- n.threads - 1
       if (n.threads < length(cl)) { 
         warning("Too many cluster nodes to use all efficiently")
       }
     } else n.threads <- 1

     G$coefficients <- rep(0,ncol(G$X))
     class(G) <- "gam"

     if (n.threads>1) { ## set up thread argument lists
       ## number of obs per thread
       nt <- rep(ceiling(n/n.threads),n.threads)
       nt[n.threads] <- n - sum(nt[-n.threads])
       arg <- list()
       n1 <- 0
       for (i in 1:n.threads) { 
         n0 <- n1+1;n1 <- n1+nt[i]
         if (i>1&&rho!=0) { ## need to start from end of last block if rho!=0
           n0 <- n0-1;nt[i] <- nt[i]+1 
         }   
         ind <- n0:n1 ## this thread's data block from mf
         n.block <- nt[i]%/%chunk.size ## number of full sized blocks
         stub <- nt[i]%%chunk.size ## size of end block
         if (n.block>0) { 
           ## each block is of size 
           start <- (0:(n.block-1))*chunk.size+1
           end <- start + chunk.size - 1
           if (stub>0) {
             start[n.block+1] <- end[n.block]+1
             end[n.block+1] <- nt[i]
             n.block <- n.block+1
           } 
           if (rho!=0) { ## then blocks must overlap
             ns <- length(start)
             if (ns>1) start[2:ns] <- start[2:ns]-1 
           }
         } else {
           n.block <- 1
           start <- 1
           end <- nt[i]
         }
         arg[[i]] <- list(nobs= nt[i],start=start,end=end,n.block=n.block,
                         rho=rho,mf = mf[ind,],gc.level=gc.level,
                         offset = G$offset[ind],G = G,response=gp$response,
                         first=FALSE,last=FALSE,use.chol=use.chol)
         if (i==1) arg[[1]]$first <- TRUE
         if (i==n.threads) arg[[i]]$last <- TRUE 
         arg[[i]]$G$w <- G$w[ind];arg[[i]]$G$model <- NULL
       }
     } else { ## single thread, requires single indices 
       n.block <- n%/%chunk.size ## number of full sized blocks
       stub <- n%%chunk.size ## size of end block
       if (stub>0) n.block <- n.block + 1
       start <- 0:(n.block-1)*chunk.size    ## block starts
       end <- start + chunk.size;           ## block ends
       end[n.block] <- n
       if (rho==0) start <- start + 1  ## otherwise most blocks go to 1 before block start
       start[1] <- 1  
     } 
    
     if (n.threads==1) { ## use original single thread method...
       qrx <- list(R=NULL,f=array(0,0),y.norm2=0) ## initial empty qr object
       for (i in 1:n.block) {
         ind <- start[i]:end[i] 
         if (rho!=0) {
           N <- end[i]-start[i]+1

           row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)])
           weight <- c(1,rep(c(sd,ld),N-1))
           stop <- c(1,1:(N-1)*2+1) 
           if (!is.null(mf$"(AR.start)")) { ## need to correct the start of new AR sections...
             ii <- which(mf$"(AR.start)"[ind]==TRUE)
             if (length(ii)>0) {
               if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
               weight[ii*2-2] <- 0 ## zero sub diagonal
               weight[ii*2-1] <- 1 ## set leading diagonal to 1
             }
           }
         } 
         w <- sqrt(G$w[ind])
         X <- w*predict(G,newdata=mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))
         y <- w*(mf[ind,gp$response]-G$offset[ind])  ## w*(G$model[[gp$response]] - G$offset[ind])
         if (rho!=0) {
           ## Apply transform...
           if (end[i]==n) yX.last <- c(y[nrow(X)],X[nrow(X),]) ## store final row, in case of update
           if (i==1) {
             X <- rwMatrix(stop,row,weight,X)
             y <- rwMatrix(stop,row,weight,y)
           } else {
             X <- rwMatrix(stop,row,weight,X)[-1,]
             y <- rwMatrix(stop,row,weight,y)[-1]
           } 
         }      

         qrx <- qr_update(X,y,qrx$R,qrx$f,qrx$y.norm2,use.chol=use.chol,nt=npt)
         rm(X)
         if (gc.level>1) {gc()} ## X can be large: remove and reclaim
       } ## end of single thread block loop
       if (use.chol) { ## post proc to get R and f...
          y.norm2 <- qrx$y.norm2 
          qrx <- chol2qr(qrx$R,qrx$f,nt=npt)
          qrx$y.norm2 <- y.norm2
        }
     } else { ## use parallel accumulation
     
       res <- parallel::parLapply(cl,arg,ar.qr_up)

       ## now consolidate the results from the parallel threads...
       R <- res[[1]]$R;f <- res[[1]]$f; ## dev <- res[[1]]$dev
       y.norm2 <- res[[1]]$y.norm2
       for (i in 2:n.threads) {
         if (use.chol) {
           R <- R + res[[i]]$R; f <- f + res[[i]]$f
         } else {
           R <- rbind(R,res[[i]]$R); f <- c(f,res[[i]]$f)
         }
         y.norm2 <- y.norm2 + res[[i]]$y.norm2
       } 
       if (use.chol) {
         qrx <- chol2qr(R,f,nt=npt)
         qrx$y.norm2 <- y.norm2
       } else { ## proper QR         
         ## use parallel QR if npt>1...
         qrx <- if (npt>1) pqr2(R,npt) else qr(R,tol=0,LAPACK=TRUE) 
         f <- qr.qty(qrx,f)[1:ncol(R)]
         rp <- qrx$pivot;rp[rp] <- 1:ncol(R) # reverse pivot
         qrx <- list(R=qr.R(qrx)[,rp],f=f,y.norm2=y.norm2)
       }
       yX.last <- res[[n.threads]]$yX.last
     } 
     G$n <- n
   
   } else { ## n <= chunk.size
     if (rho==0) qrx <- qr_update(sqrt(G$w)*G$X,sqrt(G$w)*(G$y-G$offset),use.chol=use.chol,nt=npt) else {
       row <- c(1,rep(1:n,rep(2,n))[-c(1,2*n)])
       weight <- c(1,rep(c(sd,ld),n-1))
       stop <- c(1,1:(n-1)*2+1)
       if (!is.null(mf$"(AR.start)")) { ## need to correct the start of new AR sections...
         ii <- which(mf$"(AR.start)"==TRUE)
         if (length(ii)>0) {
           if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
           weight[ii*2-2] <- 0 ## zero sub diagonal
           weight[ii*2-1] <- 1 ## set leading diagonal to 1
         }
       }
       yX.last <- c(G$y[n],G$X[n,])  ## store final row, in case of update
       X <- rwMatrix(stop,row,weight,sqrt(G$w)*G$X)
       y <- rwMatrix(stop,row,weight,sqrt(G$w)*G$y)
       qrx <- qr_update(X,y,use.chol=use.chol,nt=npt)
   
       rm(X); if (gc.level>1) gc() ## X can be large: remove and reclaim
     } 
     if (use.chol) { ## post proc to get R and f...
        y.norm2 <- qrx$y.norm2 
        qrx <- chol2qr(qrx$R,qrx$f,nt=npt)
        qrx$y.norm2 <- y.norm2
     }
   }

   rss.extra <- qrx$y.norm2 - sum(qrx$f^2)
 
   if (method=="GCV.Cp") {
     fit <- magic(qrx$f,qrx$R,G$sp,G$S,G$off,L=G$L,lsp0=G$lsp0,rank=G$rank,
                H=G$H,C=matrix(0,0,ncol(qrx$R)),  ##C=G$C,
                gamma=gamma,scale=scale,gcv=(scale<=0),
                extra.rss=rss.extra,n.score=n)
 
     post <- magic.post.proc(qrx$R,fit,qrx$f*0+1) 
   } else if (method=="fREML"){ ## use fast REML code
     Sl <- Sl.setup(G) ## setup block diagonal penalty object
     um <- Sl.Xprep(Sl,qrx$R,nt=npt)
     lambda.0 <- if (is.null(in.out)) initial.sp(qrx$R,G$S,G$off) else in.out$sp
     lsp0 <- log(lambda.0) ## initial s.p.
     log.phi <- if (scale<=0) { if (is.null(in.out)) log(var(as.numeric(G$y))*.05) else log(in.out$scale) } else log(scale)
     fit <- fast.REML.fit(um$Sl,um$X,qrx$f,rho=lsp0,L=G$L,rho.0=G$lsp0,
            log.phi=log.phi,phi.fixed=scale>0,rss.extra=rss.extra,
            nobs =n,Mp=um$Mp,nt=npt,gamma=gamma)
     res <- Sl.postproc(Sl,fit,um$undrop,qrx$R,cov=TRUE,scale=scale,L=G$L,nt=npt)
     object <- list(coefficients=res$beta,edf=res$edf,edf1=res$edf1,edf2=res$edf2,##F=res$F,
                    db.drho=fit$d1b,
                    gcv.ubre=fit$reml,hat=res$hat,mgcv.conv=list(iter=fit$iter,
                    message=fit$conv),rank=ncol(um$X),
                    Ve=res$Ve,Vp=res$Vp,Vc=res$Vc,V.sp=res$V.sp,
                    scale.estimated = scale<=0,outer.info=fit$outer.info,
                    optimizer=c("perf","newton"))
     if (scale<=0) { ## get sp's and scale estimate
       nsp <- length(fit$rho)
       object$sig2 <- object$scale <- exp(fit$rho[nsp])
       object$sp <- exp(fit$rho[-nsp])
       nsp <- length(fit$rho.full)
       object$full.sp <- exp(fit$rho.full[-nsp])
     } else { ## get sp's
       object$sig2 <- object$scale <- scale  
       object$sp <- exp(fit$rho)
       object$full.sp <- exp(fit$rho.full)
     }
     
     if (rho!=0) { ## correct RE/ML score for AR1 transform
       df <- if (is.null(mf$"(AR.start)")) 1 else sum(mf$"(AR.start)")
       object$gcv.ubre <- object$gcv.ubre - (n-df)*log(ld)
     }

     G$X <- qrx$R;G$dev.extra <- rss.extra
     G$pearson.extra <- rss.extra;G$n.true <- n
     object$Sl <- Sl ## to allow for efficient update
     class(object)<-c("gam")
   } else { ## method is "ML", "P-REML" or similar
     y <- G$y; w <- G$w; n <- G$n;offset <- G$offset
     G$y <- qrx$f
     G$w <- G$y*0+1
     G$X <- qrx$R
     G$n <- length(G$y)
     G$offset <- G$y*0
     G$dev.extra <- rss.extra
     G$pearson.extra <- rss.extra
     G$n.true <- n
     object <- gam(G=G,method=method,gamma=gamma,scale=scale,control=gam.control(nthreads=npt))
     object$null.deviance <- object$fitted.values <- NULL
     y -> G$y; w -> G$w; n -> G$n;offset -> G$offset
     if (rho!=0) { ## correct RE/ML score for AR1 transform 
       df <- if (is.null(mf$"(AR.start)")) 1 else sum(mf$"(AR.start)")
       object$gcv.ubre <- object$gcv.ubre - (n-df)*log(ld)
     }
   }
   if (method=="GCV.Cp") { 
     object <- list()
     object$coefficients <- fit$b
     object$edf <- post$edf
     object$edf1 <- post$edf1
     object$full.sp <- fit$sp.full
     object$gcv.ubre <- fit$score
     object$hat <- post$hat
     object$mgcv.conv <- fit$gcv.info 
     object$optimizer="magic"
     object$rank <- fit$gcv.info$rank
     object$Ve <- post$Ve
     object$Vp <- post$Vb
     object$sig2 <- object$scale <- fit$scale
     object$sp <- fit$sp
     class(object)<-c("gam")
   } else {
    
   }
   G$smooth <- G$X <- NULL
   object$prior.weights <- G$w
   object$AR1.rho <- rho
   if (rho!=0) { ## need to store last model matrix row, to allow update
     object$yX.last <- yX.last
   }
  
   object$R <- qrx$R
   object$gamma <- gamma;object$G <- G;object$qrx <- qrx ## to allow updating of the model
   object$y <- mf[[gp$response]]
   object$iter <- 1
   object
} # end of bam.fit

predict.bamd <- function(object,newdata,type="link",se.fit=FALSE,terms=NULL,exclude=NULL,
                         block.size=50000,newdata.guaranteed=FALSE,na.action=na.pass,
			 n.threads=1,gc.level=0,...) {
## function for prediction from a bam object, by discrete methods
## Note that behaviour with factor levels not in fit differs from predict.gam - the dummies
## for the factor all get set to zero (there is a warning). 

  ## remove some un-needed stuff from object
  object$Sl <- object$qrx <- object$R <- object$F <- object$Ve <-
  object$Vc <- object$G <- object$residuals <- object$fitted.values <-
  object$linear.predictors <- NULL
  if (gc.level>0) gc()
  if (missing(newdata)) newdata <- object$model
   
  #convert2mf <- is.null(attr(newdata,"terms")) ### was needed to get variables post any transform for discretization

  if (type=="iterms") {
    type <- "terms"
    warning("iterms reset to terms")
  }

  lpi <- attr(object$formula,"lpi") ## lpi[[i]] indexes coefs for ith linear predoctor
  nlp <- if (is.null(lpi)) 1 else length(lpi) ## number of linear predictors
  lpid <- if (nlp>1) object$dinfo$lpid else NULL ## index of discrete terms involved in each linear predictor
 
  ## newdata has to be processed first to avoid, e.g. dropping different subsets of data
  ## for parametric and smooth components....

  newterms <- attr(newdata,"terms") ## non NULL for model frame

  ## for discretization to work, we need the model frame names, not data frame ones. i.e. variables post any transform
  object$pred.formula <- reformulate(object$dinfo$gp$fake.names)

  newd <- predict.gam(object,newdata=newdata,type="newdata",se.fit=se.fit,terms=terms,exclude=exclude,
            block.size=block.size,newdata.guaranteed=newdata.guaranteed,
            na.action=na.action,...) 
  if (nrow(newd)>0) {
    newdata <- newd;rm(newd)
    ## could use xlev indicator of extra factor levels to attach xlev attributes to Xd below,
    ## then modify Xbd to use these to insert NA at extra level positions - seems convoluted
    ## for the tiny gain (NAs when term not excluded. )
    #for (i in which(sapply(newdata,function(x) !is.null(attr(x,"xlev"))))) {
      ## EXPERIMENTAL: problem is that "xlev" lost at model frame stage...
      #newdata[[i]][attr(newdata[[i]],"xlev")] <- levels(newdata[[i]])[1] ## reset extra levels to first level
    #}
  } else { ## no non NA data, might as well call predict.gam
    return(predict.gam(object,newdata=newdata,type=type,se.fit=se.fit,terms=terms,exclude=exclude,
            block.size=block.size,newdata.guaranteed=newdata.guaranteed,
            na.action=na.action,...))
  }

  ## Next line needed to avoid treating newdata as a model frame if it was supplied not as a model frame.
  ## Otherwise names of e.g. offset are messed up (as they will also be if it was supplied as a model frame
  ## or was set to object$model and we set terms to NULL)
  
  if (is.null(newterms)) attr(newdata,"terms") <- NULL
  
  na.act <- attr(newdata,"na.action") ## save the NA action for later

  yname <- attr(attr(object$terms,"dataClasses"),"names")[attr(object$terms,"response")]
  response <- newdata[[yname]]
  
  ## Parametric terms have to be dealt with safely, but without forming all terms 
  ## or a full model matrix. Strategy here is to use predict.gam, having removed
  ## key smooth related components from model object, so that it appears to be
  ## a parametric model... 
  offset <- if (nlp>1) as.list(rep(0,nlp)) else 0

  pp <- if (se.fit) list(fit=rep(0,0),se.fit=rep(0,0)) else rep(0,0) ## note: needed in terms branch below

  ## now discretize covariates...
  
  #if (convert2mf) newdata <- model.frame(object$dinfo$gp$fake.formula[-2],newdata) ## Why??
  
  dk <- discrete.mf(object$dinfo$gp,mf=newdata,names.pmf=object$dinfo$pmf.names,full=TRUE)
  Xd <- list() ### list of discrete model matrices...
  n <- if (is.matrix(newdata[[1]])) nrow(newdata[[1]]) else length(newdata[[1]])
  kb <- k <- 1;kd <- dk$k
  ks <- matrix(0,0,2) ## NOTE: slightly more efficient not to repeatedly extend
  for (j in 1:nlp) { ## loop over parametric components of linear predictors
      ## get discretized marginal matrices for the parametric model (if contrast warnings see predict.gam fix)
      ptens <- if (nlp==1&&!is.list(object$pterms)) terms2tensor(object$pterms,dk$mf,object$contrasts) else
	             terms2tensor(object$pterms[[j]],dk$mf,object$contrasts)
      offi <- if (nlp==1&&!is.list(object$pterms)) attr(object$pterms,"offset") else attr(object$pterms[[j]],"offset")
      if (!is.null(offi)) {
        offname <- if (nlp==1&&!is.list(object$pterms)) names(attr(object$pterms,"dataClasses"))[offi] else
	                                                names(attr(object$pterms[[j]],"dataClasses"))[offi]
	if (j==1&&nlp>1) offset <- list()
	if (nlp>1) offset[[j]] <- newdata[[offname]] else  offset <- newdata[[offname]]
      }
      if (!is.null(ptens)) {
        np <-length(ptens$X);
        for (i in 1:np) {
	  ks <-  rbind(ks,dk$ks[ptens$xname[i],])
	  Xd[[k]] <- ptens$X[[i]][1:dk$nr[ptens$xname[i]],,drop=FALSE]
	  k <- k + 1
        }	
      }
      kb <- kb + length(ptens$ts)
  } ## parametric lp loop
 
  ts <- object$dinfo$ts
  dt <- object$dinfo$dt

  ## remove padding from the discretized data...
  class(dk$mf) <- "list" 
  for (i in 1:length(dk$mf)) dk$mf[[i]] <- if (is.matrix(dk$mf[[i]])) dk$mf[[i]][1:dk$nr[i],] else dk$mf[[i]][1:dk$nr[i]]

  if (length(object$smooth)) for (i in 1:length(object$smooth)) { ## work through the smooth list
    ## potentially each smoother model matrix can be made up of a sequence
    ## of row-tensor products, nead to loop over such sub blocks...
    nsub <- if (!is.null(object$smooth[[i]]$ts)) length(object$smooth[[i]]$ts) else 1
    ## first deal with any by variable (as first marginal of tensor)...
    for (sb in 1:nsub) { ## loop over sub-blocks
      if (object$smooth[[i]]$by!="NA") {
        termk <- object$smooth[[i]]$by
        by.var <- dk$mf[[object$smooth[[i]]$by]][1:dk$nr[termk]]
        if (is.factor(by.var)) { 
          ## create dummy by variable...
          by.var <- as.numeric(by.var==object$smooth[[i]]$by.level)  
        }
        Xd[[k]] <- matrix(by.var,dk$nr[termk],1)
        ks <- rbind(ks,dk$ks[termk,])
        k <- k + 1
        by.present <- 1
      } else by.present <- 0
      ## ... by done
      if (inherits(object$smooth[[i]],"tensor.smooth")) { 
        nmar <- if (is.null(object$smooth[[i]]$dt)) length(object$smooth[[i]]$margin) else object$smooth[[i]]$dt[sb]
        XP <- object$smooth[[i]]$XP
	jind <- if (sb>1) object$smooth[[i]]$ts[sb] + 1:object$smooth[[i]]$dt[sb] - 1 else 1:nmar       
        for (j in jind) {
          object$smooth[[i]]$margin[[j]]$by<- "NA" ## should be no by's here (any by dealt with above)
	  termk <- object$smooth[[i]]$margin[[j]]$term[1]
          Xd[[k]] <-if (termk=="(Intercept)") matrix(1,dk$nr[termk],1) else PredictMat(object$smooth[[i]]$margin[[j]],dk$mf,n=dk$nr[termk])
	  ks <- rbind(ks,dk$ks[termk,])
          if (!is.null(XP)&&(j<=length(XP))&&!is.null(XP[[j]])) Xd[[k]] <- Xd[[k]]%*%XP[[j]]
          k <- k + 1 
        }
      } else { ## not a tensor smooth
        object$smooth[[i]]$by <- "NA" ## have to ensure by not applied here (it's dealt with as a tensor marginal)!
        termk <- object$smooth[[i]]$term[1]
        ks <- rbind(ks,dk$ks[termk,])
        Xd[[k]] <- PredictMat(object$smooth[[i]],dk$mf,n=dk$nr[termk])
        k <- k + 1
      }
      kb <- kb + 1
    } ## sub block loop  
  }

  attr(Xd,"lpip") <- object$dinfo$lpip ## list of coef indices for each term

  ## end of discrete set up
  se <- se.fit
  if (type=="terms") {
    term.lab <- names(object$dinfo$ts) ## names of all terms
    termi <- rep(TRUE,length(term.lab)) ## do we want term included?
    termi[term.lab=="(Intercept)"] <- FALSE ## not if it's intercept
    if (!is.null(terms)) termi <- termi & term.lab %in% terms ## and only if it's in terms
    if (!is.null(exclude)) termi <- termi & !(term.lab %in% exclude) ## and not in exclude
    uterms <- unique(term.lab[termi])
    nterms <- length(uterms) ## sum(termi) 
    fit <- matrix(0,nrow(kd),nterms)
    if (se) se.fit <- matrix(0,nrow(kd),nterms)
    for (k in 1:nterms) {
      lt <- which(term.lab%in%uterms[k]) ## discrete terms involved in this smooth (can be more than one)
      if (termi[lt]) {
        fit[,k] <- Xbd(Xd,object$coefficients,kd,ks,                          
                    ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop,lt=lt)
        if (se) se.fit[,k] <- diagXVXd(Xd,object$Vp,kd,ks,                
                     ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop,nthreads=n.threads,lt=lt,rt=lt)^.5
      }	
    }
    colnames(fit) <- uterms
    if (se) {
      colnames(se.fit) <- uterms
      fit <- list(fit=fit,se.fit=se.fit)
    } 
  } else if (type=="lpmatrix") {
    lt <- 1:length(ts);names(lt) <- names(ts)
    if (!is.null(terms)) lt <- lt[names(lt)%in%terms]
    if (!is.null(exclude)) lt <- lt[!names(lt)%in%exclude]
    fit <- Xbd(Xd,diag(length(object$coefficients)),kd,ks,ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop,lt=lt)
    if (nlp>1) attr(fit,"lpi") <- lpi
  } else { ## link or response 
    if (is.null(object$family$predict)||type=="link") {
      lt <- 1:length(ts);names(lt) <- names(ts)
      if (nlp>1) {
        fit <- matrix(0,nrow(kd),nlp)
	if (se) se.fit <- matrix(0,nrow(kd),nlp)
        for (i in 1:nlp) {
          lt <- 1:length(ts);names(lt) <- names(ts)
	  lt <- lt[lpid[[i]]]
          if (!is.null(terms)) lt <- lt[names(lt)%in%terms]
          if (!is.null(exclude)) lt <- lt[!names(lt)%in%exclude]
          fit[,i] <- Xbd(Xd,object$coefficients,kd,ks,ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop,lt=lt) + offset[[i]]
	  if (se) se.fit[,i] <- diagXVXd(Xd,object$Vp,kd,ks,ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop,
	                                        lt=lt,rt=lt,nthreads=n.threads)^.5
	}  
      } else {
        lt <- 1:length(ts);names(lt) <- names(ts)
        if (!is.null(terms)) lt <- lt[names(lt)%in%terms]
        if (!is.null(exclude)) lt <- lt[!names(lt)%in%exclude]
        fit <- Xbd(Xd,object$coefficients,kd,ks,ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop,lt=lt) + offset
	if (se) se.fit <- diagXVXd(Xd,object$Vp,kd,ks,ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop,lt=lt,rt=lt,nthreads=n.threads)^.5
      }
      if (type=="response") {
        linkinv <- object$family$linkinv
        dmu.deta <- object$family$mu.eta
      } else linkinv <- dmu.deta <- NULL
      if (se==TRUE) {
        if (type=="response") {
	  if (nlp>1) for (i in 1:nlp) {
	    se.fit[,i] <- se.fit[,i] * abs(object$family$linfo[[i]]$mu.eta[fit[,i]])
	    fit[,i] <- object$family$linfo[[i]]$linkinv[fit[,i]]
          } else {
            se.fit <- se.fit * abs(object$family$mu.eta(fit))
            fit <- object$family$linkinv(fit)
	  }  
        }
        fit <- list(fit=fit,se.fit=se.fit)
      } else if (type=="response") fit <- object$family$linkinv(fit)
    } else { ## family has its own response prediction code
      X <- list(Xd=Xd,kd=kd,ks=ks,ts=ts,dt=dt,v=object$dinfo$v,qc=object$dinfo$qc,drop=object$dinfo$drop,lpid=lpid)
      if (nlp>1) attr(X,"lpi") <- lpi
      ## NOTE: not set up for families needing response for prediction (e.g. cox.ph)
      fampred <- object$family$predict ## just eases debugging
      offs <- if (length(offset)==1) rep(offset,nrow(kd)) else offset ## avoid subsetting failure
      ffv <- fampred(object$family,se=se,y=response,X=X,beta=object$coefficients,
                             off=offs,Vb=object$Vp)  ## NOTE: offsets not handled
      fit <- ffv[[1]]
      if (se) fit <- list(fit=fit,se.fit =ffv[[2]])
    }
  }
  rn <- rownames(newdata)
  if (type=="lpmatrix") {
    colnames(fit) <- names(object$coefficients)
    rownames(fit) <- rn
    attr(fit,"model.offset") <- offset
    fit <- napredict(na.act,fit)
  } else {
     if (se) { 
      if (is.null(nrow(fit$fit))) {
        names(fit$fit) <- rn
        names(fit$se.fit) <- rn
        fit$fit <- napredict(na.act,fit$fit)
        fit$se.fit <- napredict(na.act,fit$se.fit) 
      } else { 
        rownames(fit$fit) <- rn
        rownames(fit$se.fit) <- rn
        fit$fit <- napredict(na.act,fit$fit)
        fit$se.fit <- napredict(na.act,fit$se.fit)
      }
    } else { 
      if (is.null(nrow(fit))) names(fit) <- rn else
      rownames(fit) <- rn
      fit <- napredict(na.act,fit)
    }
  }
  fit
} ## end predict.bamd 



tero <- function(sm) {
## te smooth spec re-order so that largest marginal is last.
  maxd <- 0
  ns <- length(sm$margin)
  for (i in 1:ns) if (sm$margin[[i]]$bs.dim>=maxd) {
    maxi <- i;maxd <- sm$margin[[i]]$bs.dim
  }
  if (maxi<ns) { ## re-ordering required
    ind <- 1:ns;ind[maxi] <- ns;ind[ns] <- maxi
    sm$margin <- sm$margin[ind]
    sm$fix <- sm$fix[ind]
    if (!is.null(sm$mc)) sm$mc <- sm$mc[ind]
    sm$term <- rep("",0)
    for (i in 1:ns) sm$term <- c(sm$term,sm$margin[[i]]$term)
    sm$label <- paste0(substr(sm$label,1,3),paste0(sm$term,collapse=","),")",collapse="")
  }
  sm
} ## tero

AR.resid <- function(rsd,rho=0,AR.start=NULL) {
## standardised residuals for AR1 model
  if (rho==0) return(rsd)
  ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
  sd <- -rho*ld         ## sub diagonal
  N <- length(rsd)    
  ## see rwMatrix() for how following are used...
  ar.row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)]) ## index of rows to reweight
  ar.weight <- c(1,rep(c(sd,ld),N-1))     ## row weights
  ar.stop <- c(1,1:(N-1)*2+1)    ## (stop[i-1]+1):stop[i] are the rows to reweight to get ith row
  if (!is.null(AR.start)) { ## need to correct the start of new AR sections...
    ii <- which(AR.start==TRUE)
    if (length(ii)>0) {
          if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
          ar.weight[ii*2-2] <- 0 ## zero sub diagonal
          ar.weight[ii*2-1] <- 1 ## set leading diagonal to 1
    }
  }
  rwMatrix(ar.stop,ar.row,ar.weight,rsd)
} ## AR.resid

tens2matrix <- function(X,ts,dt) {
## converts result of terms2tensor to a full model matrix. Should
## match direct model.matrix call on the original terms object.
## Purely for checking purposes (insufficiently efficient for
## general use)
  Xt <- matrix(0,nrow(X[[1]]),0)
  for (i in 1:length(ts)) {
    ind <- ts[i]:(ts[i]+dt[i]-1)
    Xt <- if (dt[i]==1) cbind(Xt,X[[ind]]) else cbind(Xt,tensor.prod.model.matrix(X[ind]))
  }
  Xt
} ## tens2matrix


terms2tensor <- function(terms,data=NULL,contrasts.arg=NULL,drop.intercept=FALSE,identify=TRUE,sparse=FALSE) {
## takes a terms object or formula and converts it into a sequence of marginal model matrices, X
## and associated indices ts and dt, using the data in 'data'. drops the intercept if needed.
## If 'identify' then identifiability constraints/contrasts imposed, otherwise not. 
## The ith block of the model matrix is the row tensor product
## of the dt[i] elements of X starting at the ts[i].
## i.e. tensor.prod.model.matrix(X[ts[i]:(ts[i]+dt[i]-1)])
## xname[i] is the name of the variable in data used to create the marginal matrix X[[i]], with
## `(Intercept)' returned for the intercept term
## If data=NULL does a dummy run, not returning X or p
  if (!inherits(terms,"formula")) stop("requires a terms or formula object as first argument")
  if (!inherits(terms,"terms")) terms <- terms(terms)
  fac <- attr(terms,"factors")
  intercept <- attr(terms,"intercept")==1 ## is there an intercept?
  nt <- if (is.matrix(fac)) ncol(fac) else 0
  nt <- nt + intercept - drop.intercept ## number of terms.
  if (nt==0) return(NULL)
  p <- ts <- dt <- rep(1,nt)
  X <- list() ## marginal model matrix list
  form <- list() ## marginal formulae list
  k <- 1
  xname <- rep("",nt)
  dummy <- is.null(data)
  if (!dummy) n <- if (is.matrix(data[[1]])) nrow(data[[1]]) else length(data[[1]])
  if (intercept && !drop.intercept) {
    ts[k] <- 1;dt[k] <- 1
    if (!dummy) X[[k]] <- matrix(1,n,1)
    xname[k] <- "(Intercept)"
    form[[k]] <- ~1
    environment(form[[k]]) <- NULL 
    k <- k + 1
    term.labels <- c("(Intercept)",attr(terms,"term.labels"))
  } else term.labels <- attr(terms,"term.labels")
  varn <- rownames(fac)
  ord <- attr(terms,"order")
  no.int <- !intercept ## indicator of whether missing intercept still unhandled
  if (drop.intercept) intercept <- 0;
  if (nt-intercept>0) for (i in 1:(nt-intercept)) { ## work through the terms
    ts[i+intercept] <- k;dt[i+intercept] <- ord[i]
    if (ord[i]==1) { ## special case due to odd intercept handling
      vn <- varn[as.logical(fac[,i])]
      if (no.int||!identify) {
        fm <- as.formula(paste("~",vn,"-1"))
	## unfortunately a warning has been introduced in model.matrix that contradicts the documentation and will complain about any
	## element of contrast.arg that does not match a variable name in the formula. Rather than introduce redundant work around
	## code it seems better to just suppress the warning...
        if (!dummy) X[[k]] <- if (sparse) sparse.model.matrix(fm,data,contrasts.arg) else suppressWarnings(model.matrix(fm,data,contrasts.arg))
        no.int <- FALSE
      } else {
        fm <- as.formula(paste("~",vn))
        if (!dummy) X[[k]] <-
	if (sparse) sparse.model.matrix(fm,data,contrasts.arg)[,-1,drop=FALSE] else suppressWarnings(model.matrix(fm,data,contrasts.arg)[,-1,drop=FALSE])
      }
      xname[k] <- if (!dummy && vn %in% names(data)) vn else all.vars(fm)
      form[[k]] <- fm; environment(form[[k]]) <- NULL
      if (!dummy) p[i+intercept] = ncol(X[[k]])
      k <- k + 1
    } else { ## interaction term
      m <- which(fac[,i]!=0) ## marginal terms involved
      for (j in ord[i]:1) { ## reverse ordering is to conform with model.matrix(terms,data) called directly
        vn <- varn[m[j]]
        if (fac[m[j],i]==2||!identify) { ## no contrast
	  fm <- as.formula(paste("~",vn,"-1"))
	  if (!dummy) X[[k]] <- if (sparse) sparse.model.matrix(fm,data,contrasts.arg) else suppressWarnings(model.matrix(fm,data,contrasts.arg))
	} else { ## with contrast
          fm <- as.formula(paste("~",vn))
	  if (!dummy) X[[k]] <-
	  if (sparse) sparse.model.matrix(fm,data,contrasts.arg)[,-1,drop=FALSE] else suppressWarnings(model.matrix(fm,data,contrasts.arg)[,-1,drop=FALSE])
        }
	xname[k] <- if (!dummy && vn %in% names(data)) vn else all.vars(fm)
	form[[k]] <- fm; environment(form[[k]]) <- NULL
	if (!dummy) p[i+intercept] <- p[i+intercept] * ncol(X[[k]])
	k <- k + 1
      } ## j marginal loop end
    } ## finished interaction term
  } ## i term loop end
  list(X=X, ## list of marginal model matrices
       ts=ts, ## ts[i] is starting marginal matrix for ith term 
       dt=dt, ## dt[i] is number of marginal matrices for ith term
       xname=xname, ## xname[k] is name of covariate associated with X[[k]]
       form=form, ## form[[k]] is formula used to produce X[[k]]
       term.labels=term.labels, ## term.labels[i] is label of ith term
       p=p) ## p[i] is total coef count for ith term  
} # terms2tensor


bam <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action=na.omit,
                offset=NULL,method="fREML",control=list(),select=FALSE,scale=0,gamma=1,knots=NULL,sp=NULL,
                min.sp=NULL,paraPen=NULL,chunk.size=10000,rho=0,AR.start=NULL,discrete=FALSE,
                cluster=NULL,nthreads=1,gc.level=0,use.chol=FALSE,samfrac=1,coef=NULL,
                drop.unused.levels=TRUE,G=NULL,fit=TRUE,drop.intercept=NULL,in.out=NULL,nei=NULL,...) {

## Routine to fit an additive model to a large dataset. The model is stated in the formula, 
## which is then interpreted to figure out which bits relate to smooth terms and which to 
## parametric terms.
## This is a modification of `gam' designed to build the QR decomposition of the model matrix 
## up in chunks, to keep memory costs down.
## If cluster is a parallel package cluster uses parallel QR build on cluster. 
## 'n.threads' is number of threads to use for non-cluster computation (e.g. combining 
## results from cluster nodes). If 'NA' then is set to max(1,length(cluster)).
  control <- do.call("gam.control",control)
 
  if (control$trace) t3 <- t2 <- t1 <- t0 <- proc.time()
  if (length(nthreads)==1) nthreads <- rep(nthreads,2)
  if (is.null(G)) { ## need to set up model!
    weights.name <- substitute(weights) ## needed in case of update
    AR.sname <- substitute(AR.start) ## needed in case of update
    if (is.character(family))
            family <- eval(parse(text = family))
    if (is.function(family))
            family <- family()
    if (is.null(family$family))
            stop("family not recognized")
    family <- fix.family(family) ## apply any general family patches (e.g. gaussian log link initialization)  
    if (family$family=="gaussian"&&family$link=="identity") am <- TRUE else am <- FALSE
    if (scale==0) { if (family$family %in% c("poisson","binomial")) scale <- 1 else scale <- -1} 
    if (!method%in%c("fREML","GACV.Cp","GCV.Cp","REML",
                    "ML","P-REML","P-ML","NCV")) stop("un-supported smoothness selection method")
    if (is.logical(discrete)) {
      discretize <- discrete
      discrete <- NULL ## use default discretization, if any
    } else {
      discretize <- if (is.numeric(discrete)) TRUE else FALSE
    }
    if (method=="NCV") {
      if (!discretize) {
        discretize <- TRUE
	warning("NCV only available with discrete=TRUE - resetting")
      }
      if (rho!=0) {
        rho <- 0
	warning("AR1 residuals not available with NCV, resetting rho=0")
      }
    } ## NCV pre-processing
    if (discretize) { 
      if (!method %in% c("fREML","NCV")) { 
        discretize <- FALSE
        warning("discretization only available with fREML or NCV")
      } else {
        if (!is.null(cluster)) warning("discrete method does not use parallel cluster - use nthreads instead")
	if (all(is.finite(nthreads)) && any(nthreads>1) && !mgcv.omp()) warning("openMP not available: single threaded computation only")
      }
    }  
    if (inherits(family,"extended.family")) {
      family <- fix.family.link(family); efam <- TRUE
    } else efam <- FALSE
    
    if (method%in%c("fREML","NCV")&&!is.null(min.sp)) {
      min.sp <- NULL
      warning("min.sp not supported with fast REML and NCV computation, and ignored.")
    }
   
    gp <- interpret.gam(formula) # interpret the formula
    if (discretize && length(gp$smooth.spec)==0) {
      ok <- TRUE
      ## check it's not a list formula
      if (!is.null(gp$nlp)) for (i in 1:gp$nlp) if (length(gp[[i]]$smooth.spec)>0) ok <- FALSE
      if (ok) {
        warning("no smooths, ignoring `discrete=TRUE'")
        discretize <- FALSE
      }	
    }
    if (discretize) { 
      ## re-order the tensor terms for maximum efficiency, and 
      ## signal that "re"/"fs" terms should be constructed with marginals
      ## also for efficiency
      
      if (is.null(gp$nlp)) for (i in 1:length(gp$smooth.spec)) { 
        if (inherits(gp$smooth.spec[[i]],"tensor.smooth.spec")) gp$smooth.spec[[i]] <- tero(gp$smooth.spec[[i]])
        #if (inherits(gp$smooth.spec[[i]],c("re.smooth.spec","fs.smooth.spec"))&&gp$smooth.spec[[i]]$dim>1)
	 if (!is.null(gp$smooth.spec[[i]]$tensor.possible)&&gp$smooth.spec[[i]]$dim>1){
          class(gp$smooth.spec[[i]]) <- c(class(gp$smooth.spec[[i]]),"tensor.smooth.spec")
	  if (is.null(gp$smooth.spec[[i]]$margin)) {
            gp$smooth.spec[[i]]$margin <- list()
            ## only ok for 'fs' with univariate metric variable (caught in 'fs' construcor)...
            for (j in 1:gp$smooth.spec[[i]]$dim) gp$smooth.spec[[i]]$margin[[j]] <- list(term=gp$smooth.spec[[i]]$term[j])
	  }  
        }
      } else for (j in 1:length(formula)) if (length(gp[[j]]$smooth.spec)>0) for (i in 1:length(gp[[j]]$smooth.spec)) {
        if (inherits(gp[[j]]$smooth.spec[[i]],"tensor.smooth.spec")) gp[[j]]$smooth.spec[[i]] <- tero(gp[[j]]$smooth.spec[[i]])
        #if (inherits(gp[[j]]$smooth.spec[[i]],c("re.smooth.spec","fs.smooth.spec"))&&gp[[j]]$smooth.spec[[i]]$dim>1)
	if (!is.null(gp[[j]]$smooth.spec[[i]]$tensor.possible)&&gp[[j]]$smooth.spec[[i]]$dim>1) {
          class(gp[[j]]$smooth.spec[[i]]) <- c(class(gp[[j]]$smooth.spec[[i]]),"tensor.smooth.spec")
          if (is.null(gp[[j]]$smooth.spec[[i]]$margin)) {
	    gp[[j]]$smooth.spec[[i]]$margin <- list()
            ## only ok for 'fs' with univariate metric variable (caught in 'fs' construcor)...
            for (k in 1:gp[[j]]$smooth.spec[[i]]$dim) gp[[j]]$smooth.spec[[i]]$margin[[k]] <- list(term=gp[[j]]$smooth.spec[[i]]$term[k])
          }
        }
      }
    } ## if (discretize)
    cl <- match.call() # call needed in gam object for update to work
    mf <- match.call(expand.dots=FALSE)
    mf$formula <- gp$fake.formula 
    mf$method <-  mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp <- mf$gc.level <-
    mf$gamma <- mf$paraPen<- mf$chunk.size <- mf$rho  <- mf$cluster <- mf$discrete <-
    mf$use.chol <- mf$samfrac <- mf$nthreads <- mf$G <- mf$fit <- mf$select <- mf$drop.intercept <-
    mf$coef <- mf$in.out <- mf$nei <- mf$... <-NULL
    mf$drop.unused.levels <- drop.unused.levels
    mf[[1]] <- quote(stats::model.frame) ## as.name("model.frame")

    if (is.list(formula)) { ## then there are several linear predictors
      environment(formula) <- environment(formula[[1]]) ## e.g. termplots needs this
      pterms <- list()
      tlab <- rep("",0)
      pmf.names <- rep("",0)
      for (i in 1:length(formula)) {
        pmf <- mf
        pmf$formula <- gp[[i]]$pf
	pmf <- eval(pmf, parent.frame())
	pmf.names <- c(pmf.names,names(pmf))
        pterms[[i]] <- attr(pmf,"terms")
        tlabi <- attr(pterms[[i]],"term.labels")
        if (i>1&&length(tlabi)>0) tlabi <- paste(tlabi,i-1,sep=".")
        tlab <- c(tlab,tlabi)
      }
      pmf.names <- unique(pmf.names)
      attr(pterms,"term.labels") <- tlab ## labels for all parametric terms, distinguished by predictor
      nlp <- gp$nlp
      lpid <- list() ## list of terms for each lp
      lpid[[nlp]] <- rep(0,0)
    } else { ## single linear predictor case
      nlp <- 1
      pmf <- mf
      pmf$formula <- gp$pf
      pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
      pterms <- attr(pmf,"terms") ## pmf only used for this and discretization, if selected.
      pmf.names <- names(pmf)
    }
    
    if (gc.level>0) gc()

    mf <- eval(mf, parent.frame()) # the model frame now contains all the data

    if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf,"terms")
    if (gc.level>0) gc()  
    if (rho!=0&&!is.null(mf$"(AR.start)")) if (!is.logical(mf$"(AR.start)")) stop("AR.start must be logical")

    if (!is.null(nei)) { ## check if data dropped
      k <- attr(mf,"na.action")
      if (!is.null(k)) { ## need to adjust nei for dropped data
        nei <- nanei(nei,as.numeric(k))
      }
    }

    if (method=="NCV") { ## pre-process the neighbourhood structure 'nei'
      n <- nrow(mf)
      if (is.null(nei)||is.null(nei$a)||is.null(nei$ma)) {
        nei <- list(a=1:n,ma=1:n,d=1:n,md=1:n) ## LOOCV  
      } else if (is.null(nei$d)||is.null(nei$md)) {
        if (length(nei$ma)!=n) stop("nei$d and nei$md must be supplied if number of neighbourhoods is not number of data")
	nei$d <- 1:n; nei$md <- 1:n
      }
      if (max(nei$ma)>length(nei$a)||min(nei$ma)<1) stop("nei$ma does not match nei$a")
      if (max(nei$md)>length(nei$d)||min(nei$md)<1) stop("nei$md does not match nei$d")
      if (max(nei$d)>n||min(nei$d)<1||max(nei$a)>n||min(nei$a)<1) stop("supplied nei index out of range")
      nei <- onei(nei,TRUE) ## order indices within neighbourhoods and chnage to C indexing - needed for discrete NCV C code
    } ## nei pre-processing
    
    ## summarize the *raw* input variables
    ## note can't use get_all_vars here -- buggy with matrices
    vars <- all_vars1(gp$fake.formula[-2]) ## drop response here
    inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

    ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
    if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data) 

    dl <- eval(inp, data, parent.frame())
    if (!control$keepData) { rm(data);if (gc.level>0) gc()} ## save space
    names(dl) <- vars ## list of all variables needed
    var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
    rm(dl); if (gc.level>0) gc() ## save space    

    ## should we force the intercept to be dropped, meaning that the constant is removed
    ## from the span of the parametric effects?
    if (is.null(family$drop.intercept)) { ## family does not provide information
      lengthf <- if (is.list(formula)) length(formula) else 1
      if (is.null(drop.intercept)) drop.intercept <- rep(FALSE,lengthf) else {
        drop.intercept <- rep(drop.intercept,length=lengthf) ## force drop.intercept to correct length
	if (sum(drop.intercept)) family$drop.intercept <- drop.intercept ## ensure prediction works
      }
    } else drop.intercept <- as.logical(family$drop.intercept) ## family overrides argument
 

    ## need mini.mf for basis setup, then accumulate full X, y, w and offset
    if (discretize) {
      ## discretize the data, creating list mf0 with discrete values
      ## and indices giving the discretized value for each element of model frame.
      ## 'discrete' can be null, or contain a discretization size, or
      ## a discretization size per smooth term.   
      dk <- discrete.mf(gp,mf,pmf.names,m=discrete)
      mf0 <- dk$mf ## padded discretized model frame
      sparse.cons <- 0 ## default constraints required for tensor terms

    } else { 
      mf0 <- mini.mf(mf,chunk.size)
      sparse.cons <- -1
    }
    rm(pmf); ## no further use

    ## allow bam to set up general families, even if it can not (yet) process them...
    if (inherits(family,"general.family")&&!is.null(family$presetup)) eval(family$presetup)
    gsname <- if (is.list(formula)) "gam.setup.list" else "gam.setup"
    
    if (control$trace) t1 <- proc.time()
    reset <- TRUE
    while (reset) {
     G <- do.call(gsname,list(formula=gp,pterms=pterms,
                  data=mf0,knots=knots,sp=sp,min.sp=min.sp,
                  H=NULL,absorb.cons=TRUE,sparse.cons=sparse.cons,select=select,
                  idLinksBases=!discretize,scale.penalty=control$scalePenalty,
                  paraPen=paraPen,apply.by=!discretize,drop.intercept=drop.intercept,modCon=2))

    if (!discretize&&ncol(G$X)>=chunk.size) { ## no point having chunk.size < p
        chunk.size <- 4*ncol(G$X)
        warning(gettextf("chunk.size < number of coefficients. Reset to %d",chunk.size))
        if (chunk.size>=nrow(mf)) { ## no sense splitting up computation
          mf0 <- mf ## just use full dataset
        } else reset <- FALSE
      } else reset <- FALSE
    }
    if (control$trace) t2 <- proc.time()
    if (discretize) {
      ks <- matrix(0,0,2) ## NOTE: slightly more efficient not to repeatedly extend
      if (nlp>1) lpi <- attr(G$X,"lpi")
      v <- G$Xd <- list()
      kb <- k <- 1 ## kb indexes blocks, k indexes marginal matrices
      G$kd <- dk$k
      qc <- dt <- ts <- rep(0,length(G$smooth))
      ## have to extract full parametric model matrix from pterms and mf
      npt <- if (nlp==1) 1 else length(G$pterms)
      lpip <- list() ## record coef indices for each discretized term
      for (j in 1:npt) { ## loop over parametric terms in each formula
        paratens <- TRUE
        ## get the parametric model split into tensor components...
        ptens <- if (nlp==1&&!is.list(G$pterms)) terms2tensor(G$pterms,mf0,drop.intercept=drop.intercept[j]) else
	             terms2tensor(G$pterms[[j]],mf0,drop.intercept=drop.intercept[j])
        if (j>1) ptens$term.labels <- paste(ptens$term.labels,".",j-1,sep="")		     
        ## now locate the index vectors for each parametric marginal discrete matrix ptens$X		     
        if (!is.null(ptens)) {
          np <-length(ptens$X);n <- nrow(mf)
	  qc <- c(rep(0,length(ptens$ts)),qc) ## extend (empty) constraint indicator
	  coef.ind <-1
	  kk <- 1
          for (i in 1:length(ptens$ts)) {
            jj <- 1
	    ts[kb] = k;names(ts)[kb] <- ptens$term.labels[i]
	    dt[kb] = ptens$dt[i]
            for (ii in 1:dt[kb]) {
	      ks <- rbind(ks,dk$ks[ptens$xname[kk],])
	      G$Xd[[k]] <- ptens$X[[kk]][1:dk$nr[ptens$xname[kk]],,drop=FALSE];
	      kk <- kk + 1
	      jj <- jj * ncol(G$Xd[[k]]) ## number of coeffs for this term
	      k <- k + 1 ## update matrix counter
	    }  
            lpip[[kb]] <- 1:jj - 1 + if (nlp==1) coef.ind else attr(G$nsdf,"pstart")[j] - 1 + coef.ind
	    coef.ind <- coef.ind + jj
	    if (nlp>1) for (ii in 1:length(lpi)) if (any(lpip[[kb]]%in%lpi[[ii]])) lpid[[ii]] <- c(lpid[[ii]],kb)	 
            kb <- kb + 1 ## update block counter
	  }
	} ## is.null(ptens)  
      } ## loop over parametric terms in each formula	
      ## k is marginal counter, kb is block counter
      ## G$kd[,ks[j,1]:ks[j,2]] (dk$k) gives index columns for term j, thereby allowing 
      ## summation over matrix covariates....
      #G$ks <- cbind(dk$k.start[-length(dk$k.start)],dk$k.start[-1])
     
      
      drop <- rep(0,0) ## index of te related columns to drop
      if (length(G$smooth)>0) for (i in 1:length(G$smooth)) { ## loop over smooths
        ## potentially each smoother model matrix can be made up of a sequence
	## of row-tensor products, nead to loop over such sub blocks...
        nsub <- if (!is.null(G$smooth[[i]]$ts)) length(G$smooth[[i]]$ts) else 1
	lp0 <- G$smooth[[i]]$first.para -1  ## offset for start of coeffs for this sub block 
 	for (sb in 1:nsub) { ## loop over sub-blocks
	  np <- 1 ## compute number of sub-block coeffs 
          ts[kb] <- k;names(ts)[kb] <- G$smooth[[i]]$label
          ## first deal with any by variable (as first marginal of tensor)...
          if (G$smooth[[i]]$by!="NA") {
            dt[kb] <- 1
	    termk <- G$smooth[[i]]$by
            by.var <- dk$mf[[termk]][1:dk$nr[termk]]
            if (is.factor(by.var)) { 
              ## create dummy by variable...
              by.var <- as.numeric(by.var==G$smooth[[i]]$by.level)  
            }
            G$Xd[[k]] <- matrix(by.var,dk$nr[termk],1)
	    np <- ncol(G$Xd[[k]])
	    ks <- rbind(ks,dk$ks[termk,])
            k <- k + 1
	    by.present <- 1
          } else by.present <- dt[kb] <- 0
          ## ... by done
          if (inherits(G$smooth[[i]],"tensor.smooth")) { 
            nmar <- if (is.null(G$smooth[[i]]$dt)) length(G$smooth[[i]]$margin) else G$smooth[[i]]$dt[sb] 
            dt[kb] <- dt[kb] + nmar
 
	    jind <- if (sb>1) G$smooth[[i]]$ts[sb] + 1:G$smooth[[i]]$dt[sb] - 1 else 1:nmar
            for (j in jind) {
	      termk <- G$smooth[[i]]$margin[[j]]$term[1]
              G$Xd[[k]] <- G$smooth[[i]]$margin[[j]]$X[1:dk$nr[termk],,drop=FALSE]
	      np <- np * ncol(G$Xd[[k]])
	      ks <- rbind(ks,dk$ks[termk,])
              k <- k + 1 
            }
            ## deal with any side constraints on tensor terms
	    if (sb==1) { ## only once per smooth!
              di <- attr(G$smooth[[i]],"del.index")
              if (!is.null(di)&&length(di>0)) {
                di <- di + G$smooth[[i]]$first.para + length(drop)  - 1 
                drop <- c(drop,di)
              }
	      
              ## deal with tensor smooth constraint
              qrc <- attr(G$smooth[[i]],"qrc")
              ## compute v such that Q = I-vv' and Q[,-1] is constraint null space basis
              if (inherits(qrc,"qr")) {
                v[[kb]] <- qrc$qr/sqrt(qrc$qraux);v[[kb]][1] <- sqrt(qrc$qraux)
                qc[kb] <- 1 ## indicate a constraint
              } else if (length(qrc)>1) { ## Kronecker product of set to zero constraints
	        ## on entry qrc is [unused.index, dim1, dim2,..., total number of constraints]
                v[[kb]] <- c(length(qrc)-2,qrc[-1]) ## number of sum-to-zero contrasts, their dimensions, number of constraints
		qc[kb] <- -1 
              } else { 
                v[[kb]] <- rep(0,0) ##
                if (!inherits(qrc,"character")||qrc!="no constraints") warning("unknown tensor constraint type")
              }
	    } else { ## sb==1 once per smooth stuff
              qc <- c(qc,0) ## extend
	      v[[kb]] <- rep(0,0)
            }
          } else { ## not a tensor smooth
            v[[kb]] <- rep(0,0)
            dt[kb] <- dt[kb] + 1
	    termk <- G$smooth[[i]]$term[1]
            G$Xd[[k]] <- G$X[1:dk$nr[termk],G$smooth[[i]]$first.para:G$smooth[[i]]$last.para,drop=FALSE]
	    np <- np * ncol(G$Xd[[k]])
            ks <- rbind(ks,dk$ks[termk,])
            k <- k + 1
          }  
	  #jj <- G$smooth[[i]]$first.para:G$smooth[[i]]$last.para;
	  if (sb==1&&qc[kb]) {
	    ncon <- if (qc[kb]<0) v[[kb]][length(v[[kb]])] else 1 
            jj <- 1:(np-ncon) + lp0; lp0 <- lp0 + np - ncon 
	    ## Hard to think of an application requiring constraint when nsub>1, hence not 
	    ## worked out yet. Add warning to make sure this is flagged if attempt made
	    ## to do this in future....
	    if (nsub>1) warning("constraints for sub blocked tensor products un-tested")
          } else {
	    jj <- 1:np + lp0; lp0 <- lp0 + np
	  }  
	  lpip[[kb]] <- jj
	  if (nlp>1) { ## record which lp each discrete term belongs to (can be more than one)
            for (j in 1:nlp) if (any(jj %in% lpi[[j]])) lpid[[j]] <- c(lpid[[j]],kb)
          }
          kb <- kb + 1
	} ## sub block loop  
      } ## looping over smooths
      ## put lpid indices into coefficient index order...
      if (nlp>1) {
        for (j in 1:nlp) lpid[[j]] <- lpid[[j]][order(unlist(lapply(lpip[lpid[[j]]],max)))]
	G$lpid <- lpid
      }	
      if (length(drop>0)) G$drop <- drop ## index of terms to drop as a result of side cons on tensor terms
      attr(G$Xd,"lpip") <- lpip ## index of coefs by term
      ## ... Xd is the list of discretized model matrices, or marginal model matrices
      ## kd contains indexing vectors, so the ith model matrix or margin is Xd[[i]][kd[i,],]
      ## ts[i] is the starting matrix in Xd for the ith model matrix, while dt[i] is the number 
      ## of elements of Xd that make it up (1 for a singleton, more for a tensor). 
      ## v is list of Householder vectors encoding constraints and qc the constraint indicator.
      G$v <- v;G$ts <- ts;G$dt <- dt;G$qc <- qc
    
      G$ks <- ks
      jj <- max(G$ks)-1
      if (ncol(G$kd) > jj) G$kd <- G$kd[,1:jj]
      ## bundle up discretization information needed for discrete prediction...
      G$dinfo <- list(gp=gp, v = G$v, ts = G$ts, dt = G$dt, qc = G$qc, drop = G$drop, pmf.names=pmf.names,lpip=lpip)
      if (nlp>1) G$dinfo$lpid <- lpid
      if (paratens) G$dinfo$para.discrete <- TRUE 
    } ## if (discretize)

    if (control$trace) t3 <- proc.time()

    ## no advantage to "fREML" with no free smooths...
    if (((!is.null(G$L)&&ncol(G$L) < 1)||(length(G$sp)==0))&&method=="fREML") method <- "REML"

    G$var.summary <- var.summary 
    G$family <- family
    G$terms<-terms;
    G$pred.formula <- gp$pred.formula

    n <- nrow(mf)
  
    if (is.null(mf$"(weights)")) G$w<-rep(1,n)
    else G$w<-mf$"(weights)"    

    G$y <- mf[[gp$response]]
    ## now get offset, dealing with possibility of multiple predictors (see gam.setup)
    ## the point is that G$offset relates to the compressed or discretized model frame,
    ## so we need to correct it to the full data version...
    if (discretize) {
      if (is.list(pterms)) { ## multiple predictors
        for (i in 1:length(pterms)) {
	  offi <- attr(pterms[[i]],"offset")
	  if (is.null(offi)) G$offset[[i]] <- rep(0,n) else {
	    G$offset[[i]] <- mf[[names(attr(pterms[[i]],"dataClasses"))[offi]]]
	    if (is.null(G$offset[[i]])) G$offset[[i]] <- rep(0,n)
	  }
	}
      } else { ## single predictor, handle as non-discrete
        G$offset <- model.offset(mf)
	if (is.null(G$offset)) G$offset <- rep(0,n)
      }
    } else { ## non-discrete
      G$offset <- model.offset(mf)
      if (is.null(G$offset)) G$offset <- rep(0,n)
    }
     
    if (ncol(G$X) > chunk.size && !discretize) { ## no sense having chunk.size < p
      chunk.size <- 4*ncol(G$X)
      warning(gettextf("chunk.size < number of coefficients. Reset to %d",chunk.size))    
    }

    G$cl <- cl
    G$am <- am
     
    G$min.edf<-G$nsdf #-dim(G$C)[1]
    if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim
    G$discretize <- discretize
    G$formula<-formula
    G$weights.name <- weights.name; G$AR.sname <- AR.sname
    ## environment(G$formula)<-environment(formula)
    environment(G$pterms) <- environment(G$terms) <- environment(G$pred.formula) <- 
    environment(G$formula) <- .BaseNamespaceEnv

  } else { ## G supplied
    if (scale<=0) scale <- G$scale
    efam <- G$efam
    mf <- G$mf; G$mf <- NULL
    gp <- G$gp; G$gp <- NULL
    na.action <- G$na.action; G$na.action <- NULL
    if (!is.null(sp)&&any(sp>=0)) { ## request to modify smoothing parameters
      if (is.null(G$L)) G$L <- diag(length(G$sp))
      if (length(sp)!=ncol(G$L)) stop('length of sp must be number of free smoothing parameters in original model')
      ind <- sp>=0 ## which smoothing parameters are now fixed
      spind <- log(sp[ind]); 
      spind[!is.finite(spind)] <- -30 ## set any zero parameters to effective zero
      G$lsp0 <- G$lsp0 + drop(G$L[,ind,drop=FALSE] %*% spind) ## add fix to lsp0
      G$L <- G$L[,!ind,drop=FALSE] ## drop the cols of G
      G$sp <- rep(-1,ncol(G$L))
    }
  } ## end of G setup 

  if (!fit) {
    G$efam <- efam
    G$scale <- scale
    G$mf <- mf;G$na.action <- na.action;G$gp <- gp
    class(G) <- "bam.prefit"
    return(G)
  }

  if (inherits(G$family,"general.family")) stop("general families not supported by bam")

  ## number of threads to use for non-cluster node computation
  if (!is.finite(nthreads[1])||nthreads[1]<1) nthreads[1] <- max(1,length(cluster))

  G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
  G$max.half<-control$mgcv.half     # max step halving in bfgs optimization


  ## now build up proper model matrix, and deal with y, w, and offset...

  if (control$trace) cat("Setup complete. Calling fit\n")
  
  colnamesX <- colnames(G$X)  

  if (G$am&&!G$discretize) {
    if (nrow(mf)>chunk.size) G$X <- matrix(0,0,ncol(G$X)); if (gc.level>1) gc() 
    object <- bam.fit(G,mf,chunk.size,gp,scale,gamma,method,rho=rho,cl=cluster,
                      gc.level=gc.level,use.chol=use.chol,in.out=in.out,npt=nthreads[1])
  } else if (G$discretize) {
    object <- bgam.fitd(G, mf, gp ,scale ,nobs.extra=0,rho=rho,coef=coef,
                       control = control,npt=nthreads,gc.level=gc.level,
		       gamma=gamma,in.out=in.out,method=method,nei=nei,...)
                       
  } else {
    G$X  <- matrix(0,0,ncol(G$X)); if (gc.level>1) gc()
    if (rho!=0) warning("AR1 parameter rho unused with generalized model")
    if (samfrac<1 && samfrac>0) { ## sub-sample first to get close to right answer...
      ind <- sample(1:nrow(mf),ceiling(nrow(mf)*samfrac))
      if (length(ind)<2*ncol(G$X)) warning("samfrac too small - ignored") else {
        Gw <- G$w;Goffset <- G$offset
        G$w <- G$w[ind];G$offset <- G$offset[ind]
        control1 <- control
        control1$epsilon <- 1e-2
        object <- bgam.fit(G, mf[ind,], chunk.size, gp ,scale ,gamma,method=method,nobs.extra=0,
                       control = control1,cl=cluster,npt=nthreads[1],gc.level=gc.level,coef=coef,
                       use.chol=use.chol,samfrac=1,in.out=in.out,...)
        G$w <- Gw;G$offset <- Goffset
        coef <- object$coefficients
      }
    }
    ## fit full dataset
    object <- bgam.fit(G, mf, chunk.size, gp ,scale ,gamma,method=method,coef=coef,
                       control = control,cl=cluster,npt=nthreads[1],gc.level=gc.level,
                       use.chol=use.chol,in.out=in.out,...)
  }

  if (gc.level>0) gc()

  if (control$trace) t4 <- proc.time()

  if (control$trace) cat("Fit complete. Finishing gam object.\n")

  if (scale < 0) { object$scale.estimated <- TRUE;object$scale <- object$scale.est} else {
    object$scale.estimated <- FALSE; object$scale <- scale
  }

  object$assign <- G$assign # applies only to pterms  
  object$boundary <- FALSE  # always FALSE for this case
  object$call<-G$cl # needed for update() to work 
  object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
  object$weights.name <- G$weights.name
  object$AR.sname <- G$AR.sname

  object$contrasts <- G$contrasts
  object$control <- control
  object$converged <- TRUE ## no iteration
  object$data <- NA ## not saving it in this case
  object$df.null <- nrow(mf)
  object$df.residual <- object$df.null - sum(object$edf) 
 
  if (is.null(object$family)) object$family <- family
  object$formula <- G$formula 
 
  if (method=="GCV.Cp") {
    if (scale<=0) object$method <- "GCV" else object$method <- "UBRE"
  } else {
    object$method <- method
  }
  object$min.edf<-G$min.edf
  object$model <- mf;rm(mf);if (gc.level>0) gc()
  object$na.action <- attr(object$model,"na.action") # how to deal with NA's
  object$nsdf <- G$nsdf
  if (G$nsdf>0) names(object$coefficients)[1:G$nsdf] <- colnamesX[1:G$nsdf]
  object$offset <- G$offset
  ##object$prior.weights <- G$w
  object$pterms <- G$pterms
  object$pred.formula <- G$pred.formula 
  object$smooth <- G$smooth

  object$terms <- G$terms
  object$var.summary <- G$var.summary 
  if (is.null(object$wt)) object$weights <- object$prior.weights else
  object$weights <- object$wt
  object$xlevels <- G$xlevels
  #object$y <- object$model[[gp$response]]
  object$NA.action <- na.action ## version to use in bam.update
  names(object$sp) <- names(G$sp)
  if (!is.null(object$full.sp)) names(object$full.sp) <- names(G$lsp0)

  names(object$coefficients) <- G$term.names
  names(object$edf) <- G$term.names

  ## note that predict.gam assumes that it must be ok not to split the 
  ## model frame, if no new data supplied, so need to supply explicitly
  class(object) <- c("bam","gam","glm","lm")
  if (!G$discretize) { object$linear.predictors <- 
          as.numeric(predict.bam(object,newdata=object$model,block.size=chunk.size,cluster=cluster))
    dr <- dim(object$R) ## R can have fewer rows than columns if model rank def without penalization 
    if (dr[1]<dr[2]) object$R <- rbind(object$R,matrix(0,dr[2]-dr[1],dr[2]))
  } else { ## store discretization specific information to help with discrete prediction
    object$dinfo <- G$dinfo
  } 
  rm(G);if (gc.level>0) gc()

  if (is.null(object$fitted.values)) object$fitted.values <- family$linkinv(object$linear.predictors)
   
  object$residuals <- if (is.null(family$residuals)) sqrt(family$dev.resids(object$y,object$fitted.values,object$prior.weights)) * 
                      sign(object$y-object$fitted.values) else residuals(object)
  if (rho!=0) object$std.rsd <- AR.resid(object$residuals,rho,object$model$"(AR.start)")

  if (!efam || is.null(object$deviance)) object$deviance <- sum(object$residuals^2)
  ## 'dev' is used in family$aic to estimate scale. That's standard and fine for Gaussian data, but
  ## can lead to badly biased estimates for e.g. low count data with the Tweedie (see Fletcher Biometrika paper)
  ## So set dev to give object $sig2 estimate when divided by sum(prior.weights)... 
  dev <- if (family$family!="gaussian"&&!is.null(object$sig2)) object$sig2*sum(object$prior.weights) else object$deviance ## used to give scale in family$aic
  if (rho!=0&&family$family=="gaussian") dev <- sum(object$std.rsd^2)
  object$aic <- if (efam) family$aic(object$y,object$fitted.values,family$getTheta(),object$prior.weights,dev) else
                family$aic(object$y,1,object$fitted.values,object$prior.weights,dev)
  object$aic <- object$aic -
                2 * (length(object$y) - sum(sum(object$model[["(AR.start)"]])))*log(1/sqrt(1-rho^2)) + ## correction for AR
                2*sum(object$edf)
  if (!is.null(object$edf2)&&sum(object$edf2)>sum(object$edf1)) object$edf2 <- object$edf1
  if (is.null(object$null.deviance)) object$null.deviance <- sum(family$dev.resids(object$y,weighted.mean(object$y,object$prior.weights),object$prior.weights))
  if (!is.null(object$full.sp)) {
    if (length(object$full.sp)==length(object$sp)&&
        all.equal(object$sp,object$full.sp)==TRUE) object$full.sp <- NULL
  }
  environment(object$formula) <- environment(object$pred.formula) <-
  environment(object$terms) <- environment(object$pterms) <- 
  environment(attr(object$model,"terms"))  <- .GlobalEnv  
  if (control$trace) { 
    t5 <- proc.time()
    t5 <- rbind(t1-t0,t2-t1,t3-t2,t4-t3,t5-t4)[,1:3]
    row.names(t5) <- c("initial","gam.setup","pre-fit","fit","finalise")
    print(t5)
  }
  names(object$gcv.ubre) <- method
  object
} ## end of bam


bam.update <- function(b,data,chunk.size=10000) {
## update the strictly additive gaussian model `b' in the light of new data in `data'
## Need to update modelframe (b$model) 
  if (is.null(b$qrx)) { 
    stop("Model can not be updated")
  }
  gp <- interpret.gam(b$formula) # interpret the formula 

  ## next 2 lines problematic if there are missings in the response, so now constructed from mf below...
  ## X <- predict(b,newdata=data,type="lpmatrix",na.action=b$NA.action) ## extra part of model matrix
  ## rownames(X) <- NULL
  cnames <- names(b$coefficients)

  AR.start <- NULL ## keep R checks happy 

  ## now get the new data in model frame form...
  getw <- "(weights)"%in%names(b$model)
  getARs <- "(AR.start)"%in%names(b$model)
 
  if (getw||getARs) {
    mfstr <- "mf <- model.frame(gp$fake.formula,data,"
    if (getw) mfstr <- paste(mfstr,"weights=",b$weights.name,",",sep="")
    if (getARs) mfstr <- paste(mfstr,"AR.start=",b$AR.sname,",",sep="")
    mfstr <- paste(mfstr,"xlev=b$xlev,na.action=b$NA.action)",sep="")
    eval(parse(text=mfstr))
    w <- if (getw) mf[["(weights)"]] else  rep(1,nrow(mf))
  } else {
    mf <- model.frame(gp$fake.formula,data,xlev=b$xlev,na.action=b$NA.action)
    w <- rep(1,nrow(mf))
  }

  X <- predict(b,newdata=mf,type="lpmatrix",na.action=b$NA.action) ## extra part of model matrix
  rownames(X) <- NULL

  b$model <- rbind(b$model,mf) ## complete model frame --- old + new

  ## get response and offset...

  off.col <- attr(attr(b$model,"terms"),"offset")
  if (is.null(off.col)) offset <- rep(0,nrow(mf)) else offset <-  mf[,off.col]
  y <-  mf[,attr(attr(b$model,"terms"),"response")] - offset
  
  ## update G
  b$G$y <- c(b$G$y,y)
  b$G$offset <- c(b$G$offset,offset)
  b$G$w <- c(b$G$w,w)
  b$G$n <- length(b$G$y) ## nrow(b$model) - allow user to drop b$model...
  n <- b$G$n;
  ## update the qr decomposition...

  w <- sqrt(w)

  if (b$AR1.rho!=0) { ## original model had AR1 error structure...
    rho <- b$AR1.rho
    ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
    sd <- -rho*ld         ## sub diagonal
    ## append the final row of weighted X and y from original fit, first
    wy <- c(b$yX.last[1],w*y)
    wX <- rbind(b$yX.last[-1],w*X)
    m <- nrow(wX)
    b$yX.last <- c(wy[m],wX[m,])

    row <- c(1,rep(1:m,rep(2,m))[-c(1,2*m)])
    weight <- c(1,rep(c(sd,ld),m-1))
    stop <- c(1,1:(m-1)*2+1)
    if (!is.null(mf$"(AR.start)")) { ## need to correct the start of new AR sections...
         ii <- which(mf$"(AR.start)"==TRUE)
         if (length(ii)>0) {
           if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
           weight[ii*2-2] <- 0 ## zero sub diagonal
           weight[ii*2-1] <- 1 ## set leading diagonal to 1
         }
    }
   
    ## re-weight to independence....
    wX <- rwMatrix(stop,row,weight,wX)[-1,]
    wy <- rwMatrix(stop,row,weight,wy)[-1]    

    ## update
    b$qrx <- qr_update(wX,wy,b$qrx$R,b$qrx$f,b$qrx$y.norm2)
  } else {
    b$qrx <- qr_update(w*X,w*y,b$qrx$R,b$qrx$f,b$qrx$y.norm2)
  }

  ## now do the refit...
  rss.extra <- b$qrx$y.norm2 - sum(b$qrx$f^2)

  if (b$method=="GCV"||b$method=="UBRE") method <- "GCV.Cp" else method <- b$method

 
  if (method=="GCV.Cp") {
    if (b$method=="GCV") scale <- -1 else scale = b$sig2
   
    fit <- magic(b$qrx$f,b$qrx$R,b$sp,b$G$S,b$G$off,L=b$G$L,lsp0=b$G$lsp0,rank=b$G$rank,
               H=b$G$H,C= matrix(0,0,ncol(b$qrx$R)),##C=b$G$C,
               gamma=b$gamma,scale=scale,gcv=(scale<=0),
               extra.rss=rss.extra,n.score=n)
 
    post <- magic.post.proc(b$qrx$R,fit,b$qrx$f*0+1) 
    b$y <- b$G$y;b$offset <- b$G$offset; b$G$w -> b$weights -> b$prior.weights;
    
  } else if (method=="fREML") { ## fast REML

     um <- Sl.Xprep(b$Sl,b$qrx$R)
     lsp0 <- log(b$sp) ## initial s.p.
     log.phi <- log(b$sig2) ## initial or fixed scale
     fit <- fast.REML.fit(um$Sl,um$X,b$qrx$f,rho=lsp0,L=b$G$L,rho.0=b$G$lsp0,
            log.phi=log.phi,phi.fixed = !b$scale.estimated,rss.extra=rss.extra,
            nobs =n,Mp=um$Mp,nt=1,gamma=b$gamma)
     if (b$scale.estimated) scale <- -1 else scale=b$sig2
     res <- Sl.postproc(b$Sl,fit,um$undrop,b$qrx$R,cov=TRUE,scale=scale,L=b$g$L)


     object <- list(coefficients=res$beta,edf=res$edf,edf1=res$edf1,edf2=res$edf2,##F=res$F,
                    gcv.ubre=fit$reml,hat=res$hat,outer.info=list(iter=fit$iter,
                    message=fit$conv),optimizer="fast-REML",rank=ncol(um$X),
                    Ve=res$Ve,Vp=res$Vp,Vc=res$Vc,db.drho=fit$d1b,scale.estimated = scale<=0)
     if (scale<=0) { ## get sp's and scale estimate
       nsp <- length(fit$rho)
       object$sig2 <- object$scale <- exp(fit$rho[nsp])
       object$sp <- exp(fit$rho[-nsp]) 
       nsp <- length(fit$rho.full)
       object$full.sp <- exp(fit$rho.full[-nsp])
     } else { ## get sp's
       object$sig2 <- object$scale <- scale  
       object$sp <- exp(fit$rho)
       object$full.sp <- exp(fit$rho.full)
     }
     
     if (b$AR1.rho!=0) { ## correct RE/ML score for AR1 transform
       df <- if (getARs) sum(b$model$"(AR.start)") else 1
       object$gcv.ubre <- object$gcv.ubre - (n-df)*log(ld)
     }

     b$G$X <- b$qrx$R;b$G$dev.extra <- rss.extra
     b$G$pearson.extra <- rss.extra;b$G$n.true <- n
     b$y <- b$G$y;b$offset <- b$G$offset; b$G$w -> b$weights -> b$prior.weights;

  } else { ## method is "REML" or "ML"
    y <- b$G$y; w <- b$G$w;offset <- b$G$offset
    b$G$y <- b$qrx$f
    b$G$w <- b$G$y*0+1
    b$G$X <- b$qrx$R
    b$G$n <- length(b$G$y)
    b$G$offset <- b$G$y*0
    b$G$dev.extra <- rss.extra
    b$G$pearson.extra <- rss.extra
    b$G$n.true <- n
    if (b$scale.estimated) scale <- -1 else scale = b$sig2
    in.out <- list(sp=b$sp,scale=b$reml.scale)
    object <- gam(G=b$G,method=method,gamma=b$gamma,scale=scale,in.out=in.out) 
    if (b$AR1.rho!=0) { ## correct RE/ML score for AR1 transform
       df <- if (getARs) sum(b$model$"(AR.start)") else 1
       object$gcv.ubre <- object$gcv.ubre - (n-df)*log(ld)
    }
    offset -> b$G$offset -> b$offset
    w -> b$G$w -> b$weights -> b$prior.weights; n -> b$G$n
    y -> b$G$y -> b$y;
  }
 
  if (method=="GCV.Cp") { 

    b$coefficients <- fit$b
    b$edf <- post$edf
    b$edf1 <- post$edf1
    ##b$F <- post$F
    b$full.sp <- fit$sp.full
    b$gcv.ubre <- fit$score
    b$hat <- post$hat
    b$mgcv.conv <- fit$gcv.info 
    b$optimizer="magic"
    b$rank <- fit$gcv.info$rank
    b$Ve <- post$Ve
    b$Vp <- post$Vb
    b$sig2 <- b$scale <- fit$scale
    b$sp <- fit$sp

  } else { ## REML or ML
    b$coefficients <- object$coefficients
    b$edf <- object$edf
    b$edf1 <- object$edf1
    ##b$F <- object$F
    b$full.sp <- object$sp.full
    b$gcv.ubre <- object$gcv.ubre
    b$hat <- object$hat
    b$outer.info <- object$outer.info 
    b$rank <- object$rank
    b$Ve <- object$Ve
    b$Vp <- object$Vp
    b$sig2 <- b$scale <- object$sig2
    b$sp <- object$sp
    if (b$AR1.rho!=0) { ## correct RE/ML score for AR1 transform
      b$gcv.ubre <- b$gcv.ubre - (n-1)*log(ld)
    }
  }

  b$R <- b$qrx$R
  b$G$X <- NULL
  b$linear.predictors <- as.numeric(predict.gam(b,newdata=b$model,block.size=chunk.size))
  b$fitted.values <- b$linear.predictor ## strictly additive only!
  
  b$residuals <- sqrt(b$family$dev.resids(b$y,b$fitted.values,b$prior.weights)) * 
                      sign(b$y-b$fitted.values)
  b$deviance <- sum(b$residuals^2)
  b$aic <- b$family$aic(b$y,1,b$fitted.values,b$prior.weights,b$deviance) +
           2 * sum(b$edf) 
  if (b$AR1.rho!=0) { ## correct aic for AR1 transform
    df <- if (getARs) sum(b$model$"(AR.start)") else 1
    b$aic <- b$aic + 2*(n-df)*log(ld)
  }
  b$null.deviance <- sum(b$family$dev.resids(b$y,mean(b$y),b$prior.weights))
  names(b$coefficients) <- names(b$edf) <- cnames
  b
} ## end of bam.update


#### ISSUES:   
## ? negative binomial support --- docs say it's there...
## offset unused in bam/bgam.fit, also gp only needed for "response",
## so could efficiently be replaced

## (c) Simon Wood 2018. Released under GPL2.
## Implements the version of INLA (Rue et al. 2009, JRSSB) described in
## Wood "Simplified Integrated Nested Laplace Approximation" (submitted, 2018)


FFdes <- function (size=5,ccd=FALSE) {
## creates level 5 fractional factorial designs, up to size=120
## according to Sanchez and Sanchez (2005)
## ACM Transactions on Modeling and Computer Simulation 15(4), 362-377
## If ccd==TRUE, appends this to make the outer points of a
## Central Composite design (origin is added for full design).

  fwt <- function(x) {
  ## fast Walsh transform
    lag <- 1
    while (lag < length(x)) {
      offset <-  lag * 2
      ngroups <- length(x)/offset
      for (group in 0:(ngroups-1)) { ## vectorized
        j <- 1:lag + group*offset
        k <- j + lag
	xj <- x[j];xk <- x[k]
	x[j] <- xj + xk
	x[k] <- xj - xk
      }
      lag <- offset
    } ## while lag
    x
  } ## fwt

  index <- c(1, 2, 4, 8, 15, 16, 32, 51, 64, 85, 106, 128,
  150, 171, 219, 237, 247, 256, 279, 297, 455, 512, 537,
  557, 594, 643, 803, 863, 998, 1024, 1051, 1070, 1112,
  1169, 1333, 1345, 1620, 1866, 2048, 2076, 2085, 2185,
  2372, 2456, 2618, 2800, 2873, 3127, 3284, 3483, 3557,
  3763, 4096, 4125, 4135, 4174, 4435, 4459, 4469, 4497,
  4752, 5255, 5732, 5804, 5915, 6100, 6369, 6907, 7069,
  8192, 8263, 8351, 8422, 8458, 8571, 8750, 8858, 9124,
  9314, 9500, 10026, 10455, 10556, 11778, 11885, 11984,
  13548, 14007, 14514, 14965, 15125, 15554, 16384, 16457,
  16517, 16609, 16771, 16853, 17022, 17453, 17891, 18073,
  18562, 18980, 19030, 19932, 20075, 20745, 21544, 22633,
  23200, 24167, 25700, 26360, 26591, 26776, 28443, 28905,
  29577, 32705)

  power <- index;p <- 1
  for (i in 1:length(index)) {
    if (index[i]>=p) p <- p * 2
    power[i] <- p
  }

  if (size > 120||size<1) stop("size must be in [1,120]")
  design <- matrix(0,power[size],size)
  for (i in 1:size) {
    design[index[i]+1,i] <- 1
    design[,i] <- fwt(design[,i])
  }
  if (ccd&&size>1) {
    design <- rbind(design,diag(size)*sqrt(size),-diag(size)*sqrt(size))

  }
  design
} ## FFdes


logf <- function(beta,b,Bi=NULL,Xm=NULL,deriv=0) {
## get log joint density and first deriv w.r.t. coefs for a gam...
## Bi is a matrix mapping from interesting parameters to model parameters
  ## first deal with the log likelihood...
  if (is.null(Xm)) Xm <- if (is.null(b$X)) model.matrix(b) else b$X
  dd <- NULL
  if (!is.null(Bi)) beta <- drop(Bi %*% beta)
  if (inherits(b$family,"general.family")) {
    foo <- b$family$ll(b$y,Xm,beta,b$prior.weights,b$family,offset=b$offset,deriv=1)
    sll <- -2*foo$l
    dd <- -foo$lb
  } else if (inherits(b$family,"extended.family")) {
    theta <- b$family$getTheta()
    eta <- if (is.null(b$Xd)) as.numeric(Xm%*%beta + b$offset) else Xbd(b$Xd,beta,b$kd,b$ks,b$ts,b$dt,b$v,b$qc,b$drop) + b$offset
    mu <- b$family$linkinv(eta)
    sll <- sum(b$family$dev.resids(b$y,mu,b$prior.weights,theta)) ## deviance
    if (deriv) {
      #dd <- colSums((b$family$mu.eta(eta)*b$family$Dd(b$y,mu,theta,b$prior.weights)$Dmu)*Xm)/2
      dd <- if (is.null(b$Xd)) drop((b$family$mu.eta(eta)*b$family$Dd(b$y,mu,theta,b$prior.weights)$Dmu) %*% Xm)/2 else
      XWyd(b$Xd,b$family$mu.eta(eta),b$family$Dd(b$y,mu,theta,b$prior.weights)$Dmu,b$kd,b$ks,b$ts,b$dt,b$v,b$qc,b$drop)/2
    }
  } else { ## regular exponential family
    eta <- if (is.null(b$Xd)) as.numeric(Xm%*%beta + b$offset) else Xbd(b$Xd,beta,b$kd,b$ks,b$ts,b$dt,b$v,b$qc,b$drop) + b$offset
    mu <- b$family$linkinv(eta)
    sll <- sum(b$family$dev.resids(b$y,mu,b$prior.weights)) ## deviance
    if (deriv) {
      ##dd <- -colSums(b$prior.weights*(b$family$mu.eta(eta)*(b$y-mu)/b$family$variance(mu))*Xm)
      dd <- if (is.null(b$Xd)) -drop((b$prior.weights*(b$family$mu.eta(eta)*(b$y-mu)/b$family$variance(mu))) %*% Xm) else
      -XWyd(b$Xd,b$prior.weights,b$family$mu.eta(eta)*(b$y-mu)/b$family$variance(mu),b$kd,b$ks,b$ts,b$dt,b$v,b$qc,b$drop)
    }
  } ## deviance done
  ## now the smoothing prior/penalty
  ## NOTE: id's, fixed sp ???
  if (length(b$smooth)) {
    k <- 1;pen <- 0
    for (i in 1:length(b$smooth)) for (j in 1:length(b$smooth[[i]]$S)) {
      ind <- b$smooth[[i]]$first.para:b$smooth[[i]]$last.para
      b0 <- beta[ind]
      Sb <- b$smooth[[i]]$S[[j]] %*% b0 * b$sp[k]
      pen <- pen + sum(b0*Sb)
      if (deriv) dd[ind] <- dd[ind] + Sb
      k <- k + 1 
    }
  }
  if (!is.null(Bi)) dd <- drop(t(dd) %*% Bi)
  list(ll =(sll + pen)/(2*b$sig2),dd=dd) ## on neg log lik scale
} ## logf

glogf <- function(beta,b,Bi=NULL,Xm=NULL) {
  logf(beta,b,Bi=Bi,Xm=Xm,deriv=1)$dd
}

flogf <- function(beta,b,Bi=NULL,Xm=NULL) {
  logf(beta,b,Bi=Bi,Xm=Xm,deriv=0)$ll
}

Acomp <- function(A,ortho=TRUE) {
## simple matrix completion, if A is p by n, p <= n
## then returns full rank n by n matrix, B, whose last n-p
## rows are orthogonal to A, and its inverse Bi.
  p <- nrow(A);n<- ncol(A)
  if (ortho) { ## use orthogonal methods - very stable
    qra <- qr(t(A))  
    R <- qr.R(qra)
    if (Rrank(R)<p) stop("rank deficient re-parameterization")
    Q <- qr.Q(qra,complete=TRUE)
    B <- if (p<n) rbind(A,t(Q[,(p+1):n])) else A
    Bi <- t(backsolve(R,t(Q[,1:p])))
    if (p<n) Bi <- cbind(Bi,Q[,(p+1):n])
  } else { ## cross-products and one LU based solve - BLAS friendly
    C <- A[,-(1:p)]
    D <- A[,1:p] - C %*% t(C)
    Di <- try(solve(D),silent=TRUE)
    if (inherits(Di,"try-error")) stop("rank deficient re-parameterization")
    B <- rbind(A,cbind(t(C),diag(n-p)))
    tCDi <- t(C)%*%Di 
    Bi <- rbind(cbind(Di,-Di %*% C),cbind(-tCDi,diag(n-p)+tCDi%*%C))
  }
  list(B=B,Bi=Bi)
} ## Acomp


dg <- function(m,f0=1.5) {
## inla hyperparameter design point generator...
  D <- FFdes(m,ccd=TRUE)*f0
  if (f0<sqrt(2)) warning("modal weight <=0 in integration step!!")
  ## the Rue et al (2009) weights
  delta <- 1/((f0^2-1)*(1 + exp(-m*f0^2/2)))
  k0 <- 1 - delta
  k1 <- delta/nrow(D)
  k1 <- rep(k1,nrow(D))
  list(D=D,k0=k0,k1=k1)
} ## dg

cholinv <- function(A) {
## invert a +ve definite matrix via pivoted Cholesky with diagonal
## pre-conditioning. 
  d <- 1/sqrt(diag(A))
  R <- chol(d*t(A*d),pivot=TRUE)
  ipiv <- piv <- attr(R,"pivot")
  ipiv[piv] <- 1:ncol(A)
  d*t(chol2inv(R)[ipiv,ipiv]*d)
} # cholinv

Rsolve <- function(R,b) {
## solves R'Ra=b, where A = R'R, possibly with pivoting and
## possibly with diagonal pre-conditioning. If diagonal pre-conditioning
## has been used then attribute "dpc" of R should contain the 1/sqrt(diag(A))
## in unpivoted order, and R = chol(dpc*t(A*dpc))
  piv <- attr(R,"pivot")
  d <- attr(R,"dpc")
  if (!is.null(d)) b <- b * d
  if (is.null(piv)) {
    a <- backsolve(R,forwardsolve(R,b,upper.tri=TRUE,transpose=TRUE))
  } else {
    ipiv <- piv; ipiv[piv] <- 1:ncol(R)
    a <- if (is.matrix(b)) backsolve(R,forwardsolve(R,b[piv,],upper.tri=TRUE,transpose=TRUE))[ipiv,] else
    backsolve(R,forwardsolve(R,b[piv],upper.tri=TRUE,transpose=TRUE))[ipiv]
  }
  if (!is.null(d)) a <- a * d
  a
} ## Rsolve

ginla <- function(G,A=NULL,nk=16,nb=100,J=1,interactive=FALSE,int=0,approx=0) {
## apply inla to a gam post fit
## A is matrix or vector of linear transforms of interest, or an indices
## of the coefficients of interest (only if length!=p).
## TODO:
##      * bam and gam option to return Hessian (Sl.fitChol computes hessian,
##        but in gam its buried deeper)?
##      * handling of rank deficiency?
  prog <- interactive()&&interactive<2
  if (!inherits(G,"gam.prefit")&&!inherits(G,"bam.prefit")) stop("Requires a gam or bam prefit object")
  if (int !=0 ) G0 <- G ## need un-manipulated copy for calling gam
  if (inherits(G,"gam.prefit")) b <- gam(G=G,method="REML") else {
    if (!inherits(G,"bam.prefit")) stop("Requires a gam or bam prefit object")
    if (is.null(G$Xd)) stop("bam fits only supported with discrete==TRUE")
    b <- bam(G=G)
  }
  V <- sp.vcov(b,reg=.01)
  lsp <- log(b$sp) ## smoothing parameters
  scale.estimated <- b$scale.estimated
  if (scale.estimated) lsp <- c(lsp,log(b$sig2)) ## scale parameter
  ## other hyper-parameters...
  if (inherits(b$family,"extended.family")) {
    n.theta <- if (is.null(b$family$n.theta)) 0 else b$family$n.theta
    if (n.theta > 0) lsp <- c(b$family$getTheta(),lsp)
  }  else n.theta <- 0

  ## check that family supports enough derivatives to allow integration step,
  ## otherwise just use empirical Bayes.
  if (!is.null(G$family$available.derivs)&&G$family$available.derivs==0&&int>0) {
    int <- 0
    warning("integration not available with this family - insufficient derivatives")
  }
  
  ## Gaussian approximation is that log(sp) ~ N(lsp,V)
  if (int>0) { ## integration requested
    rV <- chol(V) ## Rv'z + lsp gives trial sp
    ip <- dg(ncol(rV)) ## integration points
    nip <- nrow(ip$D)
  } else nip <- 0

  if (!is.null(G$family$preinitialize)) {
    if (inherits(G$family,"general.family")) {
      Gmod <- G$family$preinitialize(G)
      for (gnam in names(Gmod)) G[[gnam]] <- Gmod[[gnam]] ## copy these into G 
    } else {
      ## extended family - just initializes theta and possibly y
      pini <- G$family$preinitialize(G$y,G$family)
      #if (!is.null(pini$Theta)) G$family$putTheta(pini$Theta) ## DON'T - will reset stored theta!
      if (!is.null(pini$y)) G$y <- pini$y
    }
  }

  X <- G$X
  G$prior.weights <- G$w
  G$sp <- b$sp
  reml <- rep(0,nip)
  dens <- list() ## list of lists of splines representing target density
  beta.lim <- matrix(0,ncol(X),2) ## evaluation limits for each A beta
  p <- ncol(b$Vp)
  if (!is.null(A)) { ## a transformation of original parameters is required
    if (is.matrix(A)||length(A)==p) { ## linear transforms needed
      B <- Acomp(A,is.null(G$Xd)) ## use orthogonal method only with gam fitting
      pa <- nrow(A)
      kind <- 1:pa
    } else { ## just a list of elements of beta
      A <- round(A)
      pa <- length(A)
      if (max(A)>p||min(A)<1||pa>p) stop("something wrong with A index vector")
      kind <- A
      B <- list()
    }
  } else { pa=p;kind <- 1:pa;B <- list()}
  iprog <- 0
  if (prog) prg <- txtProgressBar(min = 0, max = (nip+1)*pa, initial = 0,
              char = "=",width = NA, title="Progress", style = 3)
  for (qq in 0:nip) { ## integration loop
    dens[[qq+1]] <- list() ## elements will be log densities for each target parameter
    if (qq>0) { ## then a new fit is needed
       sp <- drop(t(rV) %*% ip$D[qq,]) + lsp
       sp <- pmin(pmax(sp,lsp-10),lsp+10) ## avoid +/- inf
       wprior <- -sum((sp-lsp)^2/(200)) ## include the prior assumed in sp.vcov
       if (n.theta>0) { ## set family hyper-parameters
         ii <- 1:n.theta
         G$family$putTheta(sp[ii])
	 G$n.theta <- 0 ## fixed not estimated now.
	 sp <- sp[-ii]
       }
       if (scale.estimated) {
         scale <- exp(sp[length(sp)])
	 sp <- sp[-length(sp)]
       } else scale <- 1	 
       sp <- exp(sp)
       if (inherits(G,"gam.prefit")) b <- gam(G=G0,method="REML",sp=sp,scale=scale) else
       b <- bam(G=G0,sp=sp,scale=scale)
       G$sp <- sp
       reml[qq] <- -b$gcv.ubre + wprior
    } else max.reml <- -b$gcv.ubre ## maximum LAML
    G$family <- b$family
    G$sig2 <- b$sig2
    
    beta <- coef(b)
    if (!is.null(B$B)) { ## a transformation of original parameters is required
      b$Vp <- B$B%*%b$Vp%*%t(B$B)
      beta <- drop(B$B%*%beta)
    }
    if (approx<2) {
      H <- cholinv(b$Vp) ## get Hessian - would be better returned directly by gam/bam
      dpc <- 1/sqrt(diag(H)) ## diagonal pre-conditioning
      R1 <- chol(dpc*t(H*dpc),pivot=TRUE)
      piv <- attr(R1,"pivot")
    }  
    sd <- diag(b$Vp)^.5 ## standard dev of Gaussian approximation
    BM <- matrix(0,p,nk) ## storage for beta[k] conditional posterior modes
    inla <- list(density=matrix(0,pa,nb),beta=matrix(0,pa,nb)) ## storage for marginals
    ldet <- dens0 <- rep(0,nk)
    eps <- .0001
    qn <- qnorm(seq(eps,1-eps,length=nk))
    kk <- 0
    for (k in kind) {
      kk <- kk + 1 ## counter for output arrays
      if (approx<2) { ## need R'R = H[-k,-k] (pivoted & pre-conditioned)
        kd <- which(piv==k) ## identify column of pivoted R1 corresponding to col k in H
        R <- choldrop(R1,kd) ## update R
        pivk <- piv[-kd]; pivk[pivk>k] <- pivk[pivk>k]-1 ## pivots updated
        attr(R,"pivot") <- pivk ## pivots of updated R
        attr(R,"dpc") <- dpc[-k] ## diagonal pre-conditioner
        ldetH <- 2*(sum(log(diag(R)))-sum(log(dpc[-k]))) ## log det of H[-k,-k]
      }
      bg <- qn*sd[k]+beta[k]
      BM[k,] <- bg
      BM[-k,] <- beta[-k] + b$Vp[-k,k]%*%((t(bg)-beta[k])/b$Vp[k,k]) ## Gaussian approx.
      if (approx==0) { ## get actual modes
        db <- db0 <- beta*0
        for (i in c((nk/2):1,(nk/2):nk)) {
          beta0 <- BM[,i] + db0
          nn <- logf(beta0,G,B$Bi,X,deriv=1)
          if (is.finite(nn$ll)) for (j in 1:20) { ## newton loop
	    if (max(abs(nn$dd[-k]))<1e-4*abs(nn$ll)) break
	    # db[-k] <- -backsolve(R,forwardsolve(Rt,nn$dd[-k]))
	    db[-k] <- -Rsolve(R,nn$dd[-k])
            beta1 <- beta0 + db
            nn1 <- logf(beta1,G,B$Bi,X,deriv=1)
	    get.deriv <- FALSE
	    hstep <- 0 
	    while (!is.finite(nn1$ll) || nn1$ll>nn$ll) {
	      db <- db/2;
	      hstep <- hstep+1
	      beta1 <- beta0 + db
	      nn1 <- logf(beta1,G,B$Bi,X,deriv=0)
	      get.deriv <- TRUE
	    }  
	    if (get.deriv) nn1 <- logf(beta1,G,B$Bi,X,deriv=1)
	    nn <- nn1
	    beta0 <- beta1
          } ## newton loop
          db0 <- if (i==1) 0 else beta0 - BM[,i]
          BM[,i] <- beta0
          dens0[i] <- nn$ll
        }
      } else for (i in 1:nk) dens0[i] <- logf(BM[,i],G,B$Bi,X,deriv=0)$ll	
      ## now get the log determinant correction...
      if (approx<2) {
        if (J>1) {
          vb <- apply(BM,1,var);vb[k] <- 0
          j <- length(vb)
          del <- which(rank(vb)%in%(j:(j-J+2)))
        }	
        step.length <- mean(colSums((BM - beta)^2)^.5)/20
        D <- rep(c(-1,1),J)
        ## create matrix of steps and matrix of evaluated gradient at steps
        for (i in 1:nk) if (is.finite(dens0[i])) {
          bm <- BM[,i]
          db <- beta - bm;db[k] <- 0
          db <- db/sqrt(sum(db^2))*step.length
          for (j in 1:J) {
            h = H[-k,-k] %*% db[-k] + if (j>1)  u%*%(D[1:(2*(j-1))]*(t(u)%*%db[-k])) else 0
            g1 <-  glogf(bm+db/2,G,B$Bi,X) - glogf(bm-db/2,G,B$Bi,X)
	    v <- cbind(h/sqrt(sum(db[-k]*h)),g1[-k]/sqrt(sum(db*g1)))
            u <- if (j>1) cbind(v,u) else v
	    db <- -db; if (j<J) db[del[j]] <- 0
          }
          #Hu <- backsolve(R,forwardsolve(Rt,u)) %*% diag(D)
	  Hu <- Rsolve(R,u) %*% diag(D)
          ldet[i] <- (ldetH + as.numeric(determinant(diag(2*J)+t(u)%*%Hu)$modulus))
        }
      } else ldet <-  0 ## constant	
      dens0 <- -dens0 - ldet/2
      dens0[!is.finite(dens0)] <- min(dens0,na.rm=TRUE) - 10
      dens0 <- dens0 - max(dens0) ## overflow proof
      din <- interpSpline(bg,dens0) ## interpolant of log density
      ok <- FALSE
      bg0 <- min(bg)
      bg1 <- max(bg)
      while (!ok) {
        ## normalize the pdf
        inla$beta[kk,] <- seq(bg0,bg1,length=nb) ## output beta sequence
        inla$density[kk,] <- exp(predict(din,inla$beta[kk,])$y) ## un-normalized density
        n.const <- sum(inla$density[kk,])*(inla$beta[kk,2]-inla$beta[kk,1])
        inla$density[kk,] <- inla$density[kk,]/n.const
        ## check that evaluation range is wide enough to cover most of density
        maxd <- max(inla$density[kk,])
        ok <- TRUE
        if (inla$density[kk,1]>maxd*5e-3) {
          ok <- FALSE
	  bg0 <- bg0 - sd[k]
        }
        if (inla$density[kk,nb]>maxd*5e-3) {
          ok <- FALSE
	  bg1 <- bg1 + sd[k]
        }
      }
      ## normalizing interpolant as well...
      din$coefficients[,1] <- din$coefficients[,1] - log(n.const)
      dens[[qq+1]][[kk]] <- din
      if (qq==0||bg0<beta.lim[k,1]) beta.lim[k,1] <- bg0
      if (qq==0||bg1>beta.lim[k,2]) beta.lim[k,2] <- bg1
      ## normalize
      iprog <- iprog + 1
      if (prog) setTxtProgressBar(prg, iprog)
      if (interactive)
      { plot(inla$beta[kk,],inla$density[kk,],ylim=range(inla$density[kk,])*1.2,type="l",col=2,
           xlab=bquote(beta[.(k)]),ylab="density")
        lines(inla$beta[kk,],dnorm(inla$beta[kk,],mean=beta[k],sd=sd[k]))
        if (interactive==2) readline()
      }  
    } ## beta loop
  } ## integration loop
  ## now the actual integration, if nip>0, otherwise we are done
  if (nip) {
    ## normalize the integration weights
    reml <- c(ip$k0,ip$k1*exp(reml-max.reml))
    reml <- reml/sum(reml)
    for (k in 1:pa) { ## parameter loop
      inla$beta[k,] <- seq(beta.lim[k,1],beta.lim[k,2],length=nb) ## output beta sequence
      inla$density[k,] <- exp(predict(dens[[1]][[k]],inla$beta[k,])$y)*reml[1] ## normalized density  
    }
    for (qq in 2:(nip+1)) {
      for (k in 1:pa) { ## parameter loop
        inla$density[k,] <- inla$density[k,] +
	    exp(predict(dens[[qq]][[k]],inla$beta[k,])$y)*reml[qq] ## normalized density
      }
    }
    inla$reml <- reml
  } ## if nip
  if (prog) cat("\n")
  inla
} ## ginla or gam inla newton enhanced (ginlane)

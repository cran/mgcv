## grouped family constructor (c) Simon N Wood (2023)

gfam <- function(fl) {
## gfam implements a grouped family, made up of the families in the list 'fl'.
## The response is supplied as a two column matrix, where the first column is
## the actual response, and the second is the index (starting at one) indicating
## the family in 'fl' to which each observation belongs. The elements of 'fl'
## can be ordinary exponential families or extended families. Extended families
## that expect a matrix response are not usable. 'gfam' has class 'extended.family'
## and a fixed scale parameter of 1. Scale parameters of component families
## are estimated using REML/NCV.
## fl has an attribut "fi" indicating which element of fl each response observation
## comes from. 
  nf <- length(fl)
  ## evaluate and check type of families, adding extra derivatives as
  ## required...
  fam_class <- "extended.family" ## overall class of gfam
  fam <- "gfam{"; link <- "{"
  n.theta <- 0 ## number of theta parameters including scale parameters
  Theta <- rep(0,0)
  need.rsd <- FALSE
  for (i in 1:nf) {
    if (is.character(fl[[i]])) fl[[i]] <- eval(parse(text=fl[[i]]))
    if (is.function(fl[[i]])) fl[[i]] <- fl[[i]]()
    if (is.null(fl[[i]]$family)) stop("family not recognized")
    if (fam_class != "general.family" && inherits(fl[[i]],"general.family")) 
      fam_class <- "general.family"
    if (!inherits(fl[[i]],"extended.family")) { ## regular exponential family
      fl[[i]] <- fix.family.ls(fix.family.var(fl[[i]]))
      fl[[i]]$scale <- if (fl[[i]]$family %in% c("poisson","binomial")) 1 else -1 ## fixed scale?
    } else if (is.null(fl[[i]]$scale)) fl[[i]]$scale <- 1
    fl[[i]] <- fix.family.link.extended.family(fl[[i]])
    sep <- if (i==1) "" else ","
    fam <- paste(fam,fl[[i]]$family,sep=sep)
    link <- paste(link,fl[[i]]$link,sep=sep)
    if (inherits(fl[[i]],"extended.family")) {
      if (fl[[i]]$scale < 0) { 
        n.theta <- n.theta + fl[[i]]$n.theta + 1 
        theta <- c(fl[[i]]$getTheta(),0)
      } else {
        n.theta <- n.theta + fl[[i]]$n.theta
        theta <- fl[[i]]$getTheta()
      }
    } else { ## exponential family
      if (fl[[i]]$scale<0) {
        n.theta <- n.theta + 1
	theta <- 0
      } else theta <- rep(0,0)
    }
    Theta <- c(Theta,theta)
    if (!is.null(fl[[i]]$residuals)) need.rsd <- TRUE
  } ## family list setup loop
  fam <- paste(fam,"}",sep=""); link <- paste(link,"}",sep="")

  if (fam_class!="extended.family") stop("general familes not implemented so far")
  ## stash the family list...
  attr(fl,"n.theta") <- n.theta
  env <- new.env(parent = environment(gfam)) 
  assign(".fl",fl, envir = env)
  getfl <- function( ) get(".fl")
  putfl <- function(fl) assign(".fl", fl, envir=environment(sys.function()))
  assign(".Theta", Theta, envir = env)
  getTheta <- function() get(".Theta") 
  putTheta <- function(theta) {
    assign(".Theta", theta,envir=environment(sys.function()))
    fl <- get(".fl");i0 <- 1 
    for (i in 1:length(fl)) if (inherits(fl[[i]],"extended.family")) { ## putTheta to components as well
      nth <- fl[[i]]$n.theta + if (fl[[i]]$scale<0) 1 else 0 # number of params including any scale 
      if (fl[[i]]$n.theta>0) { ## only store non-scale params
        th <- if (fl[[i]]$scale<0) theta[i0:(i0+nth-2)] else theta[i0:(i0+nth-1)]
	fl[[i]]$putTheta(th)
      }
      i0 <- i0 + nth
    }
  } ## putTheta gfam
  
  ## function needed by bam and predict.gam to subset data...
  setInd <- function(ind) {
  ## enable the fi index to be subsetted by ind, and restored is ind==NULL
    fl <- get(".fl")
    if (is.null(attr(fl,"fifull"))) {
      if (is.null(ind)) return() else attr(fl,"fifull") <- attr(fl,"fi") ## store full fi
    }  
    attr(fl,"fi") <- if (is.null(ind)) attr(fl,"fifull") else attr(fl,"fifull")[ind]
    assign(".fl", fl, envir=environment(sys.function()))
  }  

  dev.resids <- function(y,mu,wt,theta=NULL) {
    if (is.null(theta)) theta <- get(".Theta")
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    r <- mu;i0 <- 1 ## i0 is current location in theta vector
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      if (inherits(fl[[i]],"extended.family")) {
        nth <- fl[[i]]$n.theta
	th <- if (nth) theta[i0:(i0+nth-1)] else NULL
	r[ii] <- fl[[i]]$dev.resids(y[ii],mu[ii],wt[ii],th) ## dev resids
      } else { ## standard exponential family 
        nth <- 0
        r[ii] <- fl[[i]]$dev.resids(y[ii],mu[ii],wt[ii]) ## dev resids
      }
      if (fl[[i]]$scale<0) {
        r[ii] <- r[ii]/exp(theta[i0+nth]) ## scaled deviance
	i0 <- i0 + nth + 1
      }	else i0 <- i0 + nth
    }
    r
  } ## dev.resids gfam

  Dd <- function(y, mu, theta, wt, level=0) {
  ## the derivatives of the scaled deviance for grouped families, treating
  ## any scale parameters as part of the family parameter vector theta, so
  ## that the overall scale parameter is set to 1.
    kj <- function(j,n) (2*n-j+1)*j/2 
    filth <- function(n,j,nth) {
      ## create index in total 2nd deriv wrt theta vector
      ## for the 2nd derivs for nth theta values of a family,
      ## starting at jth element of total theta vector
      A <- matrix(1:nth,nth,nth);a <- A[A<=t(A)[,nth:1]]
      rep(kj(1:nth-1,n-j+1),nth:1) + a + kj(j-1,n)
    } ## filth
    filsc <- function(n,j,nth) {
      ## create index in total 2nd deriv wrt theta vector
      ## for 2nd derivs w.r.t. scale for a family which has 
      ## nth theta values. j is where all derivs for
      ## family start in 2nd deriv vector
      kj(1:(nth+1)-1,n-j+1) +(nth+1):1+ kj(j-1,n)
    } ## filth
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    r <- list(Dmu=y,Dmu2=y,EDmu2=y)
    n.theta <- length(theta)
    n <- length(mu)
    if (level>0) {
      r$EDmu2th <- r$Dmu2th <- r$Dth <- r$Dmuth <- matrix(0,n,n.theta)
      r$Dmu3 <- y
    }
    if (level>1) {
      r$Dth2 <- r$Dmuth2 <- r$Dmu2th2 <- matrix(0,n,n.theta*(n.theta+1)/2)
      r$Dmu3th <- matrix(0,n,n.theta)
      r$Dmu4 <- y
    }
    i0 <- 1 ##  position in theta
    kk <- 1 ## position in second derivative arrays 
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      if (inherits(fl[[i]],"extended.family")) {
        nth <- fl[[i]]$n.theta
	ith <- i0:(i0+nth-1) ## index of theta params for this family in theta
	th <- if (nth) theta[ith] else NULL
      } else nth <- 0
      if (fl[[i]]$scale<0) {
        rho <- theta[i0+nth] ## log scale parameter
	isc <- i0 + nth ## index of log scale param
	nth1 <- nth + 1
      }	else { nth1 <- nth; rho=0}
      scale <- exp(rho)
      if (inherits(fl[[i]],"extended.family")) {
        ri <- fl[[i]]$Dd(y[ii],mu[ii],th,wt[ii],level=level)
	r$Dmu[ii] <- ri$Dmu/scale
	r$Dmu2[ii] <- ri$Dmu2/scale
	r$EDmu2[ii] <- ri$EDmu2/scale
	if (level>0) { ## quantities needed for first derivatives
          if (nth) {
	    r$Dth[ii,ith] <- ri$Dth/scale
	    r$Dmuth[ii,ith] <- ri$Dmuth/scale
	    r$Dmu2th[ii,ith] <- ri$Dmu2th/scale
	    r$EDmu2th[ii,ith] <- ri$EDmu2th/scale
	  }
	  r$Dmu3[ii] <- ri$Dmu3/scale
	  if (fl[[i]]$scale<0) { ## need to append derivs w.r.t. log scale as well
            D <- fl[[i]]$dev.resids(y[ii],mu[ii],wt[ii],th)
            r$Dth[ii,isc] <- -D/scale 
	    r$Dmuth[ii,isc] <- -ri$Dmu/scale
	    r$Dmu2th[ii,isc] <- -ri$Dmu2/scale
	    r$EDmu2th[ii,isc] <- -ri$EDmu2/scale
          }
        }
	if (level>1) { ## whole damn lot - packing convention here is a pain
           r$Dmu4[ii] <- ri$Dmu4/scale
	   if (nth>0) {
	     ijth <- filth(n.theta,i0,nth)
             r$Dmu3th[ii,ith] <- ri$Dmu3th/scale
	     r$Dth2[ii,ijth] <- ri$Dth2/scale
	     r$Dmuth2[ii,ijth] <- ri$Dmuth2/scale
	     r$Dmu2th2[ii,ijth] <- ri$Dmu2th2/scale
           }
           if (fl[[i]]$scale<0) { ## 2nd derivs w.r.t log scale
             ijsc <- filsc(n.theta,i0,nth)
	     its <- if (nth) c(ith,isc) else isc 
	     r$Dmu3th[ii,isc] <- -ri$Dmu3/scale
	     r$Dth2[ii,ijsc] <- -r$Dth[ii,its]
	     r$Dmuth2[ii,ijsc] <- -r$Dmuth[ii,its]
	     r$Dmu2th2[ii,ijsc] <- -r$Dmu2th[ii,its]
           }
        }
      } else { ## exponential families
        vi <- fl[[i]]$variance(mu[ii])
	dv <- fl[[i]]$dvar(mu[ii])
	ri <- y[ii] - mu[ii]
        r$Dmu[ii] <- -2*ri/(vi*scale)
	r$Dmu2[ii] <- 2*(1 + ri*dv/vi)/(vi*scale)
	r$EDmu2[ii] <- 2/(vi*scale)
	if (level>0) {
	  d2v <- fl[[i]]$d2var(mu[ii]) 
          r$Dmu3[ii] <- -r$Dmu2[ii]*dv/vi + 2*(ri*(d2v/vi-(dv/vi)^2) -dv/vi)/(vi*scale)
	  if (fl[[i]]$scale<0) { ## need to append derivs w.r.t. log scale as well
            D <- fl[[i]]$dev.resids(y[ii],mu[ii],wt[ii])
	    r$Dth[ii,isc] <- -D/scale
	    r$Dmuth[ii,isc] <- -r$Dmu[ii]
	    r$Dmu2th[ii,isc] <- -r$Dmu2[ii]
	    r$EDmu2th[ii,isc] <- -r$EDmu2[ii]
	  }  
        }
	if (level>1) { ## whole damn lot - packing convention here is a pain
	  d3v <- fl[[i]]$d3var(mu[ii])
          r$Dmu4[ii] <- -r$Dmu2[ii]*d2v/vi -2*r$Dmu3[ii]*dv/vi +
	  2*(2*((dv/vi)^2-d2v/vi)+ri*(d3v/vi-3*dv*d2v/vi^2+2*(dv/vi)^3))/(vi*scale)
          if (fl[[i]]$scale<0) { ## 2nd derivs w.r.t log scale
            ijsc <- filsc(n.theta,i0,nth)
	    r$Dmu3th[ii,isc] <- -r$Dmu3[ii]
	    r$Dth2[ii,ijsc] <- -r$Dth[ii,isc]
	    r$Dmuth2[ii,ijsc] <- -r$Dmuth[ii,isc]
	    r$Dmu2th2[ii,ijsc] <- -r$Dmu2th[ii,isc]
          }
        }
      }
      i0 <- i0 + nth + if (fl[[i]]$scale<0) 1 else 0
    } ## family loop
    r
  } ## Dd gfam


  linkfun <- function(mu) {
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    eta <- mu
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      eta[ii] <- fl[[i]]$linkfun(mu[ii])  
    }
    eta
  } ## linkfun gfam

  linkinv <- function(eta) {
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    eta -> mu
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      mu[ii] <- fl[[i]]$linkinv(eta[ii])  
    }
    mu
  } ## linkinv gfam

  mu.eta <- function(eta) {
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    eta -> z
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      z[ii] <- fl[[i]]$mu.eta(eta[ii])  
    }
    z
  } ## mu.eta gfam

  validmu <- function(mu) {
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      if (!fl[[i]]$validmu(mu[ii])) return(FALSE)  
    }
    TRUE
  } ## validmu gfam

  valideta <- function(eta) {
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      if (!fl[[i]]$valideta(eta[ii])) return(FALSE)  
    }
    TRUE
  } ## validmu gfam


  g2g <- function(mu) {
    fl <- get(".fl");fi <- attr(fl,"fi")
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      mu[ii] <- fl[[i]]$g2g(mu[ii]) ## dev resids 
    }
    mu
  } ## d2link gfam

  g3g <- function(mu) {
    fl <- get(".fl");fi <- attr(fl,"fi")
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      mu[ii] <- fl[[i]]$g3g(mu[ii]) ## dev resids 
    }
    mu
  } ## d3link gfam

  g4g <- function(mu) {
    fl <- get(".fl");fi <- attr(fl,"fi")
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      mu[ii] <- fl[[i]]$g4g(mu[ii]) ## dev resids 
    }
    mu
  } ## d4link gfam


  ls <- function(y,w,theta,scale) {
  ## log saturated likelihood
    fl <- get(".fl");fi <- attr(fl,"fi")
    n.theta <- attr(fl,"n.theta")
    n <- rep(1,length(y)) ## basically dummy here - will always be 1
    ls <- 0; lsth1 <- rep(0,n.theta); LSTH1 <- matrix(0,length(y),n.theta)
    lsth2 <- matrix(0,n.theta,n.theta)
    i0 <- 1
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      if (inherits(fl[[i]],"extended.family")) {
        nth <- fl[[i]]$n.theta 
	ith <- i0 + 1:nth - 1
	th <- if (nth>0) theta[ith] else 0
	sca <- if (fl[[i]]$scale<0) exp(theta[i0+nth]) else 1
        li <- fl[[i]]$ls(y[ii],w[ii],th,sca)
	ls <- ls + li$ls
	if (fl[[i]]$scale<0) {
	  nth <- nth + 1; ith <- c(ith,i0+nth-1)
	}  
        if (nth>0) { 
	  lsth1[ith] <- li$lsth1
	  LSTH1[ii,ith] <- li$LSTH1;lsth2[ith,ith] <- li$lsth2
	}  
      } else { ## exponential family
        if (fl[[i]]$scale<0) {
	  nth <- 1; sca <- exp(theta[i0]) 
	} else { nth <- 0; sca <- 1 }
        li <- fl[[i]]$ls(y[ii],w[ii],n[ii],scale=sca)
	ls <- ls + li[1]
	if (nth) { ## derivs are w.r.t. log scale - note ls returns dervis w.r.t. scale.
          lsth1[i0] <- li[2]*sca
	  lsth2[i0,i0] <- li[3]*sca^2 + li[2]*sca
	  LSTH1[ii,i0] <- sca*as.numeric(w[ii]>0)*li[2]/sum(w[ii]>0)
        }
      }
      i0 <- i0 + nth 
    }
    list(ls=ls,lsth1=lsth1,LSTH1=LSTH1,lsth2=lsth2)
  } ## ls gfam

  aic <- function(y, mu, theta=NULL, wt, dev) { ## (y, n, mu, wt, dev) {
  ## note dev has to be ignored and re-computed component wise here
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    n <- rep(1,length(y)) ## dummy here 
    aic <- 0;i0 <- 1
    for (i in 1:length(fl)) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      if (inherits(fl[[i]],"extended.family")) {
        nth <- fl[[i]]$n.theta; ith <- 1:nth-1 + i0
        dev <- sum(fl[[i]]$dev.resids(y[ii],mu[ii],wt[ii],theta[ith]))
        aic <- aic + fl[[i]]$aic(y[ii],mu[ii],theta[ith],wt[ii],dev) ## dev resids
	i0 <- i0 + nth + if (fl[[i]]$scale<0) 1 else 0
      } else {
        dev <- sum(fl[[i]]$dev.resids(y[ii],mu[ii],wt[ii])) 
        aic <- aic + fl[[i]]$aic(y[ii],n[ii],mu[ii],wt[ii],dev) ## dev resids
	if (fl[[i]]$scale<0) i0 <- i0 + 1
      }
    }
    aic
  } ## aic gfam

  preinitialize <- function(y,family) {
    ## may return Theta and modified y.
    if (is.matrix(y)) {
      if (ncol(y)!=2) stop("gfam response must have 2 columns")
      fi <- y[,2]
      y <- y[,1]
    }
    fl <- get(".fl")
    nf <- length(fl)
    ui <- unique(fi) ## check index vector 
    if (any(!(ui%in%1:nf))||any(!(1:nf%in%ui))) stop("family index does not match family list")
    fi -> attr(fl,"fi") ## add fi as attributed to fl
    family$putfl(fl)    ## store modified fl 
    n.theta <- attr(fl,"n.theta")
    Theta <- rep(0,n.theta)
    theta.mod <- FALSE
    i0 <- 1 ## current start in theta vector
    for (i in 1:length(fl)) if (inherits(fl[[i]],"extended.family")) { ## family loop
      ii <- which(fi==i) ## index of data to which this family applies
      nth <- fl[[i]]$n.theta
      if (!is.null(fl[[i]]$preinitialize)) {
        ith <- i0 + 1:nth - 1
        pri <- fl[[i]]$preinitialize(y[ii],fl[[i]])
        if (!is.null(pri$y)) {
	  y[ii] <- pri$y
        }
        if (is.null(pri$Theta)) {
          theta[ith] <- fl[[i]]$getTheta()
        } else {
          theta.mod <- TRUE
	  Theta[ith] <- pri$Theta
        }
      }	
      i0 <- i0 + nth + if (fl[[i]]$scale<0) 1 else 0
    } else if (fl[[i]]$scale<0) i0 <- i0 + 1
    ret <- list(y=y)
    if (theta.mod) ret$Theta <- Theta
    ret
  }  ## preinitialize gfam

  postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept){
    ## deviance (sometimes), null.deviance (always) and family$family changed (always).
    ## NOTE: null.deviance is assuming a different intercept for each family -
    ##       see get.null.coef, perhaps
    fl <- get(".fl")
    fi <- attr(fl,"fi")
    #ext <- FALSE
    #for (f in fl) if (inherits(f,"extended.family")) ext <- TRUE
    #if (!ext) return(list())
    dev.mod <- FALSE; dev <- 0
    nulldev <- 0
    fam <- "gfam{"
    for (i in 1:length(fl)) {
      ii <- which(fi==i) ## index of data to which this family applies
      if (inherits(fl[[i]],"extended.family")) {
        pp <- fl[[i]]$postproc(fl[[i]],y[ii],prior.weights[ii],fitted[ii],
	                       linear.predictors[ii],offset[ii],intercept)
	nulldev <- nulldev + pp$null.deviance
        sep <- if (i==1) "" else ","
        fam <- paste(fam,pp$family,sep=sep)
	if (is.null(pp$deviance)) { ## still need to compute deviance for current subset in case modified later
	  dev <- dev + sum(fl[[i]]$dev.resids(y[ii],fitted[ii],prior.weights[ii]))
        } else {
          dev.mod <- TRUE
	  dev <- dev + pp$deviance
        }
      } else { ## regular exponential family
        sep <- if (i==1) "" else ","
        fam <- paste(fam,fl[[i]]$family,sep=sep)
	dev <- dev + sum(fl[[i]]$dev.resids(y[ii],fitted[ii],prior.weights[ii]))
	wtdmu <- if (intercept) sum(prior.weights[ii] * y[ii])/sum(prior.weights[ii]) else linkinv(offset[ii])
        nulldev <- nulldev + sum(fl[[i]]$dev.resids(y[ii], rep(wtdmu,length(ii)), prior.weights[ii]))
      }
    }
    fam <- paste(fam,"}",sep="")
    posr <- list(family=fam,null.deviance = nulldev)
    if (dev.mod) posr$deviance = dev
    posr
  }  ## postproc gfam
  
  residuals <- if (need.rsd) function(object,type=c("deviance","working","response","pearson")) {
  ## NOTE: NA action may require work
    if (type == "working") { 
      rsd <- object$residuals 
    } else {
      rsd <- object$residuals
      fl <- get(".fl")
      fi <- attr(fl,"fi")
      for (i in 1:length(fl)) {
        ii <- which(fi==i) ## index of data to which this family applies
        ## A bit of a hack...
        ob <- object
	ob$family <- fl[[i]]
	ob$y <- object$y[ii]
	ob$linear.predictors <- object$linear.predictors[ii]
	ob$fitted.values <- object$fitted.values[ii]
	ob$prior.weights <- object$prior.weights[ii]
	ob$model <- object$model[ii,]
	rsd[ii] <- residuals.gam(ob,type)
      }
    }
    rsd
  } else NULL ## residuals gfam
  
  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## returns a vector of predictions (response scale), or a list with elements fit and se.fit
  ## bam prediction not set up to use this!
    fl <- get(".fl")
    n <- if (is.null(eta)) { if (is.list(X)) NROW(X$kd) else nrow(X) } else length(eta)
    if (is.null(y)) {
      fi <- attr(fl,"fi")
      if (length(fi)!=n) stop("no family index")
    } else {
      if (is.matrix(y)) {
        if (ncol(y)!=2) stop("if response is a matrix it must have 2 columns")
	fi <- y[,2]
      } else fi <- y
    }
    nf <- length(fl)
    if (!all(unique(fi)%in%1:nf)) stop("family index does not match list of families")

    fit <- rep(0,n)
    se.fit <- if (se) fit else NULL
    if (is.null(eta)) {
      for (i in 1:nf) {
        ii <- which(fi==i) ## index of data to which this family applies
        if (is.null(fl[[i]]$predict)) {
	  fit[ii] <- off[ii] + if (is.list(X)) Xbd(X$Xd,beta,k=X$kd[ii,],ks=X$ks,ts=X$ts,
	                      dt=X$dt,v=X$v,qc=X$qc,drop=X$drop) else drop(X[ii,]%*%beta)
	  if (se) {
            se.fit[ii] <- if (is.list(X)) sqrt(pmax(0,diagXVXd(X$Xd,Vb,k=X$kd[ii,],ks=X$ks,ts=X$ts,
		                               dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1))) else
	                  sqrt(pmax(0,rowSums((X[ii,]%*%Vb)*X[ii,])))
	    se.fit[ii] <- se.fit[ii]*abs(fl[[i]]$mu.eta(fit[ii]))
          }
          fit[ii] <- fl[[i]]$linkinv(fit[ii])
	} else {
	  if (se) {
	    if (is.list(X)) {
              Xi <- X;Xi$kd <- X$kd[ii,]
            } else Xi <- X[ii,]
            pr <- fl[[i]]$predict(fl[[i]],se=se,eta=eta,y=y[ii],X=Xi,
                  beta=beta,off=off[ii],Vb=Vb)
	    if (is.matrix(pr[[1]])) stop("gfam mixed response scale prediction not possible here")	
	    se.fit[ii] <- pr$se.fit
	    fit[ii] <- pr$fit
	  } else {
            fit[ii] <- off[ii] + if (is.list(X)) Xbd(X$Xd,beta,k=X$kd[ii,],ks=X$ks,ts=X$ts,
	                         dt=X$dt,v=X$v,qc=X$qc,drop=X$drop) else drop(X[ii,]%*%beta)
            pr <- fl[[i]]$predict(fl[[i]],se=se,eta=fit[ii])[[1]]
            if (is.matrix(pr)) stop("gfam mixed response scale prediction not possible here")
            fit[ii] <- pr
	  }  
        } 
      }
      if (se) return(list(fit=fit,se.fit=se.fit)) else return(list(fit))
    } else {
      for (i in 1:length(fl)) {
        ii <- which(fi==i)
        fit[ii] <- fl[[i]]$linkinv(eta[ii])
      }
      return(list(fit))
    }
  } ## predict gfam


  ## initialize sets up mustart and n
  initialize <- expression({
    #if (is.null(mustart)) mustart <- y
    .family <- family
    .mustart <- y
    .weights <- weights
    .nobs <- nobs ## store nobs
    .yini <- y ## store y for later restoration
    fl <- family$getfl()
    fi <- attr(fl,"fi")
    for (kk in 1:length(fl)) {
      family <- fl[[kk]]
      ii <- which(kk == fi) 
      y <- .yini[ii]; nobs <- length(ii)
      weights <- .weights[ii]
      eval(family$initialize)
      .mustart[ii] <- mustart 
    }
    y <- .yini; nobs <- .nobs;weights <- .weights;
    mustart <- .mustart;family <- .family #n <- .n
  }) ## initialize gfam

  get.null.coef <- function(G,...) { 
     y <- G$y
     mustart <- etastart <- start <- NULL
     weights <- G$w
     nobs <- G$n
     family <- G$family
     eval(family$initialize)
     fl <- family$getfl()
     mum <- etam <- fi <- attr(fl,"fi")
     for (i in 1:length(fl)) {
       ii <- which(fi==i)
       mum[ii] <- mean(y[ii]) + 0 * ii
       etam[ii] <- fl[[i]]$linkfun(mum[ii])
     }	
     null.coef <- qr.coef(qr(G$X), etam)
     null.coef[is.na(null.coef)] <- 0
     null.scale <- sum(family$dev.resids(y, mum, weights))/nrow(G$X)
     list(null.coef = null.coef, null.scale = null.scale)
  } ## get.null.coef

  environment(setInd) <-
  environment(get.null.coef) <- environment(aic) <- environment(getfl) <- environment(ls) <-
  environment(getTheta) <- environment(putTheta) <- environment(putfl) <- environment(preinitialize) <-
  environment(Dd) <- environment(linkfun) <- environment(linkinv) <- environment(valideta) <- 
  environment(mu.eta) <- environment(g2g) <- environment(g3g) <- environment(postproc) <-
  environment(g4g) <- environment(dev.resids) <- environment(validmu) <- environment(predict) <- env
  if (need.rsd) environment(residuals) <- env

  structure(list(family = fam, dev.resids = dev.resids,aic = aic,
            link = link, linkfun = linkfun, linkinv = linkinv,mu.eta = mu.eta,
	    initialize = initialize, validmu = validmu,putTheta=putTheta,getTheta=getTheta,
	    g2g=g2g,g3g=g3g,g4g=g4g,Dd=Dd,preinitialize=preinitialize,valideta=valideta,
	    postproc=postproc,predict=predict,residuals=residuals,n.theta=n.theta,
	    ls=ls,getfl=getfl,putfl=putfl,setInd=setInd,
	    get.null.coef=get.null.coef,canonical="none"), class = c("extended.family","family"))
} ## gfam


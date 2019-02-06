## (c) Simon N. Wood (2013, 2014) coxph model general family. 
## Released under GPL2 ...

cox.ph <- function (link = "identity") { 
## Extended family object for Cox PH.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("identity")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for coxph family; available link is \"identity\" ")
  }
  env <- new.env(parent = .GlobalEnv)
  validmu <- function(mu) all(is.finite(mu))

    

    ## initialization is tough here... need data frame in reverse time order,
    ## and intercept removed from X...
  
#    preinitialize <- expression({
    ## code to evaluate in estimate.gam...
      ## sort y (time) into decending order, and
      ## re-order weights and rows of X accordingly
      ## matrix y has strat as second column
#      G$family$data <- list()
#      y.order <- if (is.matrix(G$y)) order(G$y[,2],G$y[,1],decreasing=TRUE) else
#                 order(G$y,decreasing=TRUE)
#      G$family$data$y.order <- y.order
#      G$y <- if (is.matrix(G$y)) G$y[y.order,] else G$y[y.order]
#      attrX <- attributes(G$X)
#      G$X <- G$X[y.order,,drop=FALSE]
#      attributes(G$X) <- attrX
#      G$w <- G$w[y.order]
#      G$offset <- G$offset[y.order]
#    })

    preinitialize <- function(G) {
      ## G is a gam pre-fit object. Pre-initialize can manipulate some of its
      ## elements, returning a named list of the modified ones.
      ## sort y (time) into decending order, and
      ## re-order weights and rows of X accordingly
      ## matrix y has strat as second column
      G$family$data <- list()
      y.order <- if (is.matrix(G$y)) order(G$y[,2],G$y[,1],decreasing=TRUE) else
                 order(G$y,decreasing=TRUE)
      G$family$data$y.order <- y.order
      G$y <- if (is.matrix(G$y)) G$y[y.order,] else G$y[y.order]
      attrX <- attributes(G$X)
      G$X <- G$X[y.order,,drop=FALSE]
      attributes(G$X) <- attrX
      G$w <- G$w[y.order]
      G$offset <- G$offset[y.order]
      list(family=G$family,y=G$y,X=G$X,w=G$w,offset=G$offset)
    } ##  preinitialize
    
    postproc <- expression({
    ## code to evaluate in estimate.gam, to do with data ordering and 
    ## baseline hazard estimation...
      ## first get the estimated hazard and prediction information...
      G$X <- Sl.initial.repara(G$Sl,G$X,inverse=TRUE,cov=FALSE,both.sides=FALSE)
      object$family$data <- G$family$hazard(G$y,G$X,object$coefficients,G$w,G$offset)
      rumblefish <- G$family$hazard(G$y,matrix(0,nrow(G$X),0),object$coefficients,G$w)
      s0.base <- exp(-rumblefish$h[rumblefish$r]) ## no model baseline survival 
      s0.base[s0.base >= 1] <- 1 - 2*.Machine$double.eps ## avoid NA later
      ## now put the survivor function in object$fitted
      object$fitted.values <- exp(-object$family$data$h[object$family$data$r]*exp(object$linear.predictors))
      ## compute the null deviance...
      s.base <- exp(-object$family$data$h[object$family$data$r]) ## baseline survival
      s.base[s.base >= 1] <- 1 - 2*.Machine$double.eps ## avoid NA later
      object$null.deviance <- ## sum of squares of null deviance residuals
      2*sum(abs((object$prior.weights + log(s0.base) + object$prior.weights*(log(-log(s0.base)))))) 
      ## and undo the re-ordering...
      y.order <- G$family$data$y.order
      object$linear.predictors[y.order] <- object$linear.predictors
      object$fitted.values[y.order] <- object$fitted.values
      if (is.matrix(object$y)) object$y[y.order,] <- object$y else object$y[y.order] <- object$y  
      object$prior.weights[y.order] <- object$prior.weights
    })
    
    initialize <- expression({
        n <- rep(1, nobs)
        if (is.null(start)) start <- rep(0,ncol(x))
    })

    hazard <- function(y, X,beta,wt,offset=0) {
    ## get the baseline hazard function information, given times in descending order in y
    ## model matrix (same ordering) in X, coefs in beta and censoring in wt (1 = death, 0
    ## = censoring)
    
      if (is.matrix(y)) { 
        ## first column is time, second is *numeric* code indicating strata
        strat <- y[,2] ## stratification variable
        y <- y[,1] ## event times 
	strat.lev <- unique(strat) 
	ns <- length(strat.lev) ## number of strata
	nt <- 0;for (i in 1:ns) nt <- nt + length(unique(y[strat==strat.lev[i]]))
        tr.strat <- tr <- rep(0,nt)
        k <- 1
        for (i in 1:ns) {
	  tr0 <- unique(y[strat==strat.lev[i]])
          ind <- k:(k+length(tr0)-1);k <- k + length(tr0)
          tr[ind] <- tr0 ## unique times at this stratification level
	  tr.strat[ind] <- strat.lev[i] ## stratication index for tr,h,q,km and a
        }
      }	else {
        ns <- 1
	tr <- unique(y)
	nt <- length(tr)
      }
      #tr <- unique(y);nt <- length(tr)
      p <- ncol(X)
      eta <- as.double(X%*%beta) + offset
      if (ns==1) {
        r <- match(y,tr)
        oo <- .C("coxpp",eta,A=as.double(X),as.integer(r),d=as.integer(wt),
               h=as.double(rep(0,nt)),q=as.double(rep(0,nt)),km=as.double(rep(0,nt)),
               n=as.integer(nrow(X)),p=as.integer(ncol(X)),
               nt=as.integer(nt),PACKAGE="mgcv")
	return(list(tr=tr,h=oo$h,q=oo$q,a=matrix(oo$A[1:(p*nt)],p,nt),nt=nt,r=r,km=oo$km))
      } else {
        r <- y*0;a <- matrix(0,p,nt)
	h <- q <- km <- rep(0,nt)
        for (i in 1:ns) { ## loop over strata
          ind <- which(strat==strat.lev[i])
	  trind <- which(tr.strat==strat.lev[i])
	  r0 <- match(y[ind],tr[trind])
	  nti <- length(trind)
	  etai <- if (p>0) eta[ind] else eta
          oo <- .C("coxpp",etai,A=as.double(X[ind,]),as.integer(r0),d=as.integer(wt[ind]),
               h=as.double(rep(0,nti)),q=as.double(rep(0,nti)),km=as.double(rep(0,nti)),
               n=as.integer(length(ind)),p=as.integer(p),
               nt=as.integer(nti),PACKAGE="mgcv")
	  ## now paste all this into return fields
	  h[trind] <- oo$h
	  q[trind] <- oo$q
	  km[trind] <- oo$km
	  r[ind] <- r0 ## note that indexing is to subsetted tr
	  a[,trind] <- matrix(oo$A[1:(p*nti)],p,nti)
        }
	return(list(tr=tr,h=h,q=q,a=a,nt=nt,r=r,km=km,strat=strat,tr.strat=tr.strat))
      }
      #p <- ncol(X)
      #list(tr=tr,h=oo$h,q=oo$q,a=matrix(oo$A[1:(p*nt)],p,nt),nt=nt,r=r,km=oo$km,strat=strat,tr.strat=tr.strat)
    }

    residuals <- function(object,type=c("deviance","martingale")) {
      type <- match.arg(type)
      w <- object$prior.weights;log.s <- log(object$fitted.values)
      res <- w + log.s ## martingale residuals
      if (type=="deviance") { 
        log.s[log.s>-1e-50] <- -1e-50
        res <- sign(res)*sqrt(-2*(res + w * log(-log.s)))
      }
      res 
    }


    predict <- function(family,se=FALSE,eta=NULL,y=NULL,
               X=NULL,beta=NULL,off=NULL,Vb=NULL) {
      ## prediction function.
      if (is.matrix(y)) { 
	strat <- y[,2];y <- y[,1]
	if (is.null(family$data$strat)) stop("something wrong with stratified prediction")
	ii <- order(strat,y,decreasing=TRUE) ## C code expects non-increasing
	strat <- strat[ii]
      } else {
        ii <- order(y,decreasing=TRUE) ## C code expects non-increasing
        strat <- NULL
      }
      if (sum(is.na(y))>0) stop("NA times supplied for cox.ph prediction")
      X <- X[ii,,drop=FALSE];y <- y[ii];
      n <- nrow(X)
      if (is.null(off)) off <- rep(0,n)
      if (is.null(strat)) {
        oo <- .C("coxpred",as.double(X),t=as.double(y),as.double(beta),as.double(off),as.double(Vb),
                a=as.double(family$data$a),h=as.double(family$data$h),q=as.double(family$data$q),
                tr = as.double(family$data$tr),
                n=as.integer(n),p=as.integer(ncol(X)),nt = as.integer(family$data$nt),
                s=as.double(rep(0,n)),se=as.double(rep(0,n)),PACKAGE="mgcv")
        s <- sef <- oo$s
        s[ii] <- oo$s
        sef[ii] <- oo$se 
      } else { ## stratified fit, need to unravel everything by strata
        pstrata <- unique(strat)
	ns <- length(pstrata)
	p <- ncol(X)
	s <- sef <- rep(0,length(y)) 
	for (i in 1:ns) {
          ind <- which(strat==pstrata[i]) ## prediction data index
	  trind <- which(family$data$tr.strat == pstrata[i])
	  n <- length(ind)
	  oo <- .C("coxpred",as.double(X[ind,]),t=as.double(y[ind]),as.double(beta),as.double(off),as.double(Vb),
                a=as.double(family$data$a[,trind]),h=as.double(family$data$h[trind]),q=as.double(family$data$q[trind]),
                tr = as.double(family$data$tr[trind]),
                n=as.integer(n),p=as.integer(p),nt = as.integer(length(trind)),
                s=as.double(rep(0,n)),se=as.double(rep(0,n)),PACKAGE="mgcv")
	  s[ind] <- oo$s
	  sef[ind] <- oo$se
        } ## strata loop
	s[ii] <- s
        sef[ii] <- sef   
      }    
      if (se) return(list(fit=s,se.fit=sef)) else return(list(fit=s))
    } ## predict

    rd <- qf <- NULL ## these functions currently undefined for Cox PH

    ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## function defining the cox model log lik.
    ## Calls C code "coxlpl"
    ## deriv codes: 0   - evaluate the log likelihood
    ##              1   - evaluate the grad and Hessian, H, of log lik w.r.t. coefs (beta)
    ##              2/3 - evaluate d1H =dH/drho given db/drho in d1b 
    ##                    (2 is for evaluation of diagonal only)
    ##              4 -  given d1b and d2b evaluate trHid2H= tr(Hp^{-1}d2H/drhodrho')
    ## Hp is the preconditioned penalized hessian of the log lik
    ##    which is of rank 'rank'.
    ## fh is a factorization of Hp - either its eigen decomposition 
    ##    or its Choleski factor
    ## D is the diagonal pre-conditioning matrix used to obtain Hp
    ##   if Hr is the raw Hp then Hp = D*t(D*Hr)

#!! if (!is.null(offset)&&sum(offset!=0)) stop("cox.ph does not yet handle offsets")
      
      if (is.matrix(y)) {
        ## first column is time, second is *numeric* code indicating strata
        strat <- y[,2] ## stratification variable
        y <- y[,1] ## event times 
	strat.lev <- unique(strat) 
	ns <- length(strat.lev) ## number of strata
      }	else ns <- 1
     
      p <- ncol(X)
      deriv <- deriv - 1
      mu <- X%*%coef + offset
      g <- rep(0,p);H <- rep(0,p*p)
      if (deriv > 0) {
        M <- ncol(d1b)
        d1H <- if (deriv==1) rep(0,p*M) else rep(0,p*p*M)
      } else M <- d1Ho <- d1H <- 0
      if (deriv > 2) {
        d2H <- rep(0,p*M*(M+1)/2)
        if (is.list(fh)) {
          ev <- fh
        } else  { ## need to compute eigen here
          ev <- eigen(Hp,symmetric=TRUE)
          if (rank < p) ev$values[(rank+1):p] <- 0
        } 
        X <- X%*%(ev$vectors*D)
        d1b <- t(ev$vectors)%*%(d1b/D); d2b <- t(ev$vectors)%*%(d2b/D)
      } else d2Ho <- trHid2H <- d2H <- 0

      for (j in 1:ns) { ## loop over strata
        ind <- if (ns==1) 1:length(y) else which(strat==strat.lev[j]) ## index for points in this strata
        tr <- unique(y[ind])
        r <- match(y[ind],tr)
        ## note that the following call can not use .C(C_coxlpl,...) since the ll
        ## function is not in the mgcv namespace.
        oo <- .C("coxlpl",as.double(mu[ind]),as.double(X[ind,]),as.integer(r),as.integer(wt[ind]),
            as.double(tr),n=as.integer(length(y[ind])),p=as.integer(p),nt=as.integer(length(tr)),
            lp=as.double(0),g=as.double(g),H=as.double(H),
            d1b=as.double(d1b),d1H=as.double(d1H),d2b=as.double(d2b),d2H=as.double(d2H),
            n.sp=as.integer(M),deriv=as.integer(deriv),PACKAGE="mgcv");
	if (j==1) {
          lp <- oo$lp
	  lb <- oo$g
	  lbb <- matrix(oo$H,p,p)
        } else { ## accumulating over strata
          lp <- oo$lp + lp
	  lb <- oo$g + lb
	  lbb <- matrix(oo$H,p,p) + lbb
	}
        if (deriv==1) { d1Ho <- matrix(oo$d1H,p,M) + if (j==1) 0 else d1Ho } else
        if (deriv>1) {
          ind <- 1:(p^2)
          if (j==1) d1Ho <- list()
          for (i in 1:M) { 
            d1Ho[[i]] <- if (j==1) matrix(oo$d1H[ind],p,p) else matrix(oo$d1H[ind],p,p) + d1Ho[[i]]
            ind <- ind + p^2
          }
        } 
        if (deriv>2) { 
          d2Ho <- if (j==1) matrix(oo$d2H,p,M*(M+1)/2) else d2Ho + matrix(oo$d2H,p,M*(M+1)/2)
          #d <- ev$values
          #d[d>0] <- 1/d[d>0];d[d<=0] <- 0
          #trHid2H <- colSums(d2H*d)
        }
      } ## strata loop
      if (deriv>2) {
        d <- ev$values
        d[d>0] <- 1/d[d>0];d[d<=0] <- 0
        trHid2H <- colSums(d2Ho*d)
      }
      assign(".log.partial.likelihood", oo$lp, envir=environment(sys.function()))
      list(l=lp,lb=lb,lbb=lbb,d1H=d1Ho,d2H=d2Ho,trHid2H=trHid2H)
    } ## ll

    # environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    # environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) 
    #environment(aic) <- 
    environment(ll) <- env
    structure(list(family = "Cox PH", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, ll=ll, ## aic = aic, 
        mu.eta = stats$mu.eta, 
        initialize = initialize,preinitialize=preinitialize,postproc=postproc,
        hazard=hazard,predict=predict,residuals=residuals,
        validmu = validmu, valideta = stats$valideta, 
        rd=rd,qf=qf,drop.intercept = TRUE,
        ls=1, ## signal ls not needed
        available.derivs = 2 ## can use full Newton here
        ),
        class = c("general.family","extended.family","family"))
} ## cox.ph


## (c) Simon N. Wood (2013, 2014) mvn model extended family. 
## Released under GPL2 ...

mvn <- function(d=2) { 
## Extended family object for multivariate normal additive model.
  if (d<2) stop("mvn requires 2 or more dimensional data")
  stats <- list()
  for (i in 1:d) {
    stats[[i]] <- make.link("identity") 
  }
  
  env <- new.env(parent = .GlobalEnv)
  validmu <- function(mu) all(is.finite(mu))

  ## initialization has to add in the extra parameters of 
  ## the cov matrix...
  
    preinitialize <- expression({
    ## code to evaluate in estimate.gam...
    ## extends model matrix with dummy columns and 
    ## finds initial coefficients
      ydim <- ncol(G$y) ## dimension of response
      nbeta <- ncol(G$X)
      ntheta <- ydim*(ydim+1)/2 ## numer of cov matrix factor params
      lpi <- attr(G$X,"lpi")
      XX <- crossprod(G$X)
      G$X <- cbind(G$X,matrix(0,nrow(G$X),ntheta)) ## add dummy columns to G$X
      #G$cmX <- c(G$cmX,rep(0,ntheta)) ## and corresponding column means
      G$term.names <- c(G$term.names,paste("R",1:ntheta,sep="."))
      attr(G$X,"lpi") <- lpi
      attr(G$X,"XX") <- XX
      ## pad out sqrt of balanced penalty matrix to account for extra params
      attr(G$Sl,"E") <- cbind(attr(G$Sl,"E"),matrix(0,nbeta,ntheta))
      G$family$data <- list(ydim = ydim,nbeta=nbeta)
      G$family$ibeta = rep(0,ncol(G$X))
      ## now get initial parameters and store in family...
      for (k in 1:ydim) {
        sin <- G$off %in% lpi[[k]]
        Sk <- G$S[sin]
        um <- magic(G$y[,k],G$X[,lpi[[k]]],rep(-1,sum(sin)),G$S[sin],G$off[sin]-lpi[[k]][1]+1,nt=control$nthreads)
        G$family$ibeta[lpi[[k]]] <- um$b
        G$family$ibeta[nbeta+1] <- -.5*log(um$scale) ## initial log root precision
        nbeta <- nbeta + ydim - k + 1
      }
    })
    
    postproc <- expression({
    ## code to evaluate in estimate.gam, to do with estimated factor of
    ## precision matrix, etc...
      ydim <- G$family$data$ydim
      R <- matrix(0,ydim,ydim)
      ind <- G$family$data$nbeta + 1:(ydim*(ydim+1)/2);
      theta <- object$coefficients[ind]
      k <- 1;for (i in 1:ydim) for (j in i:ydim) {
        if (i==j) R[i,j] <- exp(theta[k]) else R[i,j] <- theta[k]
        k <- k + 1
      }
      object$family$data <- list(R=R) 
      rsd <- R%*%t(object$y-object$fitted.values)
      object$deviance <- sum(rsd^2)
      rsd <- R%*%(t(object$y)-colMeans(object$y))
      object$null.deviance <- sum(rsd^2)
    })
    
    initialize <- expression({
      ## called in gam.fit5 and initial.spg
      ## Ideally fit separate models to each component and
      ## extract initial coefs, s.p.s and variances this way 
        n <- rep(1, nobs)
        if (is.null(start)) start <- family$ibeta
        ## need to re-parameterize XX is non-standard
        if (exists("rp",inherits=FALSE)&&length(rp$rp)>0) 
           attr(x,"XX") <- Sl.repara(rp$rp,t(Sl.repara(rp$rp,attr(x,"XX"))))
    })


    residuals <- function(object,type=c("response","deviance")) {
      type <- match.arg(type)
      res <- object$y - object$fitted.values
      if (type=="deviance") res <- res%*%t(object$family$data$R)
      res 
    } ## residuals


    rd <- qf <- NULL ## these functions currently undefined for Cox PH

    ll <- function(y,X,coef,wt,family,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## function defining the Multivariate Normal model log lik.
    ## Calls C code "mvn_ll"
    ## deriv codes: 0 - eval; 1 - grad and Hessian
    ##              2 - d1H (diagonal only - not implemented efficiently)
    ##              3 - d1H; 4 d2H (diag - not implemented)
    ## Hp is the preconditioned penalized hessian of the log lik
    ##    which is of rank 'rank'.
    ## fh is a factorization of Hp - either its eigen decomposition 
    ##    or its Choleski factor
    ## D is the diagonal pre-conditioning matrix used to obtain Hp
    ##   if Hr is the raw Hp then Hp = D*t(D*Hr)
      lpi <- attr(X,"lpi") ## lpi[[k]] is index of model matrix columns for dim k 
      m <- length(lpi)        ## number of dimensions of MVN
      lpstart <- rep(0,m)
      for (i in 1:(m-1)) lpstart[i] <- lpi[[i+1]][1]
      lpstart[m] <- lpi[[m]][length(lpi[[m]])]+1 
      nb <- length(coef)      ## total number of parameters
      if (deriv<2) {
        nsp = 0;d1b <- dH <- 0
      } else {
        nsp = ncol(d1b)
        dH = rep(0,nsp*nb*nb)
      }
      #cat("\nderiv=",deriv,"  lpstart=",lpstart," dim(y) = ",dim(y),
      #    "\ndim(XX)=",dim(attr(X,"XX"))," m=",m," nsp=",nsp,"\n")
      oo <- .C("mvn_ll",y=as.double(t(y)),X=as.double(X),XX=as.double(attr(X,"XX")),
               beta=as.double(coef),n=as.integer(nrow(X)),
               lpi=as.integer(lpstart-1),m=as.integer(m),ll=as.double(0),lb=as.double(coef*0),
               lbb=as.double(rep(0,nb*nb)), dbeta = as.double(d1b), dH = as.double(dH), 
               deriv = as.integer(nsp>0),nsp = as.integer(nsp),nt=as.integer(1),PACKAGE="mgcv")
      if (nsp==0) d1H <- NULL else if (deriv==2) {
        d1H <- matrix(0,nb,nsp)
        for (i in 1:nsp) { 
          d1H[,i] <- diag(matrix(oo$dH[ind],nb,nb))
          ind <- ind + nb*nb
        }
      } else { ## deriv==3
        d1H <- list();ind <- 1:(nb*nb)
        for (i in 1:nsp) { 
          d1H[[i]] <- matrix(oo$dH[ind],nb,nb)
          ind <- ind + nb*nb
        }
      }
      list(l=oo$ll,lb=oo$lb,lbb=matrix(oo$lbb,nb,nb),d1H=d1H)
    } ## ll

    # environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    # environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) 
    ##environment(aic) <- 
    environment(ll) <- env
    structure(list(family = "Multivariate normal", 
        ## link = linktemp, linkfun = stats$linkfun, linkinv = stats$linkinv, 
        ll=ll,nlp=d,
        initialize = initialize,preinitialize=preinitialize,postproc=postproc,
        residuals=residuals,
        validmu = validmu, ## valideta = stats$valideta, 
        ## rd=rd,qf=qf,  
        linfo = stats, ## link information list
        d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
        ls=1, ## signal ls not needed
        available.derivs = 1 ## signal only first derivatives available...
        ),
        class = c("general.family","extended.family","family"))
} ## mvn


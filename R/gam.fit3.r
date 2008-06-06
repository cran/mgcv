## R routines for gam fitting with calculation of derivatives w.r.t. sp.s
## (c) Simon Wood 2004-2008

## This routine is for type 3 gam fitting. The basic idea is that a P-IRLS
## is run to convergence, and only then is a scheme for evaluating the 
## derivatives iterated to convergence. The advantage is that many key
## quantities are fixed at this stage, including the key decompositions
## In addition the R side work is simplified considerably.The routine
## evaluates first and second derivatives of the deviance and tr(A).


gam.fit3 <- function (x, y, sp, S=list(),rS=list(),off, H=NULL, 
            weights = rep(1, nobs), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
            control = gam.control(), intercept = TRUE,deriv=2,use.svd=TRUE,
            gamma=1,scale=1,printWarn=TRUE,...) 
## deriv, sp, S, rS, H added to arg list. 
## need to modify family before call.
{   scale <- abs(scale)
    if (!deriv%in%c(0,1,2)) stop("unsupported order of differentiation requested of gam.fit3")
    x <- as.matrix(x)
    iter <- 0;coef <- rep(0,ncol(x))
    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    n <- nobs <- NROW(y) ## n is just to keep codetools happy
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
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

    ## Added code
    nSp <- length(S)
    if (nSp==0) deriv <- FALSE 
    St <- totalPenalty(S,H,off,sp,ncol(x))
    Sr <- mroot(St)

    ## end of added code

    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("Invalid fitted means in empty model")
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
        V <- variance(mu)
        alpha <- dev
        trA2 <- trA1 <- trA <- 0
        if (deriv) GCV2 <- GCV1<- UBRE2 <- UBRE1<-trA1 <- rep(0,nSp)
        else GCV2<-GCV1<-UBRE2<-UBRE1<-trA2<-trA1 <- NULL
        GCV <- nobs*alpha/(nobs-gamma*trA)^2
        UBRE <- alpha/nobs - scale + 2*gamma/n*trA
        scale.est <- alpha / (nobs - trA)
    } ### end if (EMPTY)
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop("Length of start should equal ", nvars, 
                  " and correspond to initial coefs for ", deparse(xnames))
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1) 
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        etaold <- eta
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("Can't find valid starting values: please specify some")
    
        boundary <- conv <- FALSE
        rV=matrix(0,ncol(x),ncol(x))   
        old.pdev <- 0     
        for (iter in 1:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("No observations informative at iteration ", 
                  iter)
                break
            }
            mevg<-mu.eta.val[good];mug<-mu[good];yg<-y[good];weg<-weights[good]
            z <- (eta - offset)[good] + (yg - mug)/mevg
            var.mug<-variance(mug)
            w <- sqrt((weg * mevg^2)/var.mug)

            ## Here a Fortran call has been replaced by update.beta call
           
            if (sum(good)<ncol(x)) stop("Not enough informative observations.")
     
            oo<-.C(C_pls_fit,y=as.double(z),as.double(x[good,]),as.double(w),as.double(Sr),as.integer(sum(good)),
            as.integer(ncol(x)),as.integer(ncol(Sr)),eta=as.double(z),penalty=as.double(1),
            as.double(.Machine$double.eps*100))
       
            start <- oo$y[1:ncol(x)];
            penalty <- oo$penalty
            eta <- drop(x%*%start)

            if (any(!is.finite(start))) {
                conv <- FALSE
                warning("Non-finite coefficients at iteration ", 
                  iter)
                break
            }

           
     
           mu <- linkinv(eta <- eta + offset)
           dev <- sum(dev.resids(y, mu, weights))
          
           if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, 
                  "\n")
            boundary <- FALSE
            
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found:please supply starting values", 
                    call. = FALSE)
                warning("Step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- (eta + etaold)/2               
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                warning("Step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- (eta + etaold)/2 
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }

            pdev <- dev + penalty  ## the penalized deviance 

            if (control$trace) 
                  cat("penalized deviance =", pdev, "\n")

            div.thresh <- 10*(.1+abs(old.pdev))*.Machine$double.eps^.5 
            ## ... threshold for judging divergence --- too tight and near
            ## perfect convergence can cause a failure here

            if (iter>1&&(pdev-old.pdev>div.thresh)) { ## solution diverging
             ii <- 1 ## step halving counter
             while (pdev -old.pdev > div.thresh)  
             { ## step halve until pdev <= old.pdev
                if (ii > 200) 
                   stop("inner loop 3; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2 
                eta <- (eta + etaold)/2               
                mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                  pdev <- dev + t(start)%*%St%*%start ## the penalized deviance
                if (control$trace) 
                  cat("Step halved: new penalized deviance =", pdev, "\n")
              }
            } 

            if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) {
                if (max(abs(start-coefold))>control$epsilon*max(abs(start+coefold))/2){
                  old.pdev <- pdev
                  coef <- coefold <- start
                  etaold <- eta 
                } else {
                  conv <- TRUE
                  coef <- start
                  break 
                }
            }
            else {  old.pdev <- pdev
                coef <- coefold <- start
                etaold <- eta 
            }
        } ### end main loop 
       
        dev <- sum(dev.resids(y, mu, weights)) 
       
        ## Now call the derivative calculation scheme. This requires the
        ## following inputs:
        ## z and w - the pseudodata and weights
        ## X the model matrix and E where EE'=S
        ## rS the single penalty square roots
        ## sp the log smoothing parameters
        ## y and mu the data and model expected values
        ## g1,g2,g3 - the first 3 derivatives of g(mu) wrt mu
        ## V,V1,V2 - V(mu) and its first two derivatives wrt mu
        ## on output it returns the gradient and hessian for
        ## the deviance and trA 

         good <- weights > 0
         varmu <- variance(mu)[good]
         if (any(is.na(varmu))) stop("NAs in V(mu)")
         if (any(varmu == 0)) stop("0s in V(mu)")
         mu.eta.val <- mu.eta(eta)
         if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
         good <- (weights > 0) & (mu.eta.val != 0)
   
         mevg <- mu.eta.val[good];mug <- mu[good];yg <- y[good]
         weg <- weights[good];etag <- eta[good]
         z <- (eta - offset)[good] + (yg - mug)/mevg
         var.mug<-variance(mug)
         w <- sqrt((weg * mevg^2)/var.mug)
        
         g1 <- 1/mevg
         g2 <- family$d2link(mug)
         g3 <- family$d3link(mug)

         V <- family$variance(mug)
         V1 <- family$dvar(mug)
         V2 <- family$d2var(mug)      

         P1 <- D1 <- array(0,nSp);P2 <- D2 <- matrix(0,nSp,nSp) # for derivs of deviance/ Pearson
         trA1 <- array(0,nSp);trA2 <- matrix(0,nSp,nSp) # for derivs of tr(A)
         rV=matrix(0,ncol(x),ncol(x));
         dum <- 1
         if (control$trace) cat("calling gdi...")

       oo <- .C(C_gdi,X=as.double(x[good,]),E=as.double(Sr),rS = as.double(unlist(rS)),
           sp=as.double(exp(sp)),z=as.double(z),w=as.double(w),mu=as.double(mug),eta=as.double(etag),y=as.double(yg),
           p.weights=as.double(weg),g1=as.double(g1),g2=as.double(g2),g3=as.double(g3),V0=as.double(V),
           V1=as.double(V1),V2=as.double(V2),beta=as.double(coef),D1=as.double(D1),D2=as.double(D2),
           P=as.double(dum),P1=as.double(P1),P2=as.double(P2),trA=as.double(dum),
           trA1=as.double(trA1),trA2=as.double(trA2),rV=as.double(rV),rank.tol=as.double(.Machine$double.eps*100),
           conv.tol=as.double(control$epsilon),rank.est=as.integer(1),n=as.integer(length(z)),
           p=as.integer(ncol(x)),M=as.integer(nSp),Encol = as.integer(ncol(Sr)),
           rSncol=as.integer(unlist(lapply(rS,ncol))),deriv=as.integer(deriv),use.svd=as.integer(use.svd))      
       
         if (control$trace) cat("done!\n")
 
         rV <- matrix(oo$rV,ncol(x),ncol(x))
         coef <- oo$beta;
         trA <- oo$trA;
         P <- oo$P
         
         delta <- nobs - gamma * trA
         delta.2 <- delta*delta           
  
         GCV <- nobs*dev/delta.2
         GACV <- dev/nobs + P * 2*gamma*trA/(delta * nobs) 

         UBRE <- dev/nobs - 2*delta*scale/nobs + scale
         scale.est <- dev/(nobs-trA)

         if (deriv) {
           trA1 <- oo$trA1
           
           D1 <- oo$D1
           P1 <- oo$P1
          
           if (sum(!is.finite(D1))||sum(!is.finite(P1))||sum(!is.finite(trA1))) { 
             stop("Smoothing parameter derivate iteration diverging. Decrease fit tolerance! See `epsilon' in `gam.contol'")}
         
           delta.3 <- delta*delta.2

           GCV1 <- nobs*D1/delta.2 + 2*nobs*dev*trA1*gamma/delta.3
           GACV1 <- D1/nobs + 2*P/delta.2 * trA1 + 2*gamma*trA*P1/(delta*nobs)

           UBRE1 <- D1/nobs + gamma * trA1 *2*scale/nobs
           if (deriv==2) {
             trA2 <- matrix(oo$trA2,nSp,nSp) 
             D2 <- matrix(oo$D2,nSp,nSp)
             P2 <- matrix(oo$P2,nSp,nSp)
              
             if (sum(!is.finite(D2))||sum(!is.finite(P2))||sum(!is.finite(trA2))) { 
               stop("Smoothing parameter derivate iteration diverging. Decrease fit tolerance! See `epsilon' in `gam.contol'")}
             
             GCV2 <- outer(trA1,D1)
             GCV2 <- (GCV2 + t(GCV2))*gamma*2*nobs/delta.3 +
                      6*nobs*dev*outer(trA1,trA1)*gamma*gamma/(delta.2*delta.2) + 
                      nobs*D2/delta.2 + 2*nobs*dev*gamma*trA2/delta.3  
             GACV2 <- D2/nobs + outer(trA1,trA1)*4*P/(delta.3) +
                      2 * P * trA2 / delta.2 + 2 * outer(trA1,P1)/delta.2 +
                      2 * outer(P1,trA1) *(1/(delta * nobs) + trA/(nobs*delta.2)) +
                      2 * trA * P2 /(delta * nobs) 
             GACV2 <- (GACV2 + t(GACV2))*.5
             UBRE2 <- D2/nobs +2*gamma * trA2 * scale / nobs
           } else {GACV2<-P2<-trA2<-D2<-UBRE2<-GCV2 <- NULL}
         } else {
           GACV1<-P1<-P2<-trA1<-trA2<-D1<-D2<-GACV2<-UBRE2<-GCV2<-UBRE1<-GCV1<-NULL
         }
         
        # end of inserted code
        if (!conv&&printWarn) 
            warning("Algorithm did not converge")
        if (printWarn&&boundary) 
            warning("Algorithm stopped at boundary value")
        eps <- 10 * .Machine$double.eps
        if (printWarn&&family$family[1] == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("fitted probabilities numerically 0 or 1 occurred")
        }
        if (printWarn&&family$family[1] == "poisson") {
            if (any(mu < eps)) 
                warning("fitted rates numerically 0 occurred")
        }
 
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
          
        names(coef) <- xnames 
    } ### end if (!EMPTY)
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
   
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
   
    aic.model <- aic(y, n, mu, weights, dev) # note: incomplete 2*edf needs to be added

    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
         family = family, linear.predictors = eta, deviance = dev, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.null = nulldf, y = y, converged = conv,
        boundary = boundary,D1=D1,D2=D2,P=P,P1=P1,P2=P2,trA=trA,trA1=trA1,trA2=trA2,
        GCV=GCV,GCV1=GCV1,GCV2=GCV2,GACV=GACV,GACV1=GACV1,GACV2=GACV2,UBRE=UBRE,UBRE1=UBRE1,UBRE2=UBRE2,rV=rV,
        scale.est=scale.est,aic=aic.model,rank=oo$rank.est)
}

newton <- function(lsp,X,y,S,rS,off,L,lsp0,H,offset,family,weights,
                   control,gamma,scale,conv.tol=1e-6,maxNstep=5,maxSstep=2,
                   maxHalf=30,printWarn=FALSE,scoreType="deviance",use.svd=TRUE,
                   mustart = NULL,...)
## Newton optimizer for GAM gcv/aic optimization that can cope with an 
## indefinite Hessian! Main enhancements are: i) always peturbs the Hessian
## to +ve definite ii) step halves on step 
## failure, without obtaining derivatives until success; (iii) carries start
## values forward from one evaluation to next to speed convergence.    
## L is the matrix such that L%*%lsp + lsp0 gives the logs of the smoothing 
## parameters actually multiplying the S[[i]]'s
{ ## sanity check L
  if (is.null(L)) L <- diag(length(S)) else {
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=length(S)||ncol(L)!=length(lsp)) stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) lsp0 <- rep(0,ncol(L))
  ## initial fit
  b<-gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0, S=S,rS=rS,off=off, H=H,
     offset = offset,family = family,weights=weights,deriv=2,
     control=control,gamma=gamma,scale=scale,
     printWarn=FALSE,use.svd=use.svd,mustart=mustart,...)

  mustart<-b$fitted.values

  if (scoreType=="GACV") {
    old.score <- score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2 
  } else if (scoreType=="UBRE"){
    old.score <- score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
  } else { ## default to deviance based GCV
    old.score <- score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2
  }
  
  grad <- t(L)%*%grad
  hess <- t(L)%*%hess%*%L

  Slength <- maxSstep 
  score.scale <- b$scale.est + score;    
  uconv.ind <- abs(grad) > score.scale*conv.tol
  ## check for all converged too soon, and undo !
  if (!sum(uconv.ind)) uconv.ind <- uconv.ind | TRUE
  for (i in 1:200) {
    ## exclude apparently converged gradients from computation
    hess1 <- hess[uconv.ind,uconv.ind] 
    grad1 <- grad[uconv.ind]
    ## get the trial step ...
    eh <- eigen(hess1,symmetric=TRUE)
    d <- eh$values;U <- eh$vectors
    ind <- d < 0
    d[ind] <- -d[ind] ## see Gill Murray and Wright p107/8
    d <- 1/d
    
    Nstep <- 0 * grad
    Nstep[uconv.ind] <- -drop(U%*%(d*(t(U)%*%grad1))) # (modified) Newton direction
   
    Sstep <- -Slength * grad/max(abs(grad)) # steepest descent direction 
    
    ms <- max(abs(Nstep))
    if (ms>maxNstep) Nstep <- maxNstep * Nstep/ms

    ## try the step ...
    lsp1 <- lsp + Nstep
    b<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
       offset = offset,family = family,weights=weights,deriv=2,
       control=control,gamma=gamma,scale=scale,
       printWarn=FALSE,mustart=mustart,use.svd=use.svd,...)
    if (scoreType=="GACV") {
      score1 <- b$GACV
    } else if (scoreType=="UBRE") {
      score1 <- b$UBRE
    } else score1 <- b$GCV
    ## accept if improvement, else step halve
    ii <- 0 ## step halving counter
    if (score1<score) { ## accept
      old.score <- score 
      mustart <- b$fitted.values
      lsp <- lsp1
      if (scoreType=="GACV") {
          score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2
      } else if (scoreType=="UBRE") {
          score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
      } else { score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2} 
      grad <- t(L)%*%grad
      hess <- t(L)%*%hess%*%L
    } else { ## step halving ...
      step <- Nstep ## start with the (pseudo) Newton direction
      while (score1>score && ii < maxHalf) {
        if (ii==3) { ## Newton really not working
          step <- Sstep ## use steepest descent direction
        } else step <- step/2
        if (ii>3) Slength <- Slength/2 ## keep track of SD step length
        lsp1 <- lsp + step
        b1<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
           offset = offset,family = family,weights=weights,deriv=0,
           control=control,gamma=gamma,scale=scale,
           printWarn=FALSE,mustart=mustart,use.svd=use.svd,...)
       
        if (scoreType=="GACV") {
          score1 <- b1$GACV
        } else if (scoreType=="UBRE") {
          score1 <- b1$UBRE
        } else score1 <- b1$GCV

        if (score1 <= score) { ## accept
          b<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
             offset = offset,family = family,weights=weights,deriv=2,
             control=control,gamma=gamma,scale=scale,
             printWarn=FALSE,mustart=mustart,use.svd=use.svd,...)
          mustart <- b$fitted.values
          old.score <- score;lsp <- lsp1
          if (scoreType=="GACV") {
            score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2
          } else if (scoreType=="UBRE") {
            score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
          } else { score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2}
          grad <- t(L)%*%grad
          hess <- t(L)%*%hess%*%L
          if (ii>3) Slength <- min(Slength*2,maxSstep) ## try increasing SD step length
        }  # end of if (score1<= score )
        ii <- ii + 1
      } # end of step halving
    }
    ## test for convergence
    converged <- TRUE
    score.scale <- b$scale.est + score;    
    uconv.ind <- abs(grad) > score.scale*conv.tol
    if (sum(uconv.ind)) converged <- FALSE
    if (abs(old.score-score)>score.scale*conv.tol) { 
      if (converged) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
      converged <- FALSE      
    }
    if (ii==maxHalf) converged <- TRUE ## step failure
    if (converged) break
  } ## end of iteration loop
  if (ii==maxHalf) ct <- "step failed"
  else if (i==200) ct <- "iteration limit reached" 
  else ct <- "full convergence"
  list(score=score,lsp=lsp,lsp.full=L%*%lsp,grad=grad,hess=hess,iter=i,conv =ct,object=b)
}

bfgs <- function(lsp,X,y,S,rS,off,L,lsp0,H,offset,family,weights,
                   control,gamma,scale,conv.tol=1e-6,maxNstep=5,maxSstep=2,
                   maxHalf=30,printWarn=FALSE,scoreType="deviance",use.svd=TRUE,
                   mustart = NULL,...)
## This optimizer is experimental... The main feature is to alternate infrequent 
## Newton steps with BFGS Quasi-Newton steps. In theory this should be faster 
## than Newton, because of the cost of full Hessian calculation, but
## in practice the extra steps required by QN tends to mean that the advantage
## is not realized...
## Newton optimizer for GAM gcv/aic optimization that can cope with an 
## indefinite Hessian, and alternates BFGS and Newton steps for speed reasons
## Main enhancements are: i) always peturbs the Hessian
## to +ve definite ii) step halves on step 
## failure, without obtaining derivatives until success; (iii) carries start
## values forward from one evaluation to next to speed convergence.    
## L is the matrix such that L%*%lsp + lsp0 gives the logs of the smoothing 
## parameters actually multiplying the S[[i]]'s
{ ## sanity check L
  if (is.null(L)) L <- diag(length(S)) else {
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=length(S)||ncol(L)!=length(lsp)) stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) lsp0 <- rep(0,ncol(L))
  ## initial fit
#  ptm <- proc.time()
  b<-gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0, S=S,rS=rS,off=off, H=H,
     offset = offset,family = family,weights=weights,deriv=2,
     control=control,gamma=gamma,scale=scale,
     printWarn=FALSE,use.svd=use.svd,mustart=mustart,...)
#  ptm <- proc.time()-ptm
#  cat("deriv=2 ",ptm,"\n")

  mustart<-b$fitted.values

  QNsteps <- floor(length(S)/2) ## how often to Newton should depend on cost...

  if (scoreType=="GACV") {
    old.score <- score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2 
  } else if (scoreType=="UBRE"){
    old.score <- score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
  } else { ## default to deviance based GCV
    old.score <- score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2
  }
  
  grad <- t(L)%*%grad
  hess <- t(L)%*%hess%*%L

  score.scale <- b$scale.est + score;    
  uconv.ind <- abs(grad) > score.scale*conv.tol
  ## check for all converged too soon, and undo !
  if (!sum(uconv.ind)) uconv.ind <- uconv.ind | TRUE
  kk <- 0 ## counter for QN steps between Newton steps
  for (i in 1:200) {
   
    if (kk==0) { ## time to reset B
      eh <- eigen(hess,symmetric=TRUE)
      d <- eh$values;U <- eh$vectors
      ind <- d < 0
      d[ind] <- -d[ind] ## see Gill Murray and Wright p107/8
      d <- 1/d
      d[d==0] <- min(d)*.Machine$double.eps^.5
      B <- U%*%(d*t(U)) ## Newton based inverse Hessian
    }
     
    kk <- kk + 1
    if (kk > QNsteps) kk <- 0 
 
    ## get the trial step ...
    
    Nstep <- 0 * grad
    Nstep[uconv.ind] <- -drop(B[uconv.ind,uconv.ind]%*%grad[uconv.ind]) # (modified) Newton direction
    
    ms <- max(abs(Nstep))
    if (ms>maxNstep) Nstep <- maxNstep * Nstep/ms

    ## try the step ...
    sc.extra <- 1e-4*sum(grad*Nstep) ## -ve sufficient decrease 
    ii <- 0 ## step halving counter
    step <- Nstep*2
    score1 <- abs(score)*2
    while (score1>score+sc.extra && ii < maxHalf) { ## reject and step halve
      ii <- ii + 1
      step <- step/2
      lsp1 <- lsp + step
  
#      ptm <- proc.time()
      if (kk!=0||ii==1) deriv <- 1 else deriv <- 0
      b1<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
          offset = offset,family = family,weights=weights,deriv=deriv,
          control=control,gamma=gamma,scale=scale,
          printWarn=FALSE,mustart=mustart,use.svd=use.svd,...)
#       ptm <- proc.time()-ptm
#       cat("deriv= ",deriv,"  ",ptm,"\n")
      
 
      if (scoreType=="GACV") {
          score1 <- b1$GACV
      } else if (scoreType=="UBRE") {
          score1 <- b1$UBRE
      } else score1 <- b1$GCV
    } ## accepted step or step failed to lead to decrease

    if (ii < maxHalf) { ## step succeeded 
      mustart <- b1$fitted.values
      if (kk==0) { ## time for a full Newton step ...
#        ptm <- proc.time()
        b<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
               offset = offset,family = family,weights=weights,deriv=2,
               control=control,gamma=gamma,scale=scale,
               printWarn=FALSE,mustart=mustart,use.svd=use.svd,...)
#         ptm <- proc.time()-ptm
#         cat("deriv=2 ",ptm,"\n")

        mustart <- b$fitted.values
        old.score <- score;lsp <- lsp1
        if (scoreType=="GACV") {
          score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2
        } else if (scoreType=="UBRE") {
          score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
        } else { score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2}
        grad <- t(L)%*%grad
        hess <- t(L)%*%hess%*%L
      } else { ## just a BFGS update
        ## first derivatives only.... 
#        ptm <- proc.time()
         if (ii==1) b <- b1 else  
         b<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
               offset = offset,family = family,weights=weights,deriv=1,
               control=control,gamma=gamma,scale=scale,
               printWarn=FALSE,mustart=mustart,use.svd=use.svd,...)
#         ptm <- proc.time()-ptm
#         cat("deriv=1 ",ptm,"\n")

        mustart <- b$fitted.values
        old.score <- score;lsp <- lsp1
        old.grad <- grad
        if (scoreType=="GACV") {
          score <- b$GACV;grad <- b$GACV1
        } else if (scoreType=="UBRE") {
          score <- b$UBRE;grad <- b$UBRE1
        } else { score <- b$GCV;grad <- b$GCV1}
        grad <- t(L)%*%grad
        ## BFGS update of the inverse Hessian...
        yg <- grad-old.grad
        rho <- 1/sum(yg*step)
        B <- B - rho*step%*%(t(yg)%*%B)
        B <- B - rho*(B%*%yg)%*%t(step) + rho*step%*%t(step)
      } ## end of BFGS
    } ## end of successful step updating
    ## test for convergence
    converged <- TRUE
    score.scale <- b$scale.est + score;    
    uconv.ind <- abs(grad) > score.scale*conv.tol
    if (sum(uconv.ind)) converged <- FALSE
    if (abs(old.score-score)>score.scale*conv.tol) { 
      if (converged) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
      converged <- FALSE      
    }
    if (ii==maxHalf) converged <- TRUE ## step failure
    if (converged) break
  } ## end of iteration loop
  if (ii==maxHalf) ct <- "step failed"
  else if (i==200) ct <- "iteration limit reached" 
  else ct <- "full convergence"
  list(score=score,lsp=lsp,lsp.full=L%*%lsp,grad=grad,hess=hess,iter=i,conv =ct,object=b)
}



gam2derivative <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the derivatives of the GCV or UBRE score w.r.t the 
## smoothing parameters for the model.
## args is a list containing the arguments for gam.fit3
## For use as optim() objective gradient
{ if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=1,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     use.svd=FALSE,...)
  if (args$scoreType == "deviance") ret <- b$GCV1 else ret <- b$UBRE1
  if (!is.null(args$L)) ret <- t(args$L)%*%ret
  ret
}


gam2objective <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the GCV or UBRE score for the model.
## args is a list containing the arguments for gam.fit3
## For use as optim() objective
{ if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=0,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     use.svd=FALSE,...)
  if (args$scoreType == "deviance") ret <- b$GCV else ret <- b$UBRE
  attr(ret,"full.fit") <- b
  ret
}



gam4objective <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the GCV or UBRE score for the model.
## args is a list containing the arguments for gam.fit3
## For use as nlm() objective
{ if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=1,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     use.svd=FALSE,...)
  if (args$scoreType == "deviance") ret <- b$GCV else ret <- b$UBRE
  attr(ret,"full.fit") <- b
  if (args$scoreType == "deviance") at <- b$GCV1 else at <- b$UBRE1
  if (!is.null(args$L)) at <- t(args$L)%*%at
  attr(ret,"gradient") <- at
  ret
}

##
## The following fix up family objects for use with gam.fit3
##


fix.family.link<-function(fam)
# adds d2link the second derivative of the link function w.r.t. mu
# to the family supplied, as well as a 3rd derivative function 
# d3link...
# All d2link and d3link functions have been checked numerically. 
{ if (!inherits(fam,"family")) stop("fam not a family object")
  if (!is.null(fam$d2link)&&!is.null(fam$d3link)) return(fam) 
  link <- fam$link
  if (length(link)>1) if (fam$family=="quasi") # then it's a power link
  { lambda <- log(fam$linkfun(exp(1))) ## the power, if > 0
    if (lambda<=0) { fam$d2link <- function(mu) -1/mu^2
      fam$d3link <- function(mu) 2/mu^3
    }
    else { fam$d2link <- function(mu) lambda*(lambda-1)*mu^(lambda-2)
      fam$d3link <- function(mu) (lambda-2)*(lambda-1)*lambda*mu^(lambda-3)
    }
    return(fam)
  } else stop("unrecognized (vector?) link")

  if (link=="identity") {
    fam$d3link <- fam$d2link <- function(mu) rep.int(0,length(mu))
    return(fam)
  } 
  if (link == "log") {
    fam$d2link <- function(mu) -1/mu^2
    fam$d3link <- function(mu) 2/mu^3
    return(fam)
  }
  if (link == "inverse") {
    fam$d2link <- function(mu) 2/mu^3
    fam$d3link <- function(mu) {mu <- mu*mu;-6/(mu*mu)}
    return(fam)
  }
  if (link == "logit") {
    fam$d2link <- function(mu) 1/(1 - mu)^2 - 1/mu^2
    fam$d3link <- function(mu) 2/(1 - mu)^3 + 2/mu^3
    return(fam)
  }
  if (link == "probit") {
    fam$d2link <- function(mu) { 
      eta <- fam$linkfun(mu)
      eta/fam$mu.eta(eta)^2
    }
    fam$d3link <- function(mu) {
      eta <-  fam$linkfun(mu)
      (1 + 2*eta^2)/fam$mu.eta(eta)^3
    }
    return(fam)
  }
  if (link == "cloglog") {
    fam$d2link <- function(mu) { l1m <- log(1-mu)
      -1/((1 - mu)^2*l1m) *(1+ 1/l1m)
    }
    fam$d3link <- function(mu) { l1m <- log(1-mu)
       mu3 <- (1-mu)^3
      -1/(mu3 * l1m^3) -(1 + 2 * l1m)/
       (mu3 * l1m^2) * (1 + 1/l1m)
    }
    return(fam)
  }
  if (link == "sqrt") {
    fam$d2link <- function(mu) -.25 * mu^-1.5
    fam$d3link <- function(mu) .375 * mu^-2.5
    return(fam)
  }
  if (link == "cauchit") {
    fam$d2link <- function(mu) { 
     eta <- fam$linkfun(mu)
     2*pi*pi*eta*(1+eta*eta)
    }
    fam$d3link <- function(mu) { 
     eta <- fam$linkfun(mu)
     eta2 <- eta*eta
     2*pi*pi*pi*(1+3*eta2)*(1+eta2)
    }
    return(fam)
  }
  if (link == "1/mu^2") {
    fam$d2link <- function(mu) 6 * mu^-4
    fam$d3link <- function(mu) -24* mu^-5
    return(fam)
  }
  stop("link not recognised")
}


fix.family.var<-function(fam)
# adds dvar the derivative of the variance function w.r.t. mu
# to the family supplied, as well as d2var the 2nd derivative of 
# the variance function w.r.t. the mean. (All checked numerically). 
{ if (!inherits(fam,"family")) stop("fam not a family object")
  if (!is.null(fam$dvar)&&!is.null(fam$d2var)) return(fam) 
  family <- fam$family
  if (family=="gaussian") {
    fam$d2var <- fam$dvar <- function(mu) rep.int(0,length(mu))
    return(fam)
  } 
  if (family=="poisson"||family=="quasipoisson") {
    fam$dvar <- function(mu) rep.int(1,length(mu))
    fam$d2var <- function(mu) rep.int(0,length(mu))
    return(fam)
  } 
  if (family=="binomial"||family=="quasibinomial") {
    fam$dvar <- function(mu) 1-2*mu
    fam$d2var <- function(mu) rep.int(-2,length(mu))
    return(fam)
  }
  if (family=="Gamma") {
    fam$dvar <- function(mu) 2*mu
    fam$d2var <- function(mu) rep.int(2,length(mu))
    return(fam)
  }
  if (family=="quasi") {
    fam$dvar <- switch(fam$varfun,
       constant = function(mu) rep.int(0,length(mu)),
       "mu(1-mu)" = function(mu) 1-2*mu,
       mu = function(mu) rep.int(1,length(mu)),
       "mu^2" = function(mu) 2*mu,
       "mu^3" = function(mu) 3*mu^2           
    )
    if (is.null(fam$dvar)) stop("variance function not recognized for quasi")
    fam$d2var <- switch(fam$varfun,
       constant = function(mu) rep.int(0,length(mu)),
       "mu(1-mu)" = function(mu) rep.int(-2,length(mu)),
       mu = function(mu) rep.int(0,length(mu)),
       "mu^2" = function(mu) rep.int(2,length(mu)),
       "mu^3" = function(mu) 6*mu           
    )
    return(fam)
  }
  if (family=="inverse.gaussian") {
    fam$dvar <- function(mu) 3*mu^2
    fam$d2var <- function(mu) 6*mu
    return(fam)
  }
  stop("family not recognised")
}


negbin <- function (theta = stop("'theta' must be specified"), link = "log") { 
## modified from Venables and Ripley's MASS library to work with gam.fit3,
## and to allow a range of `theta' values to be specified
## single `theta' to specify fixed value; 2 theta values (first smaller that second)
## are limits within which to search for theta; otherwise supplied values make up 
## search set.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"")
    }
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", theta, envir = env)
    variance <- function(mu) mu + mu^2/get(".Theta")
    ## dvaraince/dmu needed as well
    dvar <- function(mu) 1 + 2*mu/get(".Theta")
    ## d2variance/dmu...
    d2var <- function(mu) rep(2/get(".Theta"),length(mu))
    getTheta <- function() get(".Theta")
    validmu <- function(mu) all(mu > 0)

    dev.resids <- function(y, mu, wt) { Theta <- get(".Theta")
      2 * wt * (y * log(pmax(1, y)/mu) - 
        (y + Theta) * log((y + Theta)/(mu + Theta))) 
    }
    aic <- function(y, n, mu, wt, dev) {
        Theta <- get(".Theta")
        term <- (y + Theta) * log(mu + Theta) - y * log(mu) +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
    environment(dvar) <- environment(d2var) <- environment(variance) <- environment(validmu) <- 
                         environment(dev.resids) <- environment(aic) <- environment(getTheta) <- env
    famname <- paste("Negative Binomial(", format(round(theta,3)), ")", sep = "")
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance,dvar=dvar,d2var=d2var, dev.resids = dev.resids,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
        validmu = validmu, valideta = stats$valideta,getTheta = getTheta), class = "family")
}





totalPenalty <- function(S,H,off,theta,p)
{ if (is.null(H)) St <- matrix(0,p,p)
  else { St <- H; 
    if (ncol(H)!=p||nrow(H)!=p) stop("H has wrong dimension")
  }
  theta <- exp(theta)
  m <- length(theta)
  if (m>0) for (i in 1:m) {
    k0 <- off[i]
    k1 <- k0 + nrow(S[[i]]) - 1
    St[k0:k1,k0:k1] <- St[k0:k1,k0:k1] + S[[i]] * theta[i]
  }
  St
}

mini.roots <- function(S,off,np)
# function to obtain square roots, B[[i]], of S[[i]]'s having as few
# columns as possible. S[[i]]=B[[i]]%*%t(B[[i]]). np is the total number
# of parameters. S is in packed form. 
{ m<-length(S)
  B<-S
  for (i in 1:m)
  { b<-mroot(S[[i]])
    B[[i]]<-matrix(0,np,ncol(b))
    B[[i]][off[i]:(off[i]+nrow(b)-1),]<-b
  }
  B
}



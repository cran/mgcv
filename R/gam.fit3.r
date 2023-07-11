## R routines for gam fitting with calculation of derivatives w.r.t. sp.s
## (c) Simon Wood 2004-2022

## These routines are for type 3 gam fitting. The basic idea is that a P-IRLS
## is run to convergence, and only then is a scheme for evaluating the 
## derivatives via the implicit function theorem used. 


gam.reparam <- function(rS,lsp,deriv) 
## Finds an orthogonal reparameterization which avoids `dominant machine zero leakage' between 
## components of the square root penalty.
## rS is the list of the square root penalties: last entry is root of fixed. 
##    penalty, if fixed.penalty=TRUE (i.e. length(rS)>length(sp))
## lsp is the vector of log smoothing parameters.
## *Assumption* here is that rS[[i]] are in a null space of total penalty already;
## see e.g. totalPenaltySpace & mini.roots
## Ouputs:
## S -- the total penalty matrix similarity transformed for stability
## rS -- the component square roots, transformed in the same way
##       - tcrossprod(rS[[i]]) = rS[[i]] %*% t(rS[[i]]) gives the matrix penalty component.
## Qs -- the orthogonal transformation matrix S = t(Qs)%*%S0%*%Qs, where S0 is the 
##       untransformed total penalty implied by sp and rS on input
## E -- the square root of the transformed S (obtained in a stable way by pre-conditioning)
## det -- log |S|
## det1 -- dlog|S|/dlog(sp) if deriv >0
## det2 -- hessian of log|S| wrt log(sp) if deriv>1  
{ q <- nrow(rS[[1]])
  rSncol <- unlist(lapply(rS,ncol))
  M <- length(lsp) 
  if (length(rS)>M) fixed.penalty <- TRUE else fixed.penalty <- FALSE
 
  d.tol <- .Machine$double.eps^.3 ## group `similar sized' penalties, to save work

  r.tol <- .Machine$double.eps^.75 ## This is a bit delicate -- too large and penalty range space can be supressed.

  oo <- .C(C_get_stableS,S=as.double(matrix(0,q,q)),Qs=as.double(matrix(0,q,q)),sp=as.double(exp(lsp)),
                  rS=as.double(unlist(rS)), rSncol = as.integer(rSncol), q = as.integer(q),
                  M = as.integer(M), deriv=as.integer(deriv), det = as.double(0), 
                  det1 = as.double(rep(0,M)),det2 = as.double(matrix(0,M,M)), 
                  d.tol = as.double(d.tol),
                  r.tol = as.double(r.tol),
                  fixed_penalty = as.integer(fixed.penalty))
  S <- matrix(oo$S,q,q)  
  S <- (S+t(S))*.5
  p <- abs(diag(S))^.5            ## by Choleski, p can not be zero if S +ve def
  p[p==0] <- 1                    ## but it's possible to make a mistake!!
  ##E <-  t(t(chol(t(t(S/p)/p)))*p) 
  St <- t(t(S/p)/p)
  St <- (St + t(St))*.5 ## force exact symmetry -- avoids very rare mroot fails 
  E <- t(mroot(St,rank=q)*p) ## the square root S, with column separation
  Qs <- matrix(oo$Qs,q,q)         ## the reparameterization matrix t(Qs)%*%S%*%Qs -> S
  k0 <- 1
  for (i in 1:length(rS)) { ## unpack the rS in the new space
    crs <- ncol(rS[[i]]);
    k1 <- k0 + crs * q - 1 
    rS[[i]] <- matrix(oo$rS[k0:k1],q,crs)
    k0 <- k1 + 1
  }
  ## now get determinant + derivatives, if required...
  if (deriv > 0) det1 <- oo$det1 else det1 <- NULL
  if (deriv > 1) det2 <- matrix(oo$det2,M,M) else det2 <- NULL  
  list(S=S,E=E,Qs=Qs,rS=rS,det=oo$det,det1=det1,det2=det2,fixed.penalty = fixed.penalty)
} ## gam.reparam



gam.fit3 <- function (x, y, sp, Eb,UrS=list(),
            weights = rep(1, nobs), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, nobs),U1=diag(ncol(x)), Mp=-1, family = gaussian(), 
            control = gam.control(), intercept = TRUE,deriv=2,
            gamma=1,scale=1,printWarn=TRUE,scoreType="REML",null.coef=rep(0,ncol(x)),
            pearson.extra=0,dev.extra=0,n.true=-1,Sl=NULL,nei=NULL,...) {
## Inputs:
## * x model matrix
## * y response
## * sp log smoothing parameters
## * Eb square root of nicely balanced total penalty matrix used for rank detection
## * UrS list of penalty square roots in range space of overall penalty. UrS[[i]]%*%t(UrS[[i]]) 
##   is penalty. See 'estimate.gam' for more.
## * weights prior weights (reciprocal variance scale)
## * start initial values for parameters. ignored if etastart or mustart present (although passed on).
## * etastart initial values for eta
## * mustart initial values for mu. discarded if etastart present.
## * control - control list.
## * intercept - indicates whether model has one.
## * deriv - order 0,1 or 2 derivatives are to be returned (lower is cheaper!)
## * gamma - multiplier for effective degrees of freedom in GCV/UBRE.
## * scale parameter. Negative signals to estimate.
## * printWarn print or supress?
## * scoreType - type of smoothness selection to use.
## * null.coef - coefficients for a null model, in order to be able to check for immediate 
##   divergence.
## * pearson.extra is an extra component to add to the pearson statistic in the P-REML/P-ML 
##   case, only.
## * dev.extra is an extra component to add to the deviance in the REML and ML cases only.
## * n.true is to be used in place of the length(y) in ML/REML calculations,
##   and the scale.est only.
## 
## Version with new reparameterization and truncation strategy. Allows iterative weights 
## to be negative. Basically the workhorse routine for Wood (2011) JRSSB.
## A much modified version of glm.fit. Purpose is to estimate regression coefficients 
## and compute a smoothness selection score along with its derivatives.
##
    if (control$trace) { t0 <- proc.time();tc <- 0} 
    warn <- list()
    if (inherits(family,"extended.family")) { ## then actually gam.fit4/5 is needed
      if (inherits(family,"general.family")) {
        return(gam.fit5(x,y,sp,Sl=Sl,weights=weights,offset=offset,deriv=deriv,
                        family=family,scoreType=scoreType,control=control,Mp=Mp,start=start,gamma=gamma,nei=nei))
      } else
      return(gam.fit4(x, y, sp, Eb,UrS=UrS,
            weights = weights, start = start, etastart = etastart, 
            mustart = mustart, offset = offset,U1=U1, Mp=Mp, family = family, 
            control = control, deriv=deriv,gamma=gamma,
            scale=scale,scoreType=scoreType,null.coef=null.coef,nei=nei,...))
    }

    if (family$link==family$canonical) fisher <- TRUE else fisher=FALSE 
    ## ... if canonical Newton = Fisher, but Fisher cheaper!
    if (scale>0) scale.known <- TRUE else scale.known <- FALSE
    if (!scale.known&&scoreType%in%c("REML","ML","EFS")) { ## the final element of sp is actually log(scale)
      nsp <- length(sp)
      scale <- exp(sp[nsp])
      sp <- sp[-nsp]
    }
    if (!deriv%in%c(0,1,2)) stop("unsupported order of differentiation requested of gam.fit3")
    x <- as.matrix(x)  
    nSp <- length(sp)  
    if (nSp==0) deriv.sp <- 0 else deriv.sp <- deriv 

    rank.tol <- .Machine$double.eps*100 ## tolerance to use for rank deficiency

    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)


    q <- ncol(x)

    if (length(UrS)) { ## find a stable reparameterization...
      
      grderiv <- if (scoreType=="EFS") 1 else deriv*as.numeric(scoreType%in%c("REML","ML","P-REML","P-ML")) 
      rp <- gam.reparam(UrS,sp,grderiv) ## note also detects fixed penalty if present
 ## Following is for debugging only...
 #     deriv.check <- FALSE
 #     if (deriv.check&&grderiv) {
 #       eps <- 1e-4
 #       fd.grad <- rp$det1
 #       for (i in 1:length(sp)) {
 #         spp <- sp; spp[i] <- spp[i] + eps/2
 #         rp1 <- gam.reparam(UrS,spp,grderiv)
 #         spp[i] <- spp[i] - eps #         rp0 <- gam.reparam(UrS,spp,grderiv)
 #         fd.grad[i] <- (rp1$det-rp0$det)/eps
 #       }
 #       print(fd.grad)
 #       print(rp$det1) 
 #     }

      T <- diag(q)
      T[1:ncol(rp$Qs),1:ncol(rp$Qs)] <- rp$Qs
      T <- U1%*%T ## new params b'=T'b old params
    
      null.coef <- t(T)%*%null.coef  
 
      if (!is.null(start)) start <- t(T)%*%start
    
      ## form x%*%T in parallel 
      x <- .Call(C_mgcv_pmmult2,x,T,0,0,control$nthreads)
      ## x <- x%*%T   ## model matrix 0(nq^2)
   
      rS <- list()
      for (i in 1:length(UrS)) {
        rS[[i]] <- rbind(rp$rS[[i]],matrix(0,Mp,ncol(rp$rS[[i]])))
      } ## square roots of penalty matrices in current parameterization
      Eb <- Eb%*%T ## balanced penalty matrix
      rows.E <- q-Mp
      Sr <- cbind(rp$E,matrix(0,nrow(rp$E),Mp))
      St <- rbind(cbind(rp$S,matrix(0,nrow(rp$S),Mp)),matrix(0,Mp,q))
    } else {
      grderiv <- 0
      T <- diag(q); 
      St <- matrix(0,q,q) 
      rSncol <- sp <- rows.E <- Eb <- Sr <- 0   
      rS <- list(0)
      rp <- list(det=0,det1 = rep(0,0),det2 = rep(0,0),fixed.penalty=FALSE)
    }
    iter <- 0;coef <- rep(0,ncol(x))

    if (scoreType=="EFS") {
      scoreType <- "REML" ## basically optimizing REML
      deriv <- 0 ## only derivatives of log|S|_+ required (see above)
    }

    conv <- FALSE
    n <- nobs <- NROW(y) ## n is just to keep codetools happy
    if (n.true <= 0) n.true <- nobs ## n.true is used in criteria in place of nobs
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
    } else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }

    ## Added code
    if (family$family=="gaussian"&&family$link=="identity") strictly.additive <- TRUE else
      strictly.additive <- FALSE

    ## end of added code

    D1 <- D2 <- P <- P1 <- P2 <- trA <- trA1 <- trA2 <- 
        GCV<- GCV1<- GCV2<- GACV<- GACV1<- GACV2<- UBRE <-
        UBRE1<- UBRE2<- REML<- REML1<- REML2 <- NCV <- NCV1 <- NULL

    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("Invalid fitted means in empty model")
        dev <- sum(dev.resids(y, mu, weights))
        w <- (weights * mu.eta(eta)^2)/variance(mu)   ### BUG: incorrect for Newton
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
        V <- variance(mu)
        alpha <- dev
        trA2 <- trA1 <- trA <- 0
        if (deriv) GCV2 <- GCV1<- UBRE2 <- UBRE1<-trA1 <- rep(0,nSp)
        GCV <- nobs*alpha/(nobs-gamma*trA)^2
        UBRE <- alpha/nobs - scale + 2*gamma/n*trA
        scale.est <- alpha / (nobs - trA)
    } ### end if (EMPTY)
    else {
       
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop(gettextf("Length of start should equal %d and correspond to initial coefs for %s", 
                     nvars, deparse(xnames)))
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1) 
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
       
        mu <- linkinv(eta)
      
        boundary <- conv <- FALSE
        rV=matrix(0,ncol(x),ncol(x))   
       
        ## need an initial `null deviance' to test for initial divergence... 
        ## Note: can be better not to shrink towards start on
        ## immediate failure, in case start is on edge of feasible space...
        ## if (!is.null(start)) null.coef <- start
        coefold <- null.coef
        etaold <- null.eta <- as.numeric(x%*%null.coef + as.numeric(offset))
        old.pdev <- sum(dev.resids(y, linkinv(null.eta), weights)) + t(null.coef)%*%St%*%null.coef 
        ## ... if the deviance exceeds this then there is an immediate problem
        ii <- 0
        while (!(validmu(mu) && valideta(eta))) { ## shrink towards null.coef if immediately invalid
          ii <- ii + 1
          if (ii>20) stop("Can't find valid starting values: please specify some")
          if (!is.null(start)) start <- start * .9 + coefold * .1
          eta <- .9 * eta + .1 * etaold  
          mu <- linkinv(eta)
        }
        zg <- rep(0,max(dim(x)))
        for (iter in 1:control$maxit) { ## start of main fitting iteration
            good <- weights > 0
            var.val <- variance(mu)
            varmu <- var.val[good]
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
                warn[[length(warn)+1]] <- gettextf("gam.fit3 no observations informative at iteration %d", iter)
                break
            }
            mevg<-mu.eta.val[good];mug<-mu[good];yg<-y[good]
            weg<-weights[good];var.mug<-var.val[good]
            if (fisher) { ## Conventional Fisher scoring
              z <- (eta - offset)[good] + (yg - mug)/mevg
              w <- (weg * mevg^2)/var.mug
            } else { ## full Newton
              c = yg - mug
              alpha <- 1 + c*(family$dvar(mug)/var.mug + family$d2link(mug)*mevg)
              alpha[alpha==0] <- .Machine$double.eps
              z <- (eta - offset)[good] + (yg-mug)/(mevg*alpha) 
              ## ... offset subtracted as eta = X%*%beta + offset
              w <- weg*alpha*mevg^2/var.mug
            }

            ##if (sum(good)<ncol(x)) stop("Not enough informative observations.")

            if (control$trace) t1 <- proc.time()

            ng <- sum(good);zg[1:ng] <- z ## ensure y dim large enough for beta in all cases (including p>n)
            oo <- .C(C_pls_fit1,y=as.double(zg),X=as.double(x[good,]),w=as.double(w),wy=as.double(w*z),
                     E=as.double(Sr),Es=as.double(Eb),n=as.integer(ng),
                     q=as.integer(ncol(x)),rE=as.integer(rows.E),eta=as.double(z),
                     penalty=as.double(1),rank.tol=as.double(rank.tol),nt=as.integer(control$nthreads),
                     use.wy=as.integer(0))
            if (control$trace) tc <- tc + sum((proc.time()-t1)[c(1,4)])

            if (!fisher&&oo$n<0) { ## likelihood indefinite - switch to Fisher for this step
              z <- (eta - offset)[good] + (yg - mug)/mevg
              w <- (weg * mevg^2)/var.mug
	      ng <- sum(good);zg[1:ng] <- z ## ensure y dim large enough for beta in all cases
              if (control$trace) t1 <- proc.time()
              oo <- .C(C_pls_fit1,y=as.double(zg),X=as.double(x[good,]),w=as.double(w),wy=as.double(w*z),
                       E=as.double(Sr),Es=as.double(Eb),n=as.integer(ng),
                       q=as.integer(ncol(x)),rE=as.integer(rows.E),eta=as.double(z),
                       penalty=as.double(1),rank.tol=as.double(rank.tol),nt=as.integer(control$nthreads),
                       use.wy=as.integer(0))
              if (control$trace) tc <- tc + sum((proc.time()-t1)[c(1,4)])
            }

            start <- oo$y[1:ncol(x)];
            penalty <- oo$penalty
            eta <- drop(x%*%start)

            if (any(!is.finite(start))) {
                conv <- FALSE
                warn[[length(warn)+1]] <-gettextf("gam.fit3 non-finite coefficients at iteration %d", 
                  iter)
                break
            }        
     
           mu <- linkinv(eta <- eta + offset)
           dev <- sum(dev.resids(y, mu, weights))
          
           if (control$trace) 
             message(gettextf("Deviance = %s Iterations - %d", dev, iter, domain = "R-mgcv"))
            boundary <- FALSE
            
            if (!is.finite(dev)) {
                if (is.null(coefold)) {
                  if (is.null(null.coef)) 
                  stop("no valid set of coefficients has been found:please supply starting values", 
                    call. = FALSE)
                  ## Try to find feasible coefficients from the null.coef and null.eta
                  coefold <- null.coef
                  etaold <- null.eta
                }
                warn[[length(warn)+1]] <- "gam.fit3 step size truncated due to divergence"
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- (eta + etaold)/2               
                  mu <- linkinv(eta)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE 
                penalty <- t(start)%*%St%*%start
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                warn[[length(warn)+1]] <- "gam.fit3 step size truncated: out of bounds"
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- (eta + etaold)/2 
                  mu <- linkinv(eta)
                }
                boundary <- TRUE
                penalty <- t(start)%*%St%*%start
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }

            pdev <- dev + penalty  ## the penalized deviance 

            if (control$trace) 
                message(gettextf("penalized deviance = %s", pdev, domain = "R-mgcv"))
               

            div.thresh <- 10*(.1+abs(old.pdev))*.Machine$double.eps^.5 
            ## ... threshold for judging divergence --- too tight and near
            ## perfect convergence can cause a failure here

            if (pdev-old.pdev>div.thresh) { ## solution diverging
             ii <- 1 ## step halving counter
             if (iter==1) { ## immediate divergence, need to shrink towards zero 
               etaold <- null.eta; coefold <- null.coef
             }
             while (pdev -old.pdev > div.thresh)  
             { ## step halve until pdev <= old.pdev
                if (ii > 100) 
                   stop("inner loop 3; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2 
                eta <- (eta + etaold)/2               
                mu <- linkinv(eta)
                  dev <- sum(dev.resids(y, mu, weights))
                  pdev <- dev + t(start)%*%St%*%start ## the penalized deviance
                if (control$trace) 
                  message(gettextf("Step halved: new penalized deviance = %g", pdev, "\n"))
              }
            } 
            
            if (strictly.additive) { conv <- TRUE;coef <- start;break;}

            if (abs(pdev - old.pdev) < control$epsilon*(abs(scale)+abs(pdev))) {
               ## Need to check coefs converged adequately, to ensure implicit differentiation
               ## ok. Testing coefs unchanged is problematic under rank deficiency (not guaranteed to
               ## drop same parameter every iteration!)       
               grad <- 2 * t(x[good,])%*%(w*((x%*%start)[good]-z))+ 2*St%*%start 
               #if (max(abs(grad)) > control$epsilon*max(abs(start+coefold))/2) {
	       if (max(abs(grad)) > control$epsilon*(abs(pdev)+abs(scale))) {
                  old.pdev <- pdev
                  coef <- coefold <- start
                  etaold <- eta 
                } else {
                  conv <- TRUE
                  coef <- start
                  etaold <- eta
                  break 
                }
            }
            else {  old.pdev <- pdev
                coef <- coefold <- start
                etaold <- eta 
            }
        } ### end main loop 
       
        wdr <- dev.resids(y, mu, weights)
        dev <- sum(wdr) 
        #wdr <- sign(y-mu)*sqrt(pmax(wdr,0)) ## used below in scale estimation 
  
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
         var.val <- variance(mu)
         varmu <- var.val[good]
         if (any(is.na(varmu))) stop("NAs in V(mu)")
         if (any(varmu == 0)) stop("0s in V(mu)")
         mu.eta.val <- mu.eta(eta)
         if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
   

         good <- (weights > 0) & (mu.eta.val != 0)
         mevg <- mu.eta.val[good];mug <- mu[good];yg <- y[good]
         weg <- weights[good];etag <- eta[good]
         var.mug<-var.val[good]

         if (fisher) { ## Conventional Fisher scoring
              z <- (eta - offset)[good] + (yg - mug)/mevg
              w <- (weg * mevg^2)/var.mug
              alpha <- wf <- 0 ## Don't need Fisher weights separately
         } else { ## full Newton
              c <- yg - mug
              alpha <- 1 + c*(family$dvar(mug)/var.mug + family$d2link(mug)*mevg)
              ### can't just drop obs when alpha==0, as they are informative, but
              ### happily using an `effective zero' is stable here, and there is 
              ### a natural effective zero, since E(alpha) = 1.
              alpha[alpha==0] <- .Machine$double.eps 
              z <- (eta - offset)[good] + (yg-mug)/(mevg*alpha) 
              ## ... offset subtracted as eta = X%*%beta + offset
              wf <- weg*mevg^2/var.mug ## Fisher weights for EDF calculation
              w <- wf * alpha   ## Full Newton weights
         }
        
         g1 <- 1/mevg
         g2 <- family$d2link(mug)
         g3 <- family$d3link(mug)

         V <- family$variance(mug)
         V1 <- family$dvar(mug)
         V2 <- family$d2var(mug)      
        
         if (fisher) {
           g4 <- V3 <- 0
         } else {
           g4 <- family$d4link(mug)
           V3 <- family$d3var(mug)
         }

         if (TRUE) { ### TEST CODE for derivative ratio based versions of code... 
           g2 <- g2/g1;g3 <- g3/g1;g4 <- g4/g1
           V1 <- V1/V;V2 <- V2/V;V3 <- V3/V
         }

         P1 <- D1 <- array(0,nSp);P2 <- D2 <- matrix(0,nSp,nSp) # for derivs of deviance/ Pearson
         trA1 <- array(0,nSp);trA2 <- matrix(0,nSp,nSp) # for derivs of tr(A)
         rV=matrix(0,ncol(x),ncol(x));
         dum <- 1
         if (control$trace) cat("calling gdi...")

       REML <- 0 ## signals GCV/AIC used
       if (scoreType%in%c("REML","P-REML","NCV")) {REML <- 1;remlInd <- 1} else 
       if (scoreType%in%c("ML","P-ML")) {REML <- -1;remlInd <- 0} 

       if (REML==0) rSncol <- unlist(lapply(rS,ncol)) else rSncol <- unlist(lapply(UrS,ncol))
       if (control$trace) t1 <- proc.time()
       oo <- .C(C_gdi1,X=as.double(x[good,]),E=as.double(Sr),Eb = as.double(Eb), 
                rS = as.double(unlist(rS)),U1=as.double(U1),sp=as.double(exp(sp)),
                z=as.double(z),w=as.double(w),wf=as.double(wf),alpha=as.double(alpha),
                mu=as.double(mug),eta=as.double(etag),y=as.double(yg),
                p.weights=as.double(weg),g1=as.double(g1),g2=as.double(g2),
                g3=as.double(g3),g4=as.double(g4),V0=as.double(V),V1=as.double(V1),
                V2=as.double(V2),V3=as.double(V3),beta=as.double(coef),b1=as.double(rep(0,nSp*ncol(x))),
                w1=as.double(rep(0,nSp*length(z))),
                D1=as.double(D1),D2=as.double(D2),P=as.double(dum),P1=as.double(P1),P2=as.double(P2),
                trA=as.double(dum),trA1=as.double(trA1),trA2=as.double(trA2),
                rV=as.double(rV),rank.tol=as.double(rank.tol),
                conv.tol=as.double(control$epsilon),rank.est=as.integer(1),n=as.integer(length(z)),
                p=as.integer(ncol(x)),M=as.integer(nSp),Mp=as.integer(Mp),Enrow = as.integer(rows.E),
                rSncol=as.integer(rSncol),deriv=as.integer(deriv.sp),
                REML = as.integer(REML),fisher=as.integer(fisher),
                fixed.penalty = as.integer(rp$fixed.penalty),nthreads=as.integer(control$nthreads),
                dVkk=as.double(rep(0,nSp*nSp)))      
         if (control$trace) { 
           tg <- sum((proc.time()-t1)[c(1,4)])
           cat("done!\n")
         }
 
         ## get dbeta/drho, original parameterization restored on return
         db.drho <- if (deriv) matrix(oo$b1,ncol(x),nSp) else NULL
         dw.drho <- if (deriv) matrix(oo$w1,length(z),nSp) else NULL ## NOTE: only deriv of Newton weights if REML=1 in gdi1 call

         rV <- matrix(oo$rV,ncol(x),ncol(x)) ## rV%*%t(rV)*scale gives covariance matrix 
         
         Kmat <- matrix(0,nrow(x),ncol(x)) 
         Kmat[good,] <- oo$X                    ## rV%*%t(K)%*%(sqrt(wf)*X) = F; diag(F) is edf array 

         coef <- oo$beta;
         eta <- drop(x%*%coef + offset)
         mu <- linkinv(eta)
         if (!(validmu(mu)&&valideta(eta))) {
           ## if iteration terminated with step halving then it can be that
           ## gdi1 returns an invalid coef, because likelihood is actually
           ## pushing coefs to invalid region. Probably not much hope in 
           ## this case, but it is worth at least returning feasible values,
           ## even though this is not quite consistent with derivs.
           coef <- start
           eta <- etaold
           mu <- linkinv(eta)
         }
         trA <- oo$trA;
                  
         if (control$scale.est%in%c("pearson","fletcher","Pearson","Fletcher")) {
            pearson <- sum(weights*(y-mu)^2/family$variance(mu))
            scale.est <- (pearson+dev.extra)/(n.true-trA)
            if (control$scale.est%in%c("fletcher","Fletcher")) { ## Apply Fletcher (2012) correction
              ## note limited to 10 times Pearson...
              #s.bar = max(-.9,mean(family$dvar(mu)*(y-mu)*sqrt(weights)/family$variance(mu)))
	      s.bar = max(-.9,mean(family$dvar(mu)*(y-mu)/family$variance(mu)))
              if (is.finite(s.bar)) scale.est <- scale.est/(1+s.bar)
            }
         } else { ## use the deviance estimator
           scale.est <- (dev+dev.extra)/(n.true-trA)
         }

        reml.scale <- NA
	Vg=NULL ## empirical cov matrix for grad (see NCV)

        if (scoreType%in%c("REML","ML")) { ## use Laplace (RE)ML
          
          ls <- family$ls(y,weights,n,scale)*n.true/nobs ## saturated likelihood and derivatives
          Dp <- dev + oo$conv.tol + dev.extra
          REML <- (Dp/(2*scale) - ls[1])/gamma + oo$rank.tol/2 - rp$det/2 -
	          remlInd*(Mp/2*(log(2*pi*scale)-log(gamma)))
          attr(REML,"Dp") <- Dp/(2*scale)
          if (deriv) {
            REML1 <- oo$D1/(2*scale*gamma) + oo$trA1/2 - rp$det1/2 
            if (deriv==2) REML2 <- (matrix(oo$D2,nSp,nSp)/(scale*gamma) + matrix(oo$trA2,nSp,nSp) - rp$det2)/2
            if (sum(!is.finite(REML2))) {
               stop("Non finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")
            }
          }
          if (!scale.known&&deriv) { ## need derivatives wrt log scale, too 
            ##ls <- family$ls(y,weights,n,scale) ## saturated likelihood and derivatives
            dlr.dlphi <- (-Dp/(2 *scale) - ls[2]*scale)/gamma - Mp/2*remlInd
            d2lr.d2lphi <- (Dp/(2*scale) - ls[3]*scale^2 - ls[2]*scale)/gamma
            d2lr.dspphi <- -oo$D1/(2*scale*gamma)
            REML1 <- c(REML1,dlr.dlphi)
            if (deriv==2) {
              REML2 <- rbind(REML2,as.numeric(d2lr.dspphi))
              REML2 <- cbind(REML2,c(as.numeric(d2lr.dspphi),d2lr.d2lphi))
            }
          }
          reml.scale <- scale
        } else if (scoreType%in%c("P-REML","P-ML")) { ## scale unknown use Pearson-Laplace REML
          reml.scale <- phi <- (oo$P*(nobs-Mp)+pearson.extra)/(n.true-Mp) ## REMLish scale estimate
          ## correct derivatives, if needed...
          oo$P1 <- oo$P1*(nobs-Mp)/(n.true-Mp)
          oo$P2 <- oo$P2*(nobs-Mp)/(n.true-Mp)

          ls <- family$ls(y,weights,n,phi)*n.true/nobs ## saturated likelihood and derivatives
        
          Dp <- dev + oo$conv.tol + dev.extra
         
          K <- oo$rank.tol/2 - rp$det/2
                 
          REML <- (Dp/(2*phi) - ls[1]) + K - Mp/2*(log(2*pi*phi))*remlInd
          attr(REML,"Dp") <- Dp/(2*phi)
          if (deriv) {
            phi1 <- oo$P1; Dp1 <- oo$D1; K1 <- oo$trA1/2 - rp$det1/2;
            REML1 <- Dp1/(2*phi) - phi1*(Dp/(2*phi^2)+Mp/(2*phi)*remlInd + ls[2]) + K1
            if (deriv==2) {
                   phi2 <- matrix(oo$P2,nSp,nSp);Dp2 <- matrix(oo$D2,nSp,nSp)
                   K2 <- matrix(oo$trA2,nSp,nSp)/2 - rp$det2/2   
                   REML2 <- 
                   Dp2/(2*phi) - (outer(Dp1,phi1)+outer(phi1,Dp1))/(2*phi^2) +
                   (Dp/phi^3 - ls[3] + Mp/(2*phi^2)*remlInd)*outer(phi1,phi1) -
                   (Dp/(2*phi^2)+ls[2]+Mp/(2*phi)*remlInd)*phi2 + K2
            }
          }
 
        } else if (scoreType %in% "NCV") { ## neighbourhood cross validation
	   ## requires that neighbours are supplied in nei (if NULL) each point is its own
	   ## neighbour recovering LOOCV. nei$k[nei$m[i-1]+1):nei$m[i]] are the indices of
	   ## neighbours of point i, where nei$m[0]=0 by convention.
	 
	   ww <- w1 <- rep(0,nobs)
	   ww[good] <- weg*(yg - mug)*mevg/var.mug
	   w1[good] <- w
	   eta.cv <- rep(0.0,length(nei$i))
	   deta.cv <- if (deriv) matrix(0.0,length(nei$i),length(rS)) else matrix(0.0,1,length(rS))
           if (nei$jackknife>2) { ## need coef changes for each NCV drop fold.
             dd <- matrix(0,ncol(x),length(nei$m))
	     if (deriv>0) stop("jackknife and derivatives requested together")
	     deriv1 <- -1 ## signal that coef changes to be returned
           } else { ## dd unused
             dd <- matrix(1.0,1,1);
	     deriv1 <- deriv
	   }  
	   R <- try(chol(crossprod(x,w1*x)+St),silent=TRUE)
	   if (inherits(R,"try-error")) { ## use CG approach...
	     Hi <- tcrossprod(rV) ## inverse of penalized Expected Hessian - inverse actual Hessian probably better
             cg.iter <- .Call(C_ncv,x,Hi,ww,w1,db.drho,dw.drho,rS,nei$i-1,nei$mi,nei$m,nei$k-1,coef,exp(sp),eta.cv, deta.cv, dd, deriv1);
	     warn[[length(warn)+1]] <- "NCV positive definite update check not possible"
           } else { ## use Cholesky update approach
	     pdef.fails <- .Call(C_Rncv,x,R,ww,w1,db.drho,dw.drho,rS,nei$i-1,nei$mi,nei$m,nei$k-1,coef,exp(sp),eta.cv, deta.cv,
	                         dd, deriv1,.Machine$double.eps,control$ncv.threads);
	     if (pdef.fails) warn[[length(warn)+1]] <- "some NCV updates not positive definite"
	   }
	  
	   if (family$qapprox) {
             NCV <- sum(wdr[nei$i]) + gamma*sum(-2*ww[nei$i]*(eta.cv-eta[nei$i]) + w1[nei$i]*(eta.cv-eta[nei$i])^2)
             if (deriv) {
	       deta <- x%*%db.drho
	       alpha1 <- if (fisher) 0 else (-(V1+g2) + (y-mu)*(V2-V1^2+g3-g2^2))/alpha
	       w3 <- w1/g1*(alpha1 - V1 - 2 * g2)
               ncv1 <- -2*ww[nei$i]*((1-gamma)*deta[nei$i,] + gamma*deta.cv) + 2*gamma*w1[nei$i]*(deta.cv*(eta.cv-eta[nei$i])) +
	                gamma*w3[nei$i]* deta[nei$i,]*(eta.cv-eta[nei$i])^2
             }
           } else { ## exact version
             if (TRUE) {
	       ## version that doesn't just drop neighbourhood, but tries equivalent (gamma==2) perturbation beyond dropping
               eta.cv <- gamma*(eta.cv) - (gamma-1)*eta[nei$i]
	       if (deriv && gamma!=1) deta.cv <- gamma*(deta.cv) - (gamma-1)*(x%*%db.drho)[nei$i,,drop=FALSE]
	       gamma <- 1
             }
     	     mu.cv <- linkinv(eta.cv)
	     NCV <- gamma*sum(dev.resids(y[nei$i],mu.cv,weights[nei$i])) - (gamma-1)*sum(wdr[nei$i]) ## the NCV score - simply LOOCV if nei(i) = i for all i
	    
             if (deriv) {
               dev1 <- if (gamma==1) 0 else -2*(ww*(x%*%db.drho))[nei$i,,drop=FALSE]
               var.mug <- variance(mu.cv)
               mevg <- mu.eta(eta.cv)
	       mug <- mu.cv
	       ww1 <- weights[nei$i]*(y[nei$i]-mug)*mevg/var.mug
	       ww1[!is.finite(ww1)] <- 0
	       ncv1 <- -2*ww1*deta.cv*gamma - (gamma-1)*dev1
             } ## if deriv
	   }

           if (deriv) {
             NCV1 <- colSums(ncv1) ## grad
	     #attr(NCV1,"gjk") <- gjk ## jackknife grad estimate - incomplete - need qapprox version as well
	     Vg <- crossprod(ncv1) ## empirical cov matrix of grad
           }

	   if (nei$jackknife>2) {
             nk <- c(nei$m[1],diff(nei$m)) ## dropped fold sizes
             jkw <- ((nobs-nk)/(nobs))^.4/sqrt(nk) ## jackknife weights
	     dd <-jkw*t(dd)%*%t(T)
	     #Vj <- crossprod(dd) ## jackknife cov matrix
             #attr(Vj,"dd") <- dd
             #attr(NCV,"Vj") <- Vj
	     attr(NCV,"dd") <- dd
	   }  
	   attr(NCV,"eta.cv") <- eta.cv
	   if (deriv) attr(NCV,"deta.cv") <- deta.cv
	} else { ## GCV/GACV etc ....

           P <- oo$P
           
           delta <- nobs - gamma * trA
           delta.2 <- delta*delta           
  
           GCV <- nobs*dev/delta.2
           GACV <- dev/nobs + P * 2*gamma*trA/(delta * nobs) 

           UBRE <- dev/nobs - 2*delta*scale/nobs + scale
        
           if (deriv) {
             trA1 <- oo$trA1
           
             D1 <- oo$D1
             P1 <- oo$P1
          
             if (sum(!is.finite(D1))||sum(!is.finite(P1))||sum(!is.finite(trA1))) { 
                 stop(
               "Non-finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")
             }
         
             delta.3 <- delta*delta.2
  
             GCV1 <- nobs*D1/delta.2 + 2*nobs*dev*trA1*gamma/delta.3
             GACV1 <- D1/nobs + 2*P/delta.2 * trA1 + 2*gamma*trA*P1/(delta*nobs)

             UBRE1 <- D1/nobs + gamma * trA1 *2*scale/nobs
             if (deriv==2) {
               trA2 <- matrix(oo$trA2,nSp,nSp) 
               D2 <- matrix(oo$D2,nSp,nSp)
               P2 <- matrix(oo$P2,nSp,nSp)
              
               if (sum(!is.finite(D2))||sum(!is.finite(P2))||sum(!is.finite(trA2))) { 
                 stop(
                 "Non-finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")
               }
             
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
             } ## end if (deriv==2)
           } ## end if (deriv)
        } ## end !REML
        # end of inserted code
        if (!conv) 
              warn[[length(warn)+1]] <- "gam.fit3 algorithm did not converge"
        if (boundary) 
              warn[[length(warn)+1]] <- "gam.fit3 algorithm stopped at boundary value"
        eps <- 10 * .Machine$double.eps
        if (printWarn&&family$family[1] == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warn[[length(warn)+1]] <- "gam.fit3 fitted probabilities numerically 0 or 1 occurred"
        }
        if (printWarn&&family$family[1] == "poisson") {
            if (any(mu < eps)) 
                warn[[length(warn)+1]] <- "gam.fit3 fitted rates numerically 0 occurred"
        }
 
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
          
        ## undo reparameterization....
        coef <- as.numeric(T %*% coef)
        rV <- T %*% rV
        names(coef) <- xnames 
    } ### end if (!EMPTY)
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    ww <- wt <- rep.int(0, nobs)
    if (fisher) { wt[good] <- w; ww <- wt} else { 
      wt[good] <- wf  ## note that Fisher weights are returned
      ww[good] <- w
    }
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (deriv) {
      db.drho <- T%*%db.drho
      if (nrow(dw.drho)!=nrow(x)) {
        w1 <- dw.drho
        dw.drho <- matrix(0,nrow(x),ncol(w1))
        dw.drho[good,] <- w1
      }
    }
    
    sumw <- sum(weights)
    wtdmu <- if (intercept) 
        sum(weights * y)/sumw
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)

    ## use exact MLE scale param for aic in gaussian case, otherwise scale.est (unless known)
    
    dev1 <- if (scale.known) scale*sumw else if (family$family=="gaussian") dev else
            if (is.na(reml.scale)) scale.est*sumw else reml.scale*sumw  

    aic.model <- aic(y, n, mu, weights, dev1) # note: incomplete 2*edf needs to be added
    if (control$trace) {
      t1 <- proc.time()
      at <- sum((t1-t0)[c(1,4)])
      cat("Proportion time in C: ",(tc+tg)/at," ls:",tc/at," gdi:",tg/at,"\n")
    } 
   
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
         family = family, linear.predictors = eta, deviance = dev, 
        null.deviance = nulldev, iter = iter, weights = wt, working.weights=ww,prior.weights = weights, z=z,
        df.null = nulldf, y = y, converged = conv,##pearson.warning = pearson.warning,
        boundary = boundary,D1=D1,D2=D2,P=P,P1=P1,P2=P2,trA=trA,trA1=trA1,trA2=trA2,NCV=NCV,NCV1=NCV1,
        GCV=GCV,GCV1=GCV1,GCV2=GCV2,GACV=GACV,GACV1=GACV1,GACV2=GACV2,UBRE=UBRE,
        UBRE1=UBRE1,UBRE2=UBRE2,REML=REML,REML1=REML1,REML2=REML2,rV=rV,Vg=Vg,db.drho=db.drho,
        dw.drho=dw.drho,dVkk = matrix(oo$dVkk,nSp,nSp),ldetS1 = if (grderiv) rp$det1 else 0,
        scale.est=scale.est,reml.scale= reml.scale,aic=aic.model,rank=oo$rank.est,K=Kmat,warn=warn)
} ## end gam.fit3

Vb.corr <- function(X,L,lsp0,S,off,dw,w,rho,Vr,nth=0,scale.est=FALSE) {
## compute higher order Vb correction...
## If w is NULL then X should be root Hessian, and 
## dw is treated as if it was 0, otherwise X should be model 
## matrix.
## dw is derivative w.r.t. all the smoothing parameters and family parameters as if these 
## were not linked, but not the scale parameter, of course. Vr includes scale uncertainty,
## if scale estimated...
## nth is the number of initial elements of rho that are not smoothing 
## parameters, scale.est is TRUE if scale estimated by REML and must be
## dropped from s.p.s
  M <- length(off) ## number of penalty terms
  if (scale.est) {
    ## drop scale param from L, rho and Vr...
    rho <- rho[-length(rho)]
    if (!is.null(L)) L <- L[-nrow(L),-ncol(L),drop=FALSE]
    Vr <- Vr[-nrow(Vr),-ncol(Vr),drop=FALSE]
  }
 
  if (is.null(lsp0)) lsp0 <- if (is.null(L)) rho*0 else rep(0,nrow(L))
  ## note that last element of lsp0 can be a scale parameter...
  lambda <- if (is.null(L)) exp(rho+lsp0[1:length(rho)]) else exp(L%*%rho + lsp0[1:nrow(L)])
  
  ## Re-create the Hessian, if is.null(w) then X assumed to be root
  ## unpenalized Hessian...
  H <- if (is.null(w)) crossprod(X) else H <- t(X)%*%(w*X)
  if (M>0) for (i in 1:M) {
      ind <- off[i] + 1:ncol(S[[i]]) - 1
      H[ind,ind] <- H[ind,ind] + lambda[i+nth] * S[[i]]
  }

  R <- try(chol(H),silent=TRUE) ## get its Choleski factor.  
  if (inherits(R,"try-error")) return(0) ## bail out as Hessian insufficiently well conditioned
  
  ## Create dH the derivatives of the hessian w.r.t. (all) the smoothing parameters...
  dH <- list()
  if (length(lambda)>0) for (i in 1:length(lambda)) {
    ## If w==NULL use constant H approx...
    dH[[i]] <- if (is.null(w)) H*0 else t(X)%*%(dw[,i]*X) 
    if (i>nth) { 
      ind <- off[i-nth] + 1:ncol(S[[i-nth]]) - 1
      dH[[i]][ind,ind] <- dH[[i]][ind,ind] + lambda[i]*S[[i-nth]]
    }
  }
  ## If L supplied then dH has to be re-weighted to give
  ## derivatives w.r.t. optimization smoothing params.
  if (!is.null(L)) {
    dH1 <- dH;dH <- list()
    if (length(rho)>0) for (j in 1:length(rho)) { 
      ok <- FALSE ## dH[[j]] not yet created
      if (nrow(L)>0) for (i in 1:nrow(L)) if (L[i,j]!=0.0) { 
        dH[[j]] <- if (ok) dH[[j]] + dH1[[i]]*L[i,j] else dH1[[i]]*L[i,j]
        ok <- TRUE
      }
    } 
    rm(dH1)
  } ## dH now w.r.t. optimization parameters 
  
  if (length(dH)==0) return(0) ## nothing to correct

  ## Get derivatives of Choleski factor w.r.t. the smoothing parameters 
  dR <- list()
  for (i in 1:length(dH)) dR[[i]] <- dchol(dH[[i]],R) 
  rm(dH)
  
  ## need to transform all dR to dR^{-1} = -R^{-1} dR R^{-1}...
  for (i in 1:length(dR)) dR[[i]] <- -t(forwardsolve(t(R),t(backsolve(R,dR[[i]]))))
 
  ## BUT: dR, now upper triangular, and it relates to RR' = Vb not R'R = Vb
  ## in consequence of which Rz is the thing with the right distribution
  ## and not R'z...
  dbg <- FALSE
  if (dbg) { ## debugging code
    n.rep <- 10000;p <- ncol(R)
    r <- rmvn(n.rep,rep(0,M),Vr)
    b <- matrix(0,n.rep,p)
    for (i in 1:n.rep) {
      z <- rnorm(p)
      if (M>0) for (j in 1:M) b[i,] <- b[i,] + dR[[j]]%*%z*(r[i,j]) 
    }
    Vfd <- crossprod(b)/n.rep
  }

  vcorr(dR,Vr,FALSE) ## NOTE: unscaled!!
} ## Vb.corr

gam.fit3.post.proc <- function(X,L,lsp0,S,off,object,gamma) {
## get edf array and covariance matrices after a gam fit. 
## X is original model matrix, L the mapping from working to full sp
  scale <- if (object$scale.estimated) object$scale.est else object$scale
  Vb <- object$rV%*%t(object$rV)*scale ## Bayesian cov.

  # PKt <- object$rV%*%t(object$K)
  PKt <- .Call(C_mgcv_pmmult2,object$rV,object$K,0,1,object$control$nthreads)
  # F <- PKt%*%(sqrt(object$weights)*X)
  F <- .Call(C_mgcv_pmmult2,PKt,sqrt(object$weights)*X,0,0,object$control$nthreads)
  edf <- diag(F) ## effective degrees of freedom
  edf1 <- 2*edf - rowSums(t(F)*F) ## alternative

  ## edf <- rowSums(PKt*t(sqrt(object$weights)*X))
  ## Ve <- PKt%*%t(PKt)*object$scale  ## frequentist cov
  Ve <- F%*%Vb ## not quite as stable as above, but quicker
  hat <- rowSums(object$K*object$K)
  ## get QR factor R of WX - more efficient to do this
  ## in gdi_1 really, but that means making QR of augmented 
  ## a two stage thing, so not clear cut...
  WX <- sqrt(object$weights)*X
  qrx <- pqr(WX,object$control$nthreads)
  R <- pqr.R(qrx);R[,qrx$pivot] <- R
  if (!is.na(object$reml.scale)&&!is.null(object$db.drho)) { ## compute sp uncertainty correction
    hess <- object$outer.info$hess
    edge.correct <- if (is.null(attr(hess,"edge.correct"))) FALSE else TRUE
    K <- if (edge.correct) 2 else 1
    for (k in 1:K) {
      if (k==1) { ## fitted model computations
        db.drho <- object$db.drho
        dw.drho <- object$dw.drho
        lsp <- log(object$sp)
      } else { ## edge corrected model computations
        db.drho <- attr(hess,"db.drho1")
        dw.drho <- attr(hess,"dw.drho1")
        lsp <- attr(hess,"lsp1")
	hess <- attr(hess,"hess1")
      }
      M <- ncol(db.drho)
      ## transform to derivs w.r.t. working, noting that an extra final row of L
      ## may be present, relating to scale parameter (for which db.drho is 0 since it's a scale parameter)  
      if (!is.null(L)) { 
        db.drho <- db.drho%*%L[1:M,,drop=FALSE] 
        M <- ncol(db.drho)
      }
      ## extract cov matrix for log smoothing parameters...
      ev <- eigen(hess,symmetric=TRUE) 
      d <- ev$values;ind <- d <= 0
      d[ind] <- 0;d[!ind] <- 1/sqrt(d[!ind])
      rV <- (d*t(ev$vectors))[,1:M] ## root of cov matrix
      Vc <- crossprod(rV%*%t(db.drho))
      ## set a prior precision on the smoothing parameters, but don't use it to 
      ## fit, only to regularize Cov matrix. exp(4*var^.5) gives approx 
      ## multiplicative range. e.g. var = 5.3 says parameter between .01 and 100 times
      ## estimate. Avoids nonsense at `infinite' smoothing parameters.   
      d <- ev$values; d[ind] <- 0;
      d <- if (k==1) 1/sqrt(d+1/10) else 1/sqrt(d+1e-7)
      Vr <- crossprod(d*t(ev$vectors))
      ## Note that db.drho and dw.drho are derivatives w.r.t. full set of smoothing 
      ## parameters excluding any scale parameter, but Vr includes info for scale parameter
      ## if it has been estimated. 
      nth <- if (is.null(object$family$n.theta)) 0 else object$family$n.theta ## any parameters of family itself
      drop.scale <- object$scale.estimated && !(object$method %in% c("P-REML","P-ML"))
      Vc2 <- scale*Vb.corr(R,L,lsp0,S,off,dw.drho,w=NULL,lsp,Vr,nth,drop.scale)
    
      Vc <- Vb + Vc + Vc2 ## Bayesian cov matrix with sp uncertainty
      ## finite sample size check on edf sanity...
      if (k==1) { ## compute edf2 only with fitted model, not edge corrected
        edf2 <- rowSums(Vc*crossprod(R))/scale
        if (sum(edf2)>sum(edf1)) { 
          edf2 <- edf1
        }
      }
    } ## k loop

    V.sp <- Vr;attr(V.sp,"L") <- L;attr(V.sp,"spind") <- spind <- (nth+1):M
    ## EXPERIMENTAL ## NOTE: no L handling - what about incomplete z/w??
    #P <- Vb %*% X/scale
    if (FALSE) {
      sp <- object$sp[spind]
      beta <- object$coefficients
      w <- object$weights
      m <- length(off)
      M <- K <- matrix(0,length(w),m)
      for (i in 1:m) {
        ii <- 1:ncol(S[[i]])+off[i]-1
        Sb <- S[[i]] %*% beta[ii]*sp[i]/object$scale
	## don't forget Vb is inverse hessian times scale!
        K[,i] <- w * drop(X %*% (Vb[,ii] %*% Sb))/object$scale
        B <- Vb%*%drop(t(w*object$z) %*% X)
        M[,i] <- X%*% (Vb[,ii] %*% (S[[i]]%*%B[ii]))*sp[i]/object$scale^2
      }
      ## Should following really just drop scale term?
      K <- K %*% Vr[spind,spind]
      edf3 <- sum(K*M)
      attr(edf2,"edf3") <- edf3
    }
    ## END EXPERIMENTAL
  } else V.sp <- edf2 <- Vc <- NULL
  ret <- list(Vp=Vb,Ve=Ve,V.sp=V.sp,edf=edf,edf1=edf1,edf2=edf2,hat=hat,F=F,R=R)
  if (is.null(object$Vc)) ret$Vc <- Vc
  ret
} ## gam.fit3.post.proc


score.transect <- function(ii, x, y, sp, Eb,UrS=list(), 
            weights = rep(1, length(y)), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, length(y)),U1,Mp,family = gaussian(), 
            control = gam.control(), intercept = TRUE,deriv=2,
            gamma=1,scale=1,printWarn=TRUE,scoreType="REML",eps=1e-7,null.coef=rep(0,ncol(x)),...) {
## plot a transect through the score for sp[ii]
  np <- 200
  if (scoreType%in%c("REML","P-REML","ML","P-ML")) reml <- TRUE else reml <- FALSE
  
  score <- spi <- seq(-30,30,length=np)
  for (i in 1:np) {

     sp[ii] <- spi[i]
     b<-gam.fit3(x=x, y=y, sp=sp,Eb=Eb,UrS=UrS,
      offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=0,
      control=control,gamma=gamma,scale=scale,
      printWarn=FALSE,mustart=mustart,scoreType=scoreType,null.coef=null.coef,...)

      if (reml) {
        score[i] <- b$REML
      } else if (scoreType=="GACV") {
        score[i] <- b$GACV
      } else if (scoreType=="UBRE"){
        score[i] <- b$UBRE 
      } else { ## default to deviance based GCV
        score[i] <- b$GCV
      }
  }
  par(mfrow=c(2,2),mar=c(4,4,1,1))
  plot(spi,score,xlab="log(sp)",ylab=scoreType,type="l")
  plot(spi[1:(np-1)],score[2:np]-score[1:(np-1)],type="l",ylab="differences")
  plot(spi,score,ylim=c(score[1]-.1,score[1]+.1),type="l")
  plot(spi,score,ylim=c(score[np]-.1,score[np]+.1),type="l")
} ## score.transect

deriv.check <- function(x, y, sp, Eb,UrS=list(), 
            weights = rep(1, length(y)), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, length(y)),U1,Mp,family = gaussian(), 
            control = gam.control(), intercept = TRUE,deriv=2,
            gamma=1,scale=1,printWarn=TRUE,scoreType="REML",eps=1e-7,
            null.coef=rep(0,ncol(x)),Sl=Sl,nei=nei,...)
## FD checking of derivatives: basically a debugging routine
{  
   if (!deriv%in%c(1,2)) stop("deriv should be 1 or 2")
   if (control$epsilon>1e-9) control$epsilon <- 1e-9 
   b<-gam.fit3(x=x, y=y, sp=sp,Eb=Eb,UrS=UrS,
      offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=deriv,
      control=control,gamma=gamma,scale=scale,printWarn=FALSE,
      start=start,etastart=etastart,mustart=mustart,scoreType=scoreType,
      null.coef=null.coef,Sl=Sl,nei=nei,...)

   P0 <- b$P;fd.P1 <- P10 <- b$P1;  if (deriv==2) fd.P2 <- P2 <- b$P2 
   trA0 <- b$trA;fd.gtrA <- gtrA0 <- b$trA1 ; if (deriv==2) fd.htrA <- htrA <- b$trA2 
   dev0 <- b$deviance;fd.D1 <- D10 <- b$D1 ; if (deriv==2) fd.D2 <- D2 <- b$D2 
   fd.db <- b$db.drho*0
   
   if (scoreType%in%c("REML","P-REML","ML","P-ML")) reml <- TRUE else reml <- FALSE
   sname <- if (reml) "REML" else scoreType
   sname1 <- paste(sname,"1",sep=""); sname2 <- paste(sname,"2",sep="")
   if (scoreType=="NCV") reml <- TRUE ## to avoid un-needed stuff

   score0 <- b[[sname]];grad0 <- b[[sname1]]; if (deriv==2) hess <- b[[sname2]] 


   fd.grad <- grad0*0
   if (deriv==2) fd.hess <- hess
   diter <- rep(20,length(sp))
   for (i in 1:length(sp)) {
     sp1 <- sp;sp1[i] <- sp[i]+eps/2
     bf<-gam.fit3(x=x, y=y, sp=sp1,Eb=Eb,UrS=UrS,
      offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=deriv,
      control=control,gamma=gamma,scale=scale,printWarn=FALSE,
      start=start,etastart=etastart,mustart=mustart,scoreType=scoreType,
      null.coef=null.coef,Sl=Sl,nei=nei,...)
      
     sp1 <- sp;sp1[i] <- sp[i]-eps/2
     bb<-gam.fit3(x=x, y=y, sp=sp1, Eb=Eb,UrS=UrS,
      offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=deriv,
      control=control,gamma=gamma,scale=scale,printWarn=FALSE,
      start=start,etastart=etastart,mustart=mustart,scoreType=scoreType,
     null.coef=null.coef,Sl=Sl,nei=nei,...)
     diter[i] <- bf$iter - bb$iter ## check iteration count same 
     if (i<=ncol(fd.db)) fd.db[,i] <- (bf$coefficients - bb$coefficients)/eps


      if (!reml) {
        Pb <- bb$P;Pf <- bf$P 
        P1b <- bb$P1;P1f <- bf$P1
        trAb <- bb$trA;trAf <- bf$trA
        gtrAb <- bb$trA1;gtrAf <- bf$trA1
        devb <- bb$deviance;devf <- bf$deviance
        D1b <- bb$D1;D1f <- bf$D1
      }

      scoreb <- bb[[sname]];scoref <- bf[[sname]];
      if (deriv==2) { gradb <- bb[[sname1]];gradf <- bf[[sname1]]}

      if (!reml) {
        fd.P1[i] <- (Pf-Pb)/eps
        fd.gtrA[i] <- (trAf-trAb)/eps
        fd.D1[i] <- (devf - devb)/eps
      }
      
     
      fd.grad[i] <- (scoref-scoreb)/eps
      if (deriv==2) { 
        fd.hess[,i] <- (gradf-gradb)/eps
        if (!reml) {
          fd.htrA[,i] <- (gtrAf-gtrAb)/eps
          fd.P2[,i] <- (P1f-P1b)/eps
          fd.D2[,i] <- (D1f-D1b)/eps
        } 
       
      }
   }
   
   if (!reml) {
     cat("\n Pearson Statistic... \n")
     cat("grad    ");print(P10)
     cat("fd.grad ");print(fd.P1)
     if (deriv==2) {
       fd.P2 <- .5*(fd.P2 + t(fd.P2))
       cat("hess\n");print(P2)
       cat("fd.hess\n");print(fd.P2)
     }

     cat("\n\n tr(A)... \n")
     cat("grad    ");print(gtrA0)
     cat("fd.grad ");print(fd.gtrA)
     if (deriv==2) {
       fd.htrA <- .5*(fd.htrA + t(fd.htrA))
       cat("hess\n");print(htrA)
       cat("fd.hess\n");print(fd.htrA)
     }
   

     cat("\n Deviance... \n")
     cat("grad    ");print(D10)
     cat("fd.grad ");print(fd.D1)
     if (deriv==2) {
       fd.D2 <- .5*(fd.D2 + t(fd.D2))
       cat("hess\n");print(D2)
       cat("fd.hess\n");print(fd.D2)
     }
   }

   plot(b$db.drho,fd.db,pch=".")
   for (i in 1:ncol(fd.db)) points(b$db.drho[,i],fd.db[,i],pch=19,cex=.3,col=i)

  
   cat("\n\n The objective...\n")
   cat("diter    ");print(diter)
   cat("grad    ");print(as.numeric(grad0))
   cat("fd.grad ");print(as.numeric(fd.grad))
   if (deriv==2) {
     fd.hess <- .5*(fd.hess + t(fd.hess))
     cat("hess\n");print(hess)
     cat("fd.hess\n");print(fd.hess)
   }
   NULL
} ## deriv.check


rt <- function(x,r1) {
## transform of x, asymptoting to values in r1
## returns derivatives wrt to x as well as transform values
## r1[i] == NA for no transform 
  x <- as.numeric(x)
  ind <- x>0 
  rho2 <- rho1 <- rho <- 0*x
  if (length(r1)==1) r1 <- x*0+r1
  h <- exp(x[ind])/(1+exp(x[ind]))
  h1 <- h*(1-h);h2 <- h1*(1-2*h)
  rho[ind] <- r1[ind]*(h-0.5)*2
  rho1[ind] <- r1[ind]*h1*2
  rho2[ind] <- r1[ind]*h2*2
  rho[!ind] <- r1[!ind]*x[!ind]/2
  rho1[!ind] <- r1[!ind]/2
  ind <- is.na(r1)
  rho[ind] <- x[ind]
  rho1[ind] <- 1
  rho2[ind] <- 0
  list(rho=rho,rho1=rho1,rho2=rho2)
} ## rt

rti <- function(r,r1) {
## inverse of rti.
  r <- as.numeric(r)
  ind <- r>0
  x <- r
  if (length(r1)==1) r1 <- x*0+r1
  r2 <- r[ind]*.5/r1[ind] + .5
  x[ind] <- log(r2/(1-r2))
  x[!ind] <- 2*r[!ind]/r1[!ind]
  ind <- is.na(r1)
  x[ind] <- r[ind]
  x
} ## rti

simplyFit <- function(lsp,X,y,Eb,UrS,L,lsp0,offset,U1,Mp,family,weights,
                   control,gamma,scale,conv.tol=1e-6,maxNstep=5,maxSstep=2,
                   maxHalf=30,printWarn=FALSE,scoreType="deviance",
                   mustart = NULL,null.coef=rep(0,ncol(X)),Sl=Sl,nei=NULL,...)
## function with same argument list as `newton' and `bfgs' which simply fits
## the model given the supplied smoothing parameters...
{ reml <- scoreType%in%c("REML","P-REML","ML","P-ML") ## REML/ML indicator
  sname <- if (reml) "REML" else scoreType
  
  ## sanity check L
  if (is.null(L)) L <- diag(length(lsp)) else {
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=length(lsp0)||ncol(L)!=length(lsp)) stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) lsp0 <- rep(0,ncol(L))
  ## initial fit

  b<-gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0, Eb=Eb,UrS=UrS,
     offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=0,
     control=control,gamma=gamma,scale=scale,
     printWarn=FALSE,mustart=mustart,scoreType=scoreType,null.coef=null.coef,Sl=Sl,nei=nei,...)

  if (!is.null(b$warn)&&length(b$warn)>0) for (i in 1:length(b$warn)) warning(b$warn[[i]])
  score <- b[[sname]]

  list(score=score,lsp=lsp,lsp.full=L%*%lsp+lsp0,grad=NULL,hess=NULL,score.hist=NULL,iter=0,conv =NULL,object=b)

} ## simplyFit


newton <- function(lsp,X,y,Eb,UrS,L,lsp0,offset,U1,Mp,family,weights,
                   control,gamma,scale,conv.tol=1e-6,maxNstep=5,maxSstep=2,
                   maxHalf=30,printWarn=FALSE,scoreType="deviance",start=NULL,
                   mustart = NULL,null.coef=rep(0,ncol(X)),pearson.extra,
                   dev.extra=0,n.true=-1,Sl=NULL,edge.correct=FALSE,nei=NULL,...)
## Newton optimizer for GAM reml/gcv/aic optimization that can cope with an 
## indefinite Hessian. Main enhancements are: 
## i) always perturbs the Hessian to +ve definite if indefinite 
## ii) step halves on step failure, without obtaining derivatives until success; 
## (iii) carries start values forward from one evaluation to next to speed convergence;
## iv) Always tries the steepest descent direction as well as the 
##     Newton direction for indefinite problems (step length on steepest trial could
##     be improved here - currently simply halves until descent achieved).    
## L is the matrix such that L%*%lsp + lsp0 gives the logs of the smoothing 
## parameters actually multiplying the S[[i]]'s
## NOTE: an obvious acceleration would use db/dsp to produce improved
##       starting values at each iteration... 
{ ##  inner iteration results need to be accurate enough for conv.tol...
  if (control$epsilon>conv.tol/100) control$epsilon <- conv.tol/100 

  reml <- scoreType%in%c("REML","P-REML","ML","P-ML") ## REML/ML indicator

  sname <- if (reml) "REML" else scoreType
  sname1 <- paste(sname,"1",sep=""); sname2 <- paste(sname,"2",sep="")

  ## sanity check L
  if (is.null(L)) L <- diag(length(lsp)) else {
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=length(lsp0)||ncol(L)!=length(lsp)) stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) lsp0 <- rep(0,nrow(L)) 

  if (reml&&FALSE) { ## NOTE: routine set up to allow upper limits on lsp, but turned off.
    frob.X <- sqrt(sum(X*X))
    lsp.max <- rep(NA,length(lsp0))
    for (i in 1:nrow(L)) { 
      lsp.max[i] <- 16 + log(frob.X/sqrt(sum(UrS[[i]]^2))) - lsp0[i]
      if (lsp.max[i]<2) lsp.max[i] <- 2
    } 
  } else lsp.max <- NULL

  if (!is.null(lsp.max)) { ## then there are upper limits on lsp's
    lsp1.max <- coef(lm(lsp.max-lsp0~L-1)) ## get upper limits on lsp1 scale
    ind <- lsp>lsp1.max
    lsp[ind] <- lsp1.max[ind]-1 ## reset lsp's already over limit
    delta <- rti(lsp,lsp1.max) ## initial optimization parameters
  } else { ## optimization parameters are just lsp
    delta <- lsp
  }

  ## code designed to be turned on during debugging...
  check.derivs <- FALSE;sp.trace <- FALSE
  if (check.derivs) {
     deriv <- 2
     eps <- 1e-4
     deriv.check(x=X, y=y, sp=L%*%lsp+lsp0, Eb=Eb,UrS=UrS,
         offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=deriv,
         control=control,gamma=gamma,scale=scale,
         printWarn=FALSE,start=start,mustart=mustart,
         scoreType=scoreType,eps=eps,null.coef=null.coef,Sl=Sl,nei=nei,...)

     
  }

  ## ... end of debugging code 


  ## initial fit
  initial.lsp <- lsp ## used if edge correcting to set direction of correction
  b<-gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
     offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=2,
     control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=start,
     mustart=mustart,scoreType=scoreType,null.coef=null.coef,pearson.extra=pearson.extra,
     dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)

  mustart <- b$fitted.values
  etastart <- b$linear.predictors
  start <- b$coefficients
  old.score <- score <- b[[sname]];grad <- b[[sname1]];hess <- b[[sname2]]
  
  grad <- t(L)%*%grad
  hess <- t(L)%*%hess%*%L

  if (!is.null(lsp.max)) { ## need to transform to delta space
    rho <- rt(delta,lsp1.max)
    nr <- length(rho$rho1)
    hess <- diag(rho$rho1,nr,nr)%*%hess%*%diag(rho$rho1,nr,nr) + diag(rho$rho2*grad)
    grad <- rho$rho1*grad
  }

  if (reml) score.scale <- abs(log(b$scale.est)) + abs(score) else 
  score.scale <- b$scale.est + abs(score)    
  uconv.ind <- abs(grad) > score.scale*conv.tol
  ## check for all converged too soon, and undo !
  if (!sum(uconv.ind)) uconv.ind <- uconv.ind | TRUE
  score.hist <- rep(NA,200)
  ################################
  ## Start of Newton iteration.... 
  ################################
  qerror.thresh <- .8 ## quadratic approx error to tolerate in a step
  for (i in 1:200) {
   if (control$trace) {
     cat("\n",i,"newton max(|grad|) =",max(abs(grad)),"\n")
   }
   ## debugging code for checking derivatives ....
   okc <- check.derivs 
   while (okc) {
     okc <- FALSE
     eps <- 1e-4
     deriv <- 2
     if (okc) { ## optional call to fitting to facilitate debugging 
       trial.der <- 2 ## can reset if derivs not wanted
       b <- gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
         offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=trial.der,
         control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=start,
         mustart=mustart,scoreType=scoreType,null.coef=null.coef,
         pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)  
     }
     deriv.check(x=X, y=y, sp=L%*%lsp+lsp0, Eb=Eb,UrS=UrS,
         offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=deriv,
         control=control,gamma=gamma,scale=scale,
         printWarn=FALSE,etastart=etastart,start=start,
         scoreType=scoreType,eps=eps,null.coef=null.coef,Sl=Sl,nei=nei,...)
     if (inherits(family,"general.family")) { ## call gam.fit5 checking
       eps <- 1e-6
       spe <- 1e-3
       er <- deriv.check5(x=X, y=y, sp=L%*%lsp+lsp0, 
            weights = weights, start = start,
            offset = offset,Mp=Mp,family = family, 
            control = control,deriv=deriv,eps=eps,spe=spe,
            Sl=Sl,nei=nei,...) ## ignore codetools warning
   
     }
   } ## end of derivative checking

    ## exclude dimensions from Newton step when the derviative is
    ## tiny relative to largest, as this space is likely to be poorly
    ## modelled on scale of Newton step...
    
    uconv.ind1 <- uconv.ind & abs(grad)>max(abs(grad))*.001 
    if (sum(uconv.ind1)==0) uconv.ind1 <- uconv.ind ## nothing left reset
    if (sum(uconv.ind)==0) uconv.ind[which(abs(grad)==max(abs(grad)))] <- TRUE ## need at least 1 to update
    
    ## exclude apparently converged gradients from computation
    hess1 <- hess[uconv.ind,uconv.ind] 
    grad1 <- grad[uconv.ind]
    ## get the trial step ...
    eh <- eigen(hess1,symmetric=TRUE)
    d <- eh$values;U <- eh$vectors
    indef <- (sum(-d > abs(d[1])*.Machine$double.eps^.5)>0) ## indefinite problem
    ## need a different test if there is only one smoothing parameter,
    ## otherwise infinite sp can count as always indefinite...
    if (indef && length(d)==1) indef <- d < -score.scale * .Machine$double.eps^.5
    ## set eigen-values to their absolute value - heuristically appealing
    ## as it avoids very long steps being proposed for indefinte components,
    ## unlike setting -ve e.v.s to very small +ve constant...
    ind <- d < 0
    pdef <- if (sum(ind)>0) FALSE else TRUE ## is it positive definite? 
    d[ind] <- -d[ind] ## see Gill Murray and Wright p107/8
    low.d <- max(d)*.Machine$double.eps^.7
    ind <- d < low.d
    if (sum(ind)>0) pdef <- FALSE ## not certain it is positive definite
    d[ind] <- low.d 
    ind <- d != 0
    d[ind] <- 1/d[ind]
    
    Nstep <- 0 * grad
    Nstep[uconv.ind] <- -drop(U%*%(d*(t(U)%*%grad1))) # (modified) Newton direction
   
    Sstep <- -grad/max(abs(grad)) # steepest descent direction 

    ms <- max(abs(Nstep)) ## note smaller permitted step if !pdef
    mns <- maxNstep

    if (ms>maxNstep) Nstep <- mns * Nstep/ms
    
    sd.unused <- TRUE ## steepest descent direction not yet tried

    ## try the step ...
    if (sp.trace) cat(lsp,"\n")

    if (!is.null(lsp.max)) { ## need to take step in delta space
      delta1 <- delta + Nstep
      lsp1 <- rt(delta1,lsp1.max)$rho ## transform to log sp space
      while (max(abs(lsp1-lsp))>maxNstep) { ## make sure step is not too long
        Nstep <- Nstep / 2 
        delta1 <- delta + Nstep
        lsp1 <- rt(delta1,lsp1.max)$rho
      } 
    } else lsp1 <- lsp + Nstep

    ## if pdef==TRUE then get grad and hess immediately, otherwise postpone as
    ## the steepest descent direction should be tried as well as Newton

    b <- gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0,Eb=Eb,UrS=UrS,
         offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=as.numeric(pdef)*2,
         control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=start,
         mustart=mustart,scoreType=scoreType,null.coef=null.coef,
         pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)    

    ## get the change predicted for this step according to the quadratic model
    pred.change <- sum(grad*Nstep) + 0.5*t(Nstep) %*% hess %*% Nstep
    score1 <- b[[sname]]
 
    ## accept if improvement, else step halve
    ii <- 0 ## step halving counter
    score.change <- score1 - score
    qerror <- abs(pred.change-score.change)/(max(abs(pred.change),abs(score.change))+score.scale*conv.tol) ## quadratic approx error
    if (is.finite(score1) && score.change<0 && pdef && qerror < qerror.thresh) { ## immediately accept step if it worked and positive definite
      old.score <- score 
      mustart <- b$fitted.values
      etastart <- b$linear.predictors
      start <- b$coefficients
      lsp <- lsp1
      score <- b[[sname]]; grad <- b[[sname1]]; hess <- b[[sname2]] 
    
      grad <- t(L)%*%grad
      hess <- t(L)%*%hess%*%L
      
      if (!is.null(lsp.max)) { ## need to transform to delta space
        delta <- delta1
        rho <- rt(delta,lsp1.max)
        nr <- length(rho$rho1)
        hess <- diag(rho$rho1,nr,nr)%*%hess%*%diag(rho$rho1,nr,nr) + diag(rho$rho2*grad)
        grad <- rho$rho1*grad
      }

    } else if (!is.finite(score1) || score1>=score||qerror >= qerror.thresh) { ## initial step failed, try step halving ...
      step <- Nstep ## start with the (pseudo) Newton direction
      while ((!is.finite(score1) || score1>=score ||qerror >= qerror.thresh) && ii < maxHalf) {
        if (ii==3&&i<10) { ## Newton really not working - switch to SD, but keeping step length 
          s.length <- min(sum(step^2)^.5,maxSstep)
          step <- Sstep*s.length/sum(Sstep^2)^.5 ## use steepest descent direction
          sd.unused <- FALSE ## signal that SD already tried
        } else step <- step/2
        if (!is.null(lsp.max)) { ## need to take step in delta space
          delta1 <- delta + step
          lsp1 <- rt(delta1,lsp1.max)$rho ## transform to log sp space
        } else lsp1 <- lsp + step
        b1<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0,Eb=Eb,UrS=UrS,offset = offset,U1=U1,Mp=Mp,
	             family = family,weights=weights,deriv=0,control=control,gamma=gamma,
		     scale=scale,printWarn=FALSE,start=start,mustart=mustart,scoreType=scoreType,
		     null.coef=null.coef,pearson.extra=pearson.extra,dev.extra=dev.extra,
		     n.true=n.true,Sl=Sl,nei=nei,...)
        pred.change <- sum(grad*step) + 0.5*t(step) %*% hess %*% step ## Taylor prediction of change 
        score1 <- b1[[sname]]

	score.change <- score1 - score
	## don't allow step to fail altogether just because of qerror
	qerror <- if (ii>min(4,maxHalf/2)) qerror.thresh/2 else
	          abs(pred.change-score.change)/(max(abs(pred.change),abs(score.change))+score.scale*conv.tol) ## quadratic approx error
        if (is.finite(score1) && score.change < 0 && qerror < qerror.thresh) { ## accept
          if (pdef||!sd.unused) { ## then accept and compute derivatives
            b <- gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0,Eb=Eb,UrS=UrS,
                 offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=2,
                 control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=start,
                 mustart=mustart,scoreType=scoreType,null.coef=null.coef,
                 pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
            mustart <- b$fitted.values 
            etastart <- b$linear.predictors
            start <- b$coefficients
            old.score <- score;lsp <- lsp1
            score <- b[[sname]];grad <- b[[sname1]];hess <- b[[sname2]] 

            grad <- t(L)%*%grad
            hess <- t(L)%*%hess%*%L
            if (!is.null(lsp.max)) { ## need to transform to delta space
               delta <- delta1
               rho <- rt(delta,lsp1.max)
               nr <- length(rho$rho1)
               hess <- diag(rho$rho1,nr,nr)%*%hess%*%diag(rho$rho1,nr,nr) + diag(rho$rho2*grad)
               grad <- rho$rho1*grad
            }
          } else { ## still need to try the steepest descent step to see if it does better
            b <- b1
            score2 <- score1 ## store temporarily and restore below
          }
          score1 <- score - abs(score) - 1 ## make sure that score1 < score (restore once out of loop)
        }  # end of if (score1<= score ) # accept
        if (!is.finite(score1) || score1>=score || qerror >= qerror.thresh) ii <- ii + 1
      } ## end while (score1>score && ii < maxHalf)
      if (!pdef&&sd.unused&&ii<maxHalf) score1 <- score2 ## restored (not needed if termination on ii==maxHalf)
    } ## end of step halving branch

    ## if the problem is not positive definite, and the sd direction has not 
    ## yet been tried then it should be tried now, and the better of the 
    ## newton and steepest steps used (before computing the derivatives)
    ## score1, lsp1 (delta1) contain the best so far...

    if (!pdef&&sd.unused) {
      step <- Sstep*2
      kk <- 0;score2 <- NA
      ok <- TRUE
      while (ok) { ## step length loop for steepest....
        step <- step/2;kk <- kk+1
        if (!is.null(lsp.max)) { ## need to take step in delta space
          delta3 <- delta + step
          lsp3 <- rt(delta3,lsp1.max)$rho ## transform to log sp space
        } else lsp3 <- lsp + step
        b1 <- gam.fit3(x=X, y=y, sp=L%*%lsp3+lsp0,Eb=Eb,UrS=UrS,
              offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=0,
              control=control,gamma=gamma,scale=scale,
              printWarn=FALSE,start=start,mustart=mustart,scoreType=scoreType,
              null.coef=null.coef,pearson.extra=pearson.extra,
              dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
	pred.change <- sum(grad*step) + 0.5*t(step) %*% hess %*% step ## Taylor prediction of change 
        score3 <- b1[[sname]]
      
	score.change <- score3 - score
        qerror <- abs(pred.change-score.change)/(max(abs(pred.change),abs(score.change))+score.scale*conv.tol) ## quadratic approx error 
        if (!is.finite(score2)||(is.finite(score3)&&score3<=score2&&qerror<qerror.thresh)) { ## record step - best SD step so far
          score2 <- score3
          lsp2 <- lsp3
          if (!is.null(lsp.max)) delta2 <- delta3
        }
        ## stop when improvement found, and shorter step is worse...
        if ((is.finite(score2)&&is.finite(score3)&&score2<score&&score3>score2)||kk==40) ok <- FALSE
      } ## while (ok) ## step length control loop

      ## now pick the step that led to the biggest decrease  

      if (is.finite(score2) && score2<score1) {
        lsp1 <- lsp2
        if (!is.null(lsp.max)) delta1 <- delta2
        score1 <- score2
      }

      ## and compute derivatives for the accepted step....

      b <- gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0,Eb=Eb,UrS=UrS,
                 offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=2,
                 control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=start,
                 mustart=mustart,scoreType=scoreType,null.coef=null.coef,
                 pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
      mustart <- b$fitted.values 
      etastart <- b$linear.predictors
      start <- b$coefficients
      old.score <- score;lsp <- lsp1
      score <- b[[sname]]; grad <- b[[sname1]]; hess <- b[[sname2]]  

      grad <- t(L)%*%grad
      hess <- t(L)%*%hess%*%L
      if (!is.null(lsp.max)) { ## need to transform to delta space
          delta <- delta1
          rho <- rt(delta,lsp1.max)
          nr <- length(rho$rho1)
          hess <- diag(rho$rho1,nr,nr)%*%hess%*%diag(rho$rho1,nr,nr) + diag(rho$rho2*grad)
          grad <- rho$rho1*grad
       }

    } ## end of steepest descent trial for indefinite problems

    ## record current score
    score.hist[i] <- score 
   
    ## test for convergence
    converged <- !indef ## not converged if indefinite
    if (reml) score.scale <- abs(log(b$scale.est)) + abs(score) else
    score.scale <- abs(b$scale.est) + abs(score)
    grad2 <- diag(hess)    
    uconv.ind <- (abs(grad) > score.scale*conv.tol*.1)|(abs(grad2)>score.scale*conv.tol*.1)
    if (sum(abs(grad)>score.scale*conv.tol*5)) converged <- FALSE
    if (abs(old.score-score)>score.scale*conv.tol) { 
      if (converged) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
      converged <- FALSE      
    }
    if (ii==maxHalf) converged <- TRUE ## step failure
    if (converged) break
  } ## end of iteration loop
  if (ii==maxHalf) { 
    ct <- "step failed"
    warning("Fitting terminated with step failure - check results carefully")
  } else if (i==200) { 
    ct <- "iteration limit reached"
    warning("Iteration limit reached without full convergence - check carefully")
  } else ct <- "full convergence"
  b$dVkk <- NULL
  
  if (as.logical(edge.correct)&&reml) {
    ## for those smoothing parameters that appear to be at working infinity
    ## reduce them until there is a detectable increase in RE/ML...
    flat <- which(abs(grad2) < abs(grad)*100) ## candidates for reduction
    REML <- b[[sname]]
    alpha <- if (is.logical(edge.correct)) .02 else abs(edge.correct) ## target RE/ML change per sp
    b1 <- b; lsp1 <- lsp
    if (length(flat)) {
      step <- as.numeric(initial.lsp > lsp)*2-1 ## could use sign here 
      for (i in flat) {
        REML <- b1$REML + alpha
        while (b1$REML < REML) {
          lsp1[i] <- lsp1[i] + step[i]
          b1 <- gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0,Eb=Eb,UrS=UrS,
              offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=0,
              control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=start,
              mustart=mustart,scoreType=scoreType,null.coef=null.coef,
              pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
        }
      }
    } ## if length(flat) 
    b1 <- gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0,Eb=Eb,UrS=UrS,
                 offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=2,
                 control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=start,
                 mustart=mustart,scoreType=scoreType,null.coef=null.coef,
                 pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
  
    score1 <- b1[[sname]];grad1 <- b1[[sname1]];hess1 <- b1[[sname2]] 
           
    grad1 <- t(L)%*%grad1
    hess1 <- t(L)%*%hess1%*%L
    if (!is.null(lsp.max)) { ## need to transform to delta space
               delta <- delta1
               rho <- rt(delta,lsp1.max)
               nr <- length(rho$rho1)
               hess1 <- diag(rho$rho1,nr,nr)%*%hess1%*%diag(rho$rho1,nr,nr) + diag(rho$rho2*grad1)
               grad1 <- rho$rho1*grad1
    }
    attr(hess,"edge.correct") <- TRUE
    attr(hess,"hess1") <- hess1
    attr(hess,"db.drho1") <- b1$db.drho
    attr(hess,"dw.drho1") <- b1$dw.drho
    attr(hess,"lsp1") <- lsp1
    attr(hess,"rp") <- b1$rp
  } ## if edge.correct
  ## report any warnings from inner loop if outer not converged...
  #if (!is.null(b$warn)&&length(b$warn)>0&&ct!="full convergence") for (i in 1:length(b$warn)) warning(b$warn[[i]])
  if (!is.null(b$warn)&&length(b$warn)>0) for (i in 1:length(b$warn)) warning(b$warn[[i]])
  list(score=score,lsp=lsp,lsp.full=L%*%lsp+lsp0,grad=grad,hess=hess,iter=i,
       conv =ct,score.hist = score.hist[!is.na(score.hist)],object=b)
} ## newton


bfgs <-  function(lsp,X,y,Eb,UrS,L,lsp0,offset,U1,Mp,family,weights,
                   control,gamma,scale,conv.tol=1e-6,maxNstep=3,maxSstep=2,
                   maxHalf=30,printWarn=FALSE,scoreType="GCV",start=NULL,
                   mustart = NULL,null.coef=rep(0,ncol(X)),pearson.extra=0,
                   dev.extra=0,n.true=-1,Sl=NULL,nei=NULL,...)

## BFGS optimizer to estimate smoothing parameters of models fitted by
## gam.fit3....
##
## L is the matrix such that L%*%lsp + lsp0 gives the logs of the smoothing 
## parameters actually multiplying the S[[i]]'s. sp's do not include the 
## log scale parameter here.
##
## BFGS is based on Nocedal & Wright (2006) Numerical Optimization, Springer.
## In particular the step lengths are chosen to meet the Wolfe conditions
## using their algorithms 3.5 (p60) and 3.6 (p61). On p143 they recommend a post step
## adjustment to the initial Hessian. I can't understand why one would do anything
## other than adjust so that the initial Hessian would give the step taken, and
## indeed the latter adjustment seems to give faster convergence than their 
## proposal, and is therefore implemented.
##
{ zoom <- function(lo,hi) {
  ## local function implementing Algorithm 3.6 of Nocedal & Wright
  ## (2006, p61) Numerical Optimization. Relies on R scoping rules. 
  ## alpha.lo and alpha.hi are the bracketing step lengths.
  ## This routine bisection searches for a step length that meets the
  ## Wolfe conditions. lo and hi are both objects containing fields
  ## `score', `alpha', `dscore', where `dscore' is the derivative of 
  ## the score in the current step direction, `grad' and `mustart'. 
  ## `dscore' will be NULL if the gradiant has yet to be evaluated.
    for (i in 1:40) {
      trial <- list(alpha = (lo$alpha+hi$alpha)/2)
      lsp <- ilsp + step * trial$alpha
      b <- gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
           offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=0,
           control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=lo$start,
           mustart=lo$mustart,scoreType=scoreType,null.coef=null.coef,
           pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)

      trial$mustart <- fitted(b)
      trial$scale.est <- b$scale.est ## previously dev, but this differs from newton
      trial$start <- coef(b)
      trial$score <- b[[sname]]
      
      rm(b)  
      if (trial$score>initial$score+trial$alpha*c1*initial$dscore||trial$score>=lo$score) {
        hi <- trial ## failed Wolfe 1 - insufficient decrease - step too long
      } else { ## met Wolfe 1 so check Wolve 2 - sufficiently positive second derivative?

        b <- gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
           offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
           control=control,gamma=gamma,scale=scale,printWarn=FALSE,
           start=trial$start,mustart=trial$mustart,
           scoreType=scoreType,null.coef=null.coef,pearson.extra=pearson.extra,
           dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)

        trial$grad <- t(L)%*%b[[sname1]];
      
        trial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0)
        trial$scale.est <- b$scale.est;rm(b);
        trial$dscore <- sum(step*trial$grad) ## directional derivative
        
        if (abs(trial$dscore) <= -c2*initial$dscore) return(trial) ## met Wolfe 2

        ## failed Wolfe 2 derivative not increased enough
        if (trial$dscore*(hi$alpha-lo$alpha)>=0) {
          hi <- lo }  
        lo <- trial 
      }  
    } ## end while(TRUE)
    return(NULL) ## failed
  } ## end zoom

  if (control$epsilon>conv.tol/100) control$epsilon <- conv.tol/100

  reml <- scoreType%in%c("REML","P-REML","ML","P-ML") ## REML/ML indicator

  sname <- if (reml) "REML" else scoreType ## name of score 
  sname1 <- paste(sname,"1",sep="")        ## names of its derivative

  ## sanity check L
  if (is.null(L)) L <- diag(length(lsp)) else {
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=length(lsp0)||ncol(L)!=length(lsp)) stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) lsp0 <- rep(0,nrow(L))
 
  ## initial fit...

  initial.lsp <- ilsp <- lsp

  b <- gam.fit3(x=X, y=y, sp=L%*%ilsp+lsp0,Eb=Eb,UrS=UrS,
               offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
               control=control,gamma=gamma,scale=scale,printWarn=FALSE,
               start=start,mustart=mustart,
               scoreType=scoreType,null.coef=null.coef,
               pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)

  initial <- list(alpha = 0,mustart=b$fitted.values,start=coef(b))
  score <- b[[sname]];grad <- t(L)%*%b[[sname1]];

  ## dVkk only refers to smoothing parameters, but sp may contain
  ## extra parameters at start and scale parameter at end. Have
  ## to reduce L accordingly... 
  if (!is.null(family$n.theta)&&family$n.theta>0) {
    ind <- 1:family$n.theta
    nind <- ncol(L) - family$n.theta - if (family$n.theta + nrow(b$dVkk)<nrow(L)) 1 else 0 
    spind <- if (nind>0) family$n.theta+1:nind else rep(0,0)
    rspind <- family$n.theta + 1:nrow(b$dVkk)
  } else {
    nind <- ncol(L) - if (nrow(b$dVkk)<nrow(L)) 1 else 0 
    spind <- if (nind>0) 1:nind else rep(0,0) ## index of smooth parameters
    rspind <- 1:nrow(b$dVkk)
  }  
  L0 <- L[rspind,spind] ##if (nrow(L)!=nrow(b$dVkk)) L[spind,spind] else L
  
  initial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0)
  initial$score <- score;initial$grad <- grad;
  initial$scale.est <- b$scale.est

  if (reml) score.scale <- 1 + abs(initial$score) 
            else score.scale <- abs(initial$scale.est) + abs(initial$score)   
  
  start0 <- coef(b)
  mustart0 <- fitted(b)
  rm(b)

  B <- diag(length(initial$grad)) ## initial Hessian
  feps <- 1e-4;fdgrad <- grad*0
  for (i in 1:length(lsp)) { ## loop to FD for Hessian
     ilsp <- lsp;ilsp[i] <- ilsp[i] + feps 
     b <- gam.fit3(x=X, y=y, sp=L%*%ilsp+lsp0,Eb=Eb,UrS=UrS,
               offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
               control=control,gamma=gamma,scale=scale,printWarn=FALSE,
               start=start0,mustart=mustart0,
               scoreType=scoreType,null.coef=null.coef,
               pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...) 
     grad1 <- t(L)%*%b[[sname1]];
     fdgrad[i] <- (b[[sname]]-score)/feps ## get FD grad for free - useful for debug checks
     B[i,] <- (grad1-grad)/feps 
     rm(b)
  } ## end of FD Hessian loop
  ## force initial Hessian to +ve def and invert... 
  B <- (B+t(B))/2
  eb <- eigen(B,symmetric=TRUE)
  eb$values <- abs(eb$values)
  thresh <- max(eb$values) * 1e-4
  eb$values[eb$values<thresh] <- thresh
  B <- eb$vectors%*%(t(eb$vectors)/eb$values)
  ilsp <- lsp
  max.step <- 200

  c1 <- 1e-4;c2 <- .9 ## Wolfe condition constants

  score.hist <- rep(NA,max.step+1)
  score.hist[1] <- initial$score

  check.derivs <- FALSE;eps <- 1e-5

  uconv.ind <- rep(TRUE,ncol(B))
  rolled.back <- FALSE 

  for (i in 1:max.step) { ## the main BFGS loop
   
    ## get the trial step ...
    step <- initial$grad*0
    step[uconv.ind] <- -B[uconv.ind,uconv.ind]%*%initial$grad[uconv.ind]

    ## following tends to have lower directional grad than above (or full version commented out below)
    #step <- -drop(B%*%initial$grad)
    ## following line would mess up conditions under which Wolfe guarantees update,
    ## *if* based only on grad and not grad and hess...  
    #step[!uconv.ind] <- 0 ## don't move if apparently converged 
    
    if (sum(step*initial$grad)>=0) { ## step not descending!
      ## Following would really be in the positive definite space... 
      ##step[uconv.ind] <- -solve(chol2inv(chol(B))[uconv.ind,uconv.ind],initial$grad[uconv.ind])
      step <- -diag(B)*initial$grad ## simple scaled steepest descent 
      step[!uconv.ind] <- 0 ## don't move if apparently converged 
    }

    ms <- max(abs(step))
    trial <- list()
    if (ms>maxNstep) { 
      trial$alpha <- maxNstep/ms
      alpha.max <- trial$alpha*1.05
      ## step <- maxNstep * step/ms
      #alpha.max <- 1 ## was 50 in place of 1 here and below
    } else {
      trial$alpha <- 1 
      alpha.max <- min(2,maxNstep/ms) ## 1*maxNstep/ms
    }
    initial$dscore <- sum(step*initial$grad)
    prev <- initial

    deriv <- 1 ## only get derivatives immediately for initial step length   
    while(TRUE) { ## step length control Alg 3.5 of N&W (2006, p60)
      lsp <- ilsp + trial$alpha*step
      b <- gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
                    offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=deriv,
                    control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=prev$start,
                    mustart=prev$mustart,scoreType=scoreType,null.coef=null.coef,
                    pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
     
     ### Derivative testing code. Not usually called and not part of BFGS...
     ok <- check.derivs
     while (ok) { ## derivative testing
       #deriv <- 1
       ok <- FALSE ## set to TRUE to re-run (e.g. with different eps)
       deriv.check(x=X, y=y, sp=L%*%lsp+lsp0, Eb=Eb,UrS=UrS,
         offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
         control=control,gamma=gamma,scale=scale,
         printWarn=FALSE,mustart=prev$mustart,start=prev$start,
         scoreType=scoreType,eps=eps,null.coef=null.coef,Sl=Sl,nei=nei,...)
       ## deal with fact that deriv might be 0...	 
       bb <- if (deriv==1) b else gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
                    offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
                    control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=prev$start,
                    mustart=prev$mustart,scoreType=scoreType,null.coef=null.coef,
                    pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
       #fdH <- bb$dH
       fdb.dr <- bb$db.drho*0
       if (!is.null(bb$NCV)) {
         deta.cv <-  attr(bb$NCV,"deta.cv")
	 fd.eta <- deta.cv*0
       }	 
       for (j in 1:ncol(fdb.dr)) { ## check dH and db.drho
         lsp1 <- lsp;lsp1[j] <- lsp[j] + eps
         ba <- gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0,Eb=Eb,UrS=UrS,
                    offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
                    control=control,gamma=gamma,scale=scale,printWarn=FALSE,start=prev$start,
                    mustart=prev$mustart,scoreType=scoreType,null.coef=null.coef,
                    pearson.extra=pearson.extra,dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
        # fdH[[j]] <- (ba$H - bb$H)/eps
         fdb.dr[,j] <- (ba$coefficients - bb$coefficients)/eps
	 if (!is.null(bb$NCV)) fd.eta[,j] <- as.numeric(attr(ba$NCV,"eta.cv")-attr(bb$NCV,"eta.cv"))/eps
       }
     } 
     ### end of derivative testing. BFGS code resumes...
      trial$score <- b[[sname]];
  

      if (deriv>0) {
        trial$grad <- t(L)%*%b[[sname1]];
       
        trial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0) ## curvature testing matrix
        trial$dscore <- sum(trial$grad*step)
        deriv <- 0 
      } else trial$grad <- trial$dscore <- NULL
      trial$mustart <- b$fitted.values
      trial$start <- b$coefficients
      trial$scale.est <- b$scale.est
      
      rm(b)
      Wolfe2 <- TRUE
      ## check the first Wolfe condition (sufficient decrease)...
      if (trial$score>initial$score+c1*trial$alpha*initial$dscore||(deriv==0&&trial$score>=prev$score)) {
         trial <- zoom(prev,trial) ## Wolfe 1 not met so backtracking
         break
      } 

      if (is.null(trial$dscore)) { ## getting gradients
        b <- gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
                      offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
                      control=control,gamma=gamma,scale=scale,printWarn=FALSE,
                      start=trial$start,mustart=trial$mustart,
                      scoreType=scoreType,null.coef=null.coef,pearson.extra=pearson.extra,
                      dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
        trial$grad <- t(L)%*%b[[sname1]];
       
        trial$dscore <- sum(trial$grad*step)
        trial$scale.est <- b$scale.est
        trial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0) ## curvature testing matrix
        rm(b)
      }
      
      ## Note that written this way so that we can pass on to next test when appropriate...
     
      if (abs(trial$dscore) <= -c2*initial$dscore) break; ## `trial' is ok. (2nd Wolfe condition met).
      Wolfe2 <- FALSE

      if (trial$dscore>=0) { ## increase at end of trial step
        trial <- zoom(trial,prev)
        Wolfe2 <- if (is.null(trial)) FALSE else TRUE
        break
      }
      
      prev <- trial
      if (trial$alpha == alpha.max) break ## { trial <- NULL;break;} ## step failed
      trial <- list(alpha = min(prev$alpha*1.3, alpha.max)) ## increase trial step to try to meet Wolfe 2
    } ## end of while(TRUE)

    ## Now `trial' contains a suitable step, or is NULL on complete failure to meet Wolfe,
    ## or contains a step that fails to meet Wolfe2, so that B can not be updated  
    if (is.null(trial)) { ## step failed
      lsp <- ilsp
      if (rolled.back) break ## failed to move, so nothing more can be done.
      ## check for working infinite smoothing params... 
      uconv.ind <- abs(initial$grad) > score.scale*conv.tol*.1 
      uconv.ind[spind] <- uconv.ind[spind] | abs(initial$dVkk) > score.scale * conv.tol*.1
      if (sum(!uconv.ind)==0) break ## nothing to move back so nothing more can be done.
      trial <- initial ## reset to allow roll back
      converged <- TRUE ## only to signal that roll back should be tried
    } else { ## update the Hessian etc...
     
      yg <- trial$grad-initial$grad
      step <- step*trial$alpha
      rho <- sum(yg*step)
      if (rho>0) { #Wolfe2) { ## only update if Wolfe2 is met, otherwise B can fail to be +ve def.
        if (i==1) { ## initial step --- adjust Hessian as p143 of N&W
          B <- B * trial$alpha ## this is my version 
          ## B <- B * sum(yg*step)/sum(yg*yg) ## this is N&W
        }
        rho <- 1/rho # sum(yg*step)
        B <- B - rho*step%*%(t(yg)%*%B)

        ## Note that Wolfe 2 guarantees that rho>0 and updated B is 
        ## +ve definite (left as an exercise for the reader)...
        B <- B - rho*(B%*%yg)%*%t(step) + rho*step%*%t(step)
      }

      score.hist[i+1] <- trial$score

      lsp <- ilsp <- ilsp + step 

      ## test for convergence
      converged <- TRUE
      if (reml) score.scale <- 1 + abs(trial$score) ## abs(log(trial$dev/nrow(X))) + abs(trial$score)
      else score.scale <- abs(trial$scale.est) + abs(trial$score)  ##trial$dev/nrow(X) + abs(trial$score)    
      uconv.ind <- abs(trial$grad) > score.scale*conv.tol 
      if (sum(uconv.ind)) converged <- FALSE
      #if (length(uconv.ind)>length(trial$dVkk)) trial$dVkk <- c(trial$dVkk,score.scale)
      ## following must be tighter than convergence...
      uconv.ind <- abs(trial$grad) > score.scale*conv.tol*.1 
      uconv.ind[spind] <- uconv.ind[spind] | abs(trial$dVkk) > score.scale * conv.tol*.1 
      if (abs(initial$score-trial$score) > score.scale*conv.tol) { 
        if (!sum(uconv.ind)) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
        converged <- FALSE      
      }
    }
    ## roll back any `infinite' smoothing parameters to the point at
    ## which score carries some information about them and continue 
    ## optimization. Guards against early long steps missing shallow minimum. 
    if (converged) { ## try roll back for `working inf' sps...
        if (sum(!uconv.ind)==0||rolled.back) break
        rolled.back <- TRUE
        counter <- 0
        uconv.ind0 <- uconv.ind 
        while (sum(!uconv.ind0)>0&&counter<5) {
          ## shrink towards initial values...
          lsp[!uconv.ind0] <- lsp[!uconv.ind0]*.8 + initial.lsp[!uconv.ind0]*.2
          b <- gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
                      offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
                      control=control,gamma=gamma,scale=scale,printWarn=FALSE,
                      start=trial$start,mustart=trial$mustart,
                      scoreType=scoreType,null.coef=null.coef,pearson.extra=pearson.extra,
                      dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
          trial$score <- b[[sname]]
          trial$grad <- t(L)%*%b[[sname1]];

          trial$dscore <- sum(trial$grad*step)
          trial$scale.est <- b$scale.est
          trial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0) ## curvature testing matrix 
          #if (length(uconv.ind)>length(trial$dVkk)) trial$dVkk <- c(trial$dVkk,score.scale)
          rm(b);counter <- counter + 1
          ## note that following rolls back until there is clear signal in derivs...
          uconv.ind0 <- abs(trial$grad) > score.scale*conv.tol*20        
          uconv.ind0[spind] <- uconv.ind0[spind] |  abs(trial$dVkk) > score.scale * conv.tol * 20
          uconv.ind0 <- uconv.ind0 | uconv.ind ## make sure we don't start rolling back unproblematic sps 
        }
        uconv.ind <- uconv.ind | TRUE
        ## following line is tempting, but will likely reduce usefullness of B as approximtion 
        ## to inverse Hessian on return...
        ##B <- diag(diag(B),nrow=nrow(B))
        ilsp <- lsp
    }
    
    initial <- trial
    initial$alpha <- 0
      
  } ## end of iteration loop


  if (is.null(trial)) { 
    ct <- "step failed"
    lsp <- ilsp
    trial <- initial
  }
  else if (i==max.step) ct <- "iteration limit reached" 
  else ct <- "full convergence"
  ## final fit
  b <- gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0,Eb=Eb,UrS=UrS,
                offset = offset,U1=U1,Mp=Mp,family = family,weights=weights,deriv=1,
                control=control,gamma=gamma,scale=scale,printWarn=FALSE,
                start=trial$start,mustart=trial$mustart,
                scoreType=scoreType,null.coef=null.coef,pearson.extra=pearson.extra,
                dev.extra=dev.extra,n.true=n.true,Sl=Sl,nei=nei,...)
  score <- b[[sname]];grad <- t(L)%*%b[[sname1]];

  if (!is.null(b$Vg)) {
    M <- ncol(b$db.drho)
    b$Vg <- (B%*%t(L)%*%b$Vg%*%L%*%B)[1:M,1:M] ## sandwich estimate of 
    db.drho <- b$db.drho%*%L[1:M,1:M,drop=FALSE]
    b$Vc <- db.drho %*% b$Vg %*% t(db.drho) ## correction term for cov matrices
  }
  b$dVkk <- NULL
  ## get approximate Hessian...
  ev <- eigen(B,symmetric=TRUE)
  ind <- ev$values>max(ev$values)*.Machine$double.eps^.9
  ev$values[ind] <- 1/ev$values[ind]
  ev$values[!ind] <- 0
  B <- ev$vectors %*% (ev$values*t(ev$vectors))
  #if (!is.null(b$warn)&&length(b$warn)>0&&ct!="full convergence") for (j in 1:length(b$warn)) warning(b$warn[[j]])
  if (!is.null(b$warn)&&length(b$warn)>0) for (j in 1:length(b$warn)) warning(b$warn[[j]])
  list(score=score,lsp=lsp,lsp.full=L%*%lsp+lsp0,grad=grad,hess=B,iter=i,conv =ct,
       score.hist=score.hist[!is.na(score.hist)],object=b)
} ## end of bfgs



gam2derivative <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the derivatives of the GCV or UBRE score w.r.t the 
## smoothing parameters for the model.
## args is a list containing the arguments for gam.fit3
## For use as optim() objective gradient
{ reml <- args$scoreType%in%c("REML","P-REML","ML","P-ML") ## REML/ML indicator
  sname <- if (reml) "REML" else args$scoreType
  sname1 <- paste(sname,"1",sep=""); 
  if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp,Eb=args$Eb,UrS=args$UrS,
     offset = args$offset,U1=args$U1,Mp=args$Mp,family = args$family,weights=args$w,deriv=1,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     null.coef=args$null.coef,n.true=args$n.true,Sl=args$Sl,nei=args$nei,...)
  ret <- b[[sname1]]
  if (!is.null(args$L)) ret <- t(args$L)%*%ret
  ret
} ## gam2derivative

gam2objective <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the GCV or UBRE score for the model.
## args is a list containing the arguments for gam.fit3
## For use as optim() objective
{ reml <- args$scoreType%in%c("REML","P-REML","ML","P-ML") ## REML/ML indicator
  sname <- if (reml) "REML" else args$scoreType
  if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp,Eb=args$Eb,UrS=args$UrS,
     offset = args$offset,U1=args$U1,Mp=args$Mp,family = args$family,weights=args$w,deriv=0,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     null.coef=args$null.coef,n.true=args$n.true,Sl=args$Sl,start=args$start,nei=args$nei,...)
  ret <- b[[sname]]
  attr(ret,"full.fit") <- b
  ret
} ## gam2objective



gam4objective <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the GCV or UBRE score for the model.
## args is a list containing the arguments for gam.fit3
## For use as nlm() objective
{ reml <- args$scoreType%in%c("REML","P-REML","ML","P-ML") ## REML/ML indicator
  sname <- if (reml) "REML" else args$scoreType
  sname1 <- paste(sname,"1",sep=""); 
  if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp, Eb=args$Eb,UrS=args$UrS,
     offset = args$offset,U1=args$U1,Mp=args$Mp,family = args$family,weights=args$w,deriv=1,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     null.coef=args$null.coef,Sl=args$Sl,start=args$start,nei=args$nei,...)
  ret <- b[[sname]]
  at <- b[[sname1]]

  attr(ret,"full.fit") <- b

  if (!is.null(args$L)) at <- t(args$L)%*%at

  attr(ret,"gradient") <- at
  ret
} ## gam4objective

##
## The following fix up family objects for use with gam.fit3
##

fix.family.link.general.family <- function(fam) fix.family.link.family(fam)

fix.family.link.extended.family <- function(fam) {
## extended families require link derivatives in ratio form.
## g2g= g''/g'^2, g3g = g'''/g'^3, g4g = g''''/g'^4 - these quanitities are often
## less overflow prone than the raw derivatives
  if (!is.null(fam$g2g)&&!is.null(fam$g3g)&&!is.null(fam$g4g)) return(fam)  
  link <- fam$link
  if (link=="identity") {
    fam$g2g <- fam$g3g <- fam$g4g <- 
    function(mu) rep.int(0,length(mu))
  } else if (link == "log") {
    fam$g2g <- function(mu) rep(-1,length(mu))
    fam$g3g <- function(mu) rep(2,length(mu))
    fam$g4g <- function(mu) rep(-6,length(mu))
  } else if (link == "inverse") {
    ## g'(mu) = -1/mu^2
    fam$g2g <- function(mu) 2*mu     ## g'' = 2/mu^3
    fam$g3g <- function(mu) 6*mu^2   ## g''' = -6/mu^4
    fam$g4g <- function(mu) 24*mu^3     ## g'''' = 24/mu^5
  } else if (link == "logit") {
    ## g = log(mu/(1-mu)) g' = 1/(1-mu) + 1/mu = 1/(mu*(1-mu))
    fam$g2g <- function(mu) mu^2 - (1-mu)^2      ## g'' = 1/(1 - mu)^2 - 1/mu^2
    fam$g3g <- function(mu) 2*mu^3 + 2*(1-mu)^3  ## g''' = 2/(1 - mu)^3 + 2/mu^3
    fam$g4g <- function(mu) 6*mu^4 - 6*(1-mu)^4  ## g'''' = 6/(1-mu)^4 - 6/mu^4
  } else if (link == "sqrt") {
  ## g = sqrt(mu); g' = .5*mu^-.5
    fam$g2g <- function(mu) - mu^-.5  ## g'' = -.25 * mu^-1.5
    fam$g3g <- function(mu) 3 * mu^-1 ## g''' = .375 * mu^-2.5
    fam$g4g <- function(mu) -15 * mu^-1.5 ## -0.9375 * mu^-3.5
  } else if (link == "probit") {
  ## g(mu) = qnorm(mu); 1/g' = dmu/deta = 1/dnorm(eta)
    fam$g2g <- function(mu) { 
      #eta <- fam$linkfun(mu)
      eta <- qnorm(mu)
      ## g'' = eta/fam$mu.eta(eta)^2
      eta
    }
    fam$g3g <- function(mu) {
      #eta <-  fam$linkfun(mu)
      eta <- qnorm(mu)
      ## g''' = (1 + 2*eta^2)/fam$mu.eta(eta)^3
      (1 + 2*eta^2)
    }
    fam$g4g <- function(mu) {
       #eta <-  fam$linkfun(mu)
       eta <- qnorm(mu)
       ## g'''' = (7*eta + 6*eta^3)/fam$mu.eta(eta)^4
       (7*eta + 6*eta^3)
    }
  } else if (link == "cauchit") {
  ## uses general result that if link is a quantile function then 
  ## d mu / d eta = f(eta) where f is the density. Link derivative
  ## is one over this... repeated differentiation w.r.t. mu using chain
  ## rule gives results...
    fam$g2g <- function(mu) { 
     #eta <- fam$linkfun(mu)
     eta <- qcauchy(mu)
     ## g'' = 2*pi*pi*eta*(1+eta*eta)
     eta/(1+eta*eta)
    }
    fam$g3g <- function(mu) { 
     #eta <- fam$linkfun(mu)
     eta <- qcauchy(mu)
     eta2 <- eta*eta
     ## g''' = 2*pi*pi*pi*(1+3*eta2)*(1+eta2)
     (1+3*eta2)/(1+eta2)^2
    }
    fam$g4g <- function(mu) { 
     #eta <- fam$linkfun(mu)
     eta <- qcauchy(mu)
     eta2 <- eta*eta
     ## g'''' = 2*pi^4*(8*eta+12*eta2*eta)*(1+eta2)
     ((8+ 12*eta2)/(1+eta2)^2)*(eta/(1+eta2))
    }
  } else if (link == "cloglog") {
    ## g = log(-log(1-mu)), g' = -1/(log(1-mu)*(1-mu))
    fam$g2g <- function(mu) { l1m <- log1p(-mu)
      -l1m - 1
    }
    fam$g3g <- function(mu) { l1m <- log1p(-mu)
       l1m*(2*l1m + 3) + 2    
    }
    fam$g4g <- function(mu){
      l1m <- log1p(-mu)
      -l1m*(l1m*(6*l1m+11)+12)-6
    }
  } else stop("link not implemented for extended families")
  ## avoid storing the calling environment of fix.family.link... 
  environment(fam$g2g) <-  environment(fam$g3g) <-  environment(fam$g4g) <- environment(fam$linkfun)
  return(fam)
} ## fix.family.link.extended.family

fix.family.link.family <- function(fam)
# adds d2link the second derivative of the link function w.r.t. mu
# to the family supplied, as well as a 3rd derivative function 
# d3link...
# All d2link and d3link functions have been checked numerically. 
{ if (!inherits(fam,"family")) stop("fam not a family object")
  if (is.null(fam$canonical)) { ## note the canonical link - saves effort in full Newton
    if (fam$family=="gaussian") fam$canonical <- "identity" else
    if (fam$family=="poisson"||fam$family=="quasipoisson") fam$canonical <- "log" else
    if (fam$family=="binomial"||fam$family=="quasibinomial") fam$canonical <- "logit" else
    if (fam$family=="Gamma") fam$canonical <- "inverse" else
    if (fam$family=="inverse.gaussian") fam$canonical <- "1/mu^2" else
    fam$canonical <- "none"
  }
  if (!is.null(fam$d2link)&&!is.null(fam$d3link)&&!is.null(fam$d4link)) return(fam) 
  link <- fam$link
  if (length(link)>1) {
    if (fam$family=="quasi") # then it's a power link
    { lambda <- log(fam$linkfun(exp(1))) ## the power, if > 0
      if (lambda<=0) { fam$d2link <- function(mu) -1/mu^2
        fam$d3link <- function(mu) 2/mu^3
        fam$d4link <- function(mu) -6/mu^4
      } else { fam$d2link <- function(mu) lambda*(lambda-1)*mu^(lambda-2)
        fam$d3link <- function(mu) (lambda-2)*(lambda-1)*lambda*mu^(lambda-3)
        fam$d4link <- function(mu) (lambda-3)*(lambda-2)*(lambda-1)*lambda*mu^(lambda-4)
      }
    } else stop("unrecognized (vector?) link")
  } else if (link=="identity") {
    fam$d4link <- fam$d3link <- fam$d2link <- 
    function(mu) rep.int(0,length(mu))
  } else if (link == "log") {
    fam$d2link <- function(mu) -1/mu^2
    fam$d3link <- function(mu) 2/mu^3
    fam$d4link <- function(mu) -6/mu^4
  } else if (link == "inverse") {
    fam$d2link <- function(mu) 2/mu^3
    fam$d3link <- function(mu) { mu <- mu*mu;-6/(mu*mu)}
    fam$d4link <- function(mu) { mu2 <- mu*mu;24/(mu2*mu2*mu)}
  } else if (link == "logit") {
    fam$d2link <- function(mu) 1/(1 - mu)^2 - 1/mu^2
    fam$d3link <- function(mu) 2/(1 - mu)^3 + 2/mu^3
    fam$d4link <- function(mu) 6/(1-mu)^4 - 6/mu^4
  } else if (link == "probit") {
    fam$d2link <- function(mu) { 
      #eta <- fam$linkfun(mu)
      eta <- qnorm(mu)
      #eta/fam$mu.eta(eta)^2
      eta/pmax(dnorm(eta), .Machine$double.eps)^2
    }
    fam$d3link <- function(mu) {
      #eta <-  fam$linkfun(mu)
      eta <- qnorm(mu)
      #(1 + 2*eta^2)/fam$mu.eta(eta)^3
      (1 + 2*eta^2)/pmax(dnorm(eta), .Machine$double.eps)^3
    }
    fam$d4link <- function(mu) {
       #eta <-  fam$linkfun(mu)
       eta <- qnorm(mu)
       #(7*eta + 6*eta^3)/fam$mu.eta(eta)^4
       (7*eta + 6*eta^3)/pmax(dnorm(eta), .Machine$double.eps)^4
    }
  } else if (link == "cloglog") {
  ## g = log(-log(1-mu)), g' = -1/(log(1-mu)*(1-mu))
    fam$d2link <- function(mu) { l1m <- log1p(-mu)
      -1/((1 - mu)^2*l1m) *(1+ 1/l1m)
    }
    fam$d3link <- function(mu) { l1m <- log1p(-mu)
       mu3 <- (1-mu)^3
      (-2 - 3*l1m - 2*l1m^2)/mu3/l1m^3
    }
    fam$d4link <- function(mu){
      l1m <- log1p(-mu)
      mu4 <- (1-mu)^4
      ( - 12 - 11 * l1m - 6 * l1m^2 - 6/l1m )/mu4  /l1m^3
    }
  } else if (link == "sqrt") {
    fam$d2link <- function(mu) -.25 * mu^-1.5
    fam$d3link <- function(mu) .375 * mu^-2.5
    fam$d4link <- function(mu) -0.9375 * mu^-3.5
  } else if (link == "cauchit") {
  ## uses general result that if link is a quantile function then 
  ## d mu / d eta = f(eta) where f is the density. Link derivative
  ## is one over this... repeated differentiation w.r.t. mu using chain
  ## rule gives results...
    fam$d2link <- function(mu) { 
     #eta <- fam$linkfun(mu)
     eta <- qcauchy(mu)
     2*pi*pi*eta*(1+eta*eta)
    }
    fam$d3link <- function(mu) { 
     #eta <- fam$linkfun(mu)
     eta <- qcauchy(mu)
     eta2 <- eta*eta
     2*pi*pi*pi*(1+3*eta2)*(1+eta2)
    }
    fam$d4link <- function(mu) { 
     #eta <- fam$linkfun(mu)
     eta <- qcauchy(mu)
     eta2 <- eta*eta
     2*pi^4*(8*eta+12*eta2*eta)*(1+eta2)
    }
  } else if (link == "1/mu^2") {
    fam$d2link <- function(mu) 6 * mu^-4
    fam$d3link <- function(mu) -24 * mu^-5
    fam$d4link <- function(mu) 120 * mu^-6
  } else if (substr(link,1,3)=="mu^") { ## it's a power link
    ## note that lambda <=0 gives log link so don't end up here
    lambda <- get("lambda",environment(fam$linkfun))
    fam$d2link <- function(mu) (lambda*(lambda-1)) * mu^{lambda-2}
    fam$d3link <- function(mu) (lambda*(lambda-1)*(lambda-2)) * mu^{lambda-3}
    fam$d4link <- function(mu) (lambda*(lambda-1)*(lambda-2)*(lambda-3)) * mu^{lambda-4}
  } else stop("link not recognised")
  ## avoid giant environments being stored....
  environment(fam$d2link) <-  environment(fam$d3link) <-  environment(fam$d4link) <- environment(fam$linkfun)
  return(fam)
} ## fix.family.link.family


## NOTE: something horrible can happen here. The way method dispatch works, the
## environment attached to functions created in fix.family.link is the environment
## from which fix.family.link was called - and this whole environment is stored
## with the created function - in the gam context that means the model matrix is
## stored invisibly away for no useful purpose at all. pryr:::object_size will
## show the true stored size of an object with hidden environments. But environments
## of functions created in method functions should be set explicitly to something
## harmless (see ?environment for some possibilities, empty is rarely a good idea)
## 9/2017

fix.family.link <- function(fam) UseMethod("fix.family.link")

fix.family.var <- function(fam)
# adds dvar the derivative of the variance function w.r.t. mu
# to the family supplied, as well as d2var the 2nd derivative of 
# the variance function w.r.t. the mean. (All checked numerically). 
{ if (inherits(fam,"extended.family")) return(fam)
  if (!inherits(fam,"family")) stop("fam not a family object")
  if (!is.null(fam$dvar)&&!is.null(fam$d2var)&&!is.null(fam$d3var)) return(fam) 
  family <- fam$family
  fam$scale <- -1
  if (family=="gaussian") {
    fam$d3var <- fam$d2var <- fam$dvar <- function(mu) rep.int(0,length(mu))
  } else if (family=="poisson"||family=="quasipoisson") {
    fam$dvar <- function(mu) rep.int(1,length(mu))
    fam$d3var <- fam$d2var <- function(mu) rep.int(0,length(mu))
    if (family=="poisson") fam$scale <- 1
  } else if (family=="binomial"||family=="quasibinomial") {
    fam$dvar <- function(mu) 1-2*mu
    fam$d2var <- function(mu) rep.int(-2,length(mu))
    fam$d3var <- function(mu) rep.int(0,length(mu))
    if (family=="binomial") fam$scale <- 1
  } else if (family=="Gamma") {
    fam$dvar <- function(mu) 2*mu
    fam$d2var <- function(mu) rep.int(2,length(mu))
    fam$d3var <- function(mu) rep.int(0,length(mu))
  } else if (family=="quasi") {
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
    fam$d3var <- switch(fam$varfun,
       constant = function(mu) rep.int(0,length(mu)),
       "mu(1-mu)" = function(mu) rep.int(0,length(mu)),
       mu = function(mu) rep.int(0,length(mu)),
       "mu^2" = function(mu) rep.int(0,length(mu)),
       "mu^3" = function(mu) rep.int(6,length(mu))           
    )
  } else if (family=="inverse.gaussian") {
    fam$dvar <- function(mu) 3*mu^2
    fam$d2var <- function(mu) 6*mu
    fam$d3var <- function(mu) rep.int(6,length(mu)) 
  } else stop("family not recognised")
  environment(fam$dvar) <-  environment(fam$d2var) <-  environment(fam$d3var) <- environment(fam$linkfun)
  return(fam)
} ## fix.family.var


fix.family.ls <- function(fam)
# adds ls the log saturated likelihood and its derivatives
# w.r.t. the scale parameter to the family object.
{ if (!inherits(fam,"family")) stop("fam not a family object")
  if (!is.null(fam$ls)) return(fam) 
  family <- fam$family
  if (family=="gaussian") {
    fam$ls <- function(y,w,n,scale) { 
      nobs <- sum(w>0)
      c(-nobs*log(2*pi*scale)/2 + sum(log(w[w>0]))/2,-nobs/(2*scale),nobs/(2*scale*scale))
    }
  } else if (family=="poisson") {
    fam$ls <- function(y,w,n,scale) {
      res <- rep(0,3)
      res[1] <- sum(dpois(y,y,log=TRUE)*w)
      res
    }
  } else if (family=="binomial") {
    fam$ls <- function(y,w,n,scale) { 
      c(-binomial()$aic(y,n,y,w,0)/2,0,0)
    }
  } else if (family=="Gamma") {
    fam$ls <- function(y,w,n,scale) {
      res <- rep(0,3)
      y <- y[w>0];w <- w[w>0]
      scale <- scale/w
      k <- -lgamma(1/scale) - log(scale)/scale - 1/scale 
      res[1] <- sum(k-log(y))
      k <- (digamma(1/scale)+log(scale))/(scale*scale)
      res[2] <- sum(k/w)
      k <- (-trigamma(1/scale)/(scale) + (1-2*log(scale)-2*digamma(1/scale)))/(scale^3)
      res[3] <- sum(k/w^2) 
      res
    }
  } else if (family=="quasi"||family=="quasipoisson"||family=="quasibinomial") {
    ## fam$ls <- function(y,w,n,scale) rep(0,3)
    ## Uses extended quasi-likelihood form...
    fam$ls <- function(y,w,n,scale) { 
      nobs <- sum(w>0)
      c(-nobs*log(scale)/2 + sum(log(w[w>0]))/2,-nobs/(2*scale),nobs/(2*scale*scale))
    }
  } else if (family=="inverse.gaussian") {
    fam$ls <- function(y,w,n,scale) {
      ii <- w>0
      nobs <- sum(ii)
      c(-sum(log(2*pi*scale*y[ii]^3))/2 + sum(log(w[ii]))/2,-nobs/(2*scale),nobs/(2*scale*scale))
      ## c(-sum(w*log(2*pi*scale*y^3))/2,-sum(w)/(2*scale),sum(w)/(2*scale*scale))
    }
  } else stop("family not recognised")
  environment(fam$ls) <- environment(fam$linkfun)
  return(fam)
} ## fix.family.ls

fix.family <- function(fam) {
## allows families to be patched...
   if (fam$family[1]=="gaussian") { ## sensible starting values given link...
     fam$initialize <- expression({
     n <- rep.int(1, nobs)
     if (family$link == "inverse") mustart <- y + (y==0)*sd(y)*.01 else
     if (family$link == "log") mustart <- pmax(y,.01*sd(y)) else
     mustart <- y
     })
  }
  fam
} ## fix.family


negbin <- function (theta = stop("'theta' must be specified"), link = "log") { 
## modified from Venables and Ripley's MASS library to work with gam.fit3,
## and to allow a range of `theta' values to be specified
## single `theta' to specify fixed value; 2 theta values (first smaller than second)
## are limits within which to search for theta; otherwise supplied values make up 
## search set.
## Note: to avoid warnings, get(".Theta")[1] is used below. Otherwise the initialization
##       call to negbin can generate warnings since get(".Theta") returns a vector
##       during initialization (only).
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
        else stop(gettextf("%s link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"",linktemp))
    }
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", theta, envir = env)
    variance <- function(mu) mu + mu^2/get(".Theta")[1]
    ## dvaraince/dmu needed as well
    dvar <- function(mu) 1 + 2*mu/get(".Theta")[1]
    ## d2variance/dmu...
    d2var <- function(mu) rep(2/get(".Theta")[1],length(mu))
    d3var <- function(mu) rep(0,length(mu))
    getTheta <- function() get(".Theta")
    validmu <- function(mu) all(mu > 0)

    dev.resids <- function(y, mu, wt) { Theta <- get(".Theta")[1]
      2 * wt * (y * log(pmax(1, y)/mu) - 
        (y + Theta) * log((y + Theta)/(mu + Theta))) 
    }
    aic <- function(y, n, mu, wt, dev) {
        Theta <- get(".Theta")[1]
        term <- (y + Theta) * log(mu + Theta) - y * log(mu) +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
        2 * sum(term * wt)
    }
    ls <- function(y,w,n,scale) {
       Theta <- get(".Theta")[1]
       ylogy <- y;ind <- y>0;ylogy[ind] <- y[ind]*log(y[ind])
       term <- (y + Theta) * log(y + Theta) - ylogy +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
       c(-sum(term*w),0,0)
    }
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })

    rd <- function(mu,wt,scale) {
      Theta <- get(".Theta")[1]
      rnbinom(n=length(mu),size=Theta,mu=mu)
    }

    qf <- function(p,mu,wt,scale) {
      Theta <- get(".Theta")[1]
      qnbinom(p,size=Theta,mu=mu)
    }
 
    environment(qf) <- environment(rd) <- environment(dvar) <- environment(d2var) <- 
    environment(d3var) <-environment(variance) <- environment(validmu) <- 
    environment(ls) <- environment(dev.resids) <- environment(aic) <- environment(getTheta) <- env
    famname <- paste("Negative Binomial(", format(round(theta,3)), ")", sep = "")
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance,dvar=dvar,d2var=d2var,d3var=d3var, dev.resids = dev.resids,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls,
        validmu = validmu, valideta = stats$valideta,getTheta = getTheta,qf=qf,rd=rd,canonical=""), class = "family")
} ## negbin



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
} ## totalPenalty

totalPenaltySpace <- function(S,H,off,p)
{ ## function to obtain (orthogonal) basis for the null space and 
  ## range space of the penalty, and obtain actual null space dimension
  ## components are roughly rescaled to avoid any dominating

  Hscale <- sqrt(sum(H*H));
  if (Hscale==0) H <- NULL ## H was all zeroes anyway!

  if (is.null(H)) St <- matrix(0,p,p)
  else { St <- H/sqrt(sum(H*H)); 
    if (ncol(H)!=p||nrow(H)!=p) stop("H has wrong dimension")
  }
  m <- length(S)
  if (m>0) for (i in 1:m) {
    k0 <- off[i]
    k1 <- k0 + nrow(S[[i]]) - 1
    St[k0:k1,k0:k1] <- St[k0:k1,k0:k1] + S[[i]]/sqrt(sum(S[[i]]*S[[i]]))
  }
  es <- eigen(St,symmetric=TRUE)
  ind <- es$values>max(es$values)*.Machine$double.eps^.66
  Y <- es$vectors[,ind,drop=FALSE]  ## range space
  Z <- es$vectors[,!ind,drop=FALSE] ## null space - ncol(Z) is null space dimension
  E <- sqrt(as.numeric(es$values[ind]))*t(Y) ## E'E = St
  list(Y=Y,Z=Z,E=E)
} ## totalPenaltySpace



mini.roots <- function(S,off,np,rank=NULL)
# function to obtain square roots, B[[i]], of S[[i]]'s having as few
# columns as possible. S[[i]]=B[[i]]%*%t(B[[i]]). np is the total number
# of parameters. S is in packed form. rank[i] is optional supplied rank 
# of S[[i]], rank[i] < 1, or rank=NULL to estimate.
{ m<-length(S)
  if (m<=0) return(list())
  B<-S
  if (is.null(rank)) rank <- rep(-1,m)
  for (i in 1:m)
  { b <- mroot(S[[i]],rank=rank[i]) 
    B[[i]] <- matrix(0,np,ncol(b))
    B[[i]][off[i]:(off[i]+nrow(b)-1),] <- b
  }
  B
}


ldTweedie0 <- function(y,mu=y,p=1.5,phi=1,rho=NA,theta=NA,a=1.001,b=1.999) {
## evaluates log Tweedie density for 1<=p<=2, using series summation of
## Dunn & Smyth (2005) Statistics and Computing 15:267-280.
## Original fixed p and phi version.

  if (!is.na(rho)&&!is.na(theta)) { ## use rho and theta and get derivs w.r.t. these
    if (length(rho)>1||length(theta)>1) stop("only scalar `rho' and `theta' allowed.")
    if (a>=b||a<=1||b>=2) stop("1<a<b<2 (strict) required")
    work.param <- TRUE
    th <- theta;phi <- exp(rho)
    p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1) 
    dpth1 <- if (th>0) exp(-th)*(b-a)/(1+exp(-th))^2 else exp(th)*(b-a)/(exp(th)+1)^2
    dpth2 <- if (th>0) ((a-b)*exp(-th)+(b-a)*exp(-2*th))/(exp(-th)+1)^3 else
                   ((a-b)*exp(2*th)+(b-a)*exp(th))/(exp(th)+1)^3
  } else { ## still need working params for tweedious call...
    work.param <- FALSE 
    if (length(p)>1||length(phi)>1) stop("only scalar `p' and `phi' allowed.")
    rho <- log(phi)
    if (p>1&&p<2) {
      if (p <= a) a <- (1+p)/2
      if (p >= b) b <- (2+p)/2
      pabp <- (p-a)/(b-p)
      theta <- log((p-a)/(b-p))
      dthp1 <- (1+pabp)/(p-a)
      dthp2 <- (pabp+1)/((p-a)*(b-p)) -(pabp+1)/(p-a)^2
    }
  }

  if (p<1||p>2) stop("p must be in [1,2]")
  ld <- cbind(y,y,y);ld <- cbind(ld,ld*NA)
  if (p == 2) { ## It's Gamma
    if (sum(y<=0)) stop("y must be strictly positive for a Gamma density")
    ld[,1] <- dgamma(y, shape = 1/phi,rate = 1/(phi * mu),log=TRUE)
    ld[,2] <- (digamma(1/phi) + log(phi) - 1 + y/mu - log(y/mu))/(phi*phi)
    ld[,3] <- -2*ld[,2]/phi + (1-trigamma(1/phi)/phi)/(phi^3)
    return(ld)
  }  

  if (length(mu)==1) mu <- rep(mu,length(y))

  if (p == 1) { ## It's Poisson like
    ## ld[,1] <- dpois(x = y/phi, lambda = mu/phi,log=TRUE)
    if (all.equal(y/phi,round(y/phi))!=TRUE) stop("y must be an integer multiple of phi for Tweedie(p=1)")
    ind <- (y!=0)|(mu!=0) ## take care to deal with y log(mu) when y=mu=0
    bkt <- y*0
    bkt[ind] <- (y[ind]*log(mu[ind]/phi) - mu[ind])
    dig <- digamma(y/phi+1)
    trig <- trigamma(y/phi+1)
    ld[,1] <- bkt/phi - lgamma(y/phi+1)
    ld[,2] <- (-bkt - y + dig*y)/(phi*phi)
    ld[,3] <- (2*bkt + 3*y - 2*dig*y - trig *y*y/phi)/(phi^3)
    return(ld) 
  }

  ## .. otherwise need the full series thing....
  ## first deal with the zeros  
  
  ind <- y==0;ld[ind,] <- 0
  ind <- ind & mu>0 ## need mu condition otherwise may try to find log(0)
 
  ld[ind,1] <- -mu[ind]^(2-p)/(phi*(2-p))
  ld[ind,2] <- -ld[ind,1]/phi  ## dld/d phi 
  ld[ind,3] <- -2*ld[ind,2]/phi ## d2ld/dphi2
  ld[ind,4] <- -ld[ind,1] * (log(mu[ind]) - 1/(2-p)) ## dld/dp
  ld[ind,5] <- 2*ld[ind,4]/(2-p) + ld[ind,1]*log(mu[ind])^2 ## d2ld/dp2
  ld[ind,6] <- -ld[ind,4]/phi ## d2ld/dphidp

  if (sum(!ind)==0) return(ld)

  ## now the non-zeros
  ind <- y==0
  y <- y[!ind];mu <- mu[!ind]
  w <- w1 <- w2 <- y*0
  oo <- .C(C_tweedious,w=as.double(w),w1=as.double(w1),w2=as.double(w2),w1p=as.double(y*0),w2p=as.double(y*0),
           w2pp=as.double(y*0),y=as.double(y),eps=as.double(.Machine$double.eps^2),n=as.integer(length(y)),
           th=as.double(theta),rho=as.double(rho),a=as.double(a),b=as.double(b))
  
  if (!work.param) { ## transform working param derivatives to p/phi derivs...
    oo$w2 <- oo$w2/phi^2 - oo$w1/phi^2
    oo$w1 <- oo$w1/phi
    oo$w2p <- oo$w2p*dthp1^2 + dthp2 * oo$w1p
    oo$w1p <- oo$w1p*dthp1
    oo$w2pp <- oo$w2pp*dthp1/phi ## this appears to be wrong
  }


  log.mu <- log(mu)
  mu1p <- theta <- mu^(1-p)
  k.theta <- mu*theta/(2-p) ## mu^(2-p)/(2-p)
  theta <- theta/(1-p) ## mu^(1-p)/(1-p)
  l.base <-  mu1p*(y/(1-p)-mu/(2-p))/phi
  ld[!ind,1] <- l.base - log(y) ## log density
  ld[!ind,2] <- -l.base/phi  ## d log f / dphi
  ld[!ind,3] <- 2*l.base/(phi*phi)  ## d2 logf / dphi2
  x <- theta*y*(1/(1-p) - log.mu)/phi + k.theta*(log.mu-1/(2-p))/phi
  ld[!ind,4] <- x
  ld[!ind,5] <- theta * y * (log.mu^2 - 2*log.mu/(1-p) + 2/(1-p)^2)/phi -
                k.theta * (log.mu^2 - 2*log.mu/(2-p) + 2/(2-p)^2)/phi ## d2 logf / dp2
  ld[!ind,6] <- - x/phi ## d2 logf / dphi dp

  if (work.param) { ## transform derivs to derivs wrt working
    ld[,3] <- ld[,3]*phi^2 + ld[,2]*phi
    ld[,2] <- ld[,2]*phi
    ld[,5] <- ld[,5]*dpth1^2 + ld[,4]*dpth2
    ld[,4] <- ld[,4]*dpth1
    ld[,6] <- ld[,6]*dpth1*phi
  }

if (TRUE) { ## DEBUG disconnetion of a terms
  ld[!ind,1] <- ld[!ind,1] + oo$w ## log density
  ld[!ind,2] <- ld[!ind,2] + oo$w1   ## d log f / dphi
  ld[!ind,3] <- ld[!ind,3] + oo$w2 ## d2 logf / dphi2
  ld[!ind,4] <- ld[!ind,4] + oo$w1p 
  ld[!ind,5] <- ld[!ind,5] + oo$w2p  ## d2 logf / dp2
  ld[!ind,6] <- ld[!ind,6] + oo$w2pp ## d2 logf / dphi dp
} 

if (FALSE) { ## DEBUG disconnetion of density terms
  ld[!ind,1] <-  oo$w ## log density
  ld[!ind,2] <-  oo$w1   ## d log f / dphi
  ld[!ind,3] <-  oo$w2 ## d2 logf / dphi2
  ld[!ind,4] <-  oo$w1p 
  ld[!ind,5] <-  oo$w2p  ## d2 logf / dp2
  ld[!ind,6] <-  oo$w2pp ## d2 logf / dphi dp
} 

  ld
} ## ldTweedie0



ldTweedie <- function(y,mu=y,p=1.5,phi=1,rho=NA,theta=NA,a=1.001,b=1.999,all.derivs=FALSE) {
## evaluates log Tweedie density for 1<=p<=2, using series summation of
## Dunn & Smyth (2005) Statistics and Computing 15:267-280.
  n <- length(y)
  if (all(!is.na(rho))&&all(!is.na(theta))) { ## use rho and theta and get derivs w.r.t. these
    #if (length(rho)>1||length(theta)>1) stop("only scalar `rho' and `theta' allowed.")
    if (a>=b||a<=1||b>=2) stop("1<a<b<2 (strict) required")
    work.param <- TRUE
    ## should buffered code for fixed p and phi be used?
    buffer <- if (length(unique(theta))==1&&length(unique(rho))==1) TRUE else FALSE 
    theta <- th <- array(theta,dim=n);
    phi <- exp(rho)
    ind <- th > 0;dpth1 <- dpth2 <-p <- rep(0,n)
    ethi <- exp(-th[ind])
    ethni <- exp(th[!ind])
    p[ind] <- (b+a*ethi)/(1+ethi)
    p[!ind] <- (b*ethni+a)/(ethni+1)
    dpth1[ind] <- ethi*(b-a)/(1+ethi)^2
    dpth1[!ind] <- ethni*(b-a)/(ethni+1)^2
    dpth2[ind] <-((a-b)*ethi+(b-a)*ethi^2)/(ethi+1)^3
    dpth2[!ind] <- ((a-b)*ethni^2+(b-a)*ethni)/(ethni+1)^3
    #p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1) 
    #dpth1 <- if (th>0) exp(-th)*(b-a)/(1+exp(-th))^2 else exp(th)*(b-a)/(exp(th)+1)^2
    #dpth2 <- if (th>0) ((a-b)*exp(-th)+(b-a)*exp(-2*th))/(exp(-th)+1)^3 else
    #               ((a-b)*exp(2*th)+(b-a)*exp(th))/(exp(th)+1)^3
  } else { ## still need working params for tweedious call...
    work.param <- FALSE
    if (all.derivs) warning("all.derivs only available in rho, theta parameterization")
    #if (length(p)>1||length(phi)>1) stop("only scalar `p' and `phi' allowed.")
    buffer <- if (length(unique(p))==1&&length(unique(phi))==1) TRUE else FALSE 
    rho <- log(phi)
    if (min(p)>=1&&max(p)<=2) {
      ind <- p>1&p<2
      if (sum(ind)) {
        p.ind <- p[ind]
        if (min(p.ind) <= a) a <- (1+min(p.ind))/2
        if (max(p.ind) >= b) b <- (2+max(p.ind))/2
        pabp <- theta <- dthp1 <- dthp2 <- rep(0,n)
        pabp[ind] <- (p.ind-a)/(b-p.ind)
        theta[ind] <- log((p.ind-a)/(b-p.ind))
        dthp1[ind] <- (1+pabp[ind])/(p.ind-a)
        dthp2[ind] <- (pabp[ind]+1)/((p.ind-a)*(b-p.ind)) -(pabp[ind]+1)/(p.ind-a)^2
      }
    }
  }

  if (min(p)<1||max(p)>2) stop("p must be in [1,2]")
  ld <- cbind(y,y,y);ld <- cbind(ld,ld*NA)
  if (work.param&&all.derivs) ld <- cbind(ld,ld[,1:3]*0,y*0)
  if (length(p)!=n) p <- array(p,dim=n);
  if (length(phi)!=n) phi <- array(phi,dim=n)
  if (length(mu)!=n) mu <- array(mu,dim=n)
  ind <- p == 2
  if (sum(ind)) { ## It's Gamma
    if (sum(y[ind]<=0)) stop("y must be strictly positive for a Gamma density")
    ld[ind,1] <- dgamma(y[ind], shape = 1/phi[ind],rate = 1/(phi[ind] * mu[ind]),log=TRUE)
    ld[ind,2] <- (digamma(1/phi[ind]) + log(phi[ind]) - 1 + y[ind]/mu[ind] - log(y[ind]/mu[ind]))/(phi[ind]*phi[ind])
    ld[ind,3] <- -2*ld[ind,2]/phi[ind] + (1-trigamma(1/phi[ind])/phi[ind])/(phi[ind]^3)
    #return(ld)
  }  

  ind <- p == 1
  if (sum(ind)) { ## It's Poisson like
    ## ld[,1] <- dpois(x = y/phi, lambda = mu/phi,log=TRUE)
    if (all.equal(y[ind]/phi[ind],round(y[ind]/phi[ind]))!=TRUE) stop("y must be an integer multiple of phi for Tweedie(p=1)")
    indi <- (y[ind]!=0)|(mu[ind]!=0) ## take care to deal with y log(mu) when y=mu=0
    bkt <- y[ind]*0
    bkt[indi] <- ((y[ind])[indi]*log((mu[ind]/phi[ind])[indi]) - (mu[ind])[indi])
    dig <- digamma(y[ind]/phi[ind]+1)
    trig <- trigamma(y[ind]/phi[ind]+1)
    ld[ind,1] <- bkt/phi[ind] - lgamma(y[ind]/phi[ind]+1)
    ld[ind,2] <- (-bkt - y[ind] + dig[ind]*y[ind])/(phi[ind]^2)
    ld[ind,3] <- (2*bkt + 3*y[ind] - 2*dig*y[ind] - trig * y[ind]^2/phi[ind])/(phi[ind]^3)
    #return(ld) 
  }

  ## .. otherwise need the full series thing....
  ## first deal with the zeros  
  
  ind <- y==0&p>1&p<2;ld[ind,] <- 0
  ind <- ind & mu>0 ## need mu condition otherwise may try to find log(0)
  if (sum(ind)) {
    mu.ind <- mu[ind];p.ind <- p[ind];phii <- phi[ind]
    ld[ind,1] <- -mu.ind^(2-p.ind)/(phii*(2-p.ind))
    ld[ind,2] <- -ld[ind,1]/phii  ## dld/d phi 
    ld[ind,3] <- -2*ld[ind,2]/phii ## d2ld/dphi2
    ld[ind,4] <- -ld[ind,1] * (log(mu.ind) - 1/(2-p.ind)) ## dld/dp
    ld[ind,5] <- 2*ld[ind,4]/(2-p.ind) + ld[ind,1]*log(mu.ind)^2 ## d2ld/dp2
    ld[ind,6] <- -ld[ind,4]/phii ## d2ld/dphidp
    if (work.param&&all.derivs) {
      mup <- mu.ind^p.ind
      ld[ind,7] <- -mu.ind/(mup*phii)
      ld[ind,8] <- -(1-p.ind)/(mup*phii)
      ld[ind,9] <- log(mu.ind)*mu.ind/(mup*phii)
      ld[ind,10] <- -ld[ind,7]/phii
    }
  }
  if (sum(!ind)==0) return(ld)
 
  ## now the non-zeros
  ind <- which(y>0&p>1&p<2)
  y <- y[ind];mu <- mu[ind];p<- p[ind]
  w <- w1 <- w2 <- y*0
  if (length(ind)>0) {
    if (buffer) { ## use code that can buffer expensive lgamma,digamma and trigamma evaluations...
      oo <- .C(C_tweedious,w=as.double(w),w1=as.double(w1),w2=as.double(w2),w1p=as.double(y*0),w2p=as.double(y*0),
               w2pp=as.double(y*0),y=as.double(y),eps=as.double(.Machine$double.eps^2),n=as.integer(length(y)),
               th=as.double(theta[1]),rho=as.double(rho[1]),a=as.double(a),b=as.double(b))
    } else { ## use code that is not able to buffer as p and phi variable...
      if (length(theta)!=n) theta <- array(theta,dim=n)
      if (length(rho)!=n) rho <- array(rho,dim=n)
      oo <- .C(C_tweedious2,w=as.double(w),w1=as.double(w1),w2=as.double(w2),w1p=as.double(y*0),w2p=as.double(y*0),
           w2pp=as.double(y*0),y=as.double(y),eps=as.double(.Machine$double.eps^2),n=as.integer(length(y)),
           th=as.double(theta[ind]),rho=as.double(rho[ind]),a=as.double(a),b=as.double(b))
    }
    if (oo$eps < -.5) {
      if (oo$eps < -1.5) { ## failure of series in C code
        oo$w2 <- oo$w1 <- oo$w2p <- oo$w1p <- oo$w2pp <- rep(NA,length(y)) 
      }
      else warning("Tweedie density may be unreliable - series not fully converged")
    }
    phii <- phi[ind]
    if (!work.param) { ## transform working param derivatives to p/phi derivs...
      if (length(dthp1)!=n) dthp1 <- array(dthp1,dim=n)
      if (length(dthp2)!=n) dthp2 <- array(dthp2,dim=n)
      dthp1i <- dthp1[ind]
      oo$w2 <- oo$w2/phii^2 - oo$w1/phii^2
      oo$w1 <- oo$w1/phii
      oo$w2p <- oo$w2p*dthp1i^2 + dthp2[ind] * oo$w1p
      oo$w1p <- oo$w1p*dthp1i
      oo$w2pp <- oo$w2pp*dthp1i/phii 
    }


    log.mu <- log(mu)
    onep <- 1-p
    twop <- 2-p
    mu1p <- theta <- mu^onep
    k.theta <- mu*theta/twop ## mu^(2-p)/(2-p)
    theta <- theta/onep ## mu^(1-p)/(1-p)
    a1 <- (y/onep-mu/twop)
    l.base <-  mu1p*a1/phii
    ld[ind,1] <- l.base - log(y) ## log density
    ld[ind,2] <- -l.base/phii  ## d log f / dphi
    ld[ind,3] <- 2*l.base/(phii^2)  ## d2 logf / dphi2
    x <- theta*y*(1/onep - log.mu)/phii + k.theta*(log.mu-1/twop)/phii
    ld[ind,4] <- x
    ld[ind,5] <- theta * y * (log.mu^2 - 2*log.mu/onep + 2/onep^2)/phii -
                  k.theta * (log.mu^2 - 2*log.mu/twop + 2/twop^2)/phii ## d2 logf / dp2
    ld[ind,6] <- - x/phii ## d2 logf / dphi dp
  } ## length(ind)>0

  if (work.param) { ## transform derivs to derivs wrt working
    ld[,3] <- ld[,3]*phi^2 + ld[,2]*phi
    ld[,2] <- ld[,2]*phi
    ld[,5] <- ld[,5]*dpth1^2 + ld[,4]*dpth2
    ld[,4] <- ld[,4]*dpth1
    ld[,6] <- ld[,6]*dpth1*phi
    colnames(ld)[1:6] <- c("l","rho","rho.2","th","th.2","th.rho")
  }

  if (work.param&&all.derivs&&length(ind)>0) {
    #ld <- cbind(ld,ld[,1:4]*0)
    a2 <- mu1p/(mu*phii) ## 1/(mu^p*phii)
    ld[ind,7] <- a2*(onep*a1-mu/twop)   ## deriv w.r.t mu
    ld[ind,8] <- -a2*(onep*p*a1/mu+2*onep/twop) ## 2nd deriv w.r.t. mu
    ld[ind,9] <- a2*(-log.mu*onep*a1-a1 + onep*(y/onep^2-mu/twop^2)+mu*log.mu/twop-mu/twop^2) ## mu p
    ld[ind,10] <- a2*(mu/(phii*twop) - onep*a1/phii) ## mu phi
    ## transform to working...
    ld[,10] <- ld[,10]*phi
    ld[,9] <- ld[,9]*dpth1
    colnames(ld) <- c("l","rho","rho.2","th","th.2","th.rho","mu","mu.2","mu.theta","mu.rho")
  }


  if (length(ind)>0) { 
    ld[ind,1] <- ld[ind,1] + oo$w ## log density
    ld[ind,2] <- ld[ind,2] + oo$w1   ## d log f / dphi
    ld[ind,3] <- ld[ind,3] + oo$w2 ## d2 logf / dphi2
    ld[ind,4] <- ld[ind,4] + oo$w1p 
    ld[ind,5] <- ld[ind,5] + oo$w2p  ## d2 logf / dp2
    ld[ind,6] <- ld[ind,6] + oo$w2pp ## d2 logf / dphi dp
  } 

if (FALSE) { ## DEBUG disconnection of density terms
  ld[ind,1] <-  oo$w ## log density
  ld[ind,2] <-  oo$w1   ## d log f / dphi
  ld[ind,3] <-  oo$w2 ## d2 logf / dphi2
  ld[ind,4] <-  oo$w1p 
  ld[ind,5] <-  oo$w2p  ## d2 logf / dp2
  ld[ind,6] <-  oo$w2pp ## d2 logf / dphi dp
} 

  ld
} ## ldTweedie


Tweedie <- function(p=1,link=power(0)) {
## a restricted Tweedie family
  if (p<=1||p>2) stop("Only 1<p<=2 supported")
  
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt","inverse")
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp) else 
  if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
       if (!is.null(stats$name))
          linktemp <- stats$name
        } else {
            stop(gettextf("link \"%s\" not available for Tweedie family.",
                linktemp, collapse = ""),domain = NA)
        }
    }
    
    variance <- function(mu) mu^p
    dvar <- function(mu) p*mu^(p-1)
    if (p==1) d2var <- function(mu) 0*mu else
      d2var <- function(mu) p*(p-1)*mu^(p-2)
    if (p==1||p==2)  d3var <- function(mu) 0*mu else
      d3var <- function(mu) p*(p-1)*(p-2)*mu^(p-3)
    validmu <- function(mu) all(mu >= 0)

    dev.resids <- function(y, mu, wt) {
        y1 <- y + (y == 0)
        if (p == 1)
            theta <- log(y1/mu)
        else theta <- (y1^(1 - p) - mu^(1 - p))/(1 - p)
        if (p == 2)
            kappa <- log(y1/mu)
        else kappa <- (y^(2 - p) - mu^(2 - p))/(2 - p)
        pmax(2 * wt * (y * theta - kappa),0)
    }
    initialize <- expression({
        n <- rep(1, nobs)
        mustart <- y + 0.1 * (y == 0)
    })

    ls <-  function(y,w,n,scale) {
      power <- p
      colSums(w*ldTweedie(y,y,p=power,phi=scale))
    }

    aic <- function(y, n, mu, wt, dev) {
      power <- p
      scale <- dev/sum(wt)
      -2*sum(ldTweedie(y,mu,p=power,phi=scale)[,1]*wt) + 2
    }

    if (p==2) {
      rd <- function(mu,wt,scale) {
        rgamma(mu,shape=1/scale,scale=mu*scale)
      }   
    } else {
      rd <- function(mu,wt,scale) {
        rTweedie(mu,p=p,phi=scale)
      }
    }

    structure(list(family = paste("Tweedie(",p,")",sep=""), variance = variance, 
              dev.resids = dev.resids,aic = aic, link = linktemp, linkfun = stats$linkfun, linkinv = stats$linkinv,
        mu.eta = stats$mu.eta, initialize = initialize, validmu = validmu,
        valideta = stats$valideta,dvar=dvar,d2var=d2var,d3var=d3var,ls=ls,rd=rd,canonical="none"), class = "family")


} ## Tweedie



rTweedie <- function(mu,p=1.5,phi=1) {
## generate Tweedie random variables, with 1<p<2, 
## adapted from rtweedie in the tweedie package
  if (p<=1||p>=2) stop("p must be in (1,2)")
  if (sum(mu<0)) stop("mean, mu, must be non negative")
  if (phi<=0) stop("scale parameter must be positive")
  
  lambda <- mu^(2-p)/((2-p)*phi)
  shape <- (2-p)/(p-1)
  scale <- phi*(p-1)*mu^(p-1)

  n.sim <- length(mu)

  ## how many Gamma r.v.s to sum up to get Tweedie
  ## 0 => none, and a zero value

  N <- rpois(length(lambda),lambda)

  ## following is a vector of N[i] copies of each gamma.scale[i]
  ## concatonated one after the other

  gs <- rep(scale,N)

  ## simulate gamma deviates to sum to get tweedie deviates

  y <- rgamma(gs*0+1,shape=shape,scale=gs)

  ## create summation index...

  lab <- rep(1:length(N),N)

  ## sum up each gamma sharing a label. 0 deviate if label does not occur
  o <- .C(C_psum,y=as.double(rep(0,n.sim)),as.double(y),as.integer(lab),as.integer(length(lab)))  
  o$y
} ## rTweedie

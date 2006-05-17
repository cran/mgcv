## R routines for gam fitting with calculation of derivatives w.r.t. sp.s
## (c) Simon Wood 2004,2005,2006

gam.fit2 <- function (x, y, sp, S=list(),rS=list(),off, H=NULL, weights =
rep(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
    control = gam.control(), intercept =
    TRUE,deriv=TRUE,gamma=1,scale=1,pearson=FALSE,
    printWarn=TRUE) 
## deriv, sp, S, H added to arg list. 
## need to modify family before call.
{   x <- as.matrix(x)
    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
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
    z1 <- w1 <- matrix(0,nobs,nSp)
    beta1old <- matrix(0,ncol(x),nSp)
    upe <- list(beta=rep(0,ncol(x)),beta1=beta1old,trA=0,trA1=rep(0,nSp))
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
        if (pearson) {
          alpha1 <- alpha <- sum((y-mu)^2/V)
        } else {
          alpha1 <- dev
        }
        trA1 <- trA <- 0
        if (deriv) GCV1<-UBRE1<-trA1 <- alpha1 <- rep(0,nSp)
        else GCV1<-UBRE1<-trA1 <- alpha1 <- NULL
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
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("Can't find valid starting values: please specify some")
        devold <- sum(dev.resids(y, mu, weights))
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

            if (deriv&&iter>1) ## then get derivatives of z and w w.r.t. theta
            { d2g <- family$d2link(mug)
              dV <- family$dvar(mug)
              eta1 <- (x%*%upe$beta1)[good,]
              z1 <- as.vector((yg-mug)*d2g*mevg)*eta1
              w1 <- as.vector(-0.5*w^3/weg*(dV/mevg + 2*var.mug*d2g))*eta1
            }
            ngoodobs <- as.integer(nobs - sum(!good))
            ## Here a Fortran call has been replaced by update.beta call
           
            if (sum(good)<ncol(x)) stop("Not enough informative observations.")

            oo<-.C(C_update_beta,as.double(x[good,]),as.double(Sr),as.double(unlist(rS)),as.double(sp),
                   as.double(w),
                   as.double(w1),as.double(z),as.double(z1),as.integer(ncol(Sr)),
                   rSncol=as.integer(unlist(lapply(rS,ncol))),m=as.integer(length(rS)),
                   n=as.integer(sum(good)),
                   q=as.integer(ncol(x)),get.trA=as.integer(0),as.integer(deriv),
                   rank.tol= as.double(.Machine$double.eps),
                   beta=as.double(upe$beta),trA=as.double(upe$trA),beta1=as.double(upe$beta1),
                   trA1=as.double(upe$trA1),rV=as.double(rV),rank=as.integer(1))
        
            upe$beta <- oo$beta;
            upe$beta1 <- matrix(oo$beta1,oo$q,oo$m)

            if (any(!is.finite(upe$beta))) {
                conv <- FALSE
                warning("Non-finite coefficients at iteration ", 
                  iter)
                break
            }

           start <- upe$beta 
           eta <- drop(x%*%start)
     
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
                  # ... modify derivatives similarly ...
                  if (deriv) upe$beta1 <- (upe$beta1 + beta1old)/2
                 
                  eta <- drop(x %*% start)
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
                  # ... modify derivatives similarly ...
                  if (deriv) upe$beta1 <- (upe$beta1 + beta1old)/2

                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
             pdev <- dev + t(start)%*%St%*%start ## the penalized deviance 

            if (iter>1&&pdev>old.pdev) { ## solution diverging
              ii <- 1
             # while (pdev - old.pdev> (0.1 + abs(old.pdev))* control$epsilon*.9)
            while (pdev -old.pdev > (.1+abs(old.pdev))*.2)  
            {
                if (ii > 200) 
                   stop("inner loop 3; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                # ... modify derivatives similarly ...
                if (deriv) upe$beta1 <- (upe$beta1 + beta1old)/2
                 
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                  pdev <- dev + t(start)%*%St%*%start ## the penalized deviance
            
               }
            } 

            if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {  old.pdev <- pdev
                devold <- dev
                coef <- coefold <- start
                if (deriv) beta1old <- upe$beta1
            }
        } ### end main loop 
        ## Now do a final update.eta call to get trA and rV (rV%*%t(rV) is the
        ## posterior covariance matrix)....
       
        oo<-.C(C_update_beta,as.double(x[good,]),as.double(Sr),as.double(unlist(rS)),as.double(sp),
                   as.double(w),
                   as.double(w1),as.double(z),as.double(z1),as.integer(ncol(Sr)),
                   rSncol=as.integer(unlist(lapply(rS,ncol))),m=as.integer(length(rS)),
                   n=as.integer(sum(good)),
                   q=as.integer(ncol(x)),get.trA=as.integer(1),as.integer(deriv),
                   rank.tol= as.double(.Machine$double.eps),
                   beta=as.double(upe$beta),trA=as.double(upe$trA),beta1=as.double(upe$beta1),
                   trA1=as.double(upe$trA1),rV=as.double(rV),rank=as.integer(1))

        rV <- matrix(oo$rV,ncol(x),ncol(x))
        upe$beta <- oo$beta;
        upe$trA <- oo$trA;
        upe$beta1 <- matrix(oo$beta1,oo$q,oo$m)
        upe$trA1 <- oo$trA1

        V <- variance(mug)

        ### make pearson/deviance dependent
        
        if (pearson) alpha1 <- alpha <- sum(weights[good]*(yg-mug)^2/V)
        else { # devaince based GCV/UBRE
          dev <- sum(dev.resids(y, mu, weights))
          alpha1 <- alpha <- dev 
        } 
        ####

        trA1 <- trA <- upe$trA
            
        GCV <- nobs*alpha/(nobs-gamma*trA)^2        
        UBRE <- alpha/nobs + 2*gamma*trA*scale/nobs - scale
        scale.est <- alpha/(length(mug)-trA)
        if (deriv) { # need to evaluate score component derivatives
          mu.eta.val <- mu.eta(eta[good])
          if (pearson) {
            d2g <- family$d2link(mug)
            dV <- family$dvar(mug)
            eta1 <- (x%*%upe$beta1)[good,]
            temp <- as.vector(-dV/V^2*mu.eta.val*(yg-mug)^2-2/V*(yg-mug)*mu.eta.val)*eta1
            temp <- as.matrix(temp*as.vector(weights[good]))
            alpha1 <- colSums(temp) # deriv of alpha w.r.t. s.p.s
          } else { ## deviance based GCV/UBRE scores
            temp <- weights[good]*(yg-mug)*mu.eta.val/V
            temp <- t(temp)%*%x[good,] ## dl / d beta_j (unscaled)
            alpha1 <- -2 * as.numeric(temp%*% upe$beta1)
          }
          trA1 <- upe$trA1
          GCV1 <- nobs*alpha1/(nobs-gamma*trA)^2 + 
                  2*nobs*alpha*gamma*trA1/(nobs-gamma*trA)^3
          UBRE1 <- alpha1/nobs + 2*gamma*scale/nobs*trA1 
        } else UBRE1<-GCV1<-NULL
  
        # end of inserted code
        if (!conv&&printWarn) 
            warning("Algorithm did not converge")
        if (printWarn&&boundary) 
            warning("Algorithm stopped at boundary value")
        eps <- 10 * .Machine$double.eps
        if (printWarn&&family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("fitted probabilities numerically 0 or 1 occurred")
        }
        if (printWarn&&family$family == "poisson") {
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
        boundary = boundary,alpha=alpha,alpha1=alpha1,trA=trA,trA1=trA1,
        GCV=GCV,GCV1=GCV1,UBRE=UBRE,UBRE1=UBRE1,rV=rV,scale.est=scale.est,aic=aic.model,rank=oo$rank)
}



gam2derivative <- function(lsp,args)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the derivatives of the GCV or UBRE score w.r.t the 
## smoothing parameters for the model.
## args is a list containing the arguments for gam.fit2
## For use as optim() objective gradient
{ b<-gam.fit2(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=TRUE,
     control=args$control,gamma=args$gamma,scale=args$scale,pearson=args$pearson,
     printWarn=FALSE)
  if (args$scoreType == "GCV") ret <- b$GCV1 else ret <- b$UBRE1
  ret
}


gam2objective <- function(lsp,args,printWarn=FALSE)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the GCV or UBRE score for the model.
## args is a list containing the arguments for gam.fit2
## For use as optim() objective
{ 
  b<-gam.fit2(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=FALSE,
     control=args$control,gamma=args$gamma,scale=args$scale,pearson=args$pearson,
     printWarn=printWarn)
  if (args$scoreType == "GCV") ret <- b$GCV else ret <- b$UBRE
  attr(ret,"full.fit") <- b
  ret
}

gam3objective <- function(lsp,args)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the GCV or UBRE score for the model.
## args is a list containing the arguments for gam.fit2
## For use as nlm() objective
{ 
  b<-gam.fit2(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=TRUE,
     control=args$control,gamma=args$gamma,scale=args$scale,pearson=args$pearson,
     printWarn=FALSE)
  if (args$scoreType == "GCV") ret <- b$GCV else ret <- b$UBRE
  attr(ret,"full.fit") <- b
  if (args$scoreType == "GCV") at <- b$GCV1 else at <- b$UBRE1
  attr(ret,"gradient") <- at
  ret
}



fix.family.link<-function(fam)
# adds d2link the second derivative of the link function w.r.t. mu
# to the family supplied, as well as a 3rd derivative function 
# d3link...
# All d2link and d3link functions have been checked numerically. 
{ if (!inherits(fam,"family")) stop("fam not a family object")
  if (!is.null(fam$d2link)) return(fam) 
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
  if (!is.null(fam$dvar)) return(fam) 
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



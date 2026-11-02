## shape constrained additive modelling code 
## (c) Simon N. Wood 2025

ucon <- function(A,b) {
## finds unique constraints for a constraint set A x = b or Ax >= b.
  n <- length(b)
  A <- cbind(A,b)
  p <- ncol(A)
  up <- rep(TRUE,n) ## not processed yet
  ## scale constraint so that first non-zero element is 1 - this will mean that
  ## constraints identical to a multiplicative constant are correctly trweated as duplicates.
  for (i in 1:p) {
    ii <- A[,i]!=0 & up
    A[ii,] <- A[ii,]/abs(A[ii,i])
    up[ii] <- FALSE
    if (sum(up)==0) break
  }
  A <- unique(A)
  list(A=A[,1:(p-1),drop=FALSE],b = A[,p])
} ## ucon

scasm <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,offset=NULL,
                  control=list(),scale=0,select=FALSE,knots=NULL,sp=NULL,gamma=1,fit=TRUE,
                  G=NULL,drop.unused.levels=TRUE,drop.intercept=NULL,bs=0,...) {
## Shape Constrained Additive Smooth Model function. Modified from 'gam'. Smoothing parameter selection by EFS,
## constrained coefficient estimation by quadratic programming using 'pcls'.

   control <- do.call("gam.control",control)

   if (is.null(G)) {
    ## create model frame..... 
    gp <- interpret.gam(formula) # interpret the formula 
    cl <- match.call() # call needed in gam object for update to work
    mf <- match.call(expand.dots=FALSE)
    mf$formula <- gp$fake.formula 
    mf$family <- mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$select <- mf$drop.intercept  <-
                 mf$gamma<-mf$fit<-mf$G <- mf$bs <- mf$... <-NULL
    mf$drop.unused.levels <- drop.unused.levels
    mf[[1]] <- quote(stats::model.frame) ## as.name("model.frame")
    pmf <- mf
    mf <- eval(mf, parent.frame()) # the model frame now contains all the data 
    if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf,"terms")

    ## summarize the *raw* input variables
    ## note can't use get_all_vars here -- buggy with matrices
    vars <- all_vars1(gp$fake.formula[-2]) ## drop response here
    inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

    ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
    if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data) 

    dl <- eval(inp, data, parent.frame())
    names(dl) <- vars ## list of all variables needed
    var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
    rm(dl) ## save space    

    ## pterms are terms objects for the parametric model components used in 
    ## model setup - don't try obtaining by evaluating pf in mf - doesn't
    ## work in general (e.g. with offset)...

    if (is.list(formula)) { ## then there are several linear predictors
      environment(formula) <- environment(formula[[1]]) ## e.g. termplots needs this
      pterms <- list()
      tlab <- rep("",0)
      for (i in 1:length(formula)) {
        pmf$formula <- gp[[i]]$pf 
        pterms[[i]] <- attr(eval(pmf, parent.frame()),"terms")
        tlabi <- attr(pterms[[i]],"term.labels")
        if (i>1&&length(tlabi)>0) tlabi <- paste(tlabi,i-1,sep=".")
        tlab <- c(tlab,tlabi)
      }
      attr(pterms,"term.labels") <- tlab ## labels for all parametric terms, distinguished by predictor
    } else { ## single linear predictor case
      pmf$formula <- gp$pf
      pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
      pterms <- attr(pmf,"terms") ## pmf only used for this
    }

    if (is.character(family)) family <- eval(parse(text=family))
    if (is.function(family)) family <- family()
    if (is.null(family$family)) stop("family not recognized")
  
    if (family$family[1]=="gaussian" && family$link=="identity") am <- TRUE
    else am <- FALSE
    
    if (!control$keepData) rm(data) ## save space

    ## check whether family requires intercept to be dropped...
    if (is.null(family$drop.intercept)) { ## family does not provide information
      lengthf <- if (is.list(formula)) length(formula) else 1
      if (is.null(drop.intercept)) drop.intercept <- rep(FALSE, lengthf) else {
        drop.intercept <- rep(drop.intercept,length=lengthf) ## force drop.intercept to correct length
	if (sum(drop.intercept)) family$drop.intercept <- drop.intercept ## ensure prediction works
      }
    } else drop.intercept <- as.logical(family$drop.intercept) ## family overrides argument
    
    if (inherits(family,"general.family")&&!is.null(family$presetup)) eval(family$presetup)

    gsname <- if (is.list(formula)) gam.setup.list else gam.setup 

    G <- do.call(gsname,list(formula=gp,pterms=pterms,
                 data=mf,knots=knots,sp=sp,min.sp=NULL,
                 H=NULL,absorb.cons=FALSE,sparse.cons=0,select=select,
                 idLinksBases=control$idLinksBases,scale.penalty=control$scalePenalty,
                 paraPen=NULL,drop.intercept=drop.intercept))

    m <- length(G$smooth); p <- ncol(G$X)

    ## Set up constraints and feasible initial params. pcls wants initial to not exactly
    ## satisfy inequality constraints...
    Ain <- matrix(0,0,p); bi <- b <- rep(0,0)
    b0 <- numeric(p)
    for (i in 1:m) {
      ii <- G$smooth[[i]]$first.para: G$smooth[[i]]$last.para
      if (!is.null(G$C)) {
        k <- which(rowSums(G$C[,ii,drop=FALSE]!=0)!=0)
        if (length(k)&&any(G$C[k,-ii]!=0)) stop("fixed constraints span multiple smooths")
        C0 <- if (length(k)) G$C[k,ii,drop=FALSE] else matrix(0,0,length(ii))
        d0 <- if (length(k) && !is.null(G$d)) G$d[k] else numeric(0)
      } else {
        C0 <- matrix(0,0,length(ii))
	d0 <- numeric(0)
      }	
      if (!is.null(G$smooth[[i]]$Ain)) {
        con <- ucon(G$smooth[[i]]$Ain,G$smooth[[i]]$bin) ## strip any duplicate constraints (this matters!)
	G$smooth[[i]]$Ain <- con$A; G$smooth[[i]]$bin <- con$b 
        A0 <- matrix(0,nrow(G$smooth[[i]]$Ain),p)
        A0[,ii] <- G$smooth[[i]]$Ain
        Ain <- rbind(Ain,A0)
        b <- c(b,G$smooth[[i]]$bin)
      }
      if (!all(d0==0)||!is.null(G$smooth[[i]]$Ain)) {
        b0[ii] <- feasible(G$smooth[[i]]$Ain,G$smooth[[i]]$bin + .Machine$double.eps^.25,C0,d0)
	if (any(G$smooth[[i]]$Ain %*% b0[ii] <= G$smooth[[i]]$bin)) stop("initial point fails inequality constraints")
	if (nrow(C0)&&any(abs(C0 %*% b0[ii] - d0)>.Machine$double.eps^.85*max(abs(b0[ii]))*norm(C0,"M"))) stop(
	    "initial point fails equality constraints")
      }
    }
    G$Ain <- Ain; G$bin <- b; # G$C <- C;
    G$beta0 <- b0
    ## ...so beta0 satisfies Ain beta0 > bin and C beta0 = d (here zero)
    ## point is that starting from beta_0, pcls can keep constraints satisfied.
    
    G$var.summary <- var.summary
    G$family <- fix.family(family)
   
    if ((is.list(formula)&&(is.null(family$nlp)||family$nlp!=gp$nlp))||
        (!is.list(formula)&&!is.null(family$npl)&&(family$npl>1)))
	stop("incorrect number of linear predictors for family")
       
    G$terms<-terms;
    G$mf<-mf;G$cl<-cl;
    G$am <- am

    if (is.null(G$offset)) G$offset<-rep(0,G$n)
     
    G$min.edf <- G$nsdf ## -dim(G$C)[1]
    if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim

    G$formula <- formula
    G$pred.formula <- gp$pred.formula
    environment(G$formula)<-environment(formula)
  } else { ## G not null
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
  }

  if (!fit) {
    class(G) <- "gam.prefit"
    return(G)
  }  

  if (scale>0) G$family$scale <- scale

  object <- scasm.fit(G,control,gamma=gamma,bs=bs)

  object$Ain <- G$Ain; object$bin <- G$bin; object$beta0 <- G$beta0
  if (!is.null(G$L)) { 
    object$full.sp <- as.numeric(exp(G$L%*%log(object$sp)+G$lsp0))
    names(object$full.sp) <- names(G$lsp0)
  }
  names(object$sp) <- names(G$sp)
  object$paraPen <- G$pP
  object$formula <- G$formula
  ## store any lpi attribute of G$X for use in predict.gam...
  if (is.list(object$formula)) attr(object$formula,"lpi") <- attr(G$X,"lpi")
  object$var.summary <- G$var.summary 
  object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
  object$model<-G$mf # store the model frame
  object$na.action <- attr(G$mf,"na.action") # how to deal with NA's
  object$control <- control
  object$terms <- G$terms
  object$pred.formula <- G$pred.formula

  attr(object$pred.formula,"full") <- reformulate(all.vars(object$terms))
  
  object$pterms <- G$pterms
  object$assign <- G$assign # applies only to pterms
  object$contrasts <- G$contrasts
  object$xlevels <- G$xlevels
  object$offset <- G$offset
  if (!is.null(G$Xcentre)) object$Xcentre <- G$Xcentre
  if (control$keepData) object$data <- data
  object$df.residual <- nrow(G$X) - sum(object$edf)
  object$min.edf <- G$min.edf
  object$optimizer <- "EFS"
  names(object$coefficients) <- G$term.names
  object$call <- G$cl # needed for update() to work
  class(object) <- c("gam","glm","lm")
  object$y <- G$y;object$prior.weights <- G$w
  if (is.null(object$deviance)) object$deviance <- sum(residuals(object,"deviance")^2)
  object$gcv.ubre <- -object$laml
  names(object$gcv.ubre) <- "REML"
  object$method <- "REML" 
  object$smooth <- G$smooth
  object$nsdf <- G$nsdf
  object$sig2 <- object$scale
  !object$scale.known -> object$scale.estimated
  ## The following lines avoid potentially very large objects in hidden environments being stored
  ## with fitted gam objects. The downside is that functions like 'termplot' that rely on searching in
  ## the environment of the formula can fail...
  environment(object$formula) <- environment(object$pred.formula) <-
  environment(object$terms) <- environment(object$pterms) <- .GlobalEnv
  if (!is.null(object$model))  environment(attr(object$model,"terms"))  <- .GlobalEnv
  if (!is.null(attr(object$pred.formula,"full"))) environment(attr(object$pred.formula,"full")) <- .GlobalEnv
  object
} ## scasm


smooth.construct.sc.smooth.spec <- function(object,data,knots) {
## basic shape constrained smooth based on bs basis, using the
## relationship between difference constraints on B-spline basis
## coeffs and constraints on derivative of the corresponding spline. 
## Options are c+, c-, m+, m-, + for convex, concave, increasing,
## decreasing and positive. All coherent combinations are possible,
## and Ain %*% b >= bin imposes the constraints.
## For c+ with +, Ain will have more rows than columns: all other 
## combinations there are at most as many rows as columns. 
  xt <- object$xt
  sm <- smooth.construct.bs.smooth.spec(object,data,knots)
  if (is.null(xt)) return(sm)
  p <- ncol(sm$X)
  if (is.list(xt)) { 
    if (is.character(xt[1]))  con <- xt[[1]] else con <- NULL
  } else con <- xt
  Ain <- matrix(0,0,p)
  Ip <- diag(1,nrow = p)
  ## get start and end, i.e. first and last interior knots...
  xint <- c(sm$knots[sm$m[1]+1],sm$knots[length(sm$knots)-sm$m[1]])
  dp <- data.frame(xint[1]); names(dp) <- sm$term ## used for start/end constraints
  if ("c+" %in% con) { ## increasing first derivative
    Ain <- rbind(Ain,diff(Ip,diff=2))
    if ("m+" %in% con) { ## must be monotonic increasing at start
      sm$deriv <- 1
      Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
      if ("+" %in% con) { ## ... and positive at start
        sm$deriv <- NULL
	Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
      }
    } else if ("m-" %in% con) { ## must be monotonic decreasing at end
      sm$deriv <- 1; dp[[1]] <- xint[2]
      Ain <- rbind(Ain,-Predict.matrix.Bspline.smooth(sm,dp))
      if ("+" %in% con) { ## ... and positive at end
        sm$deriv <- NULL
	Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
      }
    } else  if ("+" %in% con) { ## +ve but not monotonic
      Ain <- rbind(Ain,Ip) ## all needed as minimum location not known in advance :-(
    }
  } else if ("c-" %in% con) { ## decreasing first derivative
    Ain <- rbind(Ain,-diff(Ip,diff=2))
    if ("m+" %in% con) { ## must be monotonic increasing at end
      sm$deriv <- 1; dp[[1]] <- xint[2]
      Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
      if ("+" %in% con) { ## ... and positive at start
        sm$deriv <- NULL; dp[[1]] <- xint[1]
	Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
      }
    } else if ("m-" %in% con) { ## must be monotonic decreasing at start
      sm$deriv <- 1
      Ain <- rbind(Ain,-Predict.matrix.Bspline.smooth(sm,dp))
      if ("+" %in% con) { ## ... and positive at end
        sm$deriv <- NULL; dp[[1]] <- xint[2]
	Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
      }
    } else  if ("+" %in% con) { ## +ve but not nonotonic
      dp <- data.frame(xint); names(dp) <- sm$term ## +ve both ends
      Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
    }
  } else if ("m+" %in% con) { ## not convex/cave but increasing
    Ain <- rbind(Ain,diff(Ip))
    if ("+" %in% con) { ## positive at start
      Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
    }
  } else  if ("m-" %in% con) { ## not convex/cave but decreasing
    Ain <- rbind(Ain,-diff(Ip))
    if ("+" %in% con) { ## positive at end
      dp[[1]] <- xint[2]
      Ain <- rbind(Ain,Predict.matrix.Bspline.smooth(sm,dp))
    }
  } else if ("+" %in% con) Ain <- Ip ## simply positive
  
  if ("+" %in% con) sm$C <- matrix(0,0,p) ## +ve constraint makes no sense unless centering not needed
  
  sm$Ain <- Ain; sm$bin <- rep(0,nrow(Ain))
  sm$bin0 <- rep(0.01,nrow(Ain)) ## threshold for initializing to ensure > not =
  sm$deriv <- NULL ## make sure turned off
  sm
} ## smooth.construct.sc.smooth.spec


smooth.construct.scad.smooth.spec <- function(object,data,knots) {
## basic shape constrained smooth based on ad basis...
  xt <- object$xt
  object$xt <- list(bs="ps")
  sm <- smooth.construct.ad.smooth.spec(object,data,knots)
  if (is.null(xt)) return(sm)
  if (sm$dim>1) stop("scad is for 1D smoothing only")
  p <- ncol(sm$X)
  if (is.list(xt)) { 
    if (is.character(xt[1]))  con <- xt[[1]] else con <- NULL
  } else con <- xt
  Ain <- matrix(0,0,p)
  Ip <- diag(1,nrow = p)
  if ("+" %in% con) {
    Ain <- Ip
    sm$C <- matrix(0,0,p) ## never makes any sense to use positive constraint unless centring not needed 
  }  
  if ("m+" %in% con) Ain <- rbind(Ain,diff(Ip)) else
        if ("m-" %in% con) Ain <- rbind(Ain,-diff(Ip))
  if ("c+" %in% con) Ain <- rbind(Ain,diff(Ip,diff=2)) else
        if ("c-" %in% con) Ain <- rbind(Ain,-diff(Ip,diff=2))	
  sm$Ain <- Ain; sm$bin <- rep(0,nrow(Ain))
  sm$bin0 <- rep(0.01,nrow(Ain)) ## threshold for initializing to ensure > not =
  sm
} ## smooth.construct.scad.smooth.spec



bSbfun <- function(beta,S,off) {
  bSb <- numeric(length(S))
  for (j in 1:length(S)) {
    ii <- 1:ncol(S[[j]]) + off[j]-1
    bSb[j] <- beta[ii] %*% S[[j]] %*% beta[ii]
  }
  bSb
} ## bSbfun

scasm.pirls <- function(G,lsp,control,start=NULL,etastart=NULL,mustart=NULL,
                        nstep=control$maxit+10,eps=1e-4,bs=0) {
## penalized IRLS for constrained GAM fit. Weighted least squares is performed by
## pcls. Note that 'start' is never used as the pcls initial parameter value...
## G$beta0 is designed for that purpose.

   warn <- list()
   ## initialize as glm.fit
   weights <- G$w
   offset <- G$offset
   n <- nobs <- NROW(G$y) ## really n set in initialize
   nvars <- ncol(G$X)
   family <- G$family
   variance <- family$variance
   linkinv <- family$linkinv
   if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
            call. = FALSE)
   dev.resids <- family$dev.resids
   aic <- family$aic
   mu.eta <- family$mu.eta
   valideta <- family$valideta %||% function(eta) TRUE
   validmu <- family$validmu %||% function(mu) TRUE
   y <- G$y
   if (is.null(mustart)) {
     eval(family$initialize)
   } else {
     mukeep <- mustart
     eval(family$initialize)
     mustart <- mukeep
   }
   pdevold <- coefold <- NULL
   eta <- etastart %||% {
     if (!is.null(start)) 
        if (length(start) != nvars) 
           stop(gettextf("length of 'start' should equal %d and correspond to initial coefs", 
                nvars),domain = NA)
        else {
          coefold <- start
          offset + as.vector(if (NCOL(G$X) == 1L)  G$X * start else G$X %*% start)
        }
      else G$family$linkfun(mustart)
   }
   mu <- linkinv(eta)
   etaold <- eta
   if (!(validmu(mu) && valideta(eta))) 
       stop("cannot find valid starting values: please specify some", call. = FALSE)
   boundary <- conv <- FALSE
   gausid <- family$family == "gaussian" && family$link == "id"

   if (is.null(coefold)) {
     coefold <- G$beta0 ## might as well use feasible beta as starting point
     etaold <- drop(G$X %*% coefold)+offset
   }
   pdevold <- sum(dev.resids(y,linkinv(etaold),weights))
   bSb <- bSbfun(coefold,G$S,G$off)
   pdevold <- pdevold + sum(exp(lsp)*bSb)
   if (is.null(start)||is.null(attr(start,"active"))) start <- G$beta0
   for (iter in 1L:control$maxit) {
     good <- weights > 0
     varmu <- variance(mu)[good]
     if (anyNA(varmu)) stop("NAs in V(mu)")
     if (any(varmu == 0)) stop("0s in V(mu)")
     mu.eta.val <- mu.eta(eta)
     if (any(is.na(mu.eta.val[good]))) stop("NAs in d(mu)/d(eta)")
     good <- (weights > 0) & (mu.eta.val != 0)
     if (all(!good)) {
       conv <- FALSE
       warning(gettextf("no observations informative at iteration %d", 
                  iter), domain = NA)
       break
     }
     z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
     w <- weights[good] * mu.eta.val[good]^2/variance(mu)[good]
     start <- pcls(list(y=z,w=w,X=G$X[good,],C=G$C,S=G$S,off=G$off-1L,sp=exp(lsp),
                   p=G$beta0,Ain=G$Ain,bin = G$bin,active=0*attr(start,"active")))
     bSb <- bSbfun(start,G$S,G$off)
     eta <- drop(G$X %*% start)
     mu <- linkinv(eta <- eta + offset)
     dev <- sum(dev.resids(y, mu, weights))
     pdev <- dev + sum(exp(lsp)*bSb)
     if (control$trace) cat("pdev = ", pdev, " Iterations - ", iter,"\n", sep = "")
     boundary <- FALSE
     if (!is.finite(pdev)) {
       if (is.null(coefold)) 
           stop("no valid set of coefficients has been found: please supply starting values",call. = FALSE)
       warn[length(warn)+1] <- "step size truncated due to divergence"
       ii <- 1
       while (!is.finite(pdev)) {
         if (ii > control$maxit) stop("inner loop 1; cannot correct step size", call. = FALSE)
         ii <- ii + 1
         start <- (start + coefold)/2
         eta <- drop(G$X %*% start)
         mu <- linkinv(eta <- eta + offset)
         dev <- sum(dev.resids(y, mu, weights))
	 bSb <- bSbfun(start,G$S,G$off)
	 pdev <- dev + sum(exp(lsp)*bSb)
       }
       boundary <- TRUE
       if (control$trace) cat("Step halved: new deviance = ", dev, "\n", sep = "")
     } ## finite pdev control
     
     if (!(valideta(eta) && validmu(mu))) {
       if (is.null(coefold)) 
         stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
       warn[length(warn)+1] <- "step size truncated: out of bounds"
       ii <- 1
       while (!(valideta(eta) && validmu(mu))) {
         if (ii > control$maxit) 
           stop("inner loop 2; cannot correct step size", call. = FALSE)
           ii <- ii + 1
           start <- (start + coefold)/2
           eta <- drop(G$X %*% start)
           mu <- linkinv(eta <- eta + offset)
       }
       boundary <- TRUE
       dev <- sum(dev.resids(y, mu, weights))
       bSb <- bSbfun(start,G$S,G$off)
       pdev <- dev + sum(exp(lsp)*bSb)
       if (control$trace) cat("Step halved: new pdev = ", pdev, "\n", sep = "")
    } ## valid mu/eta control
    
   
    f <- function(x) {
      p <- coefold + x*(start-coefold)
      etax <- etaold + x*(eta-etaold)
      mu <- linkinv(etax)
      dev <- sum(dev.resids(y, mu, weights))
      bSb <- bSbfun(p,G$S,G$off)
      pdev <- dev + sum(exp(lsp)*bSb)
    }

    if (pdev-pdevold > abs(pdev)*control$epsilon) { ## golden section search for best
      alpha <- tau <- 2/(1 + sqrt(5))
      a0 <- 0; a1 <- 1; a2 <- tau^2
      ft <- f(alpha); ft2 <- f(a2)

      while (alpha-a2>1e-5)
      if (ft2<ft) {
        a1 <- alpha; alpha <- a2; ft <- ft2
        a2 <- a0 + (a1-a0)*tau^2
        ft2 <- f(a2)
      } else {
        a0 <- a2; ft2 <- ft; a2 <- alpha
        alpha <- a0 + (a1-a0)*tau
        ft <- f(alpha)
      }
      pdev <- ft
      start <- coefold + alpha*(start-coefold)
      if (!is.null(G$Ain) && nrow(G$Ain)>0) {
        active <- attr(start,"active")
        attr(start,"active")  <- active[abs(G$Ain[active,,drop=FALSE]%*%start-G$bin[active])<norm(G$Ain[active,,drop=FALSE],"M")*.Machine$double.eps^.9]
      }
      eta <- etaold + alpha*(eta-etaold)
      mu <- linkinv(eta)
    }

    if (iter == nstep || abs(pdev - pdevold)/(0.1 + abs(pdev)) < eps || gausid) {
       conv <- TRUE
       coef <- start
       break
    }

    if (FALSE) {
    ii <- 1
    while (pdev-pdevold > abs(pdev)*control$epsilon) { ## step halve
      if (ii > 20) break; ii <- ii + 1
      start <- (start + coefold)/2
      eta <- (eta + etaold)/2
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
      bSb <- bSbfun(start,G$S,G$off)
      pdev <- dev + sum(exp(lsp)*bSb)
    } ## divergence control
    }
    
    pdevold <- pdev
    coef <- coefold <- start
    etaold <- eta
    
  } ## main iteration
  if (!conv) warn[length(warn)+1] <- "scasm.pirls not converged"
  bsr <- matrix(0,length(start),bs)
  ng <- length(w)
  if (bs>0) for (k in 1:bs) { ## non-parametric bootstrap   
    bsr[,k] <- pcls(list(y=z,w=w*tabulate(sample(ng,replace=TRUE),ng),X=G$X[good,],C=G$C,S=G$S,off=G$off-1L,sp=exp(lsp),
                   p=start,Ain=G$Ain,bin = G$bin,active=attr(start,"active")))
  }		   
  ## weights are iterative and w prior below (w can be changed by initialization!)
  list(coef=coef,bSb=bSb,niter=iter,pearson = sum(w*(z-eta)^2),weights=w, H = crossprod(G$X,w*G$X),
       fitted.values=mu,linear.predictors=eta,residuals=z-eta,pdev=pdev,warn=warn,n=n,y=y,w=weights,bs=bsr) 
} ## scasm.pirls

gH2ls <- function(b,g,R,eta=NULL,WX=NULL) {
## Produces a least squares problem ||z-Rb||^2 with gradient g and Hessian H at b.
## H = 2 P'R'RP, where P is a pivot matrix and R a pivoted Colesky factor. If H is
## not full rank then the least squares problem includes an extra ridge penalty on the
## unidentifiable coefficient subspace. 
  p <- ncol(R); r <- attr(R,"rank")
  if (r<p ) { ii <- (r+1):p; R[ii,ii] <- diag(1,nrow=p-r) } ## deal with rank deficiency
  pi <- attr(R,"pivot"); ip <- order(pi) ## pivot and inverse pivot sequences
  if (is.null(b)) {
    z <- forwardsolve(R,(drop(eta %*% WX)-g)[pi]/2,upper.tri=TRUE,transpose=TRUE) ## R^{-T}P(Hb-g)/2
  } else {
    z <- R %*% b[pi]-forwardsolve(R,g[pi]/2,upper.tri=TRUE,transpose=TRUE)
  }
  list(z=z,R=R[,ip]) ## R = RP here
} ## gH2ls

scasm.newton <- function(G,lsp,control,start=NULL,scale.known=TRUE,bs=0) {
## Newton sequential QP fitting of shape constrained models for extended and general families
  warn <- rep("",0)
  y <- G$y; nobs <- length(G$w); family <- G$family; weights <- G$w
  mustart <- etastart <- NULL
  if (inherits(G$family,"general.family")) { ## general family initialization
    ## get unconstrained initial coefs, from which Hessian and grad can be computed
    start0 <- start
    x <- G$X; offset <- G$offset; E <- G$Eb
    eval(G$family$initialize)
    if (!is.null(start0)) start <- start0
    efam <- FALSE; llf <- family$ll
    WX <- eta <- NULL; scale1 <- 1
    ucstart <- 2 ## number of iterations before inequality constraints used
  } else { ## extended family initialization
    ## get initial mu and hence initial Hessian and grad
    dev.resids <- G$family$dev.resids;
    eval(G$family$initialize)
    if (!is.null(start)) mustart <- G$family$linkinv(as.numeric(G$X %*% start + G$offset))
    mu <- mustart; eta <- G$family$linkfun(mu); efam <- TRUE
    scale1 <- if (is.null(G$family$scale)) 1 else G$family$scale
    theta <- G$family$getTheta()
    ucstart <- 1
  }
  if (is.null(start)) {
    feasible <- FALSE; coef <- NULL
  } else {
    coef <- as.numeric(start)
    feasible <- ((is.null(G$Ain)||!nrow(G$Ain)||min(G$Ain%*%coef - G$bin) >= -.Machine$double.eps^.5) &&
                   (is.null(G$C)||!nrow(G$C)||isTRUE(all.equal(G$C%*%coef,G$C%*%G$beta0))))	   
  }

  if (!feasible) { ## start not feasible so use G$beta0 as initial feasible point 
    if (efam) {
      mu <- G$family$linkinv(as.numeric(G$X %*% G$beta0 + G$offset))
      dev <- sum(dev.resids(y, mu, G$w,theta))
    } else { ## general family
      ll <- llf(y,G$X,G$beta0,G$w,G$family,offset=G$offset,deriv=0)
      dev <- -2*ll$l
    }
    bSb <- bSbfun(G$beta0,G$S,G$off); penalty <- sum(exp(lsp)*bSb)
    pdev.feasible <- dev + penalty ## penalized deviance according to G$beta0
    if (!is.finite(pdev.feasible))  pdev.feasible <- Inf
    pdev0 <- pdev.feasible; 
  } else ucstart <- 0 ## no unconstrained initialization if start already feasible

  w <- rep(1,ncol(G$X)); pdev <- pdev0 <- Inf;conv <- FALSE
  dth0 <- theta0 <- NULL; thmult <- 1
  tol <- .Machine$double.eps^.75

  for (iter in 1L:control$maxit) { ## Main loop
  
    if (iter>1) { ## convert to equivalent least squares problem and solve...
      suppressWarnings(R <- chol(H/2,pivot=TRUE)) ## PHP'/2 = R'R
      ls <- gH2ls(coef,g,R,eta=eta,WX=WX)
      if (iter<=ucstart) { ## run without inequality constraints to initialize
        Ain <- matrix(0,0,0);bin <- rep(0,0);feasible <- TRUE
      } else { ## inequality constraints imposed
        if (iter==ucstart+1&&is.finite(pdev.feasible)) { ## note never triggered if ucstart<=0
	  coef <- G$beta0
	  pdev0 <- pdev.feasible
	  feasible <- TRUE
	}  
        Ain <- G$Ain;bin <- G$bin
      }
      start <- pcls(list(y=ls$z,w=w,X=ls$R,C=G$C,S=G$S,off=G$off-1L,sp=exp(lsp),
                    p=G$beta0,Ain=Ain,bin=bin))
      if (efam) mu <- G$family$linkinv(as.numeric(G$X %*% start + G$offset))		    
    } ## if (iter>1)
    not.ok <- TRUE
    while (not.ok) { ## step control
      if (efam) {
        dev <- sum(dev.resids(y, mu, G$w,theta))
      } else { ## general family
        ll <- llf(y,G$X,start,G$w,G$family,offset=G$offset,deriv=1)
	dev <- -2*ll$l
      }
      if (is.null(start)) {
        if (is.finite(dev)) not.ok <- FALSE else stop("non-finite initial deviance")
      } else { 
        bSb <- bSbfun(start,G$S,G$off); penalty <- sum(exp(lsp)*bSb)
	pdev <- dev + penalty
	if (!is.null(coef) && (!is.finite(pdev) || (iter>1 && feasible && pdev - pdev0 > tol * abs(pdev)))) {
          start <- (start + coef)/2
	  if (efam) mu <- G$family$linkinv(as.numeric(G$X %*% start + G$offset))
	  if (isTRUE(all.equal(as.numeric(start),as.numeric(coef),tol=.Machine$double.eps^.75*mean(abs(coef))))) {
            warn[length(warn)+1] <- "step failure"
	    start <- coef
	    pdev <- pdev0
	    not.ok <- FALSE
          }
        } else not.ok <- FALSE
      }  
    } ## not.ok step control loop
    if (iter > ucstart+1 && feasible && abs(pdev-pdev0)<control$epsilon*abs(pdev)) { conv <- TRUE;break}
    pdev0 <- pdev
    if (!is.null(start) && efam && (G$family$n.theta>0||scale1<0)) {
      ## deal with any scale and theta parameter estimation
      theta <- estimate.theta(theta,G$family,y,mu,scale=scale1,wt=G$w,tol=1e-6)
      if (!is.null(theta0)) { ## theta zig-zag control
        dth <- theta-theta0 ## proposed change in theta
	if (is.null(dth0)) dth0 <- dth else {
          if (any(dth0*dth<0&abs(dth)>abs(dth0)*.5)) thmult <- thmult/2 else {
            thmult <- min(1,thmult*1.5)
          }
        }
	theta <- theta0 + dth*thmult
	dth0 <- dth
      }
      
      theta0 <- theta
      if (!is.null(G$family$scale) && !scale.known) {
	scale <- exp(theta[family$n.theta+1])
	theta <- theta[-(family$n.theta+1)]
      }
      if (scale1<0 && iter>4) { ## freeze scale update
        scale1 <- scale
	theta0 <- theta0[-(family$n.theta+1)]
	dth0 <- dth0[-(family$n.theta+1)]
      }
      G$family$putTheta(theta)
      ## need to recompute pdev0, with new theta...
      pdev0 <- sum(dev.resids(y, mu, G$w,theta)) + penalty
    } ## theta estimation
  
    coef <- start ## current coef acceptable
    
    if ((!feasible||iter<=ucstart) && iter > 1) { ## check whether coef vector feasible yet 
      feasible <- ((is.null(G$Ain)||!nrow(G$Ain)||min(G$Ain%*%coef - G$bin) >= -.Machine$double.eps^.5) &&
                   (is.null(G$C)||!nrow(G$C)||isTRUE(all.equal(G$C%*%coef,G$C%*%G$beta0))))
    }
    ## get the gradient and Hessian of the deviance
    if (efam) {
      dd <- dDeta(y,mu,G$w,theta,G$family,0)
      WX <- dd$Deta2*G$X
      H <- crossprod(G$X,WX)
      g <- drop(dd$Deta %*% G$X)
    } else { ## general family
      H <- -2*ll$lbb; g <- -2*ll$lb
    }
   
  } ## main loop

  if (!scale.known) { ## Get scale MLE conditional on mu and theta
    G$family$n.theta <- 0
    theta0 <- estimate.theta(theta,G$family,y,mu,scale = -1,wt=G$w,tol=1e-6)
    scale <- exp(theta0[family$n.theta+1]);
  } else scale <- scale1

  if (!conv) warn[length(warn)+1] <- "scasm.newton not converged"  
  if (efam) {
    eta <- G$family$linkfun(mu); rsd <- dev.resids(y, mu, G$w,theta)
  } else { ## general family
    lpi <- attr(G$X,"lpi")
    if (is.null(lpi)) { ## only one...
      eta <- as.numeric(G$X %*% coef) + if (is.null(G$offset)) 0 else G$offset
      mu <- G$family$linkinv(eta) 
    } else { ## multiple...
      mu <- eta <- matrix(0,nrow(G$X),length(lpi))
      if (!is.null(G$offset)) G$offset[[length(lpi)+1]] <- 0
      for (j in 1:length(lpi)) {
        eta[,j] <- as.numeric(G$X[,lpi[[j]],drop=FALSE] %*% coef[lpi[[j]]]) +
        if (is.null(G$offset[[j]])) 0 else offset[[j]]
        mu[,j] <- G$family$linfo[[j]]$linkinv(eta[,j]) 
      }
    }
    rsd <- NULL
  }
  if (!feasible) stop("failed to find feasible starting values")
  if (iter==control$maxit) warn[length(warn)+1] <- "scasm.newton failed to converge"
  aic.model <- if (inherits(G$family,"general.family")) dev else G$family$aic(y,mu,theta,G$w,scale*sum(G$w))

  bsr <- matrix(0,length(start),bs)
  if (bs>0) { ## always fixed Hessian version
    if (efam) {
      g <- matrix(0,ncol(G$X),bs)
      for (k in 1:bs) { ## bootstrap gradient vectors
         wb <- tabulate(sample(nobs,replace=TRUE),nobs) ## np bootstrap weights
         g[,k] <- drop((dd$Deta*wb) %*% G$X) ## bootstrap gradient
      }
    } else g <- -2*llf(y,G$X,start,G$w,G$family,offset=G$offset,deriv=-bs)$lb
    if (ncol(g)) for (k in 1:bs) { ## non-parametric bootstrap for coefs
      ls <- gH2ls(coef,g[,k],R,eta=eta) ## convert to LS using original Hessian
      bsr[,k] <- pcls(list(y=ls$z,w=w,X=ls$R,C=G$C,S=G$S,off=G$off-1L,sp=exp(lsp),##p=G$beta0,Ain=Ain,bin=bin)) 
                      p=coef,Ain=G$Ain,bin=G$bin,active=attr(coef,"active"))) ## seems to work just as well as above
     	    
    } else bsr <- matrix(0,length(start),0) ## for ll that can't yet do this		    
  }
  if (FALSE&&bs>0) {
    for (k in 1:bs) {
      if (efam) {
        wb <- tabulate(sample(nobs,replace=TRUE),nobs)
        g <- drop((dd$Deta*wb) %*% G$X) ## bootstrap gradient
	## NOTE: probably need to add Hessian BS
      } else {
        lf <- llf(y,G$X,start,G$w,G$family,offset=G$offset,deriv=-1)
	g <- -2*lf$lb; suppressWarnings(R <- chol(-lf$lbb,pivot=TRUE)) ## PHP'/2 = R'R
      }
      ls <- gH2ls(coef,g,R,eta=eta) ## convert to LS using original Hessian
      bsr[,k] <- pcls(list(y=ls$z,w=w,X=ls$R,C=G$C,S=G$S,off=G$off-1L,sp=exp(lsp),
                      p=coef,Ain=G$Ain,bin=G$bin,active=attr(coef,"active")))
    }		      
  }

  list(coef=coef,bSb=bSb,niter=iter,pearson = NA ,weights=G$w, H = H/2,scale=scale,pdev=pdev,aic=aic.model,
       fitted.values=mu,linear.predictors=eta,residuals=rsd,y=y,deviance= if (efam) dev else NULL,warn=warn,
       bs=bsr)
} ## scasm.newton

check.psdef <- function(A,eps=.Machine$double.eps^.75) {
## Cholesky based test for positive semi-definite
  suppressWarnings(R <- chol(A,pivot=TRUE))
  r <- attr(R,"rank")
  if (r == ncol(R)) return(TRUE) ## +def
  ipiv <- order(attr(R,"pivot"))
  ii <- (r+1):ncol(R)
  R[ii,ii] <- 0
  max(abs(crossprod(R[,ipiv])-A))<eps ## + semi def or not?
} ## check.psdef

getVb <- function(G,H,p,zob=list(),project = TRUE,HpSt=FALSE) {
## gets cov matrix from pcls fit object, to within the scale parameter
## if project==TRUE then if Z is constraint null space return
## Vb= Z(Z'HpZ)^{-1}Z', or Z'HpZ and Z'StZ if HPSt==TRUE.
## p is coefficient vector. zob is constraints as returned by
## getZ. G is pre-fit scasm list, including constraint spec.
  c2inv <- function(Hp) {
  ## returns inverse of Hp if +ve def and pseudoinverse of nearest +ve def
  ## matrix otherwise...
    suppressWarnings(R <- chol(Hp,pivot=TRUE))
    r <- attr(R,"rank"); p <- ncol(R)
    if (r<p) { ii <- (r+1):p; R[ii,ii] <- 0 }
    r <- Rrank(R)
    if (r<p) {
      ev <- eigen(Hp,symmetric=TRUE)
      ii <- ev$values > abs(ev$values[1])*.Machine$double.eps^.9
      ev$values[ii] <- 1/sqrt(ev$values[ii])
      ev$values[!ii] <- 0
      return(crossprod(ev$values*t(ev$vectors)))
    } else {
      ipiv <- order(attr(R,"pivot"))
      return(chol2inv(R)[ipiv,ipiv,drop=FALSE])
    }  
  } ## c2inv
  
  St <- H*0
  m <- length(G$S); sedf <- pSp <- numeric(m)
  if (m) {
    for (i in 1:m) { ## get total penalty matrix and p'S_jp terms
      ii <- 1:ncol(G$S[[i]]) + G$off[i] - 1
      St[ii,ii] <- St[ii,ii] + G$S[[i]]*G$sp[i]
      pSp[i] <- sum(p[ii]*drop(G$S[[i]]%*%p[ii]))
    }
    pSp <- pmax(pSp,.Machine$double.xmin)
  }  
  
  if (project&&length(zob)) {
    H <- ZtHZproj(H,zob)
    St <- ZtHZproj(St,zob)
  }
  
  Hp <- H + St ## penalized Hessian
  
  Hindef <- !check.psdef(H)
  if (Hindef) { ## H not positive semi-definite 
    eh <- eigen(H,symmetric=TRUE)
    eh$values[eh$values<0] <- 0
    H <- eh$vectors %*% (eh$values*t(eh$vectors))
    Vb1 <- H + St
  }

  if (project&&length(zob)) { ## there are constraints
   
    if (HpSt) { ## then only penalty matrix and penalized Hessian needed
      return(list(Hp=Hp,St=St))
    }
    Vb <- c2inv(Hp) ## (Z'HpZ)^{-1}
    Vb <- ZHZtproj(Vb,zob) ## Z(Z'HpZ)^{-1}Z'
    if (Hindef) {
      Vb1 <- c2inv(Vb1)
      Vb1 <- ZHZtproj(Vb1,zob)
    } else Vb1 <- Vb
    H <-  ZHZtproj(H,zob) ## project H in same way
  } else {
    if (HpSt) return(list(Hp=Hp,St=St))
    Vb <- c2inv(Hp)
    if (Hindef) {
      Vb1 <- c2inv(Vb1)
    } else Vb1 <- Vb
  }
 
  edf <- rowSums(Vb1*H) ## vector of edf per coef.
  if (m) for (i in 1:m) { ## sedf[i]*sp[i] is suppressed EDF per penalty
    ii <- 1:ncol(G$S[[i]]) + G$off[i] - 1
    sedf[i] <- sum(Vb1[ii,ii]*G$S[[i]])
  }
  fv <- G$X %*% p ## fitted values.
  list(Vb=Vb,edf = edf,fv=fv,pSp=pSp,sedf=sedf,H=H,Hp=Hp,St=St)
} ## getVb


trSiSZ <- function(G,zob) {
## need tr(S^-Sj) terms, in this case projected into the null space of any constraints
## zob is result of getZ
  sm <- ZSZproject(G,zob) 
  trS <- numeric(length(G$sp))
  k <- sum(sapply(sm,function(x) length(x$ZSZ)))
  if (k!=length(trS)) stop("sp and smooth length not matching")
  k <- 1
  for (i in 1:length(sm)) {
    m <- length(sm[[i]]$ZSZ)
    if (m==1) {
      trS[k] <- suppressWarnings(Rrank(chol(sm[[i]]$ZSZ[[1]],pivot=TRUE)))/G$sp[k]
      k <- k + 1
    } else {
      St <- sm[[i]]$ZSZ[[1]]*G$sp[k];
      for (j in 2:m) { St <- St + sm[[i]]$ZSZ[[j]]*G$sp[k+j-1] }
      R <- suppressWarnings(chol(St,pivot=TRUE)) ## otherwise chol warns not full rank!
      r <- attr(R,"rank");p <- ncol(R)
      if (r<p) { ii <- (r+1):p; R[ii,ii] <- 0}
      r <- min(Rrank(R),sm[[i]]$df-sm[[i]]$null.space.dim);
      piv <- attr(R,"pivot"); p <- ncol(St)
      if (r<p) { ii <- (r+1):p; R[ii,ii] <- diag(1,length(ii)) }
      for (j in 1:m) {
        trS[k] <- sum(diag(forwardsolve(t(R),t(forwardsolve(t(R),sm[[i]]$ZSZ[[j]][piv,piv]))))[1:r])
	k <- k + 1
      }
    }     
  }
  trS
} ## trSiSZ

lpdgDet <- function(A,rank=0) {
## Get possibly generalized determinant of +def matrix A.
## If rank>0 then it supplies the rank of A.
## S <- crossprod(matrix(runif(20)-.5,4,5))
## lpdgDet(S);lpdgDet(S,4);sum(log(eigen(S)$values[1:4]))
   R <- suppressWarnings(chol(A,pivot=TRUE))
   r <- attr(R,"rank")
   if (rank>0 && rank<r) r <- rank
   if (r < ncol(R)) R <- suppressWarnings(chol(tcrossprod(R[1:r,]),pivot=TRUE))
   return(list(ldet=2*sum(log(diag(R))),r=r))
} ## lpdgdet



ZSZproject <- function(G,zob) {
## produces a copy of G$smooth with additional elements Z'SZ
## which are the S matrices projected into the null space of
## any active constraints as stored in zob, produced by getZ
  sm <- G$smooth
  for (k in 1:length(sm)) sm[[k]]$ZSZ <- sm[[k]]$S ## default
  if (length(zob)) for (j in 1:length(zob)) { ## work through Null space object
    k <- zob[[j]]$i.smooth ## smooth this con applies to
    zq <- zob[[j]]$q
    for (i in 1:length(sm[[k]]$S)) sm[[k]]$ZSZ[[i]] <- ## Z'S_jZ
         qr.qty(zob[[j]]$qra,t(qr.qty(zob[[j]]$qra,sm[[k]]$S[[i]])[-(1:zq),,drop=FALSE]))[-(1:zq),,drop=FALSE] 
  }
  sm
} ## ZSZproject



ZtHZproj <- function(H,zob,left=TRUE) {
## Project H into null space of constraints, by forming Z'HZ.
## Here Z is in almost block diagonal form with jth block Z_j
## it's 'almost' because the Z are not square. So ijth block
## of result is Z_i'H_ijZ_j. The non-identity Z_i are stored
## in zob, produced by getZ.
## if left==FALSE, produces HZ
## Assumption is that zob is in ascending block order...
  ## work backwards so that ii still index correct rows/cols after
  ## previous dropping...
  if (length(zob)) for (j in length(zob):1) { 
    ii <- zob[[j]]$ii ## what parameters are we looking at?
    q <- zob[[j]]$q
    H[,ii] <- t(qr.qty(zob[[j]]$qra,t(H[,ii])))
    if (left) {
      H[ii,] <- qr.qty(zob[[j]]$qra,H[ii,])
      H <- H[-ii[1:q],-ii[1:q],drop=FALSE] ## now just drop the un-needed rows/cols from H
    } else H <- H[,-ii[1:q],drop=FALSE]   
  }
  H
} ## ZtHZproj

ZHZtproj <- function(H,zob) {
## counterpart of ZtHZ that forms ZHZ'. So ijth block of result
## is Z_iH_ijZ_j'. The non-identity Z_i are stored in zob,
## produced by getZ.
## Assumption is that zob is in ascending block order...
  qtot <- sum(sapply(zob,function(x) x$q))
  if (!qtot) return(H)
  H1 <- matrix(0,nrow(H)+qtot,ncol(H))
  # First pad out rows... 
  m <- length(zob)
  ## get unpenalized row/col indices...
  iup <- sort(which(!1:nrow(H1) %in% as.integer(unlist(lapply(zob,function(x) x$ii)))))
   if (length(iup)) {
    next.iup <- iup[1] ## start of next unpenalized block to process
  } else next.iup <- nrow(H1)
  j0 <- 0 ## end of previous source (H) block
  for (j in 1:m) {
    if (zob[[j]]$ii[1]>next.iup) { ## next repara is after an unconstrained block
      ii <- next.iup:(zob[[j]]$ii[1]-1) ## index of unconstrained in target
      H1[ii,] <- H[j0+1:length(ii),] ## copy unchanged from H to correct rows of H1 
      j0 <- j0 + length(ii) ## update end of processed blocks in source, H
      next.iup <- if (any(iup>max(ii))) min(iup[iup>max(ii)]) else nrow(H1) 
    }
    ## now copy required block itself...
    ii <- zob[[j]]$ii[-(1:zob[[j]]$q)] ## target (H1) rows
    jj <- 1:length(ii) + j0 ## source rows
    j0 <- max(jj) ## new source block end
    H1[ii,] <- H[jj,] 
  }
  if (next.iup < nrow(H1)) { ## deal with any trailing unconstrained elements
    ii <- iup[iup>=next.iup]
    H1[ii,] <- H[j0+1:length(ii),]
  }

  ## Now pad out columns (as above but on cols)...
  H <- matrix(0,nrow(H1),nrow(H1))
  if (length(iup)) {
    next.iup <- iup[1]
  } else next.iup <- nrow(H1)
  j0 <- 0 

  for (j in 1:m) {
    if (zob[[j]]$ii[1]>next.iup) {
      ii <- next.iup:(zob[[j]]$ii[1]-1)
      H[,ii] <- H1[,j0+1:length(ii)]
      j0 <- j0 + length(ii)
      next.iup <- if (any(iup>max(ii))) min(iup[iup>max(ii)]) else nrow(H1)
    }
    ii <- zob[[j]]$ii[-(1:zob[[j]]$q)] ## target (H) cols
    jj <- 1:length(ii) + j0 ## source cols
    j0 <- max(jj)
    H[,ii] <- H1[,jj]
  }
  if (next.iup < nrow(H1)) {
    ii <- iup[iup>=next.iup]
    H[,ii] <- H1[,j0+1:length(ii)]
  }

  ## H now contains original H, but with zero rows/cols inserted to
  ## facilitate multiplication by Z components...
  for (j in 1:m) {
    ii <- zob[[j]]$ii
    H[,ii] <- t(qr.qy(zob[[j]]$qra,t(H[,ii])))
    H[ii,] <- qr.qy(zob[[j]]$qra,H[ii,])
  }
  H
} ## ZHZtproj

Ztbproj <- function(b,zob) {
## Project b into null space of constraints, by forming Z'b.
## zob encodes Z as returned by getZ...
## Assumption is that zob is in ascending block order...
  ## work backwards so that ii still index correct rows/cols after
  ## previous dropping...
  if (length(zob)) for (j in length(zob):1) { 
    ii <- zob[[j]]$ii ## what parameters are we looking at?
    q <- zob[[j]]$q
    b[ii] <- qr.qty(zob[[j]]$qra,b[ii])
    b <- b[-ii[1:q]] ## now just drop the un-needed rows/cols from H 
  }
  b
} ## Ztbproj

Zbproj <- function(b,zob) {
## Forms Zb where b is a vector, and zob is the structure returned by getZ defining Z.
## Assumption is that zob is in ascending block order...
  qtot <- sum(sapply(zob,function(x) x$q))
  if (!qtot) return(b)
  b1 <- rep(0,length(b)+qtot)
  # First pad out rows... 
  m <- length(zob)
  ## get unpenalized row/col indices...
  iup <- sort(which(!1:length(b1) %in% as.integer(unlist(lapply(zob,function(x) x$ii)))))
  next.iup <- if (length(iup)) iup[1] else length(b1) ## start of next unpenalized block to process
  j0 <- 0 ## end of previous source (b) block
  for (j in 1:m) {
    if (zob[[j]]$ii[1]>next.iup) { ## next repara is after an unconstrained block
      ii <- next.iup:(zob[[j]]$ii[1]-1) ## index of unconstrained in target
      b1[ii] <- b[j0+1:length(ii)] ## copy unchanged from b to correct rows of b1 
      j0 <- j0 + length(ii) ## update end of processed blocks in source, b
      next.iup <- if (any(iup>max(ii))) min(iup[iup>max(ii)]) else length(b1) 
    }
    ## now copy required block itself...
    ii <- zob[[j]]$ii[-(1:zob[[j]]$q)] ## target (b1) rows
    jj <- 1:length(ii) + j0 ## source rows
    j0 <- max(jj) ## new source block end
    b1[ii] <- b[jj] 
  }
  if (next.iup < length(b1)) { ## deal with any trailing unconstrained elements
    ii <- iup[iup>=next.iup]
    b1[ii] <- b[j0+1:length(ii)]
  }
  for (j in 1:m) {
    ii <- zob[[j]]$ii
    #H[,ii] <- t(qr.qy(zob[[j]]$qra,t(H[,ii])))
    b1[ii] <- qr.qy(zob[[j]]$qra,b1[ii])
  }
  b1
} ## Zbproj

getZ <- function(G,active=rep(0,0)) {
## Given an active set gets an object representing the constraint NULL space, Z.
## this is computed smooth by smooth, to ensure that Z'StZ is still block
## diagonal. See also ZHZtproj and ZtHZproj
  zob <- list();zi <- 0;ns <- 0
  sm <- G$smooth
  nac <- length(active) ## number of active inequality constraints
  Ain <- if (nac) G$Ain[active,,drop=FALSE] else NULL ## the active inequality constraints matrix  
  for (k in 1:length(sm)) { ## work through the smooths
    ii <- sm[[k]]$first.para:sm[[k]]$last.para ## columns related to this smooth
    if (!is.null(G$C)) { ## are there fixed constraints?
      i <- rowSums(G$C[,ii,drop=FALSE]!=0)!=0 ## constraints relevant to this term
      if (any(G$C[i,-ii]!=0)) stop("constraints overlap more than one term")
      A <- G$C[i,ii,drop=FALSE] ## retain only smooth relevant part
    } else A <- matrix(0,0,length(ii))
    if (!is.null(Ain)) { ## are there active inequality constraints?
      i <- rowSums(Ain[,ii,drop=FALSE]!=0)!=0 ## constraints relevant to this term
      if (any(Ain[i,-ii]!=0)) stop("constraints overlap more than one term")
      A <- rbind(A,Ain[i,ii]) ## append smooth relevant part to any fixed
    }
    if (nrow(A)) { ## any constraints for this smooth?
      zi <- zi + 1
      ms <- length(sm[[k]]$S)
      zob[[zi]] <- list(qra=qr(t(A)), ## QR defining Z
                        q=nrow(A), ## number of constraints
                        ii=ii, ## coefs it applies to
			i.smooth=k, ## smooth it applies to
			i.S = if (ms) 1:ms + ns else numeric(0)) ## penalty matrices it applies to
      ns <- ns + ms
    }
  }
  zob
} ## getZ

Stot <- function(G,root=TRUE,sp=NULL) {
## Obtain the total penalty matrix and optionally its square root.
## If sp==NULL then a 'balanced' total penalty is produced, including
## a penalty on any fixed constraints.
  m <- length(G$S); p <- ncol(G$X)
  if (is.null(sp)) sp <- sapply(G$S,function(x) 1/norm(x))
  if (is.null(G$C)||nrow(G$C)==0||is.null(sp)) S <- matrix(0,p,p) else {
    S <- crossprod(G$C);
    nors <- norm(S)
    if (nors>0) S <- S/nors ## penalize departure from constraints
  }  
  for (i in 1:m) {
    ii <- G$off[i] + 1:ncol(G$S[[i]]) - 1
    S[ii,ii] <- S[ii,ii] + sp[i]*G$S[[i]]
  }
  if (root) {
    suppressWarnings(R <- chol(S,pivot=TRUE))
    r <- attr(R,"rank");piv <- attr(R,"pivot")
    R <- R[1:r,]; R[,piv] <- R
  } else R <- NULL
  list(S=S,R=R) ## R'R = S
} ## Stot

lsp.cov <- function(G,Vb,beta,n=1000) {
## Obtains approx cov matrix of log smoothing params by simulation.
## Based on fact that main uncertainty in EFS log sp estimation are
## the log(b'S_jb) terms...   
  m <- length(G$S); if (!m) return(NULL)
  b <- rmvn(n,as.numeric(beta),Vb) ## n by p matrix of replicate coefs
  bSb <- matrix(0,n,m)
  for (i in 1:m) {
    ii <- G$off[i] + 1:ncol(G$S[[i]]) - 1
    bSb[,i] <- rowSums((b[,ii] %*% G$S[[i]]) * b[,ii])
  }
  cov(log(bSb)) 
} ## lsp.cov

imp.diff <- function(G,Hi,sp,beta) {
## implicit differentiation to get dbeta/drho_j
## Hi is the inverse Hessian (unscaled cov matrix)
  m <- length(G$S); if (!m) return(NULL)
  b1 <- matrix(0,length(beta),m)
  for (i in 1:m) {
    ii <- G$off[i] + 1:ncol(G$S[[i]]) - 1
    b1[,i] <- Hi[,ii] %*% G$S[[i]] %*% beta[ii] * sp[i]
  }
  b1
} ## imp.diff

scasm.fit <- function(G,control=gam.control(),gamma=1,bs=0,...) {
## fits shape constrained additive model, alternating pirls and EFS,
## having first obtained initial smoothing parameters.

 G$family <- fix.family.ls(fix.family.link(G$family))
 
 if (!is.null(G$family$preinitialize)) {
    if (inherits(G$family,"general.family")) {
      Gmod <- G$family$preinitialize(G) ## modifies some elements of G
      for (gnam in names(Gmod)) G[[gnam]] <- Gmod[[gnam]] ## copy these into G  
    } else {
      ## extended family - usually just initializes theta and possibly y
      if (!is.null(attr(G$family$preinitialize,"needG"))) attr(G$family,"G") <- G ## more complicated
      pini <- G$family$preinitialize(G$y,G$family)
      attr(G$family,"G") <- NULL
      if (!is.null(pini$family)) G$family <- pini$family
      if (!is.null(pini$Theta)) G$family$putTheta(pini$Theta)
      if (!is.null(pini$y)) G$y <- pini$y
    }
  }

  G$Eb <- Stot(G,TRUE)$R ## balanced root penalty to regularize initialization

  lsp <- if (length(G$sp)>0) log(initial.spg(G$X,G$y,G$w,G$family,G$S,G$rank,G$off,
                       offset=G$offset,L=G$L,lsp0=G$lsp0,E=G$Eb,...))  else rep(0,0)
		       
  L <- if (is.null(G$L)) diag(1,length(lsp)) else G$L
  if (any(L<0)) stop("can not handle negative weighting of log smoothing parameters")
  efam <- inherits(G$family,"extended.family")
  scale.known <- (G$family$family %in% c("poisson","binomial")) ||
                 (efam && is.null(G$family$scale)) ||
		 (!is.null(G$family$scale) && G$family$scale > 0) 		 
  start <- NULL; scale <- 1;G$sp <- exp(L%*%lsp+G$lsp0)
  dlsp <- lsp*0;alpha <- dlsp+1
  eps <- .01;
  for (i in 1:200) { ## main loop
    
    fit <- if (efam) scasm.newton(G,log(G$sp),start=start,control=control,scale.known=scale.known)
           else scasm.pirls(G,log(G$sp),start=start,control=control,eps=eps)
    active <- attr(fit$coef,"active") ## active inequality constraints
    zob <- getZ(G,active)
    v <- getVb(G,H=fit$H,p=fit$coef,zob=zob,project = TRUE)
    if (efam) {
      if (!scale.known) fit$scale <- fit$scale*nrow(G$X)/(nrow(G$X)-sum(v$edf))
      G$family$scale <- fit$scale ## updated scale estimate
    }  
    if (!scale.known) scale <- if (is.null(fit$scale)) fit$pearson/(nrow(G$X)-sum(v$edf)) else fit$scale
    if (length(lsp)==0) break
    trS <- trSiSZ(G,zob)
    redf <- drop(t(L)%*%((trS-v$sedf)*G$sp)) ## residual EDF per sp
    if (any(redf<=0)) {
      if (min(redf) < -1e-4) warning(paste("min(redf)=",min(redf)))
      redf[redf<1e-7] <- 1e-7
    }
    sedf <- drop(t(L)%*%(v$sedf*G$sp)) ## suppressed EDF per sp 
    ii <- redf > .05; lsp1 <- lsp
    dlsp0 <- dlsp;
    dlsp <- log(scale*gamma) + log(redf) - log(drop(t(L)%*%(v$pSp*G$sp)))
    maxstep <- max(abs(dlsp))
    if (i==1 && maxstep > 8)  dlsp <- 8*dlsp/maxstep
    if (i>2 && maxstep > 4) dlsp <- 4*dlsp/maxstep ## avoid overly long steps
    ## don't increase if already too close to total suppression,
    ## or decrease if already too close to no suppression
    dlsp[(redf<.05 & dlsp>0)|(sedf<.05 & dlsp<0)] <- 0
    if (i>1) {
      jj <- dlsp*dlsp0 > 0 ## same sign
      alpha[!jj] <- alpha[!jj]/2
      if (i>2) {
        ii <- jj & jj.old ## same sign 3 times in a row 
        alpha[ii] <- pmin(2,alpha[ii]*1.2)
      }
      jj.old <- jj
    }
    lsp1 <- dlsp*alpha + lsp
    if (i==7) eps <- 1e-3 # pirls eps
    conv <- (i>1 && max(abs(pSpold-v$pSp)/pmax(1e-4,v$pSp))<.01)||(max(abs(lsp1-lsp))<.01)
    if (conv) { if (eps<.005) break else eps <- 1e-3}
    lsp <- lsp1
    G$sp <- exp(drop(L%*%lsp+G$lsp0))
    start <- fit$coef
    pSpold <- v$pSp
  } ## main loop
  if (length(fit$warn)) for (i in 1:length(fit$warn)) warning(fit$warn[i])
  ii <- diag(v$Vb) < 0; diag(v$Vb)[ii] <- 0 ## avoid triggering plot problem
  F = v$Vb %*% v$H;
  b1 <- imp.diff(G,Hi=v$Vb,sp=G$sp,beta=fit$coef) ## dbeta/dlsp
  
  if (!inherits(G$family,"extended.family")) { ## still need aic
    fit$aic <- sum(G$family$dev.resids(fit$y,fit$fitted.values,G$w))/scale -
               2*G$family$ls(fit$y,G$w,fit$n,scale)[1]
  }
  if (efam&&!scale.known) G$family$scale <- -1
  object <- list(coefficients = fit$coef,sp=exp(lsp),Vp=v$Vb*scale,Ve=F%*%v$Vb*scale,
       scale=scale,edf=v$edf,prior.weights=fit$weights,family=G$family,deviance=fit$deviance,
       residuals=fit$residuals,scale.known=scale.known,F = F,R=t(mroot(v$H)),
       outer.info=list(converged = i<200,niter=i),y=fit$y,aic=fit$aic+2*sum(v$edf),
       linear.predictors=fit$linear.predictors,fitted.values=fit$fitted.values)

  Vrho <- lsp.cov(G,object$Vp,fit$coef,n=1000) ## approx lsp cov matrix 
  object$Vc <- b1%*%Vrho%*%t(b1) + object$Vp
  object$edf2 <- rowSums(object$Vc*v$H)/scale
  object$edf1 <- 2*v$edf - rowSums(F*t(F))
  names(object$edf) <- names(object$edf1) <- names(object$edf2) <- names(object$coefficients)
  ## extended family may need to manipulate fit object.
  if (!is.null(G$family$postproc)) {
    if (inherits(G$family,"general.family")) eval(G$family$postproc) else {
      posr <- G$family$postproc(family=G$family,y=G$y,prior.weights=G$w,
              fitted=fit$fitted.values,linear.predictors=fit$linear.predictors,offset=G$offset,
	      intercept=G$intercept)
      if (!is.null(posr$family)) object$family$family <- posr$family
      if (!is.null(posr$deviance)) object$deviance <- posr$deviance
      if (!is.null(posr$null.deviance)) object$null.deviance <- posr$null.deviance	      
    }
  } 
  if (is.null(object$deviance)) object$deviance <- sum(residuals.gam(object,"deviance")^2)
  ## now get LAML *excluding* active constraints (removing the 'as.numeric' in next line would re-instate active)
  zob <- getZ(G) 
  z <- getVb(G,H=fit$H,p=fit$coef,zob=zob,project = TRUE,HpSt=FALSE)
  ldetS <- lpdgDet(z$St); M <- ncol(z$Hp) - ldetS$r
  laml <- (-fit$pdev/scale + ldetS$ldet -lpdgDet(z$Hp)$ldet)/2 +M*log(2*pi)/2
  if (!inherits(G$family,"general.family")) laml <- laml +
     if (inherits(G$family,"extended.family")) G$family$ls(fit$y,G$w,G$family$getTheta(),scale)$ls else
                                               G$family$ls(fit$y,fit$w,fit$n,scale)[1]
  object$laml <- laml
  if (bs) {
    object$bs <- if (inherits(G$family,"extended.family"))
               scasm.newton(G,log(G$sp),start=start,control=control,scale.known=scale.known,bs=bs)$bs else
               scasm.pirls(G,log(G$sp),start=start,control=control,eps=eps,bs=bs)$bs
    attr(object$bs,"svar.only") <- TRUE ## signal that this covers sampling var, not smoothing bias.
  }	       
  object$Vp <- z$Vb*scale; object$Vc<- b1%*%Vrho%*%t(b1) + object$Vp ## no active cons
  object$Ve = z$Vb%*%z$H%*%z$Vb*scale; attr(object$Vp,"zob") <- zob; object$St <- z$St
  object 
} ## scasm.fit









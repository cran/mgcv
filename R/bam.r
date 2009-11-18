## routines for very large dataset generalized additive modelling.
## (c) Simon N. Wood 2009 


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
}


qr.update <- function(Xn,yn,R=matrix(0,0,ncol(Xn)),f=array(0,0),y.norm2=0)
## Let X = QR and f = Q'y. This routine updates f and R
## when Xn is appended to X and yn appended to y. If R has no rows
## then initial QR of Xn is performed. ||y||^2 is accumulated as well
{ #qrx <- qr(Xn)
  p <- ncol(Xn)
  #fn <- qr.qty(qrx,yn)[1:p] 
  y.norm2 <- y.norm2+sum(yn*yn)
  if (nrow(R)) {
    Xn <- rbind(R,Xn)
    yn <- c(f,yn)
  }
  qrx <- qr(Xn,tol=0)
  fn <- qr.qty(qrx,yn)[1:p]

  list(R = qr.R(qrx),f=fn,y.norm2=y.norm2)
}

mini.mf <-function(mf,chunk.size) {
## takes a model frame and produces a representative subset of it, suitable for 
## basis setup.
  n <- nrow(mf)
  if (n<=chunk.size) return(mf)
  seed <- get(".Random.seed", envir = .GlobalEnv)
  kind <- RNGkind(NULL)
  RNGkind("default", "default")
  set.seed(66)
  ind <- sample(1:n,chunk.size)
  RNGkind(kind[1], kind[2])
  assign(".Random.seed", seed, envir = .GlobalEnv)
  mf0 <- mf[ind,] ## random sample of rows
  ## now need to ensure that max and min are in sample for each element of mf0
  ## note that min and max might occur twice, but this shouldn't matter (and
  ## is better than min overwriting max, for example)
  for (j in 1:length(mf)) if (is.numeric(mf0[[j]])) {
    if (is.matrix(mf0[[j]])) { ## find row containing minimum
      j.min <- min((1:n)[as.logical(rowSums(mf[[j]]==min(mf[[j]])))])
      j.max <- min((1:n)[as.logical(rowSums(mf[[j]]==max(mf[[j]])))])
      mf0[[j]][1,] <- mf[[j]][j.min,]
      mf0[[j]][2,] <- mf[[j]][j.max,] 
    } else { ## vector
      mf0[[j]][1] <- min(mf[[j]])
      mf0[[j]][2] <- max(mf[[j]]) 
    }
  }
  mf0
}

bgam.fit <- function (G, mf, chunk.size, gp ,scale ,gamma,method, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs),
    control = gam.control(), intercept = TRUE)
{   y <- mf[[gp$response]]
    weights <- G$w
    conv <- FALSE
    nobs <- nrow(mf)
    nvars <- ncol(G$X)
    offset <- G$offset
    family <- G$family
    G$family <- gaussian() ## needed if REML/ML used
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
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
 
    coefold <- NULL
    eta <- if (!is.null(etastart))
         etastart
    else family$linkfun(mustart)
    
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
       stop("cannot find valid starting values: please specify some")
    dev <- sum(dev.resids(y, mu, weights))*2 ## just to avoid converging at iter 1
    boundary <- conv <- FALSE
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
 
    G$coefficients <- rep(0,ncol(G$X))
    class(G) <- "gam"
    
    conv <- FALSE
    for (iter in 1L:control$maxit) { ## main fitting loop
       ## accumulate the QR decomposition of the weighted model matrix
       wt <- rep(0,0) 
       devold <- dev
       dev <- 0     
       for (b in 1:n.block) {
        
         ind <- start[b]:stop[b]
         G$model <- mf[ind,]
         X <- predict(G,type="lpmatrix")
         if (iter>1) eta1 <- drop(X%*%coef) + offset[ind] else eta1 <- eta[ind]
         mu <- linkinv(eta1) 
         y <- G$model[[gp$response]] ## - G$offset[ind]
         weights <- G$w[ind]
         mu.eta.val <- mu.eta(eta1)
         good <- (weights > 0) & (mu.eta.val != 0)
         z <- (eta1 - offset[ind])[good] + (y - mu)[good]/mu.eta.val[good]
         w <- (weights[good] * mu.eta.val[good]^2)/variance(mu)[good]
         dev <- dev + sum(dev.resids(y,mu,weights))
         wt <- c(wt,w)
         w <- sqrt(w)
         if (b == 1) qrx <- qr.update(w*X[good,],w*z) 
         else qrx <- qr.update(w*X[good,],w*z,qrx$R,qrx$f,qrx$y.norm2)
         gc()
      }
   
      G$n <- nobs
      G$y <- mf[[gp$response]]
   
      rss.extra <- qrx$y.norm2 - sum(qrx$f^2)
      
      if (control$trace)
         cat("Deviance =", dev, "Iterations -", iter,"\n")

      if (!is.finite(dev)) stop("Non-finite deviance")

      ## preparation for working model fit is ready, but need to test for convergence first
      if (iter>2 && abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
          conv <- TRUE
          coef <- start
          break
      }

      if (method=="GCV.Cp") {
         fit <- magic(qrx$f,qrx$R,G$sp,G$S,G$off,L=G$L,lsp0=G$lsp0,rank=G$rank,
                      H=G$H,C=G$C,gamma=gamma,scale=scale,gcv=(scale<=0),
                      extra.rss=rss.extra,n.score=G$n)
 
         post <- magic.post.proc(qrx$R,fit,qrx$f*0+1) 
      } else { ## method is "REML" or "ML"
        y <- G$y; w <- G$w; n <- G$n;offset <- G$offset
        G$y <- qrx$f
        G$w <- G$y*0+1
        G$X <- qrx$R
        G$n <- length(G$y)
        G$offset <- G$y*0
        G$dev.extra <- rss.extra
        G$pearson.extra <- rss.extra
        G$n.true <- n
        object <- gam(G=G,method=method,gamma=gamma,scale=scale)
        y -> G$y; w -> G$w; n -> G$n;offset -> G$offset
      }

      if (method=="GCV.Cp") { 
        object <- list()
        object$coefficients <- fit$b
        object$edf <- post$edf
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
      }

      coef <- object$coefficients
        
      if (any(!is.finite(coef))) {
          conv <- FALSE
          warning("non-finite coefficients at iteration ",
                  iter)
          break
      }
    } ## fitting iteration

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
      
  
   
  #  wtdmu <- if (intercept)
  #      sum(weights * y)/sum(weights)
  #  else linkinv(offset)
  #  nulldev <- sum(dev.resids(y, wtdmu, weights))
  object$wt <- wt
  object
} ## end bgam.fit


bam.fit <- function(G,mf,chunk.size,gp,scale,gamma,method) 
## function that does big additive model fit in strictly additive case
{  n <- nrow(mf)
   if (n>chunk.size) { 
    G$coefficients <- rep(0,ncol(G$X))
    class(G) <- "gam"
    ind <- 1:chunk.size
    G$model <- mf[ind,]
    X <- predict(G,type="lpmatrix")
    y <- G$model[[gp$response]] - G$offset[ind]
    w <- sqrt(G$w[ind])
    qrx <- qr.update(w*X,w*y)
    n.block <- n%/%chunk.size ## number of full sized blocks
    stub <- n%%chunk.size ## size of end block
    if (n.block>1) for (i in 2:n.block) {
      ind <- ind + chunk.size
      G$model <- mf[ind,]
      X <- predict(G,type="lpmatrix")
      y <- G$model[[gp$response]] - G$offset[ind]
      w <- sqrt(G$w[ind])
      qrx <- qr.update(w*X,w*y,qrx$R,qrx$f,qrx$y.norm2)
      gc()
    }
    if (stub>0) {
      ind <- (n-stub+1):n
      G$model <- mf[ind,]
      X <- predict(G,type="lpmatrix")
      y <- G$model[[gp$response]] - G$offset[ind]
      w <- sqrt(G$w[ind])
      qrx <- qr.update(w*X,w*y,qrx$R,qrx$f,qrx$y.norm2)
      gc()
    }    
    G$n <- n
    G$y <- mf[[gp$response]]
   
  } else { 
    qrx <- qr.update(sqrt(G$w)*G$X,sqrt(G$w)*G$y)
  }

  rss.extra <- qrx$y.norm2 - sum(qrx$f^2)
 
  if (method=="GCV.Cp") {
    fit <- magic(qrx$f,qrx$R,G$sp,G$S,G$off,L=G$L,lsp0=G$lsp0,rank=G$rank,
               H=G$H,C=G$C,gamma=gamma,scale=scale,gcv=(scale<=0),
               extra.rss=rss.extra,n.score=G$n)
 
    post <- magic.post.proc(qrx$R,fit,qrx$f*0+1) 
  } else { ## method is "REML" or "ML"
    y <- G$y; w <- G$w; n <- G$n;offset <- G$offset
    G$y <- qrx$f
    G$w <- G$y*0+1
    G$X <- qrx$R
    G$n <- length(G$y)
    G$offset <- G$y*0
    G$dev.extra <- rss.extra
    G$pearson.extra <- rss.extra
    G$n.true <- n
    object <- gam(G=G,method=method,gamma=gamma,scale=scale)
    y -> G$y; w -> G$w; n -> G$n;offset -> G$offset
  }

  if (method=="GCV.Cp") { 
    object <- list()
    object$coefficients <- fit$b
    object$edf <- post$edf
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
 
  object
}




bam <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,
                offset=NULL,method="REML",control=gam.control(),scale=0,gamma=1,knots=NULL,
                sp=NULL,min.sp=NULL,paraPen=NULL,chunk.size=10000,...)

## Routine to fit an additive model to a large dataset. The model is stated in the formula, 
## which is then interpreted to figure out which bits relate to smooth terms and which to 
## parametric terms.
## This is a modification of `gam' designed to build the QR decompostion of the model matrix 
## up in chunks, to keep memory costs down.
{ require(mgcv)
  if (is.character(family))
            family <- eval(parse(text = family))
  if (is.function(family))
            family <- family()
  if (is.null(family$family))
            stop("family not recognized")
  ##family = gaussian() ## no choise here
  if (family$family=="gaussian"&&family$link=="identity") am <- TRUE else am <- FALSE
  if (scale==0) { if (family$family%in%c("poisson","binomial")) scale <- 1 else scale <- -1} 
  if (!method%in%c("GCV.Cp","REML","ML","P-REML","P-ML")) stop("un-supported smoothness selection method")
  gp<-interpret.gam(formula) # interpret the formula 
  cl<-match.call() # call needed in gam object for update to work
  mf<-match.call(expand.dots=FALSE)
  mf$formula<-gp$fake.formula 
  mf$method <-  mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp <-
  mf$gamma <- mf$paraPen<- mf$chunk.size  <- mf$...<-NULL
  mf$drop.unused.levels<-TRUE
  mf[[1]]<-as.name("model.frame")
  pmf <- mf
 
  pmf$formula <- gp$pf
  pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
  pterms <- attr(pmf,"terms") ## pmf only used for this
  rm(pmf);gc()

  mf <- eval(mf, parent.frame()) # the model frame now contains all the data 
  if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
  terms <- attr(mf,"terms")
  gc()  

  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

  ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
  if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data) 

  dl <- eval(inp, data, parent.frame())
  if (!control$keepData) { rm(data);gc()} ## save space
  names(dl) <- vars ## list of all variables needed
  var.summary <- mgcv:::variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
  rm(dl);gc() ## save space    

  ## need mini.mf for basis setup, then accumulate full X, y, w and offset
  mf0 <- mini.mf(mf,chunk.size)
    
  G<-mgcv:::gam.setup(gp,pterms=pterms,data=mf0,knots=knots,sp=sp,min.sp=min.sp,
                 H=NULL,absorb.cons=TRUE,sparse.cons=0,select=FALSE,
                 idLinksBases=TRUE,scale.penalty=control$scalePenalty,
                 paraPen=paraPen)

  G$var.summary <- var.summary
  G$family <- family
  G$terms<-terms;G$pterms<-pterms
  
  n <- nrow(mf)
  
  if (is.null(mf$"(weights)")) G$w<-rep(1,n)
  else G$w<-mf$"(weights)"    
  
  G$offset <- model.offset(mf)  
  if (is.null(G$offset)) G$offset <- rep(0,n)

  if (ncol(G$X)>nrow(mf)+nrow(G$C)) stop("Model has more coefficients than data")
 
  G$cl<-cl;
  G$am <- am
     
  G$min.edf<-G$nsdf-dim(G$C)[1]
  if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim

  G$formula<-formula
  environment(G$formula)<-environment(formula)
  
  G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
  G$max.half<-control$mgcv.half     # max step halving in bfgs optimization


  ## now build up proper model matrix, and deal with y, w, and offset...
  
  if (control$trace) cat("Setup complete. Calling fit\n")

  if (am) {
    object <- bam.fit(G,mf,chunk.size,gp,scale,gamma,method)
  } else {
    object <- bgam.fit(G, mf, chunk.size, gp ,scale ,gamma,method=method,
                       control = control,...)
  }

  if (control$trace) cat("Fit complete. Finishing gam object.\n")

  if (scale < 0) { object$scale.estimated <- TRUE;object$scale <- object$scale.est} else {
    object$scale.estimated <- FALSE; object$scale <- scale
  }

  object$assign <- G$assign # applies only to pterms  
  object$boundary <- FALSE  # always FALSE for this case
  object$call<-G$cl # needed for update() to work 
  object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
 
  object$contrasts <- G$contrasts
  object$control <- control
  object$converged <- TRUE ## no iteration
  object$data <- NA ## not saving it in this case
  object$df.null <- nrow(mf)
  object$df.residual <- object$df.null - sum(object$edf) 
 
  object$family <- family
  object$formula<-G$formula 
 
  object$iter <- 1
  #object$linear.predictors <- NA
  if (method=="GCV.Cp") {
    if (scale<=0) object$method <- "GCV" else object$method <- "UBRE"
  } else {
    object$method <- method
  }
  object$min.edf<-G$min.edf
  object$model <- mf;rm(mf);gc()
  object$na.action <- attr(object$model,"na.action") # how to deal with NA's
  object$nsdf <- G$nsdf
  names(object$coefficients)[1:G$nsdf] <- colnames(G$X)[1:G$nsdf]
  object$offset <- G$offset
  object$prior.weights <- G$w
  object$pterms <- G$pterms
 
  object$smooth <- G$smooth

  object$terms <- G$terms
  object$var.summary <- G$var.summary 
 
  object$weights <- object$prior.weights
  object$xlevels <- G$xlevels
  object$y <- object$model[[gp$response]]

  gc()
  ## note that predict.gam assumes that it must be ok not to split the 
  ## model frame, if no new data supplied, so need to supply explicitly
  object$linear.predictors <- as.numeric(predict.gam(object,newdata=object$model,block.size=chunk.size))
  object$fitted.values <- family$linkinv(object$linear.predictors)
  
  object$residuals <- sqrt(family$dev.resids(object$y,object$fitted.values,object$weights)) * 
                      sign(object$y-object$fitted.values)
  object$deviance <- sum(object$residuals^2)
  object$aic <- family$aic(object$y,1,object$fitted.values,object$weights,object$deviance)
  object$null.deviance <- sum(family$dev.resids(object$y,mean(object$y),object$weights))
  class(object) <- c("gam","glm","lm")
  object
}



#### ISSUES:   
## ? negative binomial support --- docs say it's there...
## need to add `bam' examples to the test suite, and possibly update the test suite
## to work through failures, logging them. 
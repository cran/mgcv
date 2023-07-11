#  R plotting routines for the package mgcv (c) Simon Wood 2000-2017
##  With contributions from Henric Nilsson


in.out <- function(bnd,x) {
## tests whether point defined by each row of x is inside 
## or outside boundary defined by bnd. bnd my be made up of multiple 
## nested loops.
  if (!is.matrix(x)) x <- matrix(x,1,2)
  if (is.list(bnd)) { ## convert list of lists to matrix form
    b1 <- bnd[[1]][[1]]
    b2 <- bnd[[1]][[2]]
    if (length(bnd)>1) for (i in 2:length(bnd)) {
      b1 <- c(b1,NA,bnd[[i]][[1]])
      b2 <- c(b2,NA,bnd[[i]][[2]])
    }
    bnd <- cbind(b1,b2)
  }
  ## replace NA segment separators with a numeric code 
  lowLim <- min(bnd,na.rm=TRUE) - mean(abs(bnd),na.rm=TRUE)
  ind <- is.na(rowSums(bnd))
  bnd[ind,] <- lowLim
  n <- nrow(bnd)
  um <-.C(C_in_out,bx=as.double(bnd[,1]),by=as.double(bnd[,2]),break.code=as.double(lowLim),
          x=as.double(x[,1]),y=as.double(x[,2]),inside=as.integer(x[,2]*0),nb=as.integer(n),
          n=as.integer(nrow(x)))
  as.logical(um$inside)
}


fix.family.qf <- function(fam) {
## add quantile function to family object

  if (!inherits(fam, "family"))
        stop("fam not a family object")
  if (!is.null(fam$qf)) return(fam)  ## already exists
  family <- fam$family
  if (family=="poisson") {
    fam$qf <- function(p,mu,wt,scale) {
      qpois(p,mu)
    }
  } else if (family=="binomial") {
    fam$qf <- function(p,mu,wt,scale) {
      if (all.equal(wt,ceiling(wt))!=TRUE) {
        wt <- ceiling(wt)
	warning("non-integer binomial denominator: quantiles incorrect")
      }
      qbinom(p,wt,mu)/(wt + as.numeric(wt==0))
    }
  } else if (family=="Gamma") {
    fam$qf <- function(p,mu,wt,scale) {
      qgamma(p,shape=1/scale,scale=mu*scale)
    }
  } else if (family=="gaussian") {
    fam$qf <- function(p,mu,wt,scale) {
      qnorm(p,mean=mu,sd=sqrt(scale/wt))
    }
  }
  fam
}

fix.family.rd <- function(fam) {
## add random deviate function to family objet

  if (!inherits(fam, "family"))
        stop("fam not a family object")
  if (!is.null(fam$rd)) return(fam)  ## already exists
  family <- fam$family
  if (family=="poisson") {
    fam$rd <- function(mu,wt,scale) {
     rpois(length(mu),mu)
    }
  } else if (family=="binomial") {
    fam$rd <- function(mu,wt,scale) {
      rbinom(length(mu),wt,mu)/(wt + as.numeric(wt==0))
    }
  } else if (family=="Gamma") {
    fam$rd <- function(mu,wt,scale) {
      rgamma(length(mu),shape=1/scale,scale=mu*scale)
    }    
  } else if (family=="gaussian") {
    fam$rd <- function(mu,wt,scale) {
      rnorm(length(mu),mean=mu,sd=sqrt(scale/wt))
    }
  } else if (family=="inverse.gaussian") {
    fam$rd <- function(mu,wt,scale) {
      rig(length(mu),mu,scale)
    }
  } 
  fam
}


qq.gam <- function(object, rep=0, level=.9,s.rep=10,
                   type=c("deviance","pearson","response"),
                   pch=".", rl.col=2, rep.col="gray80",...) {
## get deviance residual quantiles under good fit
  type <- match.arg(type)
  ylab <- paste(type,"residuals")

  if (inherits(object,c("glm","gam"))) {
    if (is.null(object$sig2)) object$sig2 <- summary(object)$dispersion
  } else stop("object is not a glm or gam")

  ## in case of NA & na.action="na.exclude", we need the "short" residuals:
  object$na.action <- NULL
  D <- residuals(object,type=type)

  if (object$method %in% c("PQL","lme.ML","lme.REML","lmer.REML","lmer.ML","glmer.ML")) {
    ## then it's come out of a gamm fitter and qq.gam can't see the random effects
    ## that would be necessary to get quantiles. Fall back to normal QQ plot.
    qqnorm(D,ylab=ylab,pch=pch,...)
    return()
  }

  lim <- Dq <- NULL
  if (rep==0) { 
    fam <- fix.family.qf(object$family)
    if (is.null(fam$qf))
      rep <- 50 ## try simulation if quantile function not available
    level <- 0
  } 
  n <- length(D)
  if (rep > 0) { ## simulate quantiles
    fam <- fix.family.rd(object$family)
    if (!is.null(fam$rd)) {
      ##d <- rep(0,0)
      ## simulate deviates... 
      dm <- matrix(0,n,rep)
      for (i in 1:rep) { 
        yr <- fam$rd(object$fitted.values, object$prior.weights, object$sig2)
        #di <- fam$dev.resids(yr,object$fitted.values,object$prior.weights)^.5*
        #       sign(yr-object$fitted.values)
        object$y <- yr
        dm[,i] <- sort(residuals(object,type=type))
        #d <- c(d,sort(di))
      }
      # n <- length(D)
      Dq <- quantile(as.numeric(dm),(1:n - .5)/n) 
    
      ## now get simulation limits on QQ plot
      #dm <- matrix(d,length(Dq),rep)
      alpha <- (1-level)/2
      if (alpha>.5||alpha<0) alpha <- .05
      if (level>0&&level<1) lim <- apply(dm,1,FUN=quantile,p=c(alpha,1-alpha)) else
      if (level >= 1) lim <- level 
    }
  } else {
    ## ix <- sort.int(D,index.return=TRUE)$ix ## messes up under multiple ties!
    #ix <- rank(D)
    #U <- (ix-.5)/length(D) ## code used pre-randomization - not needed
    U <- (1:n-.5)/n
    if (!is.null(fam$qf)) { 
      dm <- matrix(0,n,s.rep)
      for (i in 1:s.rep) { 
        U <- sample(U,n) ## randomize uniform quantiles w.r.t. obs
        q0 <- fam$qf(U,object$fitted.values,object$prior.weights,object$sig2)
        object$y <- q0
        dm[,i] <- sort(residuals(object,type=type)) ## original proposal
      }
      Dq <- sort(rowMeans(dm))
    }
  }
 
  if (!is.null(Dq))  
  { qqplot(Dq,D,ylab=ylab,xlab="theoretical quantiles",ylim=range(c(lim,D)),
           pch=pch,...)
    abline(0,1,col=rl.col)
    if (!is.null(lim)) {
      if (level>=1) for (i in 1:rep) lines(Dq,dm[,i],col=rep.col) else {
        n <- length(Dq)
        polygon(c(Dq,Dq[n:1],Dq[1]),c(lim[1,],lim[2,n:1],lim[1,1]),col=rep.col,border=NA)
      }
      abline(0,1,col=rl.col) 
    }
    points(Dq,sort(D),pch=pch,...)
    return(invisible(Dq))
  } else qqnorm(D,ylab=ylab,pch=pch,...)
} ## qq.gam


k.check <- function(b,subsample=5000,n.rep=400) {
## function to check k in a gam fit... 
## does a randomization test looking for evidence of residual 
## pattern attributable to covariates of each smooth. 
  m <- length(b$smooth)
  if (m==0) return(NULL)
  rsd <- residuals(b)
 
  ve <- rep(0,n.rep)
  p.val<-v.obs <- kc <- edf<- rep(0,m)
  snames <- rep("",m)
  n <- nrow(b$model)
  if (n>subsample) { ## subsample to avoid excessive cost
    ind <- sample(1:n,subsample)
    modf <- b$model[ind,] 
    rsd <- rsd[ind]
  } else modf <- b$model
  nr <- length(rsd)
  for (k in 1:m) { ## work through smooths
    ok <- TRUE
    b$smooth[[k]]$by <- "NA" ## can't deal with by variables
    dat <- ExtractData(b$smooth[[k]],modf,NULL)$data
    if (!is.null(attr(dat,"index"))||!is.null(attr(dat[[1]],"matrix"))||is.matrix(dat[[1]])) ok <- FALSE
    if (ok) dat <- as.data.frame(dat)
    snames[k] <- b$smooth[[k]]$label
    ind <- b$smooth[[k]]$first.para:b$smooth[[k]]$last.para
    kc[k] <- length(ind)
    edf[k] <- sum(b$edf[ind]) 
    nc <- b$smooth[[k]]$dim
    if (ok && ncol(dat)>nc) dat <- dat[,1:nc,drop=FALSE] ## drop any by variables
    for (j in 1:nc) if (is.factor(dat[[j]])) ok <- FALSE 
    if (!ok) {
      p.val[k] <- v.obs[k] <- NA ## can't do this test with summation convention/factors
    } else { ## normal term
      if (nc==1) { ## 1-D term
        e <- diff(rsd[order(dat[,1])])
        v.obs[k] <- mean(e^2)/2
        for (i in 1:n.rep) {
          e <- diff(rsd[sample(1:nr,nr)]) ## shuffle 
          ve[i] <-  mean(e^2)/2
        }
        p.val[k] <- mean(ve<v.obs[k])
        v.obs[k] <- v.obs[k]/mean(rsd^2)
      } else { ## multi-D 
        if (!is.null(b$smooth[[k]]$margin)) { ## tensor product (have to consider scaling)
          ## get the scale factors...
          beta <- coef(b)[ind]
          f0 <- PredictMat(b$smooth[[k]],dat)%*%beta
          gr.f <- rep(0,ncol(dat))
          for (i in 1:nc) {
            datp <- dat;dx <- diff(range(dat[,i]))/1000
            datp[,i] <- datp[,i] + dx
            fp <- PredictMat(b$smooth[[k]],datp)%*%beta
            gr.f[i] <- mean(abs(fp-f0))/dx
          }
          for (i in 1:nc) { ## rescale distances
            dat[,i] <- dat[,i] - min(dat[,i])
            dat[,i] <- gr.f[i]*dat[,i]/max(dat[,i])
          }
        }
        nn <- 3
        ni <- nearest(nn,as.matrix(dat))$ni
        e <- rsd - rsd[ni[,1]]
        for (j in 2:nn) e <- c(e,rsd-rsd[ni[,j]])
        v.obs[k] <- mean(e^2)/2
        for (i in 1:n.rep) {
          rsdr <- rsd[sample(1:nr,nr)] ## shuffle
          e <- rsdr - rsdr[ni[,1]]
          for (j in 2:nn) e <- c(e,rsdr-rsdr[ni[,j]])
          ve[i] <-  mean(e^2)/2
        }
        p.val[k] <- mean(ve<v.obs[k])
        v.obs[k] <- v.obs[k]/mean(rsd^2)
      }
    }
  }
  k.table <- cbind(kc,edf,v.obs, p.val)
  dimnames(k.table) <- list(snames, c("k\'","edf","k-index", "p-value"))
  k.table
} ## end of k.check


gam.check <- function(b, old.style=FALSE,
		      type=c("deviance","pearson","response"), 
                      k.sample=5000,k.rep=200,
		      ## arguments passed to qq.gam() {w/o warnings !}:
		      rep=0, level=.9, rl.col=2, rep.col="gray80", ...)
## takes a fitted gam object and produces some standard diagnostic plots
{
  type <- match.arg(type)
  resid <- residuals(b, type=type)
  linpred <- if (is.matrix(b$linear.predictors)&&!is.matrix(resid)) 
             napredict(b$na.action, b$linear.predictors[,1]) else 
             napredict(b$na.action, b$linear.predictors)
##  if (b$method%in%c("GCV","GACV","UBRE","REML","ML","P-ML","P-REML","mle.REML","mle.ML","PQL")) { 
    if (is.null(.Platform$GUI) || .Platform$GUI != "RStudio") old.par <- par(mfrow=c(2,2))
    if (old.style)
      qqnorm(resid,...)
    else
      qq.gam(b, rep=rep, level=level, type=type, rl.col=rl.col, rep.col=rep.col, ...)
    plot(linpred, resid,main="Resids vs. linear pred.",
         xlab="linear predictor",ylab="residuals",...)
    hist(resid,xlab="Residuals",main="Histogram of residuals",...)
    fv <- if (inherits(b$family,"extended.family")) predict(b,type="response") else fitted(b)
    if (is.matrix(fv)&&!is.matrix(b$y)) fv <- fv[,1]
    plot(fv, napredict(b$na.action, b$y),
         xlab="Fitted Values",ylab="Response",main="Response vs. Fitted Values",...)

    gamm <- !(b$method%in%c("GCV","GACV","UBRE","REML","ML","P-ML","P-REML","fREML","NCV")) ## gamm `gam' object

    #if (is.null(.Platform$GUI) || .Platform$GUI != "RStudio") par(old.par)
    #   return(invisible())
    #}
    ## now summarize convergence information 
    if (gamm) {
      cat("\n\'gamm\' based fit - care required with interpretation.")
      cat("\nChecks based on working residuals may be misleading.")
    } else { 
      cat("\nMethod:",b$method,"  Optimizer:",b$optimizer)
      if (!is.null(b$outer.info)) { ## summarize convergence information
        if (b$optimizer[2]%in%c("newton","bfgs"))
        { boi <- b$outer.info
          cat("\n",boi$conv," after ",boi$iter," iteration",sep="")
          if (boi$iter==1) cat(".") else cat("s.")
          cat("\nGradient range [",min(boi$grad),",",max(boi$grad),"]",sep="")
          cat("\n(score ",b$gcv.ubre," & scale ",b$sig2,").",sep="")
          ev <- eigen(boi$hess)$values
          if (min(ev)>0) cat("\nHessian positive definite, ") else cat("\n")
          cat("eigenvalue range [",min(ev),",",max(ev),"].\n",sep="")
        } else { ## just default print of information ..
          cat("\n");print(b$outer.info)
        }
      } else { ## no sp, perf iter or AM case
        if (length(b$sp)==0) ## no sp's estimated  
          cat("\nModel required no smoothing parameter selection")
        else { 
          cat("\nSmoothing parameter selection converged after",b$mgcv.conv$iter,"iteration")       
          if (b$mgcv.conv$iter>1) cat("s")
         
          if (!b$mgcv.conv$fully.converged)
          cat(" by steepest\ndescent step failure.\n") else cat(".\n")
          cat("The RMS",b$method,"score gradient at convergence was",b$mgcv.conv$rms.grad,".\n")
          if (b$mgcv.conv$hess.pos.def)
          cat("The Hessian was positive definite.\n") else cat("The Hessian was not positive definite.\n")
          #cat("The estimated model rank was ",b$mgcv.conv$rank,
          #           " (maximum possible: ",b$mgcv.conv$full.rank,")\n",sep="")
        }
      }
      if (!is.null(b$rank)) {
        cat("Model rank = ",b$rank,"/",length(b$coefficients),"\n")
      }
    } ## if gamm
    cat("\n")
    ## now check k
    kchck <- k.check(b,subsample=k.sample,n.rep=k.rep)
    if (!is.null(kchck)) {
      cat("Basis dimension (k) checking results. Low p-value (k-index<1) may\n") 
      cat("indicate that k is too low, especially if edf is close to k\'.\n\n")
      printCoefmat(kchck,digits=3);
    }
    if (is.null(.Platform$GUI) ||.Platform$GUI != "RStudio") par(old.par)
##  } else plot(linpred,resid,xlab="linear predictor",ylab="residuals",...)
} ## end of gam.check

#############################################
## Start of plot method functions for smooths
#############################################

plot.random.effect <- function(x,P=NULL,data=NULL,label="",se1.mult=1,se2.mult=2,
                     partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,n3=3,
                     theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,by.resids=FALSE,scheme=0,...) {
## plot method for a "random.effect" smooth class
 
  if (is.null(P)) { ## get plotting information...
    if (!x$plot.me) return(NULL) else { ## shouldn't or can't plot 
      raw <- data[x$term][[1]]
      p <- x$last.para - x$first.para + 1
      X <- diag(p)   # prediction matrix for this term
      if (is.null(xlab)) xlabel<- "Gaussian quantiles" else xlabel <- xlab
      if (is.null(ylab)) ylabel <- "effects" else ylabel <- ylab
      if (!is.null(main)) label <- main
      return(list(X=X,scale=FALSE,se=FALSE,raw=raw,xlab=xlabel,ylab=ylabel,
             main=label))

    } ## end of basic plot data production 
  } else { ## produce plot
    b <- as.numeric(trans(P$fit+shift))
    qqnorm(b,main=P$main,xlab=P$xlab,ylab=P$ylab,...)
    qqline(b)
  } ## end of plot production
} ## end of plot.random.effect


repole <- function(lo,la,lop,lap) {
## painfully plodding function to get new lo, la relative to pole at
## lap,lop...
  ## x,y,z location of pole...
  yp <- sin(lap)
  xp <- cos(lap)*sin(lop)
  zp <- cos(lap)*cos(lop)
  
  ## x,y,z location of meridian point for pole - i.e. point lat pi/2
  ## from pole on pole's lon.

  ym <- sin(lap-pi/2)
  xm <- cos(lap-pi/2)*sin(lop)
  zm <- cos(lap-pi/2)*cos(lop)

  ## x,y,z locations of points in la, lo

  y <- sin(la)
  x <- cos(la)*sin(lo)
  z <- cos(la)*cos(lo)

  ## get angle between points and new equatorial plane (i.e. plane orthogonal to pole)
  d <- sqrt((x-xp)^2+(y-yp)^2+(z-zp)^2) ## distance from points to to pole 
  phi <- pi/2-2*asin(d/2)

  ## location of images of la,lo on (new) equatorial plane
  ## sin(phi) gives distance to plane, -(xp, yp, zp) is 
  ## direction... 
  x <- x - xp*sin(phi)
  y <- y - yp*sin(phi)
  z <- z - zp*sin(phi)

  ## get distances to meridian point
  d <- sqrt((x-xm)^2+(y-ym)^2+(z-zm)^2)
  ## angles to meridian plane (i.e. plane containing origin, meridian point and pole)...
  theta <- (1+cos(phi)^2-d^2)/(2*cos(phi))
  theta[theta < -1] <- -1; theta[theta > 1] <- 1
  theta <- acos(theta)
  
  ## now decide which side of meridian plane...

  ## get points at extremes of hemispheres on either side
  ## of meridian plane.... 
  y1 <- 0
  x1 <- sin(lop+pi/2)
  z1 <- cos(lop+pi/2)

  y0 <- 0
  x0 <- sin(lop-pi/2)
  z0 <- cos(lop-pi/2)

  d1 <- sqrt((x-x1)^2+(y-y1)^2+(z-z1)^2)
  d0 <- sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)

  ii <- d0 < d1 ## index -ve lon hemisphere 
  theta[ii] <- -theta[ii]
  
  list(lo=theta,la=phi)
} ## end of repole

lolaxy <- function(lo,la,theta,phi) {
## takes locations lo,la, relative to a pole at lo=theta, la=phi. 
## theta, phi are expressed relative to plotting co-ordinate system 
## with pole at top. Convert to x,y in plotting co-ordinates.
## all in radians!
  er <- repole(-lo,la,-pi,phi)
  er$lo <- er$lo - theta
  y <- sin(er$la)
  x <- cos(er$la)*sin(er$lo)
  z <- cos(er$la)*cos(er$lo)
  ind <- z<0
  list(x=x[ind],y=y[ind])
} ## end of lolaxy

plot.sos.smooth <- function(x,P=NULL,data=NULL,label="",se1.mult=2,se2.mult=1,
                     partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,n3=3,
                     theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,by.resids=FALSE,scheme=0,hcolors=heat.colors(100),
                     contour.col=4,...) {
## plot method function for sos.smooth terms
  if (scheme>1) return(plot.mgcv.smooth(x,P=P,data=data,label=label,se1.mult=se1.mult,se2.mult=se2.mult,
                     partial.resids=partial.resids,rug=rug,se=se,scale=scale,n=n,n2=n2,
                     theta=theta,phi=phi,jit=jit,xlab=xlab,ylab=ylab,main=main,
                     ylim=ylim,xlim=xlim,too.far=too.far,shade=shade,shade.col=shade.col,
                     shift=shift,trans=trans,by.resids=by.resids,scheme=scheme-2,
                     hcolors=hcolors,contour.col=contour.col,...))
  ## convert location of pole in plotting grid to radians
  phi <- phi*pi/180
  theta <- theta*pi/180
  
  ## re-map to sensible values...
  theta <- theta%%(2*pi)
  if (theta>pi) theta <- theta - 2*pi
 
  phi <- phi%%(2*pi)
  if (phi > pi) phi <- phi - 2*pi
  if (phi > pi/2) phi <- pi - phi
  if (phi < -pi/2 ) phi <- -(phi+pi)  

  if (is.null(P)) { ## get plotting information...
    if (!x$plot.me) return(NULL) ## shouldn't or can't plot
    ## get basic plot data 
    raw <- data[x$term]
    if (rug) { ## need to project data onto plotting grid...
      raw <- lolaxy(lo=raw[[2]]*pi/180,la=raw[[1]]*pi/180,theta,phi)
    }

    m <- round(n2*1.5)
    ym <- xm <- seq(-1,1,length=m)
    gr <- expand.grid(x=xm,y=ym)
    r <- z <- gr$x^2+gr$y^2
    z[z>1] <- NA
    z <- sqrt(1-z)

    ## generate la, lo in plotting grid co-ordinates...
    ind <- !is.na(z) 
    r <- r[ind]
    la <- asin(gr$y[ind])
    lo <- cos(la)
    lo <- asin(gr$x[ind]/lo)

    um <- repole(lo,la,theta,phi)
 
    dat <- data.frame(la=um$la*180/pi,lo=um$lo*180/pi)
    names(dat) <- x$term
    if (x$by!="NA") dat[[x$by]] <- la*0+1    

    X <- PredictMat(x,dat)   # prediction matrix for this term

    ## fix lo for smooth contouring
    lo <- dat[[2]]
    ii <- lo <= -177
    lo[ii] <- lo[ii] <- 360 + lo[ii]
    ii <- lo < -165 & lo > -177 
    ii <- ii | (abs(dat[[1]])>80) 
    lo[ii] <- NA

    return(list(X=X,scale=FALSE,se=FALSE,raw=raw,xlab="",ylab="",main="",
                ind=ind,xm=xm,ym=ym,lo=lo,la=dat[[1]]))
  } else { ## do plot
    op <- par(pty="s",mar=c(0,0,0,0))
    m <- length(P$xm); zz <- rep(NA,m*m)
    if (scheme == 0) {
      col <- 1# "lightgrey 
      zz[P$ind] <- trans(P$fit+shift)
      image(P$xm,P$ym,matrix(zz,m,m),col=hcolors,axes=FALSE,xlab="",ylab="",...)
      if (rug) { 
        if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
        points(P$raw$x,P$raw$y,...)
      }
      zz[P$ind] <- P$la
      contour(P$xm,P$ym,matrix(zz,m,m),add=TRUE,levels=c(-8:8*10),col=col,...)
      zz[P$ind] <- P$lo
      contour(P$xm,P$ym,matrix(zz,m,m),add=TRUE,levels=c(-8:9*20),col=col,...)
      zz[P$ind] <- P$fit
      contour(P$xm,P$ym,matrix(zz,m,m),add=TRUE,col=contour.col,...)
    } else if (scheme == 1) {
      col <- 1 
      zz[P$ind] <- trans(P$fit+shift)
      contour(P$xm,P$ym,matrix(zz,m,m),col=1,axes=FALSE,xlab="",ylab="",...) 
      if (rug) { 
        if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
        points(P$raw$x,P$raw$y,...)
      }
      zz[P$ind] <- P$la
      contour(P$xm,P$ym,matrix(zz,m,m),add=TRUE,levels=c(-8:8*10),col=col,lty=2,...)
      zz[P$ind] <- P$lo
      contour(P$xm,P$ym,matrix(zz,m,m),add=TRUE,levels=c(-8:9*20),col=col,lty=2,...)
      theta <- seq(-pi/2,pi/2,length=200)
      x <- sin(theta);y <- cos(theta)
      x <- c(x,x[200:1]);y <- c(y,-y[200:1])
      lines(x,y)
    } 
    par(op)
  }
} ## end plot.sos.smooth

poly2 <- function(x,col) {
## let x be a 2 col matrix defining some polygons. 
## Different closed loop sections are separated by
## NA rows. This routine assumes that loops nested within 
## other loops are holes (further nesting gives and island 
## in hole, etc). Holes are left unfilled.
## The first polygon should not be a hole.
  ind <- (1:nrow(x))[is.na(rowSums(x))] ## where are the splits?
  if (length(ind)==0|| ind[1]==nrow(x)) polygon(x,col=col,border="black") else {
    base <- x[1,]
    xf <- x
    xf[ind,1] <- base[1]
    xf[ind,2] <- base[2]
    if (!is.na(col)) polygon(xf,col=col,border=NA,fillOddEven=TRUE)
    polygon(x,border="black")
  }
} ## poly2

polys.plot <- function(pc,z=NULL,scheme="heat",lab="",...) { 
## pc is a list of polygons defining area boundaries
## pc[[i]] is the 2 col matrix of vertex co-ords for polygons defining 
## boundary of area i
## z gives the value associated with the area
  
  ## first find the axes ranges...

  for (i in 1:length(pc)) {
    yr <- range(pc[[i]][,2],na.rm=TRUE)
    xr <- range(pc[[i]][,1],na.rm=TRUE)

    if (i==1) {
      ylim <- yr
      xlim <- xr
    } else {
      if (yr[1]<ylim[1]) ylim[1] <- yr[1]
      if (yr[2]>ylim[2]) ylim[2] <- yr[2]
      if (xr[1]<xlim[1]) xlim[1] <- xr[1]
      if (xr[2]>xlim[2]) xlim[2] <- xr[2]
    }
  } ## end of axes range loop

  mar <- par("mar");
  oldpar <- par(mar=c(2,mar[2],2,1)) 

  if (is.null(z)) { ## no z value, no shading, no scale, just outlines...
     plot(0,0,ylim=ylim,xlim=xlim,xaxt="n",yaxt="n",type="n",bty="n",ylab=lab,xlab="",...)
     for (i in 1:length(pc)) { 
       poly2(pc[[i]],col=NA)
     }
  } else {
    
    nz <- names(z)
    npc <- names(pc)
    if (!is.null(nz)&&!is.null(npc)) { ## may have to re-order z into pc order.
      if (all.equal(sort(nz),sort(npc))!=TRUE) stop("names of z and pc must match")
      z <- z[npc]
    } 

    xmin <- xlim[1]
    xlim[1] <- xlim[1] - .1 * (xlim[2]-xlim[1]) ## allow space for scale

    n.col <- 100
    if (scheme=="heat") scheme <- heat.colors(n.col+1) else 
    scheme <- gray(0:n.col/n.col)
   
    zlim <- range(pretty(z))

    ## Now want a grey or color scale up the lhs of plot
    ## first scale the y range into the z range for plotting 

    for (i in 1:length(pc)) pc[[i]][,2] <- zlim[1] + 
         (zlim[2]-zlim[1])*(pc[[i]][,2]-ylim[1])/(ylim[2]-ylim[1])
  
    ylim <- zlim
    plot(0,0,ylim=ylim,xlim=xlim,type="n",xaxt="n",bty="n",xlab="",ylab=lab,...)
    for (i in 1:length(pc)) {
      coli <- round((z[i] - zlim[1])/(zlim[2]-zlim[1])*n.col)+1    
      poly2(pc[[i]],col=scheme[coli])
    }
  
    ## now plot the scale bar...

    xmin <- min(c(axTicks(1),xlim[1]))
    dx <- (xlim[2]-xlim[1])*.05
    x0 <- xmin-2*dx
    x1 <- xmin+dx

    dy <- (ylim[2]-ylim[1])/n.col 
    poly <- matrix(c(x0,x0,x1,x1,ylim[1],ylim[1]+dy,ylim[1]+dy,ylim[1]),4,2)
    for (i in 1:n.col) {
      polygon(poly,col=scheme[i],border=NA)
      poly[,2] <- poly[,2] + dy
    }
    poly <- matrix(c(x0,x0,x1,x1,ylim[1],ylim[2],ylim[2],ylim[1]),4,2)
    polygon(poly,border="black")
  }
  par(oldpar)
} ## polys.plot

plot.mrf.smooth <- function(x,P=NULL,data=NULL,label="",se1.mult=2,se2.mult=1,
                     partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,n3=3,
                     theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,by.resids=FALSE,scheme=0,...) {
## plot method function for mrf.smooth terms, depends heavily on polys.plot, above
  if (is.null(P)) { ## get plotting information...
    if (!x$plot.me||is.null(x$xt$polys)) return(NULL) ## shouldn't or can't plot
    ## get basic plot data 
    raw <- data[x$term][[1]]
    dat <- data.frame(x=factor(names(x$xt$polys),levels=levels(x$knots)))
    names(dat) <- x$term
    x$by <- "NA"
    X <- PredictMat(x,dat)   # prediction matrix for this term
    if (is.null(xlab)) xlabel<- "" else xlabel <- xlab
    if (is.null(ylab)) ylabel <- "" else ylabel <- ylab
    return(list(X=X,scale=FALSE,se=FALSE,raw=raw,xlab=xlabel,ylab=ylabel,
             main=label))
    } else { ## do plot
      if (scheme==0) scheme <- "heat" else scheme <- "grey"
      polys.plot(x$xt$polys,trans(P$fit+shift),scheme=scheme,lab=P$main,...)
    }

} ## end plot.mrf.smooth


plot.sz.interaction <- function(x,P=NULL,data=NULL,label="",se1.mult=2,se2.mult=1,
                     partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,n3=3,
                     theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,by.resids=FALSE,scheme=0,...) {
## plot method for factor smooth interactions designed for models such as s(x) + s(fac,x) where
## the factor level dependent smooths are strictly deviations from the main effect smooth.

  nf2i <- function(nf,i) {
    k <- length(nf)
    kk <- rep(0,k)
    i <- i-1
    for (j in k:1) {
      kk[j] <- i %% nf[j] + 1
      i <- i %/% nf[j]
    }
    kk
  } ## nf2i

  if (is.null(P)) { ## get plotting info
    if (x$base$dim!=1) return(NULL) ## no method for base smooth dim > 1
    raw <- data[x$base$term][[1]]
    xx <- seq(min(raw),max(raw),length=n) # generate x sequence for prediction
    nf <- unlist(lapply(x$flev,length)) ## levels per grouping factor
    dat <- data.frame(rep(xx,prod(nf)))
    k <- length(x$flev) ## number of factors
    for (i in 1:k) {
      re <- if (i<k) prod(nf[(i+1):k]) else 1
      rs <- if (i>1) prod(nf[1:(i-1)]) else 1
      dat[,i+1] <- factor(rep(rep(x$flev[[i]],each=re*n),rs),levels=x$flev[[i]])
    }
    names(dat) <- c(x$base$term,x$fterm)  
    if (x$by!="NA") {        # deal with any by variables
      dat[[x$by]] <- rep(1,n)
    }
    X <- PredictMat(x,dat)
    if (is.null(xlab)) xlabel <- x$base$term else xlabel <- xlab
    if (is.null(ylab)) ylabel <- label else ylabel <- ylab
    return(list(X=X,scale=TRUE,se=TRUE,se.mult=se1.mult,raw=raw,xlab=xlabel,ylab=ylabel,
             main="",x=xx,n=n,nf=nf))
  } else { ## produce the plot
    nft <- prod(P$nf) ## total number of curves
    if (scheme!=1) {
      kol <- hcl.colors(nft,palette = "viridis", alpha = .33) ## CI
      lkol <- hcl.colors(nft,palette = "viridis", alpha = .66) ## mode
      tkol <- hcl.colors(nft,palette = "viridis", alpha = 1) ## label
    }  
    xlim <- range(P$x);dx <- xlim[2]-xlim[1]
    xt <- xlim[1] + (1:nft-.5)*dx/nft ## text locations
    ind <- 1:P$n; mind <- P$n:1
    
    if(is.null(ylim)) ylim <- trans(range(c(P$fit+P$se,P$fit-P$se))+shift)

    plot(P$x[ind],trans(P$fit[ind]+shift),ylim=ylim,xlab=P$xlab,ylab=P$ylab,type="n",...)

    nfac <- length(P$nf) ## number of factors
    kk <- rep(0,nfac) ## factor level index vector
    if (scheme==1) {
      for (i in 1:nft) {
        ul <- trans(P$fit[ind] + P$se[ind]+shift)
        ll <- trans(P$fit[ind] - P$se[ind]+shift)
	lines(P$x,ul,col="grey",lty=i);lines(P$x,ll,col="grey",lty=i)
        ii <- P$x < xt[i] - dx/30
	yt <- approx(P$x,P$fit[ind],xt[i])$y
        lines(P$x[ii],(P$fit[ind])[ii],lty=i,lwd=2)
        text(xt[i],yt,paste(nf2i(P$nf,i),collapse="."))
        ii <- P$x > xt[i] + dx/30
        lines(P$x[ii],(P$fit[ind])[ii],lty=i,lwd=2)
        ind <- ind + P$n; mind <- mind + P$n
      }	
    } else {
      for (i in 1:nft) {
        ul <- trans(P$fit[ind] + P$se[ind]+shift)
        ll <- trans(P$fit[mind] - P$se[mind]+shift)
        polygon(c(P$x,P$x[P$n:1]),c(ul,ll),col=kol[i],border=kol[i])
        yt <- approx(P$x,P$fit[ind],xt[i])$y
        ii <- P$x < xt[i] - dx/30
        lines(P$x[ii],(P$fit[ind])[ii],col=lkol[i])
        text(xt[i],yt,paste(nf2i(P$nf,i),collapse="."),col=tkol[i])
        ii <- P$x > xt[i] + dx/30
        lines(P$x[ii],(P$fit[ind])[ii],col=lkol[i])
        ind <- ind + P$n; mind <- mind + P$n
      }
    }  
  }
} ## end plot.sz.interaction


plot.fs.interaction <- function(x,P=NULL,data=NULL,label="",se1.mult=2,se2.mult=1,
                     partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,n3=3,
                     theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,by.resids=FALSE,scheme=0,...) {
## plot method for simple smooth factor interactions...
  if (is.null(P)) { ## get plotting info
    if (x$dim!=1) return(NULL) ## no method for base smooth dim > 1
    raw <- data[x$base$term][[1]]
    xx <- seq(min(raw),max(raw),length=n) # generate x sequence for prediction
    nf <- length(x$flev)
    fac <- rep(x$flev,rep(n,nf))
    dat <- data.frame(fac,xx,stringsAsFactors=TRUE)
    names(dat) <- c(x$fterm,x$base$term)  
    if (x$by!="NA") {        # deal with any by variables
      dat[[x$by]] <- rep(1,n)
    }
    X <- PredictMat(x,dat)
    if (is.null(xlab)) xlabel <- x$base$term else xlabel <- xlab
    if (is.null(ylab)) ylabel <- label else ylabel <- ylab
    return(list(X=X,scale=TRUE,se=FALSE,raw=raw,xlab=xlabel,ylab=ylabel,
             main="",x=xx,n=n,nf=nf))
  } else { ## produce the plot
    ind <- 1:P$n
    if(is.null(ylim)) ylim <- trans(range(P$fit)+shift) 
    plot(P$x[ind],trans(P$fit[ind]+shift),ylim=ylim,xlab=P$xlab,ylab=P$ylab,type="l",...)
    if (P$nf>1) for (i in 2:P$nf) {
      ind <- ind + P$n
      if (scheme==0) lines(P$x,trans(P$fit[ind]+shift),lty=i,col=i) else 
      lines(P$x,trans(P$fit[ind]+shift),lty=i)
    }
  }
} ## end plot.fs.interaction

plot.mgcv.smooth <- function(x,P=NULL,data=NULL,label="",se1.mult=2,se2.mult=1,
                     partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,n3=3,
                     theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,by.resids=FALSE,scheme=0,hcolors=heat.colors(50),
                     contour.col=4,...) {
## default plot method for smooth objects `x' inheriting from "mgcv.smooth"
## `x' is a smooth object, usually part of a `gam' fit. It has an attribute
##     'coefficients' containing the coefs for the smooth, but usually these
##     are not needed.
## Usually this function is called twice. First to set up P, then to compute the
## actual plot information including standard error bands, and then to actually 
## plot... 
## `P' is a list of plot data. 
##     If `P' is NULL (first call) then the routine should compute some of this plot data
##     and return without plotting...  
##     * X the matrix mapping the smooth's coefficients to the values at
##         which the smooth must be computed for plotting.
##     * The values against which to plot.
##     * `exclude' indicates rows of X%*%p to set to NA for plotting -- NULL for none.
##     * se TRUE if plotting of the term can use standard error information.
##     * se.mult - the multiplier of the standard error used to plot confidence bands
##     * scale TRUE if the term should be considered by plot.gam if a common
##             y scale is required.
##     * any raw data information.
##     * axis labels and plot titles 
##     As an alternative, P may contain a 'fit' field directly, in which case the 
##     very little processing is done outside the routine, except for partial residual
##     computations.
##     Alternatively return P as NULL if x should not be plotted.
##     If P is not NULL (second call) it will contain the following...
##     * fit - the values for plotting 
##     * se - standard errors of fit multiplied by se.mult
##     * the values against which to plot
##     * any raw data information
##     * any partial.residuals 
## `data' is a data frame containing the raw data for the smooth, usually the 
##        model.frame of the fitted gam. Can be NULL if P is not NULL.
## `label' is the term label, usually something like e.g. `s(x,12.34)'.
## Note that if ylim is supplied it should not be transformed using trans and shift.
#############################

  sp.contour <- function(x,y,z,zse,xlab="",ylab="",zlab="",titleOnly=FALSE,
               se.plot=TRUE,se.mult=1,trans=I,shift=0,...)   
  ## function for contouring 2-d smooths with 1 s.e. limits
  { gap<-median(zse,na.rm=TRUE)  
    zr<-max(trans(z+zse+shift),na.rm=TRUE)-min(trans(z-zse+shift),na.rm=TRUE) # plotting range  
    n<-10  
    while (n>1 && zr/n<2.5*gap) n<-n-1    
    zrange<-c(min(trans(z-zse+shift),na.rm=TRUE),max(trans(z+zse+shift),na.rm=TRUE))  
    zlev<-pretty(zrange,n)  ## ignore codetools on this one  
    yrange<-range(y);yr<-yrange[2]-yrange[1]  
    xrange<-range(x);xr<-xrange[2]-xrange[1]  
    ypos<-yrange[2]+yr/10
    args <- as.list(substitute(list(...)))[-1]
    args$x <- substitute(x);args$y <- substitute(y)
    args$type="n";args$xlab<-args$ylab<-"";args$axes<-FALSE
    do.call("plot",args)

    cs<-(yr/10)/strheight(zlab);if (cs>1) cs<-1 # text scaling based on height  
  
    tl<-strwidth(zlab);  
    if (tl*cs>3*xr/10) cs<-(3*xr/10)/tl  
    args <- as.list(substitute(list(...)))[-1]
    n.args <- names(args)
    zz <- trans(z+shift) ## ignore codetools for this
    args$x<-substitute(x);args$y<-substitute(y);args$z<-substitute(zz)
    if (!"levels"%in%n.args) args$levels<-substitute(zlev)
    if (!"lwd"%in%n.args) args$lwd<-2
    if (!"labcex"%in%n.args) args$labcex<-cs*.65
    if (!"axes"%in%n.args) args$axes <- FALSE
    if (!"add"%in%n.args) args$add <- TRUE
    do.call("contour",args)
  
    if (is.null(args$cex.main)) cm <- 1 else cm <- args$cex.main
    if (titleOnly)  title(zlab,cex.main=cm) else 
    { xpos<-xrange[1]+3*xr/10  
      xl<-c(xpos,xpos+xr/10); yl<-c(ypos,ypos)   
      lines(xl,yl,xpd=TRUE,lwd=args$lwd)  
      text(xpos+xr/10,ypos,zlab,xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)  
    }
    if  (is.null(args$cex.axis)) cma <- 1 else cma <- args$cex.axis
    axis(1,cex.axis=cs*cma);axis(2,cex.axis=cs*cma);box();
    if  (is.null(args$cex.lab)) cma <- 1 else cma <- args$cex.lab  
    mtext(xlab,1,2.5,cex=cs*cma);mtext(ylab,2,2.5,cex=cs*cma)  
    if (!"lwd"%in%n.args) args$lwd<-1
    if (!"lty"%in%n.args) args$lty<-2
    if (!"col"%in%n.args) args$col<-2
    if (!"labcex"%in%n.args) args$labcex<-cs*.5
    zz <- trans(z+zse+shift)
    args$z<-substitute(zz)

    do.call("contour",args)

    if (!titleOnly) {
      xpos<-xrange[1]  
      xl<-c(xpos,xpos+xr/10)#;yl<-c(ypos,ypos)  
      lines(xl,yl,xpd=TRUE,lty=args$lty,col=args$col)  
      text(xpos+xr/10,ypos,paste("-",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)  
    }

    if (!"lty"%in%n.args) args$lty<-3
    if (!"col"%in%n.args) args$col<-3
    zz <- trans(z - zse+shift)
    args$z<-substitute(zz)
    do.call("contour",args)
    
    if (!titleOnly) {
      xpos<-xrange[2]-xr/5  
      xl<-c(xpos,xpos+xr/10);  
      lines(xl,yl,xpd=TRUE,lty=args$lty,col=args$col)  
      text(xpos+xr/10,ypos,paste("+",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)  
    }
  }  ## end of sp.contour

  if (is.null(P)) { ## get plotting information...
    if (!x$plot.me||x$dim>4) return(NULL) ## shouldn't or can't plot
    if (x$dim==1) { ## get basic plotting data for 1D terms 
      raw <- data[x$term][[1]]
      if (is.null(xlim)) xx <- seq(min(raw),max(raw),length=n) else # generate x sequence for prediction
      xx <- seq(xlim[1],xlim[2],length=n)
      if (x$by!="NA")         # deal with any by variables
      { by<-rep(1,n);dat<-data.frame(x=xx,by=by)
        names(dat)<-c(x$term,x$by)
      } else { 
        dat<-data.frame(x=xx);names(dat) <- x$term
      } ## prediction data.frame finished
      X <- PredictMat(x,dat)   # prediction matrix for this term
      if (is.null(xlab)) xlabel<- x$term else xlabel <- xlab
      if (is.null(ylab)) ylabel <- label else ylabel <- ylab
      if (is.null(xlim)) xlim <- range(xx)
      return(list(X=X,x=xx,scale=TRUE,se=TRUE,raw=raw,xlab=xlabel,ylab=ylabel,
             main=main,se.mult=se1.mult,xlim=xlim))
    } else if (x$dim==2) { ## basic plot data for 2D terms
      xterm <- x$term[1]
      if (is.null(xlab)) xlabel <- xterm else xlabel <- xlab
      yterm <- x$term[2]
      if (is.null(ylab)) ylabel <- yterm else ylabel <- ylab
      raw <- data.frame(x=as.numeric(data[xterm][[1]]),
                        y=as.numeric(data[yterm][[1]]))
      n2 <- max(10,n2)
      if (is.null(xlim)) xm <- seq(min(raw$x),max(raw$x),length=n2) else 
        xm <- seq(xlim[1],xlim[2],length=n2)
      if (is.null(ylim)) ym <- seq(min(raw$y),max(raw$y),length=n2) else
        ym <- seq(ylim[1],ylim[2],length=n2)
      xx <- rep(xm,n2)
      yy <- rep(ym,rep(n2,n2))
      if (too.far>0)
      exclude <- exclude.too.far(xx,yy,raw$x,raw$y,dist=too.far) else
      exclude <- rep(FALSE,n2*n2)
      if (x$by!="NA") {        # deal with any by variables
        by <- rep(1,n2^2);dat <- data.frame(x=xx,y=yy,by=by)
        names(dat) <- c(xterm,yterm,x$by)
      } else { 
        dat<-data.frame(x=xx,y=yy);names(dat)<-c(xterm,yterm)
      }  ## prediction data.frame complete
      X <- PredictMat(x,dat)   ## prediction matrix for this term
      if (is.null(main)) { 
        main <- label
      }
      if (is.null(ylim)) ylim <- range(ym) 
      if (is.null(xlim)) xlim <- range(xm) 
      return(list(X=X,x=xm,y=ym,scale=FALSE,se=TRUE,raw=raw,xlab=xlabel,ylab=ylabel,
             main=main,se.mult=se2.mult,ylim=ylim,xlim=xlim,exclude=exclude))
    } else { ## basic plot data for 3 or 4 d terms
      vname <- x$term
      ## if the smooth has margins and one is 2D then set that as the 
      ## term for 2D plotting, rather than conditioning....
      if (!is.null(x$margin)) {
        for (i in 1:length(x$margin)) if (x$margin[[i]]$dim==2) {
          vname <- x$margin[[i]]$term ## these are the variables to 2d plot
          vname <- c(vname,x$term[!x$term%in%vname])
          break;
        }
      }
      ## ... so first 2 terms in vname are the vars to plot in 2D. 
      ## Now get the limits for plotting...
      nv <- length(vname)
      lo <- hi <- rep(0,nv)
      for (i in 1:length(vname)) {
        xx <- data[vname[i]][[1]] 
        lo[i] <- min(xx);hi[i] <- max(xx)
      } 
      nc <- nr <- n3 ## set number cols and rows of plot
      m <- n2 ## 2d plotting grid side
      x1 <- seq(lo[1],hi[1],length=m)
      x2 <- seq(lo[2],hi[2],length=m)
      if (nv==3) {
        x3 <- seq(lo[3],hi[3],length=nr*nc)
        dat <- cbind(rep(x1,m*nr*nc),
           rep(rep(x2,each=m*nr),nc),
           x3[rep(rep((1:nr-1)*nc,each=m),m*nc) + rep(1:nc,each=m*m*nr)])
      } else {
        x3 <- seq(lo[3],hi[3],length=nr)
        x4 <- seq(lo[4],hi[4],length=nc)
        dat <- cbind(rep(x1,m*nr*nc),
             rep(rep(x2,each=m*nr),nc),
             rep(rep(x3,each=m),m*nc),
             rep(x4,each=m*m*nr))
      } ## 4D term end
      if (x$by!="NA") {
        dat <- data.frame(cbind(dat,1))
        names(dat) <- c(vname,x$by)
      } else {
        dat <- data.frame(dat)
        names(dat) <- vname
      }
      X <- PredictMat(x,dat)   ## prediction matrix for this term
      exclude <- if (too.far<=0) rep(FALSE,nrow(X)) else
      exclude.too.far(dat[,1],dat[,2],data[vname[1]][[1]],data[vname[2]][[1]],dist=too.far)
      if (is.null(main)) { 
        main <- label
      }
      return(list(X=X,scale=FALSE,se=FALSE,m=m,nc=nc,nr=nr,lo=lo,hi=hi,vname=vname,
             main=main,exclude=exclude))
    } ## end of 3/4 D case
  } else { ## produce plot
    if (se) { ## produce CI's
      if (x$dim==1) { 
        if (scheme == 1) shade <- TRUE
        ul <- P$fit + P$se ## upper CL
        ll <- P$fit - P$se ## lower CL
        if (scale==0&&is.null(ylim)) { ## get scale 
          ylimit<-c(min(ll),max(ul))
          if (partial.resids) { 
            max.r <- max(P$p.resid,na.rm=TRUE)
            if ( max.r> ylimit[2]) ylimit[2] <- max.r
            min.r <-  min(P$p.resid,na.rm=TRUE)
            if (min.r < ylimit[1]) ylimit[1] <- min.r
          }
        }
        ylimit <- if (is.null(ylim)) ylimit <- trans(ylimit + shift) else ylim
         
        ## plot the smooth... 
        if (shade) { 
          plot(P$x,trans(P$fit+shift),type="n",xlab=P$xlab,ylim=ylimit,
                 xlim=P$xlim,ylab=P$ylab,main=P$main,...)
          polygon(c(P$x,P$x[n:1],P$x[1]),
                    trans(c(ul,ll[n:1],ul[1])+shift),col = shade.col,border = NA)
          lines(P$x,trans(P$fit+shift),...)
        } else { ## ordinary plot 
          plot(P$x,trans(P$fit+shift),type="l",xlab=P$xlab,ylim=ylimit,xlim=P$xlim,
                 ylab=P$ylab,main=P$main,...)
          if (is.null(list(...)[["lty"]])) { 
            lines(P$x,trans(ul+shift),lty=2,...)
            lines(P$x,trans(ll+shift),lty=2,...)
          } else { 
            lines(P$x,trans(ul+shift),...)
            lines(P$x,trans(ll+shift),...)
          }
        } ## ... smooth plotted
       
        if (partial.resids&&(by.resids||x$by=="NA")) { ## add any partial residuals
          if (length(P$raw)==length(P$p.resid)) {
            if (is.null(list(...)[["pch"]]))
            points(P$raw,trans(P$p.resid+shift),pch=".",...) else
            points(P$raw,trans(P$p.resid+shift),...) 
          } else {
            warning("Partial residuals do not have a natural x-axis location for linear functional terms")
          }
        } ## partial residuals finished 
	 
        if (rug) { 
          if (jit) rug(jitter(as.numeric(P$raw)),...)
          else rug(as.numeric(P$raw),...)
	} ## rug plot done

      } else if (x$dim==2) { 
        P$fit[P$exclude] <- NA
        if (scheme == 1) { ## perspective plot 
          persp(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                  zlab=P$main,ylim=P$ylim,xlim=P$xlim,theta=theta,phi=phi,...)
        } else if (scheme==2||scheme==3) {
          if (scheme==3) hcolors <- grey(0:50/50)
          image(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                  main=P$main,xlim=P$xlim,ylim=P$ylim,col=hcolors,...)
          contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),add=TRUE,col=contour.col,...)
          if (rug) {  
            if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
            points(P$raw$x,P$raw$y,...)
          }
        } else { ## contour plot with error contours
          sp.contour(P$x,P$y,matrix(P$fit,n2,n2),matrix(P$se,n2,n2),
                     xlab=P$xlab,ylab=P$ylab,zlab=P$main,titleOnly=!is.null(main),
                     se.mult=1,trans=trans,shift=shift,...)
          if (rug) { 
            if (is.null(list(...)[["pch"]]))
            points(P$raw$x,P$raw$y,pch=".",...) else
            points(P$raw$x,P$raw$y,...) 
          }
        } ## contour plot done 
      } else if (x$dim<5) {
        if (scheme==1) hcolors <- grey(0:50/50)
        md.plot(P$fit,P$nr,P$nc,P$m,P$vname,P$lo,P$hi,hcolors=hcolors,scheme=scheme,P$main,...)
      } else { 
         warning("no automatic plotting for smooths of more than two variables")
      }
    } else { ## no CI's
      if (x$dim==1) { 
        if (scale==0&&is.null(ylim)) { 
          if (partial.resids) ylimit <- range(P$p.resid,na.rm=TRUE) else ylimit <-range(P$fit)
        }
        ylimit <- if (is.null(ylim)) ylimit <- trans(ylimit + shift) else ylim
        
        plot(P$x,trans(P$fit+shift),type="l",xlab=P$xlab,
             ylab=P$ylab,ylim=ylimit,xlim=P$xlim,main=P$main,...)
        if (rug) { 
          if (jit) rug(jitter(as.numeric(P$raw)),...)
          else rug(as.numeric(P$raw),...) 
        }
        if (partial.resids&&(by.resids||x$by=="NA")) { 
          if (is.null(list(...)[["pch"]]))
          points(P$raw,trans(P$p.resid+shift),pch=".",...) else
          points(P$raw,trans(P$p.resid+shift),...)
        }
      } else if (x$dim==2) { 
        P$fit[P$exclude] <- NA
        if (!is.null(main)) P$title <- main
        if (scheme==1) { 
          persp(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                          zlab=P$main,theta=theta,phi=phi,xlim=P$xlim,ylim=P$ylim,...)
        } else if (scheme==2||scheme==3) {
          if (scheme==3) hcolors <- grey(0:50/50)
          image(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                  main=P$main,xlim=P$xlim,ylim=P$ylim,col=hcolors,...)
          contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),add=TRUE,col=contour.col,...)
          if (rug) {  
            if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
            points(P$raw$x,P$raw$y,...)
          }
        } else { 
          contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                  main=P$main,xlim=P$xlim,ylim=P$ylim,...)
          if (rug) {  
            if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
            points(P$raw$x,P$raw$y,...)
          }
        }  
      } else if (x$dim<5) {
        if (scheme==1) hcolors <- grey(0:50/50)
        md.plot(P$fit,P$nr,P$nc,P$m,P$vname,P$lo,P$hi,hcolors=hcolors,scheme=scheme,P$main,...)
      } else { 
        warning("no automatic plotting for smooths of more than four variables")
      }
    } ## end of no CI code
  } ## end of plot production
} ## plot.mgcv.smooth

md.plot <- function(f,nr,nc,m,vname,lo,hi,hcolors,scheme,main,...) {
## multi-dimensional term plotter, called from plot.mgcv.smooth for
## 3 and 4 dimensional terms. 
## *f is the plot data. See `basic plot data for 3 or 4 d terms'
##   in plot.mgcv.smooth for details of the packing conventions
##   (f = X %*% coefs).
## *nr and nc the number of rows and columns of plot panels
## *m each panel is m by m
## *vname contains the variable names
## *lo and hi are the arrays of axis limits
## *hcolors is the color palette for the image plot.
## *scheme indicates b/w or color
## *main is a title.
  concol <- if (scheme==1) "white" else "black" 
  nv <- length(vname) 
  ## insert NA breaks to separate the panels within a plot...
  f1 <- matrix(NA,nr*m+nr-1,nc*m)
  ii <- rep(1:m,nr) + rep(0:(nr-1)*(m+1),each=m)
  f1[ii,] <- f
  f <- matrix(NA,nr*m+nr-1,nc*m+nc-1)
  ii <- rep(1:m,nc) + rep(0:(nc-1)*(m+1),each=m)
  f[,ii] <- f1
  xx <- seq(0,1,length=ncol(f))
  yy <- seq(0,1,length=nrow(f))
  image(xx,yy,t(f),axes=FALSE,xlab="",ylab="",col=hcolors)
  contour(xx,yy,t(f),add=TRUE,col=concol)
  dl <- list(...)
  c1 <- if (is.null(dl[["cex"]])) 1 else dl[["cex"]] 
  c2 <- if (is.null(dl[["cex.axis"]])) .6 else dl[["cex.axis"]]
  c3 <- if (is.null(dl[["cex.lab"]])) .9 else dl[["cex.lab"]]
  if (nv==4) { 
    x3 <- seq(lo[3],hi[3],length=nr)
    x4 <- seq(lo[4],hi[4],length=nc)
    mtext(vname[4],1,1.7,cex=c1*c3) ## x label
    mtext(vname[3],2,1.7,cex=c1*c3) ## y label
    at=(1:nc-.5)/nc
    lab <- format(x4,digits=2)
    for (i in 1:nc) mtext(lab[i],1,at=at[i],line=.5,cex=c1*c3)
    at=(1:nr-.5)/nr
    lab <- format(x4,digits=2)
    for (i in 1:nr) mtext(lab[i],2,at=at[i],line=.5,cex=c1*c3)
    ## now the 2d panel axes...
    xr <- axisTicks(c(lo[2],hi[2]),log=FALSE,nint=4)
    x0 <- ((nc-1)*(m+1)+1)/(nc*m+nc-1)
    xt <- (xr-lo[2])/(hi[2]-lo[2])*(1-x0)+x0
    axis(3,at=xt,labels=as.character(xr),cex.axis=c2,cex=c1)
    xr <- axisTicks(c(lo[1],hi[1]),log=FALSE,nint=4)
    x0 <- ((nr-1)*(m+1)+1)/(nr*m+nr-1)
    xt <- (xr-lo[1])/(hi[1]-lo[1])*(1-x0)+x0
    axis(4,at=xt,labels=as.character(xr),cex.axis=c2,cex=c1)
    at <- (2*nc-3)/(2*nc) 
    mtext(vname[2],3,at=at,line=.5,cex=c1*c2)
    at <- (2*nr-3)/(2*nr) 
    mtext(vname[1],4,at=at,line=.5,cex=c1*c2)
    mtext(main,3,at=0,adj=0,line=1,cex=c1*c3)
  } else {
    x3 <- seq(lo[3],hi[3],length=nr*nc)
    ## get pretty ticks
    xr <- axisTicks(c(lo[2],hi[2]),log=FALSE,nint=4)
    x0 <- (m-1)/(nc*m+nc-1)
    xt <- (xr-lo[2])/(hi[2]-lo[2])*x0
    axis(1,at=xt,labels=as.character(xr),cex.axis=c2,cex=c1)
    mtext(vname[2],1,at=x0/2,line=2,cex=c1*c2)
    xr <- axisTicks(c(lo[1],hi[1]),log=FALSE,nint=4)
    x0 <- (m-1)/(nr*m+nr-1)
    xt <- (xr-lo[1])/(hi[1]-lo[1])*x0
    axis(2,at=xt,labels=as.character(xr),cex.axis=c2,cex=c1)
    mtext(vname[1],2,at=x0/2,line=2,cex=c1*c2)
    lab <- c("",format(x3[-1],digits=2))
    at=(1:nc-.5)/nc
    for (i in 2:nc) mtext(lab[i],1,at=at[i],line=.5,cex=c1*c3)
    mtext(parse(text=paste(vname[3],"%->% \" \"")),1,at=mean(at[2:nc]),line=2,cex=c1*c3)
    ii <- ((nr-1)*nr+1):(nc*nr)
    for (i in 1:nc) mtext(lab[ii[i]],3,at=at[i],line=.5,cex=c1*c3)
    mtext(parse(text=paste(vname[3],"%->% \" \"")),3,at=mean(at),line=2,cex=c1*c3)
    mtext(main,2,at=1/nr+0.5*(nr-1)/nr,line=1,cex=c1*c3)
  }
} ## md.plot

plot.gam <- function(x,residuals=FALSE,rug=NULL,se=TRUE,pages=0,select=NULL,scale=-1,n=100,n2=40,n3=3,
                     theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,all.terms=FALSE,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,seWithMean=FALSE,unconditional=FALSE,by.resids=FALSE,scheme=0,...)

# Create an appropriate plot for each smooth term of a GAM.....
# x is a gam object
# rug determines whether a rug plot should be added to each plot
# se determines whether twice standard error bars are to be added
# pages is the number of pages over which to split output - 0 implies that 
# graphic settings should not be changed for plotting
# scale -1 for same y scale for each plot
#        0 for different y scales for each plot
# n - number of x axis points to use for plotting each term
# n2 is the square root of the number of grid points to use for contouring
# 2-d terms.

{ ######################################
  ## Local function for producing labels
  ######################################

  sub.edf <- function(lab,edf) {
    ## local function to substitute edf into brackets of label
    ## labels are e.g. smooth[[1]]$label
    pos <- regexpr(":",lab)[1]
    if (pos<0) { ## there is no by variable stuff
      pos <- nchar(lab) - 1
      lab <- paste(substr(lab,start=1,stop=pos),",",round(edf,digits=2),")",sep="")
    } else {
      lab1 <- substr(lab,start=1,stop=pos-2)
      lab2 <- substr(lab,start=pos-1,stop=nchar(lab))
      lab <- paste(lab1,",",round(edf,digits=2),lab2,sep="")
    }
    lab
  } ## end of sub.edf


  #########################
  ## start of main function
  #########################

  if (is.null(rug)) rug <- if (nrow(x$model)>10000) FALSE else TRUE

  if (unconditional) {
    if (is.null(x$Vc)) warning("Smoothness uncertainty corrected covariance not available") else 
    x$Vp <- x$Vc ## cov matrix reset to full Bayesian
  }

  w.resid <- NULL
  if (length(residuals)>1) { # residuals supplied 
    if (length(residuals)==length(x$residuals)) 
    w.resid <- residuals else
    warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  } else partial.resids <- residuals # use working residuals or none

  m <- length(x$smooth) ## number of smooth terms

  if (length(scheme)==1) scheme <- rep(scheme,m)
  if (length(scheme)!=m) { 
    warn <- paste("scheme should be a single number, or a vector with",m,"elements")
    warning(warn)
    scheme <- rep(scheme[1],m)
  }

  ## array giving order of each parametric term...
  order <- if (is.list(x$pterms))  unlist(lapply(x$pterms,attr,"order")) else attr(x$pterms,"order")

  if (all.terms) # plot parametric terms as well
  n.para <- sum(order==1) # plotable parametric terms   
  else n.para <- 0 
 
  if (se) ## sort out CI widths for 1 and 2D
  { if (is.numeric(se)) se2.mult <- se1.mult <- se else { se1.mult <- 2;se2.mult <- 1} 
    if (se1.mult<0) se1.mult<-0;if (se2.mult < 0) se2.mult <- 0
  } else se1.mult <- se2.mult <-1
  
  if (se && x$Vp[1,1] < 0) ## check that variances are actually available
  { se <- FALSE
    warning("No variance estimates available")
  }

  if (partial.resids) { ## getting information needed for partial residuals...
    if (is.null(w.resid)) { ## produce working resids if info available
      if (is.null(x$residuals)||is.null(x$weights)) partial.resids <- FALSE else {
        wr <- sqrt(x$weights)
        w.resid <- x$residuals*wr#/mean(wr) # weighted working residuals
      }
    }
    if (partial.resids) fv.terms <- predict(x,type="terms") ## get individual smooth effects
  }

  pd <- list(); ## plot data list
  i <- 1 # needs a value if no smooths, but parametric terms ...

  ##################################################
  ## First the loop to get the data for the plots...
  ##################################################

  if (m>0) for (i in 1:m) { ## work through smooth terms
    first <- x$smooth[[i]]$first.para
    last <- x$smooth[[i]]$last.para
    edf <- sum(x$edf[first:last]) ## Effective DoF for this term
    term.lab <- sub.edf(x$smooth[[i]]$label,edf)
    #P <- plot(x$smooth[[i]],P=NULL,data=x$model,n=n,n2=n2,xlab=xlab,ylab=ylab,too.far=too.far,label=term.lab,
    #          se1.mult=se1.mult,se2.mult=se2.mult,xlim=xlim,ylim=ylim,main=main,scheme=scheme[i],...)
    attr(x$smooth[[i]],"coefficients") <- x$coefficients[first:last]   ## relevent coefficients
    P <- plot(x$smooth[[i]],P=NULL,data=x$model,partial.resids=partial.resids,rug=rug,se=se,scale=scale,n=n,n2=n2,n3=n3,
                     theta=theta,phi=phi,jit=jit,xlab=xlab,ylab=ylab,main=main,label=term.lab,
                     ylim=ylim,xlim=xlim,too.far=too.far,shade=shade,shade.col=shade.col,
                     se1.mult=se1.mult,se2.mult=se2.mult,shift=shift,trans=trans,
                     by.resids=by.resids,scheme=scheme[i],...)

    if (is.null(P)) pd[[i]] <- list(plot.me=FALSE) else if (is.null(P$fit)) {
      p <- x$coefficients[first:last]   ## relevent coefficients 
      offset <- attr(P$X,"offset")      ## any term specific offset
      ## get fitted values ....
      if (is.null(offset)) P$fit <- P$X%*%p else P$fit <- P$X%*%p + offset 
      if (!is.null(P$exclude)) P$fit[P$exclude] <- NA
      if (se && P$se) { ## get standard errors for fit
        ## test whether mean variability to be added to variability (only for centred terms)
        if (seWithMean && attr(x$smooth[[i]],"nCons")>0) {
          if (length(x$cmX) < ncol(x$Vp)) x$cmX <- c(x$cmX,rep(0,ncol(x$Vp)-length(x$cmX)))
          if (seWithMean==2) x$cmX[-(1:x$nsdf)] <- 0 ## variability of fixed effects mean only
          X1 <- matrix(x$cmX,nrow(P$X),ncol(x$Vp),byrow=TRUE)
          meanL1 <- x$smooth[[i]]$meanL1
          if (!is.null(meanL1)) X1 <- X1 / meanL1
          X1[,first:last] <- P$X
          lpi <- attr(x$formula,"lpi")
          if (is.null(lpi)) se.fit <- sqrt(pmax(0,rowSums(as(X1%*%x$Vp,"matrix")*X1))) else {
            ii <- rep(0,0) ## only include constant uncertainty from relevant linear predictors  
            for (q in 1:length(lpi)) if (any(first:last%in%lpi[[q]])) ii <- c(ii,lpi[[q]])
            se.fit <- sqrt(pmax(0,rowSums(as(X1[,ii]%*%x$Vp[ii,ii],"matrix")*X1[,ii])))
          }
        } else se.fit <- ## se in centred (or anyway unconstained) space only
        sqrt(pmax(0,rowSums(as(P$X%*%x$Vp[first:last,first:last,drop=FALSE],"matrix")*P$X)))
        if (!is.null(P$exclude)) se.fit[P$exclude] <- NA
      } ## standard errors for fit completed
      if (partial.resids) { P$p.resid <- fv.terms[,length(order)+i] + w.resid }
      if (se && P$se) P$se <- se.fit*P$se.mult  # Note multiplier
      P$X <- NULL
      P$plot.me <- TRUE
      pd[[i]] <- P;rm(P) 
    } else { ## P$fit created directly
      if (partial.resids) { P$p.resid <- fv.terms[,length(order)+i] + w.resid }
      P$plot.me <- TRUE
      pd[[i]] <- P;rm(P)
    }
  } ## end of data setup loop through smooths

  
  ##############################################
  ## sort out number of pages and plots per page 
  ##############################################

  n.plots <- n.para
  if (m>0) for (i in 1:m) n.plots <- n.plots + as.numeric(pd[[i]]$plot.me) 

  if (n.plots==0) stop("No terms to plot - nothing for plot.gam() to do.")

  if (pages>n.plots) pages<-n.plots
  if (pages<0) pages<-0
  if (pages!=0)    # figure out how to display things
  { ppp<-n.plots%/%pages
    if (n.plots%%pages!=0) 
    { ppp<-ppp+1
      while (ppp*(pages-1)>=n.plots) pages<-pages-1
    } 

    # now figure out number of rows and columns
    c <- r <- trunc(sqrt(ppp))
    if (c<1) r <- c <- 1
    if (c*r < ppp) c <- c + 1
    if (c*r < ppp) r <- r + 1  
    oldpar<-par(mfrow=c(r,c))
  
  } else
  { ppp<-1;oldpar<-par()}
  

  #####################################
  ## get a common scale, if required...
  #####################################

  if (scale==-1&&is.null(ylim)) {
    k <- 0
    if (m>0) for (i in 1:m) if (pd[[i]]$plot.me&&pd[[i]]$scale) { ## loop through plot data 
      if (se&&length(pd[[i]]$se)>1) { ## require CIs on plots
        ul<-pd[[i]]$fit+pd[[i]]$se
        ll<-pd[[i]]$fit-pd[[i]]$se
        if (k==0) { 
          ylim <- c(min(ll,na.rm=TRUE),max(ul,na.rm=TRUE));k <- 1
        } else {
          if (min(ll,na.rm=TRUE)<ylim[1]) ylim[1] <- min(ll,na.rm=TRUE)
	  if (max(ul,na.rm=TRUE)>ylim[2]) ylim[2] <- max(ul,na.rm=TRUE)
        }
      } else { ## no standard errors
        if (k==0) {
          ylim <- range(pd[[i]]$fit,na.rm=TRUE);k <- 1
        } else {
          if (min(pd[[i]]$fit,na.rm=TRUE)<ylim[1]) ylim[1] <- min(pd[[i]]$fit,na.rm=TRUE)
          if (max(pd[[i]]$fit,na.rm=TRUE)>ylim[2]) ylim[2] <- max(pd[[i]]$fit,na.rm=TRUE)
        }
      }
      if (partial.resids) { 
        ul <- max(pd[[i]]$p.resid,na.rm=TRUE)
        if (ul > ylim[2]) ylim[2] <- ul
        ll <-  min(pd[[i]]$p.resid,na.rm=TRUE)
        if (ll < ylim[1]) ylim[1] <- ll
      } ## partial resids done
    } ## loop end 
    ylim <- trans(ylim+shift)
  } ## end of common scale computation
  
  ##############################################################
  ## now plot smooths, by calling plot methods with plot data...
  ##############################################################

  if ((pages==0&&prod(par("mfcol"))<n.plots&&dev.interactive())||
       pages>1&&dev.interactive()) ask <- TRUE else ask <- FALSE 
  
  if (!is.null(select)) {
    ask <- FALSE
  }
 
#  if (ask) { ## asks before plotting
#    oask <- devAskNewPage(TRUE)
#    on.exit(devAskNewPage(oask))
#  }

  if (m>0) for (i in 1:m) if (pd[[i]]$plot.me&&(is.null(select)||i==select)) {
    plot(x$smooth[[i]],P=pd[[i]],partial.resids=partial.resids,rug=rug,se=se,scale=scale,n=n,n2=n2,n3=n3,
                     theta=theta,phi=phi,jit=jit,xlab=xlab,ylab=ylab,main=main,
                     ylim=ylim,xlim=xlim,too.far=too.far,shade=shade,shade.col=shade.col,
                     shift=shift,trans=trans,by.resids=by.resids,scheme=scheme[i],...)
   if (ask) { ## this is within loop so we don't get asked before it's necessary
     oask <- devAskNewPage(TRUE)
     on.exit(devAskNewPage(oask))
     ask <- FALSE ## only need to do this once
   }

  } ## end of smooth plotting loop
  
  ####################################################
  ## Finally deal with any parametric term plotting...
  ####################################################

  if (n.para>0) # plot parameteric terms
  { class(x) <- c("gam","glm","lm") # needed to get termplot to call model.frame.glm 
    if (is.null(select)) {
      attr(x,"para.only") <- TRUE
      termplot(x,se=se,rug=rug,col.se=1,col.term=1,main=attr(x$pterms,"term.labels"),...)
    } else { # figure out which plot is required
      if (select > m) { 
        ## can't figure out how to get this to work with more than first linear predictor
        ## as termplots relies on matching terms to names in original data... 
        select <- select - m # i.e. which parametric term
        term.labels <- attr(x$pterms,"term.labels")
        term.labels <- term.labels[order==1]
        if (select <= length(term.labels)) {
          # if (interactive() && m &&i%%ppp==0) 
          termplot(x,terms=term.labels[select],se=se,rug=rug,col.se=1,col.term=1,...)
        }  
      }
    }
  }
  if (pages>0) par(oldpar)
  invisible(pd)
} ## end plot.gam



exclude.too.far <- function(g1,g2,d1,d2,dist)
# if g1 and g2 are the co-ordinates of grid modes and d1,d2 are co-ordinates of data
# then this routine returns a vector with TRUE if the grid node is too far from
# any data and FALSE otherwise. Too far is judged using dist: a positive number indicating
# distance on the unit square into which the grid is scaled prior to calculation
{ mig<-min(g1)
  d1<-d1-mig;g1<-g1-mig
  mag<-max(g1)
  d1<-d1/mag;g1<-g1/mag
  mig<-min(g2)
  d2<-d2-mig;g2<-g2-mig
  mag<-max(g2)
  d2<-d2/mag;g2<-g2/mag
  # all now in unit square
  n<-length(g1)
  m<-length(d1)
  if (length(g2)!=n) stop("grid vectors are different lengths")
  if (m!=length(d2)) stop("data vectors are of different lengths")
  if (dist<0) stop("supplied dist negative")
  distance<-array(0,n)
  o<-.C(C_MinimumSeparation,x=as.double(cbind(g1,g2)),n=as.integer(n), d=as.integer(2),
                            t=as.double(cbind(d1,d2)),m=as.integer(m),distance=as.double(distance))
 
  res <- rep(FALSE,n)
  res[o$distance > dist] <-TRUE
  res
} ## exclude.too.far



vis.gam <- function(x,view=NULL,cond=list(),n.grid=30,too.far=0,col=NA,color="heat",
           contour.col=NULL,se=-1,type="link",plot.type="persp",zlim=NULL,nCol=50,lp=1,...)
# takes a gam object and plots 2D views of it, supply ticktype="detailed" to get proper axis anotation
# (c) Simon N. Wood 23/2/03
{ fac.seq<-function(fac,n.grid)
  # generates a sequence of factor variables of length n.grid
  { fn<-length(levels(fac));gn<-n.grid;
    if (fn>gn) mf<-factor(levels(fac))[1:gn]
    else
    { ln<-floor(gn/fn) # length of runs               
      mf<-rep(levels(fac)[fn],gn)
      mf[1:(ln*fn)]<-rep(levels(fac),rep(ln,fn))
      mf<-factor(mf,levels=levels(fac))
    }
    mf
  }
  # end of local functions

  dnm <- names(list(...))

  ## basic issues in the following are that not all objects will have a useful `data'
  ## component, but they all have a `model' frame. Furthermore, `predict.gam' recognises
  ## when a model frame has been supplied

  v.names  <- names(x$var.summary) ## names of all variables

  ## Note that in what follows matrices in the parametric part of the model
  ## require special handling. Matrices arguments to smooths are different
  ## as they follow the summation convention. 
  if (is.null(view)) # get default view if none supplied
  { ## need to find first terms that can be plotted against
    k <- 0;view <- rep("",2) 
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) ok <- FALSE else
      if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]]))<=1) ok <- FALSE 
      } else {
        if (length(unique(x$var.summary[[i]]))==1) ok <- FALSE
      }
      if (ok) {
        k <- k + 1;view[k] <- v.names[i]
      }
      if (k==2) break;
    }
    if (k<2) stop("Model does not seem to have enough terms to do anything useful")
  } else { 
    if (sum(view%in%v.names)!=2) stop(gettextf("view variables must be one of %s", 
                                      paste(v.names, collapse = ", ")))

    for (i in 1:2) 
    if  (!inherits(x$var.summary[[view[i]]],c("numeric","factor"))) 
    stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }

  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]]))<=1) ok <- FALSE
  } else {
    if (length(unique(x$var.summary[[view[i]]]))<=1) ok <- FALSE 
  }
  if (!ok) stop(gettextf("View variables must contain more than one value. view = c(%s,%s).",
                view[1], view[2]))
 
  # now get the values of the variables which are not the arguments of the plotted surface

  # Make dataframe....
  if (is.factor(x$var.summary[[view[1]]]))
  m1<-fac.seq(x$var.summary[[view[1]]],n.grid)
  else { r1<-range(x$var.summary[[view[1]]]);m1<-seq(r1[1],r1[2],length=n.grid)}
  if (is.factor(x$var.summary[[view[2]]]))
  m2<-fac.seq(x$var.summary[[view[2]]],n.grid)
  else { r2<-range(x$var.summary[[view[2]]]);m2<-seq(r2[1],r2[2],length=n.grid)}
  v1<-rep(m1,n.grid);v2<-rep(m2,rep(n.grid,n.grid))
  
  newd <- data.frame(matrix(0,n.grid*n.grid,0)) ## creating prediction data frame full of conditioning values
  for (i in 1:length(x$var.summary)) { 
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) { 
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) ma <- ma[2] ## extract median
    }
    if (is.matrix(x$var.summary[[i]])) newd[[i]] <- matrix(ma,n.grid*n.grid,ncol(x$var.summary[[i]]),byrow=TRUE)
    else newd[[i]]<-rep(ma,n.grid*n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]]<-v1
  newd[[view[2]]]<-v2
  # call predict.gam to get predictions.....
  if (type=="link") zlab <- paste("linear predictor") ## ignore codetools
  else if (type=="response") zlab <- type
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x,newdata=newd,se.fit=TRUE,type=type)
  z <- fv$fit # store NA free copy now
  if (is.matrix(z)) {
    lp <- min(ncol(z),max(1,round(lp)))
    z <- z[,lp] ## retain selected linear predictor
    fv$fit <- fv$fit[,lp]
    fv$se.fit <- fv$se.fit[,lp]
  }
  if (too.far>0) { # exclude predictions too far from data
    ex.tf <- exclude.too.far(v1,v2,x$model[,view[1]],x$model[,view[2]],dist=too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  # produce a continuous scale in place of any factors
  if (is.factor(m1)) { 
    m1<-as.numeric(m1);m1<-seq(min(m1)-0.5,max(m1)+0.5,length=n.grid) 
  }
  if (is.factor(m2)) { 
    m2<-as.numeric(m2);m2<-seq(min(m1)-0.5,max(m2)+0.5,length=n.grid) 
  }
  if (se<=0) { 
    old.warn<-options(warn=-1)
    av <- matrix(c(0.5,0.5,rep(0,n.grid-1)),n.grid,n.grid-1)
    options(old.warn)
    # z is without any exclusion of gridpoints, so that averaging works nicely
    max.z <- max(z,na.rm=TRUE)
    z[is.na(z)] <- max.z*10000 # make sure NA's don't mess it up
    z <- matrix(z,n.grid,n.grid) # convert to matrix
    surf.col <- t(av)%*%z%*%av   # average over tiles  
    surf.col[surf.col>max.z*2] <- NA # restore NA's
    # use only non-NA data to set colour limits
    if (!is.null(zlim)) { 
      if (length(zlim)!=2||zlim[1]>=zlim[2]) stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    } else { 
      min.z <- min(fv$fit,na.rm=TRUE)
      max.z <- max(fv$fit,na.rm=TRUE)
    }
    if (min.z==max.z) {min.z <- min.z-1;max.z <- max.z + 1}
    surf.col <- surf.col-min.z
    surf.col <- surf.col/(max.z-min.z)  
    surf.col <- round(surf.col*nCol)
    con.col <- 1
    if (color=="heat") { pal<-heat.colors(nCol);con.col<-4;}
    else if (color=="topo") { pal<-topo.colors(nCol);con.col<-2;}
    else if (color=="cm") { pal<-cm.colors(nCol);con.col<-1;}
    else if (color=="terrain") { pal<-terrain.colors(nCol);con.col<-2;}
    else if (color=="gray"||color=="bw") {pal <- gray(seq(0.1,0.9,length=nCol));con.col<-1}
    else stop("color scheme not recognised")
    if (is.null(contour.col)) contour.col <- con.col   # default colour scheme
    surf.col[surf.col<1] <- 1; surf.col[surf.col>nCol] <- nCol # otherwise NA tiles can get e.g. -ve index
    if (is.na(col)) col <- pal[as.array(surf.col)]
    z <- matrix(fv$fit,n.grid,n.grid)
    if (plot.type=="contour")
    { stub <- paste(ifelse("xlab" %in% dnm, "" , ",xlab=view[1]"),
                    ifelse("ylab" %in% dnm, "" , ",ylab=view[2]"),
                    ifelse("main" %in% dnm, "" , ",main=zlab"),",...)",sep="")
      if (color!="bw")
      { txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)",stub,sep="") # assemble image() call
        eval(parse(text=txt))
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)",
               ifelse("add" %in% dnm, "" , ",add=TRUE"),",...)" , sep="") # assemble contour() call
         eval(parse(text=txt))       
      } else
      { txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)",stub,sep="")  # assemble contour() call
        eval(parse(text=txt))
      }
    } else
    { stub <- paste(ifelse("xlab" %in% dnm, "" , ",xlab=view[1]"),
                    ifelse("ylab" %in% dnm, "" , ",ylab=view[2]"),
                    ifelse("zlab" %in% dnm, "" , ",zlab=zlab"),",...)",sep="")
      if (color=="bw")
      { op <- par(bg="white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ",stub,sep="") # assemble persp() call
        eval(parse(text=txt))
        par(op)
      } else
      { txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)",stub,sep="")  # assemble persp() call
        eval(parse(text=txt))
      }
    }
  } else # add standard error surfaces
  { if (color=="bw"||color=="gray") 
    { subs <- paste("grey are +/-",se,"s.e.")  ## ignore codetools
      lo.col <- "gray" ## ignore codetools claims about this
      hi.col <- "gray" ## ignore codetools 
    } else
    { subs <- paste("red/green are +/-",se,"s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim)!=2||zlim[1]>=zlim[2]) stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    } else {
      max.z <- max(fv$fit+fv$se.fit*se,na.rm=TRUE)
      min.z <- min(fv$fit-fv$se.fit*se,na.rm=TRUE)
      zlim<-c(min.z,max.z)
    } 
    z <- fv$fit - fv$se.fit*se; z <- matrix(z,n.grid,n.grid)
    if (plot.type=="contour") warning("sorry no option for contouring with errors: try plot.gam")

    stub <-  paste(ifelse("xlab" %in% dnm, "" , ",xlab=view[1]"),
                   ifelse("ylab" %in% dnm, "" , ",ylab=view[2]"),
                   ifelse("zlab" %in% dnm, "" , ",zlab=zlab"),
                   ifelse("sub" %in% dnm, "" , ",sub=subs"),
                   ",...)",sep="")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim",
                 ifelse("border" %in% dnm, "" ,",border=lo.col"),
                 stub,sep="") # assemble persp() call
    eval(parse(text=txt))

    par(new=TRUE) # don't clean device
    z <- fv$fit; z <- matrix(z,n.grid,n.grid)

    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim",
                 ifelse("border" %in% dnm, "" ,",border=\"black\""),
                 stub,sep="")
    eval(parse(text=txt))

    par(new=TRUE) # don't clean device
    z <- fv$fit+se*fv$se.fit; z <- matrix(z,n.grid,n.grid)
    
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim",
                 ifelse("border" %in% dnm, "" ,",border=hi.col"),
                 stub,sep="")
    eval(parse(text=txt))
  }
} ## vis.gam


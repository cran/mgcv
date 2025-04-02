## (c) Simon N. Wood,  Natalya Pya (scat, beta), Chris Shen (clog), 
## 2013-2025. Released under GPL2.
## See gam.fit4.r for testing functions fmud.test and fetad.test.

estimate.theta <- function(theta,family,y,mu,scale=1,wt=1,tol=1e-7,attachH=FALSE) {
## Simple Newton iteration to estimate theta for an extended family,
## given y and mu. To be iterated with estimation of mu given theta.
## If used within a PIRLS loop then divergence testing of coef update
## will have to re-compute pdev after theta update.
## Not clear best way to handle scale - could optimize here as well

  if (!inherits(family,"extended.family")) stop("not an extended family")

  nlogl <- function(theta,family,y,mu,scale=1,wt=1,deriv=2) {
    ## compute the negative log likelihood and grad + hess
    nth <- length(theta) - if (scale<0) 1 else 0
    if (scale < 0) {
      scale <- exp(theta[nth+1])
      theta <- theta[-(nth+1)]
      get.scale <- TRUE
    } else get.scale <- FALSE
    dev <- sum(family$dev.resids(y, mu, wt,theta))/scale
    if (deriv>0) Dd <- family$Dd(y, mu, theta, wt=wt, level=deriv)
    ls <- family$ls(y,w=wt,theta=theta,scale=scale)
    nll <- dev/2 - ls$ls
   
    if (deriv>0) {
      g1 <- colSums(as.matrix(Dd$Dth))/(2*scale)
      g <- if (get.scale) c(g1,-dev/2) else g1
      ind <- 1:length(g)
      g <- g - ls$lsth1[ind]
    } else g <- NULL
    if (deriv>1) {
      x <- colSums(as.matrix(Dd$Dth2))/(2*scale)
      Dth2 <- matrix(0,nth,nth)
      k <- 0
      for (i in 1:nth) for (j in i:nth) {
        k <- k + 1
        Dth2[i,j] <- Dth2[j,i] <- x[k]  
      }
      if (get.scale) Dth2 <- rbind(cbind(Dth2,-g1),c(-g1,dev/2))
      H <- Dth2 - as.matrix(ls$lsth2)[ind,ind]
    } else H <- NULL
    list(nll=nll,g=g,H=H)
  } ## nlogl

  if (scale>=0&&family$n.theta==0) stop("erroneous call to estimate.theta - no free parameters")
  n.theta <- length(theta) ## dimension of theta vector (family$n.theta==0 => all fixed)
  del.ind <- 1:n.theta
  if (scale<0) theta <- c(theta,log(var(y)*.1))
  nll <- nlogl(theta,family,y,mu,scale,wt,2)
  g <- if (family$n.theta==0) nll$g[-del.ind] else nll$g
  H <- if (family$n.theta==0) nll$H[-del.ind,-del.ind,drop=FALSE] else nll$H
  step.failed <- FALSE
  uconv <- abs(g) > tol*(abs(nll$nll)+1)
  if (sum(uconv)) for (i in 1:100) { ## main Newton loop
    eh <- eigen(H[uconv,uconv,drop=FALSE],symmetric=TRUE)
    pdef <- sum(eh$values <= 0)==0
    if (!pdef) { ## Make the Hessian pos def
      eh$values <- abs(eh$values)
      thresh <- max(eh$values) * 1e-5
      eh$values[eh$values<thresh] <- thresh
    }
    ## compute step = -solve(H,g) (via eigen decomp)...
    step0 <- - drop(eh$vectors %*% ((t(eh$vectors) %*% g[uconv])/eh$values))
    if (family$n.theta==0) step0 <- c(rep(0,n.theta),step0)
    ## limit the step length...
    ms <- max(abs(step0))
    if (ms>4) step0 <- step0*4/ms
    step <- theta*0;step[uconv] <- step0
    
    nll1 <- nlogl(theta+step,family,y,mu,scale,wt,2)
    iter <- 0
    while (nll1$nll - nll$nll > .Machine$double.eps^.75*abs(nll$nll)) { ## step halving 
      step <- step/2; iter <- iter + 1
      if (sum(theta!=theta+step)==0||iter>25) {
        step.failed <- TRUE
	break
      }
      nll1 <- nlogl(theta+step,family,y,mu,scale,wt,0)
    } ## step halving
    if (step.failed) break
    theta <- theta + step ## accept updated theta
    ## accept log lik and derivs ...
    nll <- if (iter>0) nlogl(theta,family,y,mu,scale,wt,2) else nll1
    g <- if (family$n.theta==0) nll$g[-del.ind] else nll$g
    H <- if (family$n.theta==0) nll$H[-del.ind,-del.ind,drop=FALSE] else nll$H
    ## convergence checking...
    uconv <- abs(g) > tol*(abs(nll$nll)+1)
    if (sum(uconv)==0) break 
  } ## main Newton loop
  if (step.failed) warning("step failure in theta estimation")
  if (attachH) attr(theta,"H") <- H #nll$H
  theta
} ## estimate.theta

find.null.dev <- function(family,y,eta,offset,weights) {
## obtain the null deviance given y, best fit mu and
## prior weights
   fnull <- function(gamma,family,y,wt,offset) {
      ## evaluate deviance for single parameter model
      mu <- family$linkinv(gamma+offset)
      sum(family$dev.resids(y,mu, wt))
   }
   mu <- family$linkinv(eta-offset)
   mum <- mean(mu*weights)/mean(weights) ## initial value
   eta <- family$linkfun(mum) ## work on l.p. scale
   deta <- abs(eta)*.1 + 1  ## search interval half width
   ok <- FALSE
   while (!ok) {
     search.int <- c(eta-deta,eta+deta)
     op <- optimize(fnull,interval=search.int,family=family,y=y,wt = weights,offset=offset)
     if (op$minimum > search.int[1] && op$minimum < search.int[2]) ok <- TRUE else deta <- deta*2
  }
  op$objective
} ## find.null.dev


## extended families for mgcv, standard components. 
## family - name of family character string
## link - name of link character string
## linkfun - the link function
## linkinv - the inverse link function
## mu.eta - d mu/d eta function (derivative of inverse link wrt eta)
## note: for standard links this information is supplemented using 
##       function fix.family.link.extended.family with functions 
##       gkg where k is 2,3 or 4 giving the kth derivative of the 
##       link over the first derivative of the link to the power k.
##       for non standard links these functions must be supplied.
## dev.resids - function computing deviance residuals.
## Dd - function returning derivatives of deviance residuals w.r.t. mu and theta. 
## aic - function computing twice -ve log likelihood for 2df to be added to.
## initialize - expression to be evaluated in gam.fit4 and initial.spg 
##              to initialize mu or eta.
## preinitialize - optional function of y and family, returning a list with optional elements
##                 Theta - intitial Theta and y - modified y for use in fitting (see e.g. ocat and betar)
## postproc - function with arguments family, y, prior.weights, fitted, linear.predictors, offset, intercept
##            to call after fit to compute (optionally) the label for the family, deviance and null deviance.
##            See ocat for simple example and betar or ziP for complicated. Called in estimate.gam.
## ls - function to evaluated log saturated likelihood and derivatives w.r.t.
##      phi and theta for use in RE/ML optimization. If deviance used is just -2 log 
##      lik. can just return zeroes. 
## validmu, valideta - functions used to test whether mu/eta are valid.      
## n.theta - number of theta parameters.
## no.r.sq - optional TRUE/FALSE indicating whether r^2 can be computed for family
## ini.theta - function for initializing theta.
## putTheta, getTheta - functions for storing and retriving theta values in function 
##                      environment.
## rd - optional function for simulating response data from fitted model.
## residuals - optional function for computing residuals.
## predict - optional function for predicting from model, called by predict.gam.
## family$data - optional list storing any family specific data for use, e.g. in predict
##               function. - deprecated (commented out below - appears to be used nowhere)
## scale - < 0 to estimate. ignored if NULL 

#######################
## negative binomial...
#######################

nb <- function (theta = NULL, link = "log") { 
## Extended family object for negative binomial, to allow direct estimation of theta
## as part of REML optimization. Currently the template for extended family objects.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for nb family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)
  #else if (is.character(link)) {
  #  stats <- make.link(link)
  #  linktemp <- link
  #} else {
  #  if (inherits(link, "link-glm")) {
  #     stats <- link
  #          if (!is.null(stats$name))
  #              linktemp <- stats$name
  #      }
  #      else stop(linktemp, " link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"")
  #}
  ## Theta <-  NULL;
  n.theta <- 1
  if (!is.null(theta)&&theta!=0) {
      if (theta>0) { 
        iniTheta <- log(theta) ## fixed theta supplied
        n.theta <- 0 ## signal that there are no theta parameters to estimate
      } else iniTheta <- log(-theta) ## initial theta supplied
  } else iniTheta <- 0 ## inital log theta value
    
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta") # get(".Theta")
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) mu + mu^2/exp(get(".Theta")) ## Not actually needed!

    validmu <- function(mu) all(mu > 0)

    dev.resids <- function(y, mu, wt,theta=NULL) {
      if (is.null(theta)) theta <- get(".Theta")
      theta <- exp(theta) ## note log theta supplied
      mu[mu<=0] <- NA
      2 * wt * (y * log(pmax(1, y)/mu) - 
        (y + theta) * log((y + theta)/(mu + theta))) 
    } ## nb residuals
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the nb deviance...
      ##ltheta <- theta
      theta <- exp(theta)
      yth <- y + theta
      muth <- mu + theta
      r <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      r$Dmu <- 2 * wt * (yth/muth - y/mu)
      r$Dmu2 <- -2 * wt * (yth/muth^2 - y/mu^2)
      r$EDmu2 <- 2 * wt * (1/mu - 1/muth) ## exact (or estimated) expected weight
      if (level>0) { ## quantities needed for first derivatives
        r$Dth <- -2 * wt * theta * (log(yth/muth) + (1 - yth/muth) ) 
        r$Dmuth <- 2 * wt * theta * (1 - yth/muth)/muth
        r$Dmu3 <- 4 * wt * (yth/muth^3 - y/mu^3)
        r$Dmu2th <- 2 * wt * theta * (2*yth/muth - 1)/muth^2
	r$EDmu2th <- 2 * wt / muth^2
      } 
      if (level>1) { ## whole damn lot
        r$Dmu4 <- 2 * wt * (6*y/mu^4 - 6*yth/muth^4)
        r$Dth2 <- -2 * wt * theta * (log(yth/muth) +
                     theta*yth/muth^2 - yth/muth - 2*theta/muth + 1 +
                     theta /yth)
        r$Dmuth2 <- 2 * wt * theta * (2*theta*yth/muth^2 - yth/muth - 2*theta/muth + 1)/muth
        r$Dmu2th2 <- 2 * wt * theta * (- 6*yth*theta/muth^2 + 2*yth/muth + 4*theta/muth - 1) /muth^2
        r$Dmu3th <- 4 * wt * theta * (1 - 3*yth/muth)/muth^3
      }
      r
    } ## nb Dd

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        Theta <- exp(theta)
        term <- (y + Theta) * log(mu + Theta) - y * log(mu) +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
        2 * sum(term * wt)
    } ## aic nb
    
    ls <- function(y,w,theta,scale) {
       ## the log saturated likelihood function for nb
       Theta <- exp(theta)
       #vec <- !is.null(attr(theta,"vec.grad")) ## lsth by component?
       ylogy <- y;ind <- y>0;ylogy[ind] <- y[ind]*log(y[ind])
       term <- (y + Theta) * log(y + Theta) - ylogy +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
       ls <- -sum(term*w)
       ## first derivative wrt theta...
       yth <- y+Theta
       lyth <- log(yth)
       psi0.yth <- digamma(yth) 
       psi0.th <- digamma(Theta)
       term <- Theta * (lyth - psi0.yth + psi0.th-theta)
       #lsth <- if (vec) -term*w else -sum(term*w)
       LSTH <- matrix(-term*w,ncol=1)
       lsth <- sum(LSTH)
       ## second deriv wrt theta...
       psi1.yth <- trigamma(yth) 
       psi1.th <- trigamma(Theta)
       term <- Theta * (lyth - Theta*psi1.yth - psi0.yth + Theta/yth + Theta * psi1.th + psi0.th - theta -1)        
       lsth2 <- -sum(term*w)
       list(ls=ls, ## saturated log likelihood
            lsth1=lsth, ## first deriv vector w.r.t theta - last element relates to scale, if free
	    LSTH1=LSTH, ## rows are above derivs by datum
            lsth2=lsth2) ## Hessian w.r.t. theta, last row/col relates to scale, if free
    } ## ls nb

    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        ##n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
  
    postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept){
      posr <- list()
      posr$null.deviance <- find.null.dev(family,y,eta=linear.predictors,offset,prior.weights)
      posr$family <- 
      paste("Negative Binomial(",round(family$getTheta(TRUE),3),")",sep="")
      posr
    } ## postproc nb

    rd <- function(mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      rnbinom(n=length(mu),size=Theta,mu=mu)
    } ## rd nb

    qf <- function(p,mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      qnbinom(p,size=Theta,mu=mu)
    } ## qf nb
 

     environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
     environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) <- env
    structure(list(family = "negative binomial", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,rd=rd,qf=qf),
        class = c("extended.family","family"))
} ## nb

##########################
## Censored Poisson family
##########################

cpois <- function (link = "log") { 
## Extended family object for censored Poisson
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for cnorm family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)

  n.theta <- 0
 
  env <- new.env(parent = .GlobalEnv)
  getTheta <- function(trans=FALSE) {} # get(".Theta")
  putTheta <- function(theta) {}

  validmu <- if (link=="identity") function(mu) all(is.finite(mu)) else
             if (link=="log") function(mu) all(mu>0) else function(mu) all(mu>=0)

  dev.resids <- function(y, mu, wt,theta=NULL) { ## cpois
    
      yat <- attr(y,"censor") ## determines censoring type as below
      if (is.null(yat)) yat <- y 

      ii <- which(yat==y) ## uncensored observations
      d <- rep(0,length(y))
      if (length(ii)) d[ii] <- 2*(dpois(y[ii],y[ii],log=TRUE)-dpois(y[ii],mu[ii],log=TRUE))

      ii <- which(is.finite(yat)&yat!=y) ## interval censored
      if (length(ii)) {
        y1 <- pmax(yat[ii],y[ii]); y0 <- pmin(yat[ii],y[ii])
        ## get mu for saturated likelihood... 
        musat <- exp((lgamma(floor(y1)+1)-lgamma(floor(y0)+1))/(floor(y1)-floor(y0)))
        d[ii] <- 2*(log(ppois(y1,musat)-ppois(y0,musat))  - log(ppois(y1,mu[ii])-ppois(y0,mu[ii])))
      }
      
      ii <- which(yat == -Inf) ## left censored (sat lik is 1)
      if (length(ii)) d[ii] <- -2*ppois(y[ii],mu[ii],lower.tail=TRUE,log.p=TRUE)
      
      ii <- which(yat == Inf) ## right censored
      if (length(ii)) d[ii] <- -2*ppois(y[ii],mu[ii],lower.tail=FALSE,log.p=TRUE)
      d
    } ## dev.resids cpoi
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the cpois deviance...
     
      yat <- attr(y,"censor")
      if (is.null(yat)) yat <- y
      ## get case indices...
      iu <- which(yat==y) ## uncensored observations
      ii <- which(is.finite(yat*y)&yat!=y) ## interval censored
      il <- which(yat == -Inf) ## left censored
      ir <- which(yat == Inf) ## right censored
      n <- length(mu)
      f1 <- f2  <- rep(0,n)
      if (level>0) f3 <- f1
      if (level>1) f4 <- f1
      
      if (length(iu)) { ## uncensored
        yiu <- y[iu]; miu <- mu[iu]
	lf <- dpois(yiu,miu,log=TRUE)
	f1[iu] <- exp(dpois(yiu-1,miu,log=TRUE)-lf)
	f2[iu] <- exp(dpois(yiu-2,miu,log=TRUE)-lf)
	if (level>0) f3[iu] <- exp(dpois(yiu-3,miu,log=TRUE)-lf)
	if (level>1) f4[iu] <- exp(dpois(yiu-4,miu,log=TRUE)-lf)
      } ## uncensored done
      
      if (length(ii)) { ## interval censored
        y0 <- pmin(y[ii],yat[ii]); y1 <- pmax(y[ii],yat[ii])
        mii <- mu[ii]
        g <- ppois(y1,mii) - ppois(y0,mii)
	f1[ii] <- (ppois(y1-1,mii) - ppois(y0-1,mii))/g
	f2[ii] <- (ppois(y1-2,mii) - ppois(y0-2,mii))/g
	if (level>0) f3[ii] <- (ppois(y1-3,mii) - ppois(y0-3,mii))/g
	if (level>1) f4[ii] <- (ppois(y1-4,mii) - ppois(y0-4,mii))/g
      } ## interval censored done

      for (lt in c(TRUE,FALSE)) { ## do left then right censoring...
        ii <- if (lt) il else ir # left or right censoring
        if (length(ii)) {
          yil <- y[ii]; mil <- mu[ii] 
          lf <- ppois(yil,mil,lower.tail=lt,log.p=TRUE)
	  f1[ii] <- exp(ppois(yil-1,mil,lower.tail=lt,log.p=TRUE)-lf)
	  f2[ii] <- exp(ppois(yil-2,mil,lower.tail=lt,log.p=TRUE)-lf)
	  if (level>0) f3[ii] <- exp(ppois(yil-3,mil,lower.tail=lt,log.p=TRUE)-lf)
	  if (level>1) f4[ii] <- exp(ppois(yil-4,mil,lower.tail=lt,log.p=TRUE)-lf)
        } 
      } ## lt TRUE/FALSE
      
      r <- list(Dmu=-2*(f1-1),Dmu2=-2*(f2-f1^2)); r$EDmu2 = r$Dmu2
      if (level>0) r$Dmu3 <-  -2*(f3 - 3*f1*f2 + 2*f1^3)
      if (level>1) r$Dmu4 <- -2*(f4 - 4*f3*f1 + 12*f1^2*f2 - 3*f2^2 - 6*f1^4)
      r
    } ## Dd cpois

    aic <- function(y, mu, theta=NULL, wt, dev) { ## cpois AIC
      
	yat <- attr(y,"censor")
        if (is.null(yat)) yat <- y
        ii <- which(is.na(yat)|yat==y) ## uncensored observations
        d <- rep(0,length(y))
        if (length(ii)) d[ii] <- dpois(y[ii],mu[ii],log=TRUE)
	
        ii <- which(is.finite(yat)&yat!=y) ## interval censored
        if (length(ii)) {
          y1 <- pmax(yat[ii],y[ii]); y0 <- pmin(yat[ii],y[ii])
          d[ii] <- log(ppois(y1,mu[ii])-ppois(y0,mu[ii]))
        }
	
        ii <- which(yat == -Inf) ## left censored
        if (length(ii)) d[ii] <- ppois(y[ii],mu[ii],log.p=TRUE)
        ii <- which(yat == Inf) ## right censored
        if (length(ii)) d[ii] <- ppois(y[ii],mu[ii],lower.tail=FALSE,log.p=TRUE)
        -2*sum(d) ## -2*log likelihood
    } ## AIC cpois
    
    ls <- function(y,w,theta,scale) {
       ## the cpois log saturated likelihood function.
      
       yat <- attr(y,"censor")
       if (is.null(yat)) yat <- y

       ii <- which(yat==y) ## uncensored observations
       d2 <- d1 <- d <- rep(0,length(y))
       if (length(ii)) {
         d[ii] <- dpois(y[ii],y[ii],log=TRUE)
       }
       
       ii <- which(is.finite(yat)&yat!=y) ## interval censored
       if (length(ii)) {
         y1 <- pmax(yat[ii],y[ii]); y0 <- pmin(yat[ii],y[ii])
	 mus <- exp((lgamma(floor(y1)+1)-lgamma(floor(y0)+1))/(floor(y1)-floor(y0)))
         d[ii] <- log(ppois(y1,mus)-ppois(y0,mus))  ## log saturated likelihood
       }
       ## right or left censored saturated log likelihoods are zero.
       list(ls=sum(d), ## saturated log likelihood
            lsth1=0, ## first deriv vector w.r.t theta - last element relates to scale, if free
	    LSTH1=matrix(d1,ncol=1),
            lsth2=0) ## Hessian w.r.t. theta, last row/col relates to scale, if free
    } ## ls cpois

    initialize <- expression({ ## cpois
        if (is.matrix(y)) {
	 .yat <- y[,2]
	 y <- y[,1]
	 attr(y,"censor") <- .yat
	} 
        mustart <- if (family$link=="identity") y else pmax(y,min(y>0))
    })
  
    postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept){
      posr <- list()
      if (is.matrix(y)) {
	 .yat <- y[,2]
	 y <- y[,1]
	 attr(y,"censor") <- .yat
      } 
      posr$null.deviance <- find.null.dev(family,y,eta=linear.predictors,offset,prior.weights)
      posr$family <- "cpois"
      posr
    } ## postproc cpoi

    rd <- function(mu,wt,scale) { ## NOTE - not done
      
    }

    qf <- function(p,mu,wt,scale) { ## NOTE - not done
     
    }

    subsety <- function(y,ind) { ## function to subset response
      if (is.matrix(y)) return(y[ind,])
      yat <- attr(y,"censor")
      y <- y[ind]
      if (!is.null(yat)) attr(y,"censor") <- yat[ind]
      y
    } ## subsety

     environment(dev.resids) <- environment(aic) <- 
     environment(rd)<- environment(qf)<- environment(putTheta) <- env
    structure(list(family = "cnorm", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,subsety=subsety,#variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,n.theta=0,putTheta=putTheta,getTheta=getTheta),
        class = c("extended.family","family"))
} ## cpois



#########################
## Censored normal family
#########################

cnorm <- function (theta = NULL, link = "identity") { 
## Extended family object for censored Gaussian, as required for Tobit regression or log-normal
## Accelerated Failure Time models.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for cnorm family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)

  n.theta <- 1
  if (!is.null(theta)&&theta!=0) {
      if (theta>0) { 
        iniTheta <- log(theta) ## fixed theta supplied
        n.theta <- 0 ## signal that there are no theta parameters to estimate
      } else iniTheta <- log(-theta) ## initial theta supplied
  } else iniTheta <- 0 ## inital log theta value
    
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta") # get(".Theta")
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    validmu <- if (link=="identity") function(mu) all(is.finite(mu)) else function(mu) all(mu>0)

    dev.resids <- function(y, mu, wt,theta=NULL) { ## cnorm
      if (is.null(theta)) theta <- get(".Theta")
      th <- theta - log(wt)/2
      yat <- attr(y,"censor")
      if (is.null(yat)) yat <- rep(NA,length(y))
      ii <- which(yat==y) ## uncensored observations
      d <- rep(0,length(y))
      if (length(ii)) d[ii] <- (y[ii]-mu[ii])^2*exp(-2*th[ii])
      ii <- which(is.finite(yat)&yat!=y) ## interval censored
      if (length(ii)) {
        y1 <- pmax(yat[ii],y[ii]); y0 <- pmin(yat[ii],y[ii])
	y10 <- (y1-y0)*exp(-th[ii])/2
        d[ii] <- 2*log(dpnorm(-y10,y10)) - ## 2 * log saturated likelihood
	         2*log(dpnorm((y0-mu[ii])*exp(-th[ii]),(y1-mu[ii])*exp(-th[ii])))
      }
      ii <- which(yat == -Inf) ## left censored
      if (length(ii)) d[ii] <- -2*pnorm((y[ii]-mu[ii])*exp(-th[ii]),log.p=TRUE)
      ii <- which(yat == Inf) ## right censored
      if (length(ii)) d[ii] <- -2*pnorm(-(y[ii]-mu[ii])*exp(-th[ii]),log.p=TRUE)
      d
    } ## dev.resids cnorm 
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the cnorm deviance...
     
      th <- theta - log(wt)/2
      eth <- exp(-th)
      e2th <- eth*eth
      e3th <- e2th*eth
      yat <- attr(y,"censor")
      if (is.null(yat)) yat <- y
      ## get case indices...
      iu <- which(yat==y) ## uncensored observations
      ii <- which(is.finite(yat*y)&yat!=y) ## interval censored
      il <- which(yat == -Inf) ## left censored
      ir <- which(yat == Inf) ## right censored
      n <- length(mu)
      Dmu <- Dmu2 <- rep(0,n)
      if (level>0) Dth <- Dmuth <- Dmu3 <- Dmu2th <- Dmu
      if (level>1) Dmu4 <- Dth2 <- Dmuth2 <- Dmu2th2 <- Dmu3th <- Dmu
      if (length(iu)) { ## uncensored
        ethi <- eth[iu]; e2thi <- e2th[iu]
        ymeth <- (y[iu]-mu[iu])*ethi
        Dmu[iu] <- Dmui <- -2*ymeth*ethi
	Dmu2[iu] <- 2*e2thi
	if (level>0) {
          Dth[iu] <- -2*ymeth^2
	  Dmuth[iu] <- -2*Dmui
	  Dmu3[iu] <- 0
	  Dmu2th[iu] <- -4*e2thi
        }
	if (level>1) {
          Dmu4[iu] <- 0
	  Dth2[iu] <- -2*Dth[iu]
	  Dmuth2[iu] <- 4*Dmui
	  Dmu2th2[iu] <- 8*e2thi
	  Dmu3th[iu] <- 0
        }
      } ## uncensored done
      if (length(ii)) { ## interval censored
        y0 <- pmin(y[ii],yat[ii]); y1 <- pmax(y[ii],yat[ii])
	ethi <- eth[ii];e2thi <- e2th[ii];e3thi <- e3th[ii] 
	y10 <- (y1-y0)*ethi/2
	ymeth0 <- (y0-mu[ii])*ethi;ymeth1 <- (y1-mu[ii])*ethi
	D0 <- dpnorm(ymeth0,ymeth1)
	
	dnorm0 <- dnorm(ymeth0);dnorm1 <- dnorm(ymeth1)
	Dmui <- Dmu[ii] <- 2*ethi*(dnorm1-dnorm0)/D0
	Dmu2[ii] <- Dmui^2/2 + 2*e2thi*(dnorm1*ymeth1-dnorm0*ymeth0)/D0
	if (level>0) {
	  Dmu2i <- Dmu2[ii];
	  ymeth12 <- ymeth1^2; ymeth13 <- ymeth12*ymeth1
	  ymeth02 <- ymeth0^2; ymeth03 <- ymeth02*ymeth0
	  Dls <- dpnorm(-y10,y10)
          Dth[ii] <- Dthi <- 2*(dnorm1*ymeth1-dnorm0*ymeth0)/D0
	  Dmuth[ii] <- Dmuthi <- Dmui*Dthi/2 + 2*ethi*(dnorm1*(ymeth12-1)-dnorm0*(ymeth02-1))/D0
	  Dmu3[ii] <- Dmui*(3*Dmu2i/2 - Dmui^2/4 - e2thi) +
	              2*e3thi*(dnorm1*ymeth12-dnorm0*ymeth02)/D0
	  Dmu2th[ii] <- (Dmu2i*Dthi+Dmui*Dmuthi)/2 +
	  ethi*(dnorm1*(2*ymeth13*ethi+ Dmui*ymeth12-6*ymeth1*ethi - Dmui)-
                   dnorm0*(2*ymeth03*ethi+ Dmui*ymeth02-6*ymeth0*ethi - Dmui))/D0
        }
	if (level>1) {
	  ymeth14 <- ymeth13*ymeth1; ymeth15 <- ymeth12*ymeth13
	  ymeth04 <- ymeth03*ymeth0; ymeth05 <- ymeth02*ymeth03
	  Dmu3i <- Dmu3[ii]; Dmu2thi <- Dmu2th[ii]
          Dmu4[ii] <- Dmu2i*(3*Dmu2i/2-Dmui^2/4-e2thi) + Dmui*(3*Dmu3i/2-Dmui*Dmu2i/2) +
	              e3thi*(dnorm1*(2*ymeth13*ethi + Dmui*ymeth12-4*ymeth1*ethi) -
		             dnorm0*(2*ymeth03*ethi + Dmui*ymeth02-4*ymeth0*ethi))/D0
	  Dth2[ii] <- Dth2i <- Dthi^2/2 + 2*(dnorm1*(ymeth13 - ymeth1) - dnorm0*(ymeth03 - ymeth0))/D0 
	  Dmuth2[ii] <- (Dmuthi*Dthi+Dmui*Dth2i)/2 +
	             ethi*(dnorm1*(2*ymeth14 + (Dthi-8)*ymeth12 + 2-Dthi) -
	                   dnorm0*(2*ymeth04 + (Dthi-8)*ymeth02 + 2-Dthi))/D0
	  Dmu2th2[ii] <- (Dmu2thi*Dthi+Dmuthi^2 + Dmu2i*Dth2i + Dmui*Dmuth2[ii])/2 +
	  0.5*ethi*(dnorm1*(4*ymeth15*ethi+2*Dmui*ymeth14+2*(Dthi-16)*ymeth13*ethi+
	  2*(18-3*Dthi)*ymeth1*ethi+(2*Dmuthi+(Dthi-8)*Dmui)*ymeth12+(2-Dthi)*Dmui-2*Dmuthi)-
	            dnorm0*(4*ymeth05*ethi+2*Dmui*ymeth04+2*(Dthi-16)*ymeth03*ethi+
	  2*(18-3*Dthi)*ymeth0*ethi+(2*Dmuthi+(Dthi-8)*Dmui)*ymeth02+(2-Dthi)*Dmui-2*Dmuthi)
	  )/D0
	
	  Dmu3th[ii] <- Dmu3i*Dthi/2 + Dmu2i*Dmuthi + Dmui*Dmu2thi/2 +
	      0.5*ethi*(dnorm1*(4*ymeth14*e2thi + 4*Dmui*ymeth13*ethi - 12*Dmui*ymeth1*ethi +
	                        (Dmui^2+2*Dmu2i-24*e2thi)*ymeth12+ 12*e2thi -(Dmui^2+2*Dmu2i)) -
		        dnorm0*(4*ymeth04*e2thi + 4*Dmui*ymeth03*ethi - 12*Dmui*ymeth0*ethi +
			        (Dmui^2+2*Dmu2i-24*e2thi)*ymeth02+ 12*e2thi -(Dmui^2+2*Dmu2i)))/D0
        }
	if (level>0)  Dth[ii] <- Dthi -4*dnorm(y10)*y10/Dls
	if (level>1)  Dth2[ii] <- Dth2i - 8*(dnorm(y10)*y10/Dls)^2  - 4*dnorm(y10)*(y10^3 -y10)/Dls
      } ## interval censored done
      if (length(il)) { ## left censoring (y0 = -Inf, basically)
        ethi <- eth[il];e2thi <- e2th[il];e3thi <- e3th[il]
        ymeth1 <- (y[il]-mu[il])*ethi;dnorm1 <- dnorm(ymeth1)
        D0 <- pnorm(ymeth1)
	Dmui <- Dmu[il] <- 2*ethi*dnorm1/D0
	Dmu2[il] <- Dmui^2/2 + 2*e2thi*dnorm1*ymeth1/D0
	if (level>0) {
	  ymeth12 <- ymeth1^2; ymeth13 <- ymeth12*ymeth1
	  Dmu2i <- Dmu2[il]
          Dth[il] <- Dthi <- 2*dnorm1*ymeth1/D0
	  Dmuth[il] <- Dmuthi <- Dmui*Dthi/2 + 2*ethi*dnorm1*(ymeth12-1)/D0
	  Dmu3[il] <- Dmui*(3*Dmu2i/2 - Dmui^2/4 - e2thi) + 2*e3thi*dnorm1*ymeth12/D0
	  Dmu2th[il] <- (Dmu2i*Dthi+Dmui*Dmuthi)/2 +
	  ethi*dnorm1*(2*ymeth13*ethi+ Dmui*ymeth12-6*ymeth1*ethi - Dmui)/D0
        }
	if (level>1) {
	  ymeth14 <- ymeth13*ymeth1; ymeth15 <- ymeth12*ymeth13
	  Dmu3i <- Dmu3[il]; Dmu2thi <- Dmu2th[il]
          Dmu4[il] <- Dmu2i*(3*Dmu2i/2-Dmui^2/4-e2thi) + Dmui*(3*Dmu3i/2-Dmui*Dmu2i/2) +
	              e3thi*dnorm1*(2*ymeth13*ethi + Dmui*ymeth12-4*ymeth1*ethi)/D0
	  Dth2[il] <- Dth2i <- Dthi^2/2 + 2*dnorm1*(ymeth13 - ymeth1)/D0
	  Dmuth2[il] <- (Dmuthi*Dthi+Dmui*Dth2i)/2 +
	                ethi*dnorm1*(2*ymeth14 + (Dthi-8)*ymeth12 + 2-Dthi)/D0
	  Dmu2th2[il] <- (Dmu2thi*Dthi+Dmuthi^2 + Dmu2i*Dth2i + Dmui*Dmuth2[il])/2 +
	  0.5*ethi*dnorm1*(4*ymeth15*ethi+2*Dmui*ymeth14+2*(Dthi-16)*ymeth13*ethi+
	  2*(18-3*Dthi)*ymeth1*ethi+(2*Dmuthi+(Dthi-8)*Dmui)*ymeth12+(2-Dthi)*Dmui-2*Dmuthi)/D0
	
	  Dmu3th[il] <- Dmu3i*Dthi/2 + Dmu2i*Dmuthi + Dmui*Dmu2thi/2 +
	      0.5*ethi*dnorm1*(4*ymeth14*e2thi + 4*Dmui*ymeth13*ethi - 12*Dmui*ymeth1*ethi +
	      (Dmui^2+2*Dmu2i-24*e2thi)*ymeth12+ 12*e2thi -(Dmui^2+2*Dmu2i))/D0
        }
      } # left censoring done
      if (length(ir)) { ## right censoring - basically y1 = Inf
        ethi <- eth[ir];e2thi <- e2th[ir];e3thi <- e3th[ir]
	ymeth0 <- (y[ir]-mu[ir])*ethi;
	D0 <- pnorm(-ymeth0)
	dnorm0 <- dnorm(ymeth0);
	Dmu[ir] <- Dmui <- -2*ethi*dnorm0/D0
	Dmu2[ir] <- Dmui^2/2 - 2*e2thi*dnorm0*ymeth0/D0
	if (level>0) {
	  ymeth02 <- ymeth0^2; ymeth03 <- ymeth02*ymeth0
	  Dmu2i <- Dmu2[ir]
          Dth[ir] <- Dthi <- -2*dnorm0*ymeth0/D0
	  Dmuth[ir] <- Dmuthi <- Dmui*Dthi/2 - 2*ethi*dnorm0*(ymeth02-1)/D0
	  Dmu3[ir] <- Dmui*(3*Dmu2i/2 - Dmui^2/4 - e2thi) - 2*e3thi*dnorm0*ymeth02/D0
	  Dmu2th[ir] <- (Dmu2i*Dthi+Dmui*Dmuthi)/2 -
	                ethi*dnorm0*(2*ymeth0^3*ethi+ Dmui*ymeth02-6*ymeth0*ethi - Dmui)/D0
        }
	if (level>1) {
	  ymeth04 <- ymeth03*ymeth0; ymeth05 <- ymeth02*ymeth03
	  Dmu3i <- Dmu3[ir]; Dmu2thi <- Dmu2th[ir]
          Dmu4[ir] <- Dmu2i*(3*Dmu2i/2-Dmui^2/4-e2thi) + Dmui*(3*Dmu3i/2-Dmui*Dmu2i/2) -
	              e3thi*dnorm0*(2*ymeth0^3*ethi + Dmui*ymeth02-4*ymeth0*ethi)/D0
	  Dth2[ir] <- Dthi^2/2 - 2*dnorm0*(ymeth0^3 - ymeth0)/D0
	  Dmuth2[ir] <- (Dmuthi*Dthi+Dmui*Dth2[ir])/2 -
	                 ethi*dnorm0*(2*ymeth04 + (Dthi-8)*ymeth02 + 2-Dthi)/D0
	  Dmu2th2[ir] <- (Dmu2thi*Dthi+Dmuthi^2 + Dmu2i*Dth2[ir] + Dmui*Dmuth2[ir])/2 -
	  0.5*ethi*dnorm0*(4*ymeth05*ethi+2*Dmui*ymeth04+2*(Dthi-16)*ymeth0^3*ethi+
	  2*(18-3*Dthi)*ymeth0*ethi+(2*Dmuthi+(Dthi-8)*Dmui)*ymeth02+(2-Dthi)*Dmui-2*Dmuthi)/D0
	
	  Dmu3th[ir] <- Dmu3i*Dthi/2 + Dmu2i*Dmuthi + Dmui*Dmu2thi/2 -
	      0.5*ethi*(dnorm0*(4*ymeth04*e2thi + 4*Dmui*ymeth0^3*ethi -
              12*Dmui*ymeth0*ethi + (Dmui^2+2*Dmu2i-24*e2thi)*ymeth02+ 12*e2thi -(Dmui^2+2*Dmu2i)))/D0
        }
      } # right censoring done
      r <- list(Dmu=Dmu,Dmu2=Dmu2,EDmu2=Dmu2)
      if (level>0) {
        r$Dth <- Dth;r$Dmuth <- Dmuth;r$Dmu3 <- Dmu3
	r$EDmu2th <- r$Dmu2th <- Dmu2th; 
      }
      if (level>1) {
        r$Dmu4 <- Dmu4; r$Dth2 <- Dth2;  r$Dmuth2 <- Dmuth2;
	r$Dmu2th2 <- Dmu2th2; r$Dmu3th <- Dmu3th
      }
      r
    } ## Dd cnorm

    aic <- function(y, mu, theta=NULL, wt, dev) { ## cnorm AIC
        if (is.null(theta)) theta <- get(".Theta")
        th <- theta - log(wt)/2
	yat <- attr(y,"censor")
        if (is.null(yat)) yat <- y
        ii <- which(is.na(yat)|yat==y) ## uncensored observations
        d <- rep(0,length(y))
        if (length(ii)) d[ii] <- (y[ii]-mu[ii])^2*exp(-2*th[ii]) + log(2*pi) + 2*th[ii] 
        ii <- which(is.finite(yat)&yat!=y) ## interval censored
        if (length(ii)) {
          y1 <- pmax(yat[ii],y[ii]); y0 <- pmin(yat[ii],y[ii])
          d[ii] <- - 2*log(dpnorm((y0-mu[ii])*exp(-th[ii]),(y1-mu[ii])*exp(-th[ii])))
        }
        ii <- which(yat == -Inf) ## left censored
        if (length(ii)) d[ii] <- -2*pnorm((y[ii]-mu[ii])*exp(-th[ii]),log.p=TRUE)
        ii <- which(yat == Inf) ## right censored
        if (length(ii)) d[ii] <- -2*pnorm(-(y[ii]-mu[ii])*exp(-th[ii]),log.p=TRUE)
 
        sum(d) ## -2*log likelihood
    } ## AIC cnorm
    
    ls <- function(y,w,theta,scale) {
       ## the cnorm log saturated likelihood function.
       th <- theta - log(w)/2
       yat <- attr(y,"censor")
       if (is.null(yat)) yat <- y
       ii <- which(yat==y) ## uncensored observations
       d2 <- d1 <- d <- rep(0,length(y))
       if (length(ii)) {
         d[ii] <- log(2*pi)/2 - th[ii]
	 d1[ii] <- -1
       }	 
       ii <- which(is.finite(yat)&yat!=y) ## interval censored
       if (length(ii)) {
         y1 <- pmax(yat[ii],y[ii]); y0 <- pmin(yat[ii],y[ii])
	 y10 <- (y1-y0)*exp(-th[ii])/2
	 d0 <- dpnorm(-y10,y10)
         d[ii] <- log(d0)  ## log saturated likelihood
	 d1[ii] <- -2*dnorm(y10)*y10/d0
	 d2[ii] <- -d1[ii]^2 - 2*dnorm(y10)*y10^2*(y1-y0)/2
       }
       ## right or left censored saturated log likelihoods are zero.
       list(ls=sum(d), ## saturated log likelihood
            lsth1=sum(d1), ## first deriv vector w.r.t theta - last element relates to scale, if free
	    LSTH1=matrix(d1,ncol=1),
            lsth2=sum(d2)) ## Hessian w.r.t. theta, last row/col relates to scale, if free
    } ## ls cnorm

    initialize <- expression({ ## cnorm
        if (is.matrix(y)) {
	 .yat <- y[,2]
	 y <- y[,1]
	 attr(y,"censor") <- .yat
	} 
        mustart <- if (family$link=="identity") y else pmax(y,min(y>0))
    })
  
    postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept){
      posr <- list()
      if (is.matrix(y)) {
	 .yat <- y[,2]
	 y <- y[,1]
	 attr(y,"censor") <- .yat
      } 
      posr$null.deviance <- find.null.dev(family,y,eta=linear.predictors,offset,prior.weights)
      posr$family <- 
      paste("cnorm(",round(family$getTheta(TRUE),3),")",sep="")
      posr
    } ## postproc cnorm

    rd <- function(mu,wt,scale) { ## NOTE - not done
      Theta <- exp(get(".Theta"))
    }

    qf <- function(p,mu,wt,scale) { ## NOTE - not done
      Theta <- exp(get(".Theta"))
    }

    subsety <- function(y,ind) { ## function to subset response
      if (is.matrix(y)) return(y[ind,])
      yat <- attr(y,"censor")
      y <- y[ind]
      if (!is.null(yat)) attr(y,"censor") <- yat[ind]
      y
    } ## subsety

     environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
     environment(rd)<- environment(qf)<- environment(putTheta) <- env
    structure(list(family = "cnorm", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,subsety=subsety,#variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta),#,rd=rd,qf=qf),
        class = c("extended.family","family"))
} ## cnorm


#################################################
## censored logistic family (Chris Shen)
#################################################
    log1pexp <- function(x) { ## compute log(1+e^x) safely
      result <- x
      ii <- which(x<=37); result[ii] <- exp(x[ii])
      ii <- which(-37<x & x<=18); result[ii] <- log1p(exp(x[ii]))
      ii <- which(18<x & x<=33.3); result[ii] <- x[ii]+exp(-x[ii])
      return(result)
    } # log1pexp

    log1mexp <- function(a) { ## compute log(1-e^(-a)) safely
      result <- numeric(length(a))
      ii <- which(0<a & a<=log(2)); result[ii] <- log(-expm1(-a[ii]))
      ii <- which(a>log(2)); result[ii] <- log1p(-exp(-a[ii]))
      return(result)
    } # log1mexp


clog <- function(theta=NULL, link="identity") {
  # links
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for clog family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)

#  else if (is.character(link)) {
#    stats <- make.link(link)
#    linktemp <- link
#  } else {
#    if (inherits(link, "link-glm")) {
#      stats <- link
#      if (!is.null(stats$name)) linktemp <- stats$name
#    }
#    else stop(linktemp, "link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"")
#  }
  
  # theta
  n.theta <- 1
  if (!is.null(theta) && theta!=0) {
    # positive theta; implies fixed value
    if (theta>0) {
      iniTheta <- log(theta)
      n.theta <- 0
    }
    # negative theta; theta is now estimated
    else iniTheta <- log(-theta)
  } else iniTheta <- 0 # initial theta value for estimation

  env <- new.env(parent=.GlobalEnv)
  assign(".Theta", iniTheta, envir=env)
  getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta")
  putTheta <- function(theta) assign(".Theta", theta, envir=environment(sys.function()))

  validmu <- if (link=="identity") function(mu) all(is.finite(mu)) else function(mu) all(mu>0)

  # deviance residuals squared
  # basically the deviance
  dev.resids <- function(y, mu, wt, theta=NULL) {
    log1pexp <- function(x) { ## compute log(1+e^x) safely
      result <- x
      ii <- which(x<=37); result[ii] <- exp(x[ii])
      ii <- which(-37<x & x<=18); result[ii] <- log1p(exp(x[ii]))
      ii <- which(18<x & x<=33.3); result[ii] <- x[ii]+exp(-x[ii])
      return(result)
    } # log1pexp

    log1mexp <- function(a) { ## compute log(1-e^(-a)) safely
      result <- numeric(length(a))
      ii <- which(0<a & a<=log(2)); result[ii] <- log(-expm1(-a[ii]))
      ii <- which(a>log(2)); result[ii] <- log1p(-exp(-a[ii]))
      return(result)
    } # log1mexp
    
    if (is.null(theta)) theta <- get(".Theta")
    yat <- attr(y,"censor")
    if (is.null(yat)) yat <- y
    d <- rep(0,length(y))

    # get indices
    iu <- which(yat==y) # uncensored
    ii <- which(is.finite(yat*y) & yat!=y) # interval censored
    il <- which(yat==-Inf) # left censored
    ir <- which(yat==Inf) # right censored

    # uncensored
    if (length(iu)) {
      si <- exp(theta)/sqrt(wt[iu])
      mui <- mu[iu]
      yi <- y[iu]
      d[iu] <- 2*(-2*log(2)+((yi-mui)/si)+2*log1pexp(-(yi-mui)/si)) # prevents overflow
    }

    # interval censored
    if (length(ii)) {
      si <- exp(theta)/sqrt(wt[ii])
      mui <- mu[ii]
      yl <- pmin(y[ii], yat[ii])
      yr <- pmax(y[ii], yat[ii])

      lm <- ((yr-yl)/(2*si))+log1mexp((yr-yl)/(2*si))-log1pexp((yr-yl)/(2*si))
      l <- log1mexp((yr-yl)/si)-log1pexp((yl-mui)/si)-log1pexp(-(yr-mui)/si) # prevents overflow
      d[ii] <- 2*(lm-l)
    }

    # left censored
    if (length(il)) {
      si <- exp(theta)/sqrt(wt[il])
      mui <- mu[il]
      yr <- y[il]

      lm <- 0
      l <- -log1pexp(-(yr-mui)/si)
      d[il] <- 2*(lm-l)
    }

    # right censored
    if (length(ir)) {
      si <- exp(theta)/sqrt(wt[ir])
      mui <- mu[ir]
      yl <- y[ir]
      
      lm <- 0
      l <- (-(yl-mui)/si)-log1pexp(-(yl-mui)/si)
      d[ir] <- 2*(lm-l)
    }
    return(d)
  } ## dev.resids clogistic

  # deviance derivatives
  Dd <- function(y, mu, theta, wt, level=0) {
    yat <- attr(y, "censor")
    if (is.null(yat)) yat <- y

    # get indices
    iu <- which(yat==y) # uncensored
    ii <- which(is.finite(yat*y) & yat!=y) # interval censored
    il <- which(yat==-Inf) # left censored
    ir <- which(yat==Inf) # right censored

    # initialise variables
    n <- length(mu)
    Dmu <- Dmu2 <- rep(0,n)
    if (level>0) Dth <- Dmuth <- Dmu3 <- Dmu2th <- Dmu
    if (level>1) Dmu4 <- Dth2 <- Dmuth2 <- Dmu2th2 <- Dmu3th <- Dmu

    # uncensored
    if (length(iu)) {
      si <- exp(theta)/sqrt(wt[iu])
      yi <- y[iu]; mui <- mu[iu]
      alphai <- 1/(2+expm1((yi-mui)/si)) # with precision

      Dmu[iu] <- Dmui <- (2/si)*(2*alphai-1)
      Dmu2[iu] <- Dmu2i <- (4/si^2)*alphai*(1-alphai)
      if (level>0) {
        Dmuth[iu] <- Dmuthi <- -Dmui+(yi-mui)*Dmu2i
        Dth[iu] <- Dthi <- (yi-mui)*Dmui
        Dmu3[iu] <- Dmu3i <- (-1/2)*Dmui*Dmu2i
        Dmu2th[iu] <- Dmu2thi <- -2*Dmu2i+(1/2)*(mui-yi)*Dmui*Dmu2i
      }
      if (level>1) {
        Dth2[iu] <- Dth2i <- (yi-mui)*Dmuthi
        Dmuth2[iu] <- Dmuth2i <- -Dmuthi+(yi-mui)*Dmu2thi
        Dmu4[iu] <- Dmu4i <- (-1/2)*(Dmu2i^2+Dmui*Dmu3i)
        Dmu3th[iu] <- Dmu3thi <- (-1/2)*(Dmu2i*Dmuthi+Dmui*Dmu2thi)
        Dmu2th2[iu] <- Dmu2th2i <- -2*Dmu2thi+(1/2)*(mui-yi)*(Dmui*Dmu2thi+Dmu2i*Dmuthi)
      }
    }

    # interval censored
    if (length(ii)) {
      yl <- pmin(y[ii], yat[ii])
      yr <- pmax(y[ii], yat[ii])
      si <- exp(theta)/sqrt(wt[ii])
      mui <- mu[ii]

      alphai <- 1/(2+expm1(-(yr-mui)/si))
      betai <- 1/(2+expm1(-(yl-mui)/si))

      Dmu[ii] <- Dmui <- (2/si)*(1-alphai-betai)
      Dmu2[ii] <- Dmu2i <- (2/si^2)*(alphai-alphai^2+betai-betai^2)
      if (level>0) {
        Dmuth[ii] <- Dmuthi <- -Dmui+(2/si^2)*((yr-mui)*(alphai-alphai^2)+(yl-mui)*(betai-betai^2))

        lmth <- (yr-yl)*(1/si)*(1/(expm1(-(yr-yl)/(2*si))-expm1((yr-yl)/(2*si))))
        lth <- (1/si)*(-(yr-yl)*(1/(expm1((yr-yl)/si)))+
          (yl-mui)*(1/(expm1((mui-yl)/si)+2))-(yr-mui)*(1/(expm1((yr-mui)/si)+2)))
        Dth[ii] <- Dthi <- 2*(lmth-lth)

        Dmu3[ii] <- Dmu3i <- (4/si^2)*(alphai^2-alphai^3+betai^2-betai^3)-Dmu2i
        Dmu2th[ii] <- Dmu2thi <- -2*Dmu2i+(-1/si)*(Dmuthi+Dmui)+
          (4/si^3)*((yr-mui)*(alphai^2-alphai^3)+(yl-mui)*(betai^2-betai^3))
      }
      if (level>1) {
        lmth2 <- -lmth-(1/2)*(yr-yl)^2*(1/si^2)*(1/(expm1(-(yr-yl)/(2*si))+
          -expm1((yr-yl)/(2*si))))^2*(2+expm1(-(yr-yl)/(2*si))+expm1((yr-yl)/(2*si)))
        lth2 <- -lth-(1/si^2)*((yr-yl)^2*(1/(expm1((yr-yl)/si)))*(1/(-expm1(-(yr-yl)/si)))+
          (yl-mui)^2*(1/(expm1((mui-yl)/si)+2))*(1/(expm1((yl-mui)/si)+2))+
          (yr-mui)^2*(1/(expm1((yr-mui)/si)+2))*(1/(expm1((mui-yr)/si)+2)))
        Dth2[ii] <- Dth2i <- 2*(lmth2-lth2)

        Dmuth2i <- -3*Dmuthi-2*Dmui+(-2/si^3)*((yr-mui)^2*(1-2*alphai)*(alphai-alphai^2)+
          (yl-mui)^2*(1-2*betai)*(betai-betai^2))
        Dmuth2[ii] <- Dmuth2i

        Dmu4[ii] <- (-4/si^3)*((2*alphai-3*alphai^2)*(alphai-alphai^2)+
          (2*betai-3*betai^2)*(betai-betai^2))-Dmu3i
        Dmu3th[ii] <- -2*Dmu3i+2*(1/si^3)*((yr-mui)*(alphai-alphai^2)*(1-6*alphai+
          6*alphai^2)+(yl-mui)*(betai-betai^2)*(1-6*betai+6*betai^2))
        Dmu2th2[ii] <- -2*Dmu2thi+(1/si)*(Dmui-Dmuth2i)+(-12/si^3)*((yr-mui)*(alphai^2-alphai^3)+
          (yl-mui)*(betai^2-betai^3))+(-4/si^4)*((yr-mui)^2*(2*alphai-3*alphai^2)*(alphai-alphai^2)+
          (yl-mui)^2*(2*betai-3*betai^2)*(betai-betai^2))
      }
    }

    # left censored
    if (length(il)) {
      si <- exp(theta)/sqrt(wt[il]) # array to number
      yr <- y[il]
      mui <- mu[il]
      alphai <- 1/(2+expm1(-(yr-mui)/si))

      Dmu[il] <- Dmui <- (2/si)*(1-alphai)
      Dmu2[il] <- Dmu2i <- (2/si^2)*(alphai-alphai^2)
      if (level>0) {
        Dmuth[il] <- Dmuthi <- -Dmui+(yr-mui)*Dmu2i
        Dth[il] <- Dthi <- (yr-mui)*Dmui
        Dmu3[il] <- Dmu3i <- (1/si)*(2*alphai-1)*Dmu2i
        Dmu2th[il] <- Dmu2thi <- -2*Dmu2i+(1/si)*(yr-mui)*(2*alphai-1)*Dmu2i
      }
      if (level>1) {
        Dth2[il] <- Dth2i <- (yr-mui)*Dmuthi
        Dmuth2[il] <- Dmuth2i <- -Dmuthi+(yr-mui)*Dmu2thi
        Dmu4[il] <- Dmu4i <- -Dmu2i^2+(1/si)*(2*alphai-1)*Dmu3i
        Dmu3th[il] <- Dmu3thi <- -(yr-mui)*Dmu2i^2+(1/si)*(2*alphai-1)*(Dmu2thi-Dmu2i)
        Dmu2th2[il] <- -2*Dmu2thi-(yr-mui)^2*Dmu2i^2+
          (1/si)*(yr-mui)*(2*alphai-1)*(Dmu2thi-Dmu2i)
      }
    }

    # right censored
    if (length(ir)) {
      si <- exp(theta)/sqrt(wt[ir]) # array to number
      yl <- y[ir]; mui <- mu[ir]
      
      betai <- 1/(2+expm1(-(yl-mui)/si))

      #Dmu[ir] <- Dmui <- (2/si)*betai
      Dmu[ir] <- Dmui <- -(2/si)*betai
      #Dmu2[ir] <- Dmu2i <- (-2/si^2)*(betai-betai^2)
      Dmu2[ir] <- Dmu2i <- (2/si^2)*(betai-betai^2)
      if (level>0) {
        #Dmuth[ir] <- Dmuthi <- (yl-mui)*Dmu2i
        Dmuth[ir] <- Dmuthi <- -Dmui+(yl-mui)*Dmu2i
        Dth[ir] <- Dthi <- (yl-mui)*Dmui
        Dmu3[ir] <- Dmu3i <- (-1/si)*(1-2*betai)*Dmu2i
        Dmu2th[ir] <- Dmu2thi <- -(2+(1/si)*(yl-mui)*(1-2*betai))*Dmu2i
      }
      if (level>1) {
        Dth2[ir] <- Dth2i <- (yl-mui)*Dmuthi
        #Dmuth2[ir] <- Dmuth2i <- (yl-mui)*Dmu2thi
        Dmuth2[ir] <- Dmuth2i <- -Dmuthi+(yl-mui)*Dmu2thi
        #Dmu4[ir] <- Dmu4i <- Dmu2i^2+(-1/si)*(1-2*betai)*Dmu3i
        Dmu4[ir] <- Dmu4i <- -Dmu2i^2+(-1/si)*(1-2*betai)*Dmu3i
        #Dmu3th[ir] <- Dmu3thi <- (1/si)*(1-2*betai)*(Dmu2i-Dmu2thi)+(yl-mui)*Dmu2i^2
        Dmu3th[ir] <- Dmu3thi <- (1/si)*(1-2*betai)*(Dmu2i-Dmu2thi)-(yl-mui)*Dmu2i^2
        #Dmu2th2[ir] <- -2*Dmu2thi-(mui-yl)*Dmu3thi ###
        Dmu2th2[ir] <- (1/si)*(yl-mui)*(1-2*betai)*(Dmu2i-Dmu2thi)+
          -(yl-mui)^2*Dmu2i^2-2*Dmu2thi
      }
    }
    ip <- which(Dmu2<0); EDmu2t <- Dmu2
    EDmu2t[ip] <- 0

    r <- list(Dmu=Dmu, Dmu2=Dmu2, EDmu2=EDmu2t)
    if (level>0) {
      r$Dth <- Dth; r$Dmuth <- Dmuth; r$Dmu3 <- Dmu3
      r$EDmu2th <- r$Dmu2th <- Dmu2th
    }
    if (level>1) {
      r$Dmu4 <- Dmu4; r$Dth2 <- Dth2;  r$Dmuth2 <- Dmuth2;
	    r$Dmu2th2 <- Dmu2th2; r$Dmu3th <- Dmu3th
    }
    return(r)
  } ## Dd clogistic

  # akaike information criterion
  aic <- function(y, mu, theta=NULL, wt, dev) {
    log1pexp <- function(x) { ## compute log(1+e^x) safely
      result <- x
      ii <- which(x<=37); result[ii] <- exp(x[ii])
      ii <- which(-37<x & x<=18); result[ii] <- log1p(exp(x[ii]))
      ii <- which(18<x & x<=33.3); result[ii] <- x[ii]+exp(-x[ii])
      return(result)
    } # log1pexp

    log1mexp <- function(a) { ## compute log(1-e^(-a)) safely
      result <- numeric(length(a))
      ii <- which(0<a & a<=log(2)); result[ii] <- log(-expm1(-a[ii]))
      ii <- which(a>log(2)); result[ii] <- log1p(-exp(-a[ii]))
      return(result)
    } # log1mexp
    if (is.null(theta)) theta <- get(".Theta")
    th <- theta-0.5*log(wt)
    yat <- attr(y,"censor")
    if (is.null(yat)) yat <- rep(NA,length(y))
    a <- rep(0,length(y))

    # get indices
    iu <- which(yat==y) # uncensored
    ii <- which(is.finite(yat*y) & yat!=y) # interval censored
    il <- which(yat==-Inf) # left censored
    ir <- which(yat==Inf) # right censored

    # uncensored
    if (length(iu)) {
      si <- exp(theta)/sqrt(wt[iu])
      lm <- -log1p(si-1)-2*log(2)
      a[iu] <- -2*lm+2*th[iu]
    }

    # interval censored
    if (length(ii)) {
      si <- exp(theta)/sqrt(wt[ii])
      yl <- pmin(y[ii], yat[ii])
      yr <- pmax(y[ii], yat[ii])
      lm <- ((yr-yl)/(2*si))+log1mexp((yr-yl)/(2*si))-log1pexp((yr-yl)/(2*si))
      a[ii] <- -2*lm+2*th[ii]
    }

    # left censored
    if (length(il)) {
      a[il] <- 2*th[il]
    }

    # right censored
    if (length(ir)) {
      a[ir] <- 2*th[ir]
    }
    return(sum(a))
  } ## AIC clogistic

  # saturated log likelihood
  ls <- function(y, w, theta, scale) {
    log1pexp <- function(x) { ## compute log(1+e^x) safely
      result <- x
      ii <- which(x<=37); result[ii] <- exp(x[ii])
      ii <- which(-37<x & x<=18); result[ii] <- log1p(exp(x[ii]))
      ii <- which(18<x & x<=33.3); result[ii] <- x[ii]+exp(-x[ii])
      return(result)
    } # log1pexp

    log1mexp <- function(a) { ## compute log(1-e^(-a)) safely
      result <- numeric(length(a))
      ii <- which(0<a & a<=log(2)); result[ii] <- log(-expm1(-a[ii]))
      ii <- which(a>log(2)); result[ii] <- log1p(-exp(-a[ii]))
      return(result)
    } # log1mexp
    yat <- attr(y, "censor")
    if (is.null(yat)) yat <- y

    # get indices
    iu <- which(yat==y) # uncensored
    ii <- which(is.finite(yat*y) & yat!=y) # interval censored
    l2 <- l1 <- l <- rep(0, length(y))

    # uncensored
    if (length(iu)) {
      si <- exp(theta)/sqrt(w[iu])
      l[iu] <- -log1p(si-1)-2*log(2)
      l1[iu] <- -1
    }

    # interval censored
    if (length(ii)) {
      si <- exp(theta)/sqrt(w[ii])
      yl <- pmin(y[ii], yat[ii])
      yr <- pmax(y[ii], yat[ii])

      l[ii] <- ((yr-yl)/(2*si))+log1mexp((yr-yl)/(2*si))-log1pexp((yr-yl)/(2*si))
      lmth <- (yr-yl)*(1/si)*(1/(expm1(-(yr-yl)/(2*si))-expm1((yr-yl)/(2*si))))
      l1[ii] <- lmth
      lmth2 <- -lmth-(1/2)*(yr-yl)^2*(1/si^2)*(1/(expm1(-(yr-yl)/(2*si))+
        -expm1((yr-yl)/(2*si))))^2*(2+expm1(-(yr-yl)/(2*si))+expm1((yr-yl)/(2*si)))
      l2[ii] <- lmth2
    }
    # right and left censored are zero
    result <- list(ls=sum(l), lsth1=sum(l1), LSTH1=matrix(l1,ncol=1),lsth2=sum(l2))
    return(result)
  } ## ls clogistic

  initialize <- expression({
    if (is.matrix(y)) {
      .yat <- y[,2]
      y <- y[,1]
      attr(y, "censor") <- .yat
    }
    mustart <- if (family$link=="identity") y else pmax(y, min(y>0))
  })

  postproc <- function(family, y, prior.weights, fitted, linear.predictors, offset, intercept) {
    posr <- list()
    if (is.matrix(y)) {
      .yat <- y[,2]
      y <- y[,1]
      attr(y, "censor") <- .yat
    }
    posr$null.deviance <- find.null.dev(family, y, eta=linear.predictors, offset, prior.weights)
    posr$family <- paste("clogistic(",round(family$getTheta(TRUE),3),")",sep="")
    posr
  } ## postproc clogistic

  rd <- function(mu, wt, scale) { # incomplete
    Theta <- exp(get(".Theta"))
  }

  qf <- function(p ,mu, wt, scale) { # incomplete
    Theta <- exp(get(".Theta"))
  }

  subsety <- function(y, ind) {
    if (is.matrix(y)) return(y[ind,])
    yat <- attr(y, "censor")
    y <- y[ind]
    if (!is.null(yat)) attr(y, "censor") <- yat[ind]
    y
  } ## subsety

  environment(dev.resids) <- environment(aic) <- environment(getTheta) <-
    environment(rd) <- environment(qf) <- environment(putTheta) <- env

  structure(list(family="clogistic", link=linktemp, linkfun=stats$linkfun,
    linkinv=stats$linkinv, dev.resids=dev.resids, Dd=Dd, subsety=subsety,
    aic=aic, mu.eta=stats$mu.eta, initialize=initialize, postproc=postproc,
    ls=ls, validmu=validmu, valideta=stats$valideta, n.theta=n.theta, 
    ini.theta=iniTheta, putTheta=putTheta, getTheta=getTheta),
    class = c("extended.family","family"))
} ## clog

#################################################
## extended family object for ordered categorical
#################################################

ocat <- function(theta=NULL,link="identity",R=NULL) {
## extended family object for ordered categorical model.
## one of theta and R must be supplied. length(theta) == R-2.
## weights are ignored. #! is stuff removed under re-definition of ls as 0
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("identity")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for ocat family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)

#  else if (is.character(link)) {
#    stats <- make.link(link)
#    linktemp <- link
#  } else {
#    if (inherits(link, "link-glm")) {
#      stats <- link
#      if (!is.null(stats$name))
#            linktemp <- stats$name
#    } else stop(linktemp, " link not available for ordered categorical family; available links are \"identity\"")
#  }
  if (is.null(theta)&&is.null(R)) stop("Must supply theta or R to ocat")
  if (!is.null(theta)) R <- length(theta) + 2 ## number of catergories
  ## Theta <-  NULL;
  n.theta <- R-2
  ## NOTE: data based initialization is in preinitialize...
  if (!is.null(theta)&&sum(theta==0)==0) {
    if (sum(theta<0)) iniTheta <- log(abs(theta)) ## initial theta supplied
    else { 
      iniTheta <- log(theta) ## fixed theta supplied
      n.theta <- 0
    }
  } else iniTheta <- rep(-1,length=R-2) ## inital log theta value

  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", iniTheta, envir = env)

  putTheta <-function(theta) assign(".Theta", theta,envir=environment(sys.function()))

  getTheta <-function(trans=FALSE) { 
    theta <- get(".Theta")
    if (trans) { ## transform to (finite) cut points...
      R = length(theta)+2
      alpha <- rep(0,R-1) ## the thresholds
      alpha[1] <- -1
      if (R > 2) { 
        ind <- 2:(R-1)
        alpha[ind] <- alpha[1] + cumsum(exp(theta))
      } 
      theta <- alpha
    }
    theta
  } ## getTheta ocat

  postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept) {
    posr <- list()
    ## null.deviance needs to be corrected...
    posr$null.deviance <- find.null.dev(family,y,eta=linear.predictors,offset,prior.weights)
    posr$family <- 
    paste("Ordered Categorical(",paste(round(family$getTheta(TRUE),2),collapse=","),")",sep="")
    posr
  }

  validmu <- function(mu) all(is.finite(mu))

  dev.resids <- function(y, mu, wt,theta=NULL) { ## ocat dev resids
   
    Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    if (is.null(theta)) theta <- get(".Theta")
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y] ## cut points above and below y
    ## Compute sign for deviance residuals, based on sign of interval
    ## midpoint for each datum minus the fitted value of the latent 
    ## variable. This makes first and last categories like 0s and 1s in
    ## logistic regression....
    s <- sign((al1 + al0)/2-mu) ## sign for deviance residuals
    al1mu <- al1-mu;al0mu <- al0-mu
  
    f <- Fdiff(al0mu,al1mu)

    rsd <- -2*wt*log(f) 
    attr(rsd,"sign") <- s
    rsd
  } ## ocat dev.resids

  Dd <- function(y, mu, theta, wt=NULL, level=0) {
  ## derivatives of the ocat deviance...

    Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a 
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    abcd <- function(x,level=-1) {
      bj <- cj <- dj <- NULL
      ## compute f_j^2 - f_j without cancellation error
      ## x <- 10;F(x)^2-F(x);abcd(x)$aj
      h <- rep(1,length(x)); h[x>0] <- -1; ex <- exp(x*h)
      ex1 <- ex+1;ex1k <- ex1^2
      aj <- -ex/ex1k
      if (level>=0) {
        ## compute f_j - 3 f_j^2 + 2f_j^3 without cancellation error
        ex1k <- ex1k*ex1;ex2 <- ex^2
        bj <- h*(ex-ex^2)/ex1k
        if (level>0) {
          ## compute -f_j + 7 f_j^2 - 12 f_j^3 + 6 f_j^4
          ex1k <- ex1k*ex1;ex3 <- ex2*ex
          cj <- (-ex3 + 4*ex2 - ex)/ex1k
          if (level>1) {    
            ## compute d_j
            ex1k <- ex1k*ex1;ex4 <- ex3*ex
            dj <- h * (-ex4 + 11*ex3 - 11*ex2 + ex)/ex1k
          }
        }
      }        
      list(aj=aj,bj=bj,cj=cj,dj=dj)
    }
    if (is.null(wt)) wt <- 1

    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y]
    al1mu <- al1-mu;al0mu <- al0 - mu
   
    f <- pmax(Fdiff(al0mu,al1mu),.Machine$double.xmin)
    r1 <- abcd(al1mu,level); a1 <- r1$aj
    r0 <- abcd(al0mu,level); a0 <- r0$aj
    a <- a1 - a0
    
    if (level>=0) {
      b1 <- r1$bj;b0 <- r0$bj
      b <- b1 - b0
    }
    if (level>0) {
      c1 <- r1$cj;c0 <- r0$cj
      c <- c1 - c0
    }
    if (level>1) {
      d1 <- r1$dj;d0 <- r0$dj
      d <- d1 - d0
    }

    oo <- list(D=NULL,Dmu=NULL,Dmu2=NULL,Dth=NULL,Dmuth=NULL,
             Dmu2th=NULL)
    n <- length(y)
    ## deviance...
    oo$D <- -2 * wt * log(f)
    if (level >= 0) { ## get derivatives used in coefficient estimation
      oo$Dmu <- -2 * wt * a / f
      a2 <- a^2
      oo$EDmu2 <- oo$Dmu2 <- 2 * wt * (a2/f -  b)/f
    }
    if (R<3) level <- 0 ## no free parameters

    if (level > 0) { ## get first derivative related stuff
      f2 <- f^2;a3 <- a2*a
      oo$Dmu3 <- 2 * wt * (- c - 2 * a3/f2 + 3 * a * b/f)/f
      Dmua0 <- 2 *  (a0 * a/f -  b0)/f
      Dmua1 <- -2 *  (a1 * a /f - b1)/f
      Dmu2a0 <- -2 * (c0 + (a0*(2*a2/f - b)- 2*b0*a  )/f)/f
      Dmu2a1 <- 2 * (c1  + (2*(a1*a2/f - b1*a) - a1*b)/f)/f
      Da0 <- -2*a0/f
      Da1 <-  2 * a1/f 
      ## now transform to derivatives w.r.t. theta...
      oo$Dmu2th <- oo$Dmuth <- oo$Dth <- matrix(0,n,R-2)
      for (k in 1:(R-2)) { 
        etk <- exp(theta[k])
        ind <- y == k+1;wti <- wt[ind]
        oo$Dth[ind,k] <- wti * Da1[ind]*etk
        oo$Dmuth[ind,k] <- wti * Dmua1[ind]*etk
        oo$Dmu2th[ind,k] <- wti * Dmu2a1[ind]*etk
    
        if (R>k+2) { 
          ind <- y > k+1 & y < R;wti <- wt[ind]
          oo$Dth[ind,k] <- wti * (Da1[ind]+Da0[ind])*etk
          oo$Dmuth[ind,k] <- wti * (Dmua1[ind]+Dmua0[ind])*etk
          oo$Dmu2th[ind,k] <- wti * (Dmu2a1[ind]+Dmu2a0[ind])*etk
       
        }
        ind <- y == R;wti <- wt[ind]
        oo$Dth[ind,k] <- wti * Da0[ind]*etk
        oo$Dmuth[ind,k] <- wti * Dmua0[ind]*etk
        oo$Dmu2th[ind,k] <- wti * Dmu2a0[ind]*etk 
      }
      oo$EDmu2th <- oo$Dmu2th
      
    }  
    if (level >1) { ## and the second derivative components 
      oo$Dmu4 <- 2 * wt * ((3*b^2 + 4*a*c)/f + a2*(6*a2/f - 12*b)/f2 - d)/f
      Dmu3a0 <-  2 * ((a0*c  + 3*c0*a + 3*b0*b)/f - d0  + 
                      6*a*(a0*a2/f - b0*a - a0*b)/f2 )/f
      Dmu3a1 <- 2 * (d1 - (a1*c + 3*(c1*a + b1*b))/f
                     + 6*a*(b1*a - a1*a2/f  + a1*b)/f2)/f

      Dmua0a0 <- 2*(c0 + (2*a0*(b0 - a0*a/f) - b0*a)/f )/f
      Dmua1a1 <- 2*( (b1*a + 2*a1*(b1 - a1*a/f))/f - c1)/f
      Dmua0a1 <- 2*(a0*(2*a1*a/f - b1) - b0*a1)/f2 

      Dmu2a0a0 <- 2*(d0 + (b0*(2*b0 - b)  + 2*c0*(a0 - a))/f +
                     2*(b0*a2 + a0*(3*a0*a2/f  - 4*b0*a - a0*b))/f2)/f

      Dmu2a1a1 <-  2*( (2*c1*(a + a1) + b1*(2*b1 + b))/f
                + 2*(a1*(3*a1*a2/f  - a1*b) - b1*a*(a + 4*a1))/f2 - d1)/f

      Dmu2a0a1 <- 0

      Da0a0 <- 2 * (b0 + a0^2/f)/f
      Da1a1 <- -2* (b1 - a1^2/f)/f 
      Da0a1 <- -2*a0*a1/f2

      ## now transform to derivatives w.r.t. theta...
      n2d <- (R-2)*(R-1)/2
      oo$Dmu3th <- matrix(0,n,R-2)
      oo$Dmu2th2 <- oo$Dmuth2 <- oo$Dth2 <- matrix(0,n,n2d)
      i <- 0
      for (j in 1:(R-2)) for (k in j:(R-2)) { 
        i <- i + 1 ## the second deriv col
        ind <- y >= j ## rest are zero
	wti <- wt[ind]
        ar1.k <- ar.k <- rep(exp(theta[k]),n)
        ar.k[y==R | y <= k] <- 0; ar1.k[y<k+2] <- 0
        ar.j <- ar1.j <- rep(exp(theta[j]),n)
        ar.j[y==R | y <= j] <- 0; ar1.j[y<j+2] <- 0
        ar.kj <- ar1.kj <- rep(0,n)
        if (k==j) {
          ar.kj[y>k&y<R] <- exp(theta[k])
          ar1.kj[y>k+1] <- exp(theta[k])
          oo$Dmu3th[ind,k] <- wti * (Dmu3a1[ind]*ar.k[ind]  + Dmu3a0[ind]*ar1.k[ind])
        }
        oo$Dth2[,i] <- wt * (Da1a1*ar.k*ar.j + Da0a1*ar.k*ar1.j + Da1 * ar.kj +
                     Da0a0*ar1.k*ar1.j + Da0a1*ar1.k*ar.j + Da0 * ar1.kj)
        oo$Dmuth2[,i] <- wt * (Dmua1a1*ar.k*ar.j + Dmua0a1*ar.k*ar1.j + Dmua1 * ar.kj +
                     Dmua0a0*ar1.k*ar1.j + Dmua0a1*ar1.k*ar.j + Dmua0 * ar1.kj)
        oo$Dmu2th2[,i] <- wt * (Dmu2a1a1*ar.k*ar.j + Dmu2a0a1*ar.k*ar1.j + Dmu2a1 * ar.kj +
                     Dmu2a0a0*ar1.k*ar1.j + Dmu2a0a1*ar1.k*ar.j + Dmu2a0 * ar1.kj)
      } 
    }
    oo
  } ## Dd (ocat)
 
  aic <- function(y, mu, theta=NULL, wt, dev) {
  
   Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a 
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    if (is.null(theta)) theta <- get(".Theta")
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y]
    ##f1 <- F(al1-mu);f0 <- F(al0-mu);f <- f1 - f0
    f <- Fdiff(al0-mu,al1-mu)
    -2*sum(log(f)*wt)
  } ## ocat aic

  ls <- function(y,w,theta,scale) {
    ## the log saturated likelihood function. 
    return(list(ls=0,lsth1=rep(0,R-2),LSTH1=matrix(0,length(y),R-2),lsth2=matrix(0,R-2,R-2)))
  } ## ocat ls
  
  ## initialization is interesting -- needs to be with reference to initial cut-points
  ## so that mu puts each obs in correct category initially...
 
  preinitialize <- function(y,family) {
    ocat.ini <- function(R,y) {
    ## initialize theta from raw counts in each category
      if (R<3) return()
      y <- c(1:R,y) ## make sure there is *something* in each class
      p <- cumsum(tabulate(y[is.finite(y)])/length(y[is.finite(y)]))
      eta <- if (p[1]==0) 5 else -1 - log(p[1]/(1-p[1])) ## mean of latent
      theta <- rep(-1,R-1) 
      for (i in 2:(R-1)) theta[i] <- log(p[i]/(1-p[i])) + eta 
      theta <- diff(theta)
      theta[theta <= 0.01] <- 0.01
      theta <- log(theta)
    }
    R3 <- length(family$getTheta())+2
    if (!is.numeric(y)) stop("Response should be integer class labels")
    if (R3>2&&family$n.theta>0) { 
      Theta <- ocat.ini(R3,y)
      return(list(Theta=Theta))
    } 
  } ## ocat preinitialize

  initialize <- expression({ 
    R <- length(family$getTheta())+2 ## don't use n.theta as it's used to signal fixed theta
    if (any(y < 1)||any(y>R)) stop("values out of range")
    ## n <- rep(1, nobs)
    ## get the cut points... 
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -2;alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(family$getTheta()))
    } 
    alpha[R+1] <- alpha[R] + 1
    mustart <- (alpha[y+1] + alpha[y])/2
  })

  residuals <- function(object,type=c("deviance","working","response")) {
    if (type == "working") { 
      res <- object$residuals 
    } else if (type == "response") {
       theta <- object$family$getTheta()
       mu <- object$linear.predictors
       R = length(theta)+2
       alpha <- rep(0,R+1) ## the thresholds
       alpha[1] <- -Inf;alpha[R+1] <- Inf
       alpha[2] <- -1
       if (R > 2) { 
         ind <- 3:R
         alpha[ind] <- alpha[2] + cumsum(exp(theta))
       } 
       fv <- mu*NA
       for (i in 1:(R+1)) {
         ind <- mu>alpha[i] & mu<=alpha[i+1]
         fv[ind] <- i
       } 
       res <- object$y - fv
    } else if (type == "deviance") { 
      y <- object$y
      mu <- object$fitted.values
      wts <- object$prior.weights
      res <- object$family$dev.resids(y,mu,wts)
      s <- attr(res,"sign")
      if (is.null(s)) s <- sign(y-mu)
      res <- as.numeric(sqrt(pmax(res,0)) * s) 
      
    }
    res
  } ## ocat residuals

 
  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
    ocat.prob <- function(theta,lp,se=NULL) {
    ## compute probabilities for each class in ocat model
    ## theta is finite cut points, lp is linear predictor, se 
    ## is standard error on lp...
      R <- length(theta) 
      dp <- prob <- matrix(0,length(lp),R+2)
      prob[,R+2] <- 1
      for (i in 1:R) {
        x <- theta[i] - lp
        ind <- x > 0
        prob[ind,i+1] <- 1/(1+exp(-x[ind]))
        ex <- exp(x[!ind])
        prob[!ind,i+1] <- ex/(1+ex)
        dp[,i+1] <- prob[,i+1]*(prob[,i+1]-1)
      }
      prob <- t(diff(t(prob)))
      dp <- t(diff(t(dp))) ## dprob/deta
      if (!is.null(se)) se <- as.numeric(se)*abs(dp)
      list(prob,se)
    } ## ocat.prob

    theta <- family$getTheta(TRUE)
    if (is.null(eta)) { ## return probabilities
      discrete <- is.list(X) 
      mu <- off + if (discrete) Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop) else drop(X%*%beta)  
      if (se) {
        se <- if (discrete) sqrt(pmax(0,diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1))) else
	                    sqrt(pmax(0,rowSums((X%*%Vb)*X)))
      } else se <- NULL
      ##theta <- cumsum(c(-1,exp(theta)))
      p <- ocat.prob(theta,mu,se)
      if (is.null(se)) return(p) else { ## approx se on prob also returned
        names(p) <- c("fit","se.fit")
        return(p)
      } 
    } else { ## return category implied by eta (i.e mean of latent)
      R = length(theta)+2
      #alpha <- rep(0,R) ## the thresholds
      #alpha[1] <- -Inf;alpha[R] <- Inf
      alpha <- c(-Inf, theta, Inf)
      fv <- eta*NA
      for (i in 1:(R+1)) {
        ind <- eta>alpha[i] & eta<=alpha[i+1]
        fv[ind] <- i
      } 
      return(list(fv))
    }
  } ## ocat predict

  rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu 
    theta <- get(".Theta") 
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    ## ... cut points computed, now simulate latent variable, u
    y <- u <- runif(length(mu))
    u <- mu + log(u/(1-u)) 
    ## and allocate categories according to u and cut points...
    for (i in 1:R) {
      y[u > alpha[i]&u <= alpha[i+1]] <- i
    }
    y
  } ## ocat rd

  environment(dev.resids) <- environment(aic) <- environment(putTheta) <-
  environment(getTheta) <- environment(rd) <- environment(predict) <- env
  structure(list(family = "Ordered Categorical", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,
        preinitialize = preinitialize, ls=ls,rd=rd,residuals=residuals,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta,
        ini.theta = iniTheta,putTheta=putTheta,predict=predict,step = 1,
        getTheta=getTheta,no.r.sq=TRUE), class = c("extended.family","family"))
} ## end of ocat



## Tweedie....


tw <- function (theta = NULL, link = "log",a=1.01,b=1.99) { 
## Extended family object for Tweedie, to allow direct estimation of p
## as part of REML optimization. 
## p = (a+b*exp(theta))/(1+exp(theta)), i.e. a < p < b
## NOTE: The Tweedie density computation at low phi, low p is susceptible
##       to cancellation error, which seems unavoidable. Furthermore 
##       there are known problems with spurious maxima in the likelihood 
##       w.r.t. p when the data are rounded. 
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt","inverse")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for tw family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)

#  else if (is.character(link)) {
#    stats <- make.link(link)
#    linktemp <- link
#  } else {
#    if (inherits(link, "link-glm")) {
#       stats <- link
#            if (!is.null(stats$name))
#                linktemp <- stats$name
#        }
#        else  stop(gettextf("link \"%s\" not available for Tweedie family.", 
#                linktemp, collapse = ""), domain = NA)
#  }
  ## Theta <-  NULL;
  n.theta <- 1
  if (!is.null(theta)&&theta!=0) {
      if (abs(theta)<=a||abs(theta)>=b) stop("Tweedie p must be in interval (a,b)")
      if (theta>0) { ## fixed theta supplied
        iniTheta <- log((theta-a)/(b-theta)) 
        n.theta <- 0 ## so no theta to estimate
      } else iniTheta <- log((-theta-a)/(b+theta)) ## initial theta supplied
  } else iniTheta <- 0 ## inital log theta value
    
  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", iniTheta, envir = env) 
  assign(".a",a, envir = env);assign(".b",b, envir = env)
  getTheta <- function(trans=FALSE) { 
  ## trans transforms to the original scale...
    th <- get(".Theta")
    a <- get(".a");b <- get(".b")
    if (trans) th <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    th
  }
  putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

  validmu <- function(mu) all(mu > 0)
  
  variance <- function(mu) { 
    th <- get(".Theta");a <- get(".a");b <- get(".b")
    p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    mu^p
  }
  
  dev.resids <- function(y, mu, wt,theta=NULL) {
    if (is.null(theta)) theta <- get(".Theta")
    a <- get(".a");b <- get(".b")
    p <- if (theta>0) (b+a*exp(-theta))/(1+exp(-theta)) else (b*exp(theta)+a)/(exp(theta)+1)
    y1 <- y + (y == 0)
    theta <- if (p == 1) log(y1/mu) else (y1^(1 - p) - mu^(1 - p))/(1 - p)
    kappa <- if (p == 2) log(y1/mu) else (y^(2 - p) - mu^(2 - p))/(2 - p)
    pmax(2 * (y * theta - kappa) * wt,0)
  } ## tw dev.resids
    
  Dd <- function(y, mu, theta, wt, level=0) {
  ## derivatives of the tw deviance...
    a <- get(".a");b <- get(".b")
    th <- theta
    p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    dpth1 <- if (th>0) exp(-th)*(b-a)/(1+exp(-th))^2 else exp(th)*(b-a)/(exp(th)+1)^2
    dpth2 <- if (th>0) ((a-b)*exp(-th)+(b-a)*exp(-2*th))/(exp(-th)+1)^3 else
                   ((a-b)*exp(2*th)+(b-a)*exp(th))/(exp(th)+1)^3
    mu1p <- mu^(1-p)
    mup <- mu^p
    r <- list()
    ## get the quantities needed for IRLS. 
    ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
    ## Dmu is deriv w.r.t. mu once, etc...
    ymupi <- y/mup
    r$Dmu <- 2*wt*(mu1p - ymupi)
    r$Dmu2 <- 2*wt*(mu^(-1-p)*p*y + (1-p)/mup)
    r$EDmu2 <- (2*wt)/mup ## expected Dmu2 (weight)
   
    if (level>0) { ## quantities needed for first derivatives
        i1p <- 1/(1-p)
        y1 <- y + (y==0)
        ##ylogy <- y*log(y1)
        logmu <- log(mu)
        mu2p <- mu * mu1p
        r$Dth <- 2 * wt * ( (y^(2-p)*log(y1) - mu2p*logmu)/(2-p) + 
                            (y*mu1p*logmu - y^(2-p)*log(y1))/(1-p) -
                            (y^(2-p) - mu2p)/(2-p)^2 + 
                            (y^(2-p) - y*mu1p)*i1p^2) *dpth1       

        r$Dmuth <- 2 * wt * logmu * (ymupi - mu1p)*dpth1
        mup1 <-  mu^(-p-1)
        r$Dmu3 <- -2 * wt * mup1*p*(y/mu*(p+1) + 1-p)    
        r$Dmu2th <- 2 * wt  * (mup1*y*(1-p*logmu)-(logmu*(1-p)+1)/mup )*dpth1
	r$EDmu3 <- -2*wt*p*mup1
	r$EDmu2th <- -2*wt*logmu/mup*dpth1
      } 
      if (level>1) { ## whole damn lot
        mup2 <- mup1/mu
        r$Dmu4 <- 2 * wt * mup2*p*(p+1)*(y*(p+2)/mu + 1 - p)
        y2plogy <- y^(2-p)*log(y1);y2plog2y <- y2plogy*log(y1)
       
        r$Dth2 <- 2 * wt * (((mu2p*logmu^2-y2plog2y)/(2-p) + (y2plog2y - y*mu1p*logmu^2)/(1-p) +
                             2*(y2plogy-mu2p*logmu)/(2-p)^2 + 2*(y*mu1p*logmu-y2plogy)/(1-p)^2
                             + 2 * (mu2p - y^(2-p))/(2-p)^3+2*(y^(2-p)-y*mu^(1-p))/(1-p)^3)*dpth1^2) +
                              r$Dth*dpth2/dpth1
        
        r$Dmuth2 <- 2 * wt * ((mu1p * logmu^2 - logmu^2*ymupi)*dpth1^2) + r$Dmuth*dpth2/dpth1

        r$Dmu2th2 <- 2 * wt * ( (mup1 * logmu*y*(logmu*p - 2) + logmu/mup*(logmu*(1-p) + 2))*dpth1^2) +
                                r$Dmu2th * dpth2/dpth1

        r$Dmu3th <- 2 * wt * mup1*(y/mu*(logmu*(1+p)*p-p -p-1) +logmu*(1-p)*p + p - 1 + p)*dpth1
      }
      r
    } ## Dd tw

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")  
        a <- get(".a");b <- get(".b")
        p <- if (theta>0) (b+a*exp(-theta))/(1+exp(-theta)) else (b*exp(theta)+a)/(exp(theta)+1)
        scale <- dev/sum(wt)
        -2 * sum(ldTweedie(y, mu, p = p, phi = scale)[, 1] * 
            wt) + 2
    } ## aic tw

    ls <- function(y, w, theta, scale) {
        ## evaluate saturated log likelihood + derivs w.r.t. working params and log(scale) Tweedie
        a <- get(".a");b <- get(".b")
        Ls <- w * ldTweedie(y, y, rho=log(scale), theta=theta,a=a,b=b)
	LS <- colSums(Ls)
	lsth1 <- c(LS[4],LS[2]) ## deriv w.r.t. p then log scale
        lsth2 <- matrix(c(LS[5],LS[6],LS[6],LS[3]),2,2)
        list(ls=LS[1],lsth1=lsth1,LSTH1=Ls[,c(4,2)],lsth2=lsth2)
    } ## ls tw

 
    initialize <- expression({
        ##n <- rep(1, nobs)
        mustart <- y + (y == 0)*.1
    })
    
    postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept) {
      posr <- list()
      posr$null.deviance <- find.null.dev(family,y,eta=linear.predictors,offset,prior.weights)
      posr$family <- 
      paste("Tweedie(p=",round(family$getTheta(TRUE),3),")",sep="")
      posr
    } ## postproc tw
  
    rd <- function(mu,wt,scale) {
     th <- get(".Theta") 
     a <- get(".a");b <- get(".b")
     p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)

     if (p == 2) 
            rgamma(n=length(mu), shape = 1/scale, scale = mu * scale)
     else
            rTweedie(mu, p = p, phi = scale)
    } ## rd tw

     environment(Dd) <- environment(ls) <-
     environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
      environment(rd) <- environment(variance) <- environment(putTheta) <- env
    structure(list(family = "Tweedie", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,rd=rd,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,canonical="none",n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,scale = -1),
        class = c("extended.family","family"))
} ## tw

## beta regression

betar <- function (theta = NULL, link = "logit",eps=.Machine$double.eps*100) { 
## Extended family object for beta regression
## length(theta)=1; log theta supplied
## This serves as a prototype for working with -2logLik
## as deviance, and only dealing with saturated likelihood 
## at the end.
## Written by Natalya Pya. 'saturated.ll' by Simon Wood
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("logit", "probit", "cloglog", "cauchit")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for betar family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)

#  else if (is.character(link)) {
#    stats <- make.link(link)
#    linktemp <- link
#  } else {
#    if (inherits(link, "link-glm")) {
#       stats <- link
#            if (!is.null(stats$name))
#                linktemp <- stats$name
#        }
#        else stop(linktemp, " link not available for beta regression; available links are  \"logit\", \"probit\", \"cloglog\" and \"cauchit\"")
#    }
   
    n.theta <- 1
    if (!is.null(theta)&&theta!=0) {
       if (theta>0) {
           iniTheta <- log(theta) ## fixed theta supplied
           n.theta <- 0 ## signal that there are no theta parameters to estimate
       } else iniTheta <- log(-theta) ## initial theta supplied
    } else iniTheta <- 0 ##  inital log theta value
     
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    assign(".betarEps",eps, envir = env)
    getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta") 
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) { 
        th <- get(".Theta")
        mu*(1 - mu)/(1+exp(th))
    }
  
    validmu <- function(mu) all(mu > 0 & mu < 1)

    dev.resids <- function(y, mu, wt,theta=NULL) {
    ## '-2*loglik' instead of deviance in REML/ML expression
      if (is.null(theta)) theta <- get(".Theta")
      theta <- exp(theta) ## note log theta supplied
      muth <- mu*theta
      ## yth <- y*theta
      2* wt * (-lgamma(theta) +lgamma(muth) + lgamma(theta - muth) - muth*(log(y)-log1p(-y)) - 
                theta*log1p(-y) + log(y) + log1p(-y)) 
    } ## dev.resids betar
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the -2*loglik...
      ## ltheta <- theta
      theta <- exp(theta)
      onemu <- 1 - mu; ## oney <- 1 - y
      muth <- mu*theta; ## yth <- y*theta
      onemuth <- onemu*theta  ## (1-mu)*theta
      psi0.th <- digamma(theta)
      psi1.th <- trigamma(theta)
      psi0.muth <- digamma(muth) 
      psi0.onemuth <- digamma(onemuth)
      psi1.muth <- trigamma(muth)
      psi1.onemuth <- trigamma(onemuth)
      psi2.muth <- psigamma(muth,2)
      psi2.onemuth <- psigamma(onemuth,2)
      psi3.muth <- psigamma(muth,3)
      psi3.onemuth <- psigamma(onemuth,3)
      log.yoney <- log(y)-log1p(-y)
      r <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      r$Dmu <- 2 * wt * theta* (psi0.muth - psi0.onemuth - log.yoney)
      r$Dmu2 <- 2 * wt * theta^2*(psi1.muth+psi1.onemuth)
      r$EDmu2 <- r$Dmu2
      if (level>0) { ## quantities needed for first derivatives
        r$Dth <- 2 * wt *theta*(-mu*log.yoney - log1p(-y)+ mu*psi0.muth+onemu*psi0.onemuth -psi0.th) 
        r$Dmuth <- r$Dmu + 2 * wt * theta^2*(mu*psi1.muth -onemu*psi1.onemuth)
        r$Dmu3 <- 2 * wt *theta^3 * (psi2.muth - psi2.onemuth) 
        r$EDmu2th <- r$Dmu2th <- 2* r$Dmu2 + 2 * wt * theta^3* (mu*psi2.muth + onemu*psi2.onemuth)
      } 
      if (level>1) { ## whole lot
        r$Dmu4 <- 2 * wt *theta^4 * (psi3.muth+psi3.onemuth) 
        r$Dth2 <- r$Dth +2 * wt *theta^2* (mu^2*psi1.muth+ onemu^2*psi1.onemuth-psi1.th)
        r$Dmuth2 <- r$Dmuth + 2 * wt *theta^2* (mu^2*theta*psi2.muth+ 2*mu*psi1.muth -
                    theta*onemu^2*psi2.onemuth - 2*onemu*psi1.onemuth)
        r$Dmu2th2 <- 2*r$Dmu2th + 2* wt * theta^3* (mu^2*theta*psi3.muth +3*mu*psi2.muth+ 
                    onemu^2*theta*psi3.onemuth + 3*onemu*psi2.onemuth )
        r$Dmu3th <- 3*r$Dmu3 + 2 * wt *theta^4*(mu*psi3.muth-onemu*psi3.onemuth)
      }
      r
    } ## Dd betar

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        theta <- exp(theta)
        muth <- mu*theta
        term <- -lgamma(theta)+lgamma(muth)+lgamma(theta-muth)-(muth-1)*log(y)-
               (theta-muth-1)*log1p(-y) ## `-' log likelihood for each observation
        2 * sum(term * wt)
    } ## aic betar
    
    ls <- function(y,w,theta,scale) {
       ## the log saturated likelihood function for betar
       ## ls is defined as zero for REML/ML expression as deviance is defined as -2*log.lik 
       list(ls=0,## saturated log likelihood
            lsth1=0,  ## first deriv vector w.r.t theta - last element relates to scale
	    LSTH1 = matrix(0,length(y),1),
            lsth2=0) ##Hessian w.r.t. theta
     } ## ls betar

    preinitialize <- function(y,family) {
      eps <- get(".betarEps")
      y[y >= 1-eps] <- 1 - eps;y[y<= eps] <- eps
      return(list(y=y))
    }

    saturated.ll <- function(y,wt,theta=NULL){
    ## function to find the saturated loglik by Newton method,
    ## searching for the mu (on logit scale) that max loglik given theta and data...

      gbh <- function(y,eta,phi,deriv=FALSE,a=1e-8,b=1-a) {
      ## local function evaluating log likelihood (l), gradient and second deriv 
      ## vectors for beta... a and b are min and max mu values allowed. 
      ## mu = (a + b*exp(eta))/(1+exp(eta))
        ind <- eta>0
        expeta <- mu <- eta;
        expeta[ind] <- exp(-eta[ind]);expeta[!ind] <- exp(eta[!ind])
        mu[ind] <- (a*expeta[ind] + b)/(1+expeta[ind])
        mu[!ind] <- (a + b*expeta[!ind])/(1+expeta[!ind])
        l <- dbeta(y,phi*mu,phi*(1-mu),log=TRUE)
        if (deriv) {
          g <- phi * log(y) - phi * log1p(-y) - phi * digamma(mu*phi) + phi * digamma((1-mu)*phi)
          h <- -phi^2*(trigamma(mu*phi)+trigamma((1-mu)*phi))
          dmueta2 <- dmueta1 <-  eta
          dmueta1 <- expeta*(b-a)/(1+expeta)^2 
          dmueta2 <- sign(eta)* ((a-b)*expeta+(b-a)*expeta^2)/(expeta+1)^3
          h <- h * dmueta1^2 + g * dmueta2
          g <- g * dmueta1
        } else g=h=NULL
        list(l=l,g=g,h=h,mu=mu)
      } ## gbh 
      ## now Newton loop...
      eps <- get(".betarEps")
      eta <- y
      a <- eps;b <- 1 - eps
      y[y<eps] <- eps;y[y>1-eps] <- 1-eps
      eta[y<=eps*1.2] <- eps *1.2
      eta[y>=1-eps*1.2] <- 1-eps*1.2
      eta <- log((eta-a)/(b-eta)) 
      mu <- LS <- ii <- 1:length(y)
      for (i in 1:200) {
        ls <- gbh(y,eta,theta,TRUE,a=eps/10)
        conv <- abs(ls$g)<mean(abs(ls$l)+.1)*1e-8
        if (sum(conv)>0) { ## some convergences occured
          LS[ii[conv]] <- ls$l[conv] ## store converged
          mu[ii[conv]] <- ls$mu[conv] ## store mu at converged
          ii <- ii[!conv] ## drop indices
          if (length(ii)>0) { ## drop the converged
            y <- y[!conv];eta <- eta[!conv]
            ls$l <- ls$l[!conv];ls$g <- ls$g[!conv];ls$h <- ls$h[!conv]
          } else break ## nothing left to do
        }
        h <- -ls$h
        hmin <- max(h)*1e-4 
        h[h<hmin] <- hmin ## make +ve def
        delta <- ls$g/h   ## step
        ind <- abs(delta)>2
        delta[ind] <- sign(delta[ind])*2 ## step length limit
        ls1 <- gbh(y,eta+delta,theta,FALSE,a=eps/10); ## did it work?
        ind <- ls1$l<ls$l ## failure index
        k <- 0
        while (sum(ind)>0&&k<20) { ## step halve only failed steps
          k <- k + 1
          delta[ind] <- delta[ind]/2
          ls1$l[ind] <- gbh(y[ind],eta[ind]+delta[ind],theta,FALSE,a=eps/10)$l
          ind <- ls1$l<ls$l
        }
        eta <- eta + delta
      } ## end newton loop
      if (length(ii)>0) { 
        LS[ii] <- ls$l
        warning("saturated likelihood may be inaccurate")
      }

      list(f=sum(wt*LS),term=LS,mu=mu) ## fields f (sat lik) and term (individual datum sat lik) expected
    } ## saturated.ll betar


    ## computes deviance, null deviance, family label
    ## requires prior weights, family, y, fitted values, offset, intercept indicator

    postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept) {
    ## code to evaluate in estimate.gam, to find the saturated
    ## loglik by Newton method
    ## searching for the mu (on logit scale) that max loglik given theta...
      # wts <- object$prior.weights
      theta <- family$getTheta(trans=TRUE) ## exp theta
      lf <- family$saturated.ll(y, prior.weights,theta)
      ## storing the saturated loglik for each datum...
      ##object$family$data <- list(ls = lf$term,mu.ls = lf$mu)   
      l2 <- family$dev.resids(y,fitted,prior.weights)
      posr <- list()
      posr$deviance <- 2*lf$f + sum(l2)
      wtdmu <- if (intercept) sum(prior.weights * y)/sum(prior.weights) 
              else family$linkinv(offset)
      posr$null.deviance <- 2*lf$f + sum(family$dev.resids(y, wtdmu, prior.weights))
      posr$family <- 
      paste("Beta regression(",round(theta,3),")",sep="")
      posr
    } ## postproc betar

    initialize <- expression({
        ##n <- rep(1, nobs)
        mustart <- y 
    })

    residuals <- function(object,type=c("deviance","working","response","pearson")) {
      if (type == "working") { 
        res <- object$residuals 
      } else if (type == "response") {
        res <- object$y - object$fitted.values
      } else if (type == "deviance") { 
        y <- object$y
        mu <- object$fitted.values
        wts <- object$prior.weights
        lf <- object$family$saturated.ll(y, wts,object$family$getTheta(TRUE))
        #object$family$data$ls <- lf$term  
        res <- 2*lf$term + object$family$dev.resids(y,mu,wts)
        res[res<0] <- 0
        s <- sign(y-mu)
        res <- sqrt(res) * s   
      } else if (type == "pearson") {
        mu <- object$fitted.values
        res <- (object$y - mu)/object$family$variance(mu)^.5
      }
      res
     } ## residuals betar

    rd <- function(mu,wt,scale) {
     ## simulate data given fitted latent variable in mu 
      Theta <- exp(get(".Theta"))
      r <- rbeta(n=length(mu),shape1=Theta*mu,shape2=Theta*(1-mu))
      eps <- get(".betarEps")
      r[r>=1-eps] <- 1 - eps
      r[r<eps] <- eps
      r
    }

    qf <- function(p,mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      q <- qbeta(p,shape1=Theta*mu,shape2=Theta*(1-mu))
      eps <-  get(".betarEps")
      q[q>=1-eps] <- 1 - eps
      q[q<eps] <- eps
      q
    } ## rs betar

    environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    environment(rd)<- environment(qf) <- environment(variance) <- environment(putTheta) <-
    environment(saturated.ll) <- environment(preinitialize) <- env

    structure(list(family = "Beta regression", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd, variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls,
        preinitialize=preinitialize,postproc=postproc, residuals=residuals, saturated.ll=saturated.ll,
        validmu = validmu, valideta = stats$valideta, n.theta=n.theta,  
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,rd=rd,qf=qf), 
        class = c("extended.family","family"))
    
} ## betar


  
## scaled t (Natalya Pya) ...

scat <- function (theta = NULL, link = "identity",min.df = 3) { 
## Extended family object for scaled t distribution
## length(theta)=2; log theta supplied. 
## Written by Natalya Pya.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("identity", "log", "inverse")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for scat family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)

#  else if (is.character(link)) {
#    stats <- make.link(link)
#    linktemp <- link
#  } else {
#    if (inherits(link, "link-glm")) {
#       stats <- link
#            if (!is.null(stats$name))
#                linktemp <- stats$name
#        }
#        else stop(linktemp, " link not available for scaled t distribution; available links are \"identity\", \"log\",  and \"inverse\"")
#    }
    ## Theta <-  NULL;
    n.theta <- 2
    if (!is.null(theta)&&sum(theta==0)==0) {
      if (abs(theta[1])<=min.df) {
        min.df <- 0.9*abs(theta[1])
        warning("Supplied df below min.df. min.df reset")
      }	
      if (sum(theta<0)) { 
        iniTheta <- c(log(abs(theta[1])-min.df),log(abs(theta[2]))) ## initial theta supplied
      } else { ## fixed theta supplied
        iniTheta <- c(log(theta[1]-min.df),log(theta[2])) 
        n.theta <- 0 ## no thetas to estimate
      }
    } else iniTheta <- c(-2,-1) ## inital log theta value
               
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) { 
    ## trans transforms to the original scale...
      th <- get(".Theta");min.df <- get(".min.df")
      if (trans) { th <- exp(th); th[1] <- th[1] + min.df  }
      th
    }
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))
    assign(".min.df", min.df, envir = env)

    variance <- function(mu) { 
        th <- get(".Theta")
	min.df <- get(".min.df")
        nu <- exp(th[1])+min.df; sig <- exp(th[2])
        sig^2*nu/(nu-2)
    }

   validmu <- function(mu) all(is.finite(mu))

   dev.resids <- function(y, mu, wt,theta=NULL) {
      if (is.null(theta)) theta <- get(".Theta")
      min.df <- get(".min.df")
      nu <- exp(theta[1])+min.df; sig <- exp(theta[2])
      wt * (nu + 1)*log1p((1/nu)*((y-mu)/sig)^2)
    }
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the scat deviance...
      ## ltheta <- theta
      min.df <- get(".min.df")
      nu <- exp(theta[1])+min.df; sig <- exp(theta[2])
      nu1 <- nu + 1;  ym <- y - mu; nu2 <- nu - min.df;
      a <- 1 + (ym/sig)^2/nu
      oo <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      nu1ym <- nu1*ym
      sig2a <- sig^2*a
      nusig2a <- nu*sig2a
      f <- nu1ym/nusig2a
      f1 <- ym/nusig2a
      oo$Dmu <- -2 * wt * f
      oo$Dmu2 <- 2 * wt * nu1*(1/nusig2a- 2*f1^2)  # - 2*ym^2/(nu^2*sig^4*a^2)
      term <- 2*nu1/sig^2/(nu+3)
      n <- length(y) 
      oo$EDmu2 <- rep(term,n)
      
      if (level>0) { ## quantities needed for first derivatives
        nu1nusig2a <- nu1/nusig2a
        nu2nu <- nu2/nu
        fym <- f*ym; ff1 <- f*f1; f1ym <- f1*ym; fymf1 <- fym*f1
        ymsig2a <- ym/sig2a
        oo$EDmu2th <- oo$Dmu2th <- oo$Dmuth <- oo$Dth <- matrix(0,n,2)
        oo$Dth[,1] <- 1 * wt * nu2 * (log(a) - fym/nu) 
        oo$Dth[,2] <- -2 * wt * fym    
        oo$Dmuth[,1] <- 2 * wt *(f - ymsig2a - fymf1)*nu2nu
        oo$Dmuth[,2] <- 4* wt* f* (1- f1ym)
        oo$Dmu3 <- 4 * wt * f * (3/nusig2a - 4*f1^2) 
        oo$Dmu2th[,1] <- 2* wt * (-nu1nusig2a + 1/sig2a + 5*ff1- 2*f1ym/sig2a - 4*fymf1*f1)*nu2nu
        oo$Dmu2th[,2] <- 4*wt*(-nu1nusig2a + ff1*5 - 4*ff1*f1ym)
	oo$EDmu3 <- rep(0,n)
	oo$EDmu2th <- cbind(4/(sig^2*(nu+3)^2)*exp(theta[1]),-2*oo$EDmu2)
      } 
      if (level>1) { ## whole lot
        ## nu1nu2 <- nu1*nu2; 
        nu1nu <- nu1/nu
        fymf1ym <- fym*f1ym; f1ymf1 <- f1ym*f1
        oo$Dmu4 <- 12 * wt * (-nu1nusig2a/nusig2a + 8*ff1/nusig2a - 8*ff1 *f1^2) 
        n2d <- 3 # number of the 2nd order derivatives
        oo$Dmu3th <- matrix(0,n,2)
        oo$Dmu2th2 <- oo$Dmuth2 <- oo$Dth2 <- matrix(0,n,n2d)
        oo$Dmu3th[,1] <- 4*wt*(-6*f/nusig2a + 3*f1/sig2a + 18*ff1*f1 - 4*f1ymf1/sig2a - 12*nu1ym*f1^4)*nu2nu
        oo$Dmu3th[,2] <- 48*wt* f* (- 1/nusig2a + 3*f1^2 -  2*f1ymf1*f1)
        
        oo$Dth2[,1] <- 1*wt *(nu2*log(a) +nu2nu*ym^2*(-2*nu2-nu1+ 2*nu1*nu2nu - nu1*nu2nu*f1ym)/nusig2a) ## deriv of D w.r.t. theta1 theta1 
  
        oo$Dth2[,2] <- 2*wt*(fym - ym*ymsig2a - fymf1ym)*nu2nu  ## deriv of D wrt theta1 theta2
        oo$Dth2[,3] <- 4 * wt * fym *(1 - f1ym)  ## deriv of D wrt theta2 theta2

        term <- 2*nu2nu - 2*nu1nu*nu2nu -1 + nu1nu
        oo$Dmuth2[,1] <- 2*wt*f1*nu2*(term - 2*nu2nu*f1ym + 4*fym*nu2nu/nu - fym/nu - 2*fymf1ym*nu2nu/nu)
 
        oo$Dmuth2[,2] <- 4*wt* (-f + ymsig2a + 3*fymf1 - ymsig2a*f1ym - 2*fymf1*f1ym)*nu2nu

        oo$Dmuth2[,3] <- 8*wt* f * (-1 + 3*f1ym - 2*f1ym^2)

        oo$Dmu2th2[,1] <- 2*wt*nu2*(-term + 10*nu2nu*f1ym -
               16*fym*nu2nu/nu - 2*f1ym + 5*nu1nu*f1ym - 8*nu2nu*f1ym^2 + 26*fymf1ym*nu2nu/nu - 
               4*nu1nu*f1ym^2 - 12*nu1nu*nu2nu*f1ym^3)/nusig2a
       
        oo$Dmu2th2[,2] <- 4*wt*(nu1nusig2a - 1/sig2a - 11*nu1*f1^2 + 5*f1ym/sig2a + 22*nu1*f1ymf1*f1 - 
                    4*f1ym^2/sig2a - 12*nu1*f1ymf1^2)*nu2nu
        oo$Dmu2th2[,3] <- 8*wt * (nu1nusig2a - 11*nu1*f1^2 + 22*nu1*f1ymf1*f1 - 12*nu1*f1ymf1^2)

      }
      oo
    }  ## Dd scat

 
    aic <- function(y, mu, theta=NULL, wt, dev) {
        min.df <- get(".min.df")
        if (is.null(theta)) theta <- get(".Theta")
        nu <- exp(theta[1])+min.df; sig <- exp(theta[2])
        term <- -lgamma((nu+1)/2)+ lgamma(nu/2) + log(sig*(pi*nu)^.5) +
           (nu+1)*log1p(((y-mu)/sig)^2/nu)/2  ## `-'log likelihood for each observation
        2 * sum(term * wt)
    } ## aic scat
    
    ls <- function(y,w,theta,scale) {
       ## the log saturated likelihood function for scat
       ## (Note these are correct but do not correspond to NP notes)
       if (length(w)==1) w <- rep(w,length(y))
       #vec <- !is.null(attr(theta,"vec.grad"))
       min.df <- get(".min.df")
       nu <- exp(theta[1])+min.df; sig <- exp(theta[2]); nu2 <- nu-min.df;
       nu2nu <- nu2/nu; nu12 <- (nu+1)/2
       term <- lgamma(nu12) - lgamma(nu/2) - log(sig*(pi*nu)^.5)
       ls <- sum(term*w) 
       ## first derivative wrt theta...
       lsth2 <- matrix(0,2,2)  ## rep(0, 3)
       term <- nu2 * digamma(nu12)/2- nu2 * digamma(nu/2)/2 - 0.5*nu2nu
       LSTH <- cbind(w*term,-1*w)
       lsth <- colSums(LSTH)
       ## second deriv...      
       term <-  nu2^2 * trigamma(nu12)/4 + nu2 * digamma(nu12)/2 -
           nu2^2 * trigamma(nu/2)/4 - nu2 * digamma(nu/2)/2 + 0.5*(nu2nu)^2 - 0.5*nu2nu
       lsth2[1,1] <- sum(term*w)
       lsth2[1,2] <- lsth2[2,1] <- lsth2[2,2] <- 0
       list(ls=ls,## saturated log likelihood
            lsth1=lsth, ## first derivative vector wrt theta
	    LSTH1=LSTH,
            lsth2=lsth2) ## Hessian wrt theta
    } ## ls scat

    preinitialize <- function(y,family) {
      ## initialize theta from raw observations..
       if (family$n.theta>0) {
         ## low df and low variance promotes indefiniteness.
	 ## Better to start with moderate df and fairly high
	 ## variance...
         Theta <- c(1.5, log(0.8*sd(y))) 
         return(list(Theta=Theta))
       } ## otherwise fixed theta supplied
    } ## preinitialize scat

    initialize <- expression({
        if (any(is.na(y))) stop("NA values not allowed for the scaled t family")
        ##n <- rep(1, nobs)
        mustart <- y + (y == 0)*.1
    })
    
    postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept)  {
      posr <- list()
      posr$null.deviance <- find.null.dev(family,y,eta=linear.predictors,offset,prior.weights)
      th <- round(family$getTheta(TRUE),3)
      if (th[1]>999) th[1] <- Inf
      posr$family <- paste("Scaled t(",paste(th,collapse=","),")",sep="")
      posr
    } ## postproc scat
    
    rd <- function(mu,wt,scale) {
     ## simulate data given fitted latent variable in mu 
      theta <- get(".Theta");min.df <- get(".min.df")
      nu <- exp(theta[1])+min.df; sig <- exp(theta[2])
      n <- length(mu)
      stats::rt(n=n,df=nu)*sig + mu
    } ## rd scat

    environment(dev.resids) <- environment(aic) <- environment(getTheta) <- environment(Dd) <-
    environment(ls) <- environment(rd)<- environment(variance) <- environment(putTheta) <- env

    structure(list(family = "scaled t", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,postproc=postproc,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls, preinitialize=preinitialize,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta,   
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta, rd=rd),
        class = c("extended.family","family"))
} ## scat



## zero inflated Poisson (Simon Wood)...

lind <- function(l,th,deriv=0,k=0) {
## evaluate th[1] + exp(th[2])*l and some derivs
  th[2] <- exp(th[2])
  r <- list(p = th[1] + (k+th[2])*l)
  r$p.l <- k + th[2]   ## p_l
  r$p.ll <- 0 ## p_ll
  if (deriv) {
    n <- length(l);  
    r$p.lllth <- r$p.llth <- r$p.lth <- r$p.th <- matrix(0,n,2)
    r$p.th[,1] <- 1   ## dp/dth1 
    r$p.th[,2] <- th[2]*l ## dp/dth2 
    r$p.lth[,2] <- th[2] ## p_lth2
    r$p.llll <- r$p.lll <- 0   ## p_lll,p_llll
    r$p.llth2 <- r$p.lth2 <- r$p.th2 <- matrix(0,n,3) ## ordered l_th1th1,l_th1th2,l_th2th2
    r$p.th2[,3] <- l*th[2] ## p_th2th2
    r$p.lth2[,3] <- th[2]  ## p_lth2th2
  }
  r
} ## lind

logid <- function(l,th,deriv=0,a=0,trans=TRUE) {
## evaluate exp(th[1]+th[2]*l)/(1+exp(th[1]+th[2]*l))
## and some of its derivatives
## if trans==TRUE then it is assumed that the 
## transformation th[2] = exp(th[2]) is applied on input
  b <- 1-2*a
  ## x is dth[2]/dth[2]' where th[2]' is input version, xx is second deriv over first
  if (trans) { xx <- 1; x <- th[2] <- exp(th[2])} else { x <- 1;xx <- 0}
  p <- f <- th[1] + th[2] * l
  ind <- f > 0; ef <- exp(f[!ind])
  p[!ind] <- ef/(1+ef); p[ind] <- 1/(1+exp(-f[ind]))  
  r <- list(p = a + b * p)
  
  a1 <- p*(1-p); a2 <- p*(p*(2*p-3)+1)
  r$p.l <- b * th[2]*a1;   ## p_l
  r$p.ll <- b * th[2]^2*a2 ## p_ll
  if (deriv>0) { 
    n <- length(l); r$p.lth <- r$p.th <- matrix(0,n,2) 
    r$p.th[,1] <- b * a1   ## dp/dth1
    r$p.th[,2] <- b * l*a1 * x ## dp/dth2
    r$p.lth[,1] <- b * th[2]*a2  ## p_lth1
    r$p.lth[,2] <- b * (l*th[2]*a2 + a1) * x ## p_lth2

    a3 <- p*(p*(p*(-6*p + 12) -7)+1)
    r$p.lll <- b * th[2]^3*a3   ## p_lll
    r$p.llth <- matrix(0,n,2) 
    r$p.llth[,1] <- b * th[2]^2 * a3 ## p_llth1
    r$p.llth[,2] <- b * (l*th[2]^2*a3 + 2*th[2]*a2) * x ## p_ppth2

    a4 <- p*(p*(p*(p*(p*24-60)+50)-15)+1)
    r$p.llll <- b * th[2]^4*a4 ## p_llll
    r$p.lllth <- matrix(0,n,2)
    r$p.lllth[,1] <- b * th[2]^3*a4 ## p_lllth1
    r$p.lllth[,2] <- b * (th[2]^3*l*a4 + 3*th[2]^2*a3) * x ## p_lllth2

    r$p.llth2 <- r$p.lth2 <- r$p.th2 <- matrix(0,n,3) ## ordered l_th1th1,l_th1th2,l_th2th2

    r$p.th2[,1] <- b * a2   ## p_th1th1
    r$p.th2[,2] <- b * l*a2 * x ## p_th1th2  
    r$p.th2[,3] <- b * l*l*a2 * x * x + xx* r$p.th[,2] ## p_th2th2

    r$p.lth2[,1] <- b * th[2]*a3  ## p_lth1th1
    r$p.lth2[,2] <- b * (th[2]*l*a3 + a2) * x ## p_lth1th2  
    r$p.lth2[,3] <- b * (l*l*a3*th[2] + 2*l*a2) *x * x + xx*r$p.lth[,2]  ## p_lth2th2

    r$p.llth2[,1] <- b * th[2]^2*a4  ## p_llth1th1
    r$p.llth2[,2] <- b * (th[2]^2*l*a4 + 2*th[2]*a3) *x ## p_llth1th2  
    r$p.llth2[,3] <- b * (l*l*th[2]^2*a4 + 4*l*th[2]*a3 + 2*a2) *x*x + xx*r$p.llth[,2] ## p_llth2th2
  }
  r
} ## logid



ziP <- function (theta = NULL, link = "identity",b=0) { 
## zero inflated Poisson parameterized in terms of the log Poisson parameter, gamma. 
## eta = theta[1] + exp(theta[2])*gamma), and 1-p = exp(-exp(eta)) where p is 
## probability of presence.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("identity")
  if (linktemp %in% okLinks) stats <- make.link(linktemp) else 
  stop(gettextf("link \"%s\" not available for ziP family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")),domain = NA)
	
  ## Theta <-  NULL;
  n.theta <- 2
  if (!is.null(theta)) {
      ## fixed theta supplied
      iniTheta <-  c(theta[1],theta[2])
      n.theta <- 0 ## no thetas to estimate
  } else iniTheta <- c(0,0) ## inital theta value - start at Poisson
  
  env <- new.env(parent = environment(ziP))# new.env(parent = .GlobalEnv) 
  
  if (b<0) b <- 0; assign(".b", b, envir = env)
  assign(".Theta", iniTheta, envir = env)
  getTheta <- function(trans=FALSE) { 
  ## trans transforms to the original scale...
    th <- get(".Theta")
    if (trans) {
      th[2] <- get(".b") + exp(th[2])
    }
    th
  }

  putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))
  
  validmu <- function(mu) all(is.finite(mu))
  
  dev.resids <- function(y, mu, wt,theta=NULL) {
    ## this version ignores saturated likelihood
    if (is.null(theta)) theta <- get(".Theta")
    b <- get(".b")
    p <- theta[1] + (b + exp(theta[2])) * mu ## l.p. for prob present
    -2*zipll(y,mu,p,deriv=0)$l
  }
  
  Dd <- function(y, mu, theta, wt=NULL, level=0) {
    ## here mu is lin pred for Poisson mean so E(y) = exp(mu)
    ## Deviance for log lik of zero inflated Poisson. 
    ## code here is far more general than is needed - could deal 
    ## with any 2 parameter mapping of lp of mean to lp of prob presence.
    if (is.null(theta)) theta <- get(".Theta")
    deriv <- 1; if (level==1) deriv <- 2 else if (level>1) deriv <- 4 
    b <- get(".b")
    g <- lind(mu,theta,level,b) ## the derviatives of the transform mapping mu to p
    z <- zipll(y,mu,g$p,deriv)
    oo <- list();n <- length(y)
    if (is.null(wt)) wt <- rep(1,n)
    oo$Dmu <- -2*wt*(z$l1[,1] + z$l1[,2]*g$p.l)
    oo$Dmu2 <- -2*wt*(z$l2[,1] + 2*z$l2[,2]*g$p.l + z$l2[,3]*g$p.l^2 + z$l1[,2]*g$p.ll)
    ## WARNING: following requires z$El1 term to be added if transform modified so 
    ##          that g$p.ll != 0....
    oo$EDmu2 <- -2*wt*(z$El2[,1] + 2*z$El2[,2]*g$p.l + z$El2[,3]*g$p.l^2)

    if (level>0) { ## l,p - ll,lp,pp -  lll,llp,lpp,ppp - llll,lllp,llpp,lppp,pppp
      oo$Dth <- -2*wt*z$l1[,2]*g$p.th ## l_p p_th
      oo$Dmuth <- -2*wt*(z$l2[,2]*g$p.th + z$l2[,3]*g$p.l*g$p.th + z$l1[,2]*g$p.lth) 
      oo$Dmu2th <- -2*wt*(z$l3[,2]*g$p.th + 2*z$l3[,3]*g$p.l*g$p.th + 2* z$l2[,2]*g$p.lth + 
       z$l3[,4]*g$p.l^2*g$p.th + z$l2[,3]*(2*g$p.l*g$p.lth + g$p.th*g$p.ll) + z$l1[,2]*g$p.llth)
      oo$Dmu3 <- -2*wt*(z$l3[,1] + 3*z$l3[,2]*g$p.l + 3*z$l3[,3]*g$p.l^2 + 3*z$l2[,2]*g$p.ll +
       z$l3[,4]*g$p.l^3 +3*z$l2[,3]*g$p.l*g$p.ll + z$l1[,2]*g$p.lll)
    } 
    if (level>1) {
      p.thth <- matrix(0,n,3);p.thth[,1] <- g$p.th[,1]^2
      p.thth[,2] <- g$p.th[,1]*g$p.th[,2];p.thth[,3] <- g$p.th[,2]^2
      oo$Dth2 <- -2*wt*(z$l2[,3]*p.thth + z$l1[,2]*g$p.th2)
      p.lthth <- matrix(0,n,3);p.lthth[,1] <- g$p.th[,1]*g$p.lth[,1]*2
      p.lthth[,2] <- g$p.th[,1]*g$p.lth[,2] + g$p.th[,2]*g$p.lth[,1];
      p.lthth[,3] <- g$p.th[,2]*g$p.lth[,2]*2
      oo$Dmuth2 <- -2*wt*( z$l3[,3]*p.thth + z$l2[,2]*g$p.th2 + z$l3[,4]*g$p.l*p.thth + 
        z$l2[,3]*(g$p.th2*g$p.l + p.lthth) + z$l1[,2]*g$p.lth2)
      p.lthlth <- matrix(0,n,3);p.lthlth[,1] <- g$p.lth[,1]*g$p.lth[,1]*2
      p.lthlth[,2] <- g$p.lth[,1]*g$p.lth[,2] + g$p.lth[,2]*g$p.lth[,1];
      p.lthlth[,3] <- g$p.lth[,2]*g$p.lth[,2]*2
      p.llthth <- matrix(0,n,3);p.llthth[,1] <- g$p.th[,1]*g$p.llth[,1]*2
      p.llthth[,2] <- g$p.th[,1]*g$p.llth[,2] + g$p.th[,2]*g$p.llth[,1];
      p.llthth[,3] <- g$p.th[,2]*g$p.llth[,2]*2

      oo$Dmu2th2 <- -2*wt*(z$l4[,3]*p.thth + z$l3[,2]*g$p.th2 + 2*z$l4[,4] * p.thth *g$p.l + 2*z$l3[,3]*(g$p.th2*g$p.l + p.lthth) +
        2*z$l2[,2]*g$p.lth2 + z$l4[,5]*p.thth*g$p.l^2 + z$l3[,4]*(g$p.th2*g$p.l^2 + 2*p.lthth*g$p.l + p.thth*g$p.ll) +
        z$l2[,3]*(p.lthlth + 2*g$p.l*g$p.lth2 + p.llthth + g$p.th2*g$p.ll) + z$l1[,2]*g$p.llth2)

      oo$Dmu3th <- -2*wt*(z$l4[,2]*g$p.th + 3*z$l4[,3]*g$p.th*g$p.l + 3*z$l3[,2]*g$p.lth + 2*z$l4[,4]*g$p.th*g$p.l^2 + 
        z$l3[,3]*(6*g$p.lth*g$p.l + 3*g$p.th*g$p.ll) + 3*z$l2[,2]*g$p.llth + z$l4[,4]*g$p.th*g$p.l^2 + 
        z$l4[,5]*g$p.th*g$p.l^3 + 3*z$l3[,4]*(g$p.l^2*g$p.lth + g$p.th*g$p.l*g$p.ll) +
        z$l2[,3]*(3*g$p.lth*g$p.ll + 3*g$p.l*g$p.llth + g$p.th*g$p.lll) + z$l1[,2]*g$p.lllth)

      oo$Dmu4 <- -2*wt*(z$l4[,1] + 4*z$l4[,2]*g$p.l + 6*z$l4[,3]*g$p.l^2 + 6*z$l3[,2]*g$p.ll + 
        4*z$l4[,4]*g$p.l^3 + 12*z$l3[,3]*g$p.l*g$p.ll + 4*z$l2[,2]*g$p.lll + z$l4[,5] * g$p.l^4 +
        6*z$l3[,4]*g$p.l^2*g$p.ll + z$l2[,3] *(4*g$p.l*g$p.lll + 3*g$p.ll^2) + z$l1[,2]*g$p.llll)

    }
    oo
  } ## Dd ziP
  
  aic <- function(y, mu, theta=NULL, wt, dev) {
    if (is.null(theta)) theta <- get(".Theta")
    b <- get(".b")
    p <- theta[1] + (b+ exp(theta[2])) * mu ## l.p. for prob present
    sum(-2*wt*zipll(y,mu,p,0)$l)
  }

  ls <- function(y,w,theta,scale) {
       ## the log saturated likelihood function for ziP
       ## ls is defined as zero for REML/ML expression as deviance is defined as -2*log.lik
       #vec <- !is.null(attr(theta,"vec.grad"))
       #lsth1 <- if (vec) matrix(0,length(y),2) else c(0,0)
       list(ls=0,## saturated log likelihood
            lsth1=c(0,0),  ## first deriv vector w.r.t theta - last element relates to scale
            LSTH1=matrix(0,length(y),2),
            lsth2=matrix(0,2,2)) ##Hessian w.r.t. theta
  } ## ls ZiP


    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the zero inflated Poisson family")
        if (all.equal(y,round(y))!=TRUE) {
          stop("Non-integer response variables are not allowed with ziP ")
        }
        if ((min(y)==0&&max(y)==1)) stop("Using ziP for binary data makes no sense")
        ##n <- rep(1, nobs)
        mustart <- log(y + (y==0)/5) 
    })

    ## compute family label, deviance, null.deviance...
    ## requires prior weights, y, family, linear predictors
    postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept) {
      posr <- list()
      posr$family <- 
      paste("Zero inflated Poisson(",paste(round(family$getTheta(TRUE),3),collapse=","),")",sep="")
      ## need to fix deviance here!!
      ## wts <- object$prior.weights
      lf <- family$saturated.ll(y,family,prior.weights)
      ## storing the saturated loglik for each datum...
      ##object$family$data <- list(ls = lf)   
      l2 <- family$dev.resids(y,linear.predictors,prior.weights)
      posr$deviance <- sum(l2-lf)
      fnull <- function(gamma,family,y,wt) {
        ## evaluate deviance for single parameter model
        sum(family$dev.resids(y, rep(gamma,length(y)), wt))
      }
      meany <- mean(y)
      posr$null.deviance <-
      optimize(fnull,interval=c(meany/5,meany*3),family=family,y=y,wt = prior.weights)$objective - sum(lf)
 
      ## object$weights <- pmax(0,object$working.weights) ## Fisher can be too extreme
      ## E(y) = p * E(y) - but really can't mess with fitted.values if e.g. rd is to work.
      posr
    } ## postproc ZiP

    rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu 
      rzip <- function(gamma,theta) { ## generate ziP deviates according to model and lp gamma
        y <- gamma; n <- length(y)
        lambda <- exp(gamma)
        mlam <- max(c(lambda[is.finite(lambda)],.Machine$double.eps^.2))
        lambda[!is.finite(lambda)] <- mlam
        b <- get(".b")
        eta <- theta[1] + (b+exp(theta[2]))*gamma
        p <- 1- exp(-exp(eta))
        ind <- p > runif(n)
        y[!ind] <- 0
        #np <- sum(ind)
        ## generate from zero truncated Poisson, given presence...
        lami <- lambda[ind]
        yi <- p0 <- dpois(0,lami)
        nearly1 <- 1 - .Machine$double.eps*10
        ii <- p0 > nearly1 
        yi[ii] <- 1 ## lambda so low that almost certainly y=1
        yi[!ii] <- qpois(runif(sum(!ii),p0[!ii],nearly1),lami[!ii])
        y[ind] <- yi 
        y
      } 
      rzip(mu,get(".Theta"))
    } ## rd ZiP
   
   saturated.ll <- function(y,family,wt=rep(1,length(y))) {
      ## function to get saturated ll for ziP - 
      ## actually computes -2 sat ll.
      pind <- y>0 ## only these are interesting
      wt <- wt[pind]
      y <- y[pind];
      mu <- log(y)  
      keep.on <- TRUE
      theta <- family$getTheta() 
      r <- family$Dd(y,mu,theta,wt)
      l <- family$dev.resids(y,mu,wt,theta)
      lmax <- max(abs(l))
      ucov <- abs(r$Dmu) > lmax*1e-7
      k <- 0
      while (keep.on) { 
        step <- -r$Dmu/r$Dmu2
        step[!ucov] <- 0
        mu1 <- mu + step
        l1 <- family$dev.resids(y,mu1,wt,theta)
        ind <- l1>l & ucov
        kk <- 0
        while (sum(ind)>0&&kk<50) {
            step[ind] <- step[ind]/2
            mu1 <- mu + step
            l1 <- family$dev.resids(y,mu1,wt,theta) 
            ind <- l1>l & ucov
            kk <- kk + 1
        }       
        mu <- mu1;l <- l1
        r <- family$Dd(y,mu,theta,wt)
        ucov <- abs(r$Dmu) > lmax*1e-7
        k <- k + 1
        if (all(!ucov)||k==100) keep.on <- FALSE
      }
      l1 <- rep(0,length(pind));l1[pind] <- l
      l1
    } ## saturated.ll ZiP

  residuals <- function(object,type=c("deviance","working","response")) {
    if (type == "working") { 
      res <- object$residuals 
    } else if (type == "response") {
      res <- object$y - predict.gam(object,type="response")
    } else if (type == "deviance") { 
      y <- object$y
      mu <- object$linear.predictors
      wts <- object$prior.weights
      res <- object$family$dev.resids(y,mu,wts)
      ## next line is correct as function returns -2*saturated.log.lik
      res <- res - object$family$saturated.ll(y,object$family,wts)
      fv <- predict.gam(object,type="response")
      s <- attr(res,"sign")
      if (is.null(s)) s <- sign(y-fv)
      res <- as.numeric(sqrt(pmax(res,0)) * s) 
    }
    res
  } ## residuals (ziP)

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
 
    theta <- family$getTheta()

    if (is.null(eta)) { ## return probabilities
      discrete <- is.list(X)
      ## linear predictor for poisson parameter... 
      gamma <- off + if (discrete) Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop) else drop(X%*%beta)  
      if (se) {
        se <- if (discrete) sqrt(pmax(0,diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1))) else
	                    sqrt(pmax(0,rowSums((X%*%Vb)*X))) ## se of lin pred
      } else se <- NULL
      #gamma <- drop(X%*%beta + off) ## linear predictor for poisson parameter 
      #se <- if (se) drop(sqrt(pmax(0,rowSums((X%*%Vb)*X)))) else NULL ## se of lin pred
    } else { se <- NULL; gamma <- eta}
    ## now compute linear predictor for probability of presence...
    b <- get(".b")
    eta <- theta[1] + (b+exp(theta[2]))*gamma
    et <- exp(eta)
    mu <- p <- 1 - exp(-et)
    fv <- lambda <- exp(gamma)  
    ind <- gamma < log(.Machine$double.eps)/2
    mu[!ind] <- lambda[!ind]/(1-exp(-lambda[!ind]))
    mu[ind] <- 1
    fv <- list(p*mu)    ## E(y)    
    if (is.null(se)) return(fv) else {
      dp.dg <- p  
      # ind <- eta < log(.Machine$double.xmax)/2
      # dp.dg[!ind] <- 0
      dp.dg <- exp(-et)*et*(b + exp(theta[2]))
      dmu.dg <- (lambda + 1)*mu - mu^2
      fv[[2]] <- abs(dp.dg*mu+dmu.dg*p)*se   
      names(fv) <- c("fit","se.fit")
      return(fv)
    }
  } ## predict ZiP


   
  environment(saturated.ll) <- environment(dev.resids) <- environment(Dd) <-
  environment(aic) <- environment(getTheta) <- environment(rd) <- environment(predict) <-
  environment(putTheta) <- env

  structure(list(family = "zero inflated Poisson", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd, rd=rd,residuals=residuals,
        aic = aic, mu.eta = stats$mu.eta, g2g = stats$g2g,g3g=stats$g3g, g4g=stats$g4g, 
        #preinitialize=preinitialize,
        initialize = initialize,postproc=postproc,ls=ls,no.r.sq=TRUE,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta,predict=predict,
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,saturated.ll = saturated.ll),
        class = c("extended.family","family"))

} ## ziP




##  R routines for the package mgcv (c) Simon Wood 2000-2008
##  With contributions from Henric Nilsson



pcls <- function(M)
# Function to perform penalized constrained least squares.
# Problem to be solved is:
#
#  minimise      ||W^0.5 (y - Xp)||^2 + p'Bp
#  subject to    Ain p >= b  & C p = "constant"
#
# where B = \sum_{i=1}^m \theta_i S_i and W=diag(w)
# on entry this routine requires a list M, with the following elements:
# M$X - the design matrix for this problem.
# M$p - a feasible initial parameter vector - note that this should not be set up to
#       lie exactly on all the inequality constraints - which can easily happen if M$p=0!
# M$y - response variable
# M$w - weight vector: W= diag(M$w)
# M$Ain - matrix of inequality constraints
# M$bin - b above
# M$C  - fixed constraint matrix
# M$S  - List of (minimal) penalty matrices
# M$off - used for unpacking M$S
# M$sp - array of theta_i's 
# Ain, bin and p are not in the object needed to call mgcv....
#
{ nar<-c(length(M$y),length(M$p),dim(M$Ain)[1],dim(M$C)[1],0)
  H<-0
  ## sanity checking ...
  if (nrow(M$X)!=nar[1]) stop("nrow(M$X) != length(M$y)") 
  if (ncol(M$X)!=nar[2]) stop("ncol(M$X) != length(M$p)")
  if (length(M$w)!=nar[1]) stop("length(M$w) != length(M$y)")
  if (nar[3]!=length(M$bin)) stop("nrow(M$Ain) != length(M$bin)")
  if (nrow(M$Ain)>0)
  { if (ncol(M$Ain)!=nar[2]) stop("nrow(M$Ain) != length(M$p)") 
    res <- as.numeric(M$Ain%*%M$p) - as.numeric(M$bin)
    res <- mean(abs(res))
    if (res<.Machine$double.eps^.5) 
    warning("initial parameters very close to inequality constraints")
  }
  
  if (nrow(M$C)>0) if (ncol(M$C)!=nar[2]) stop("ncol(M$C) != length(M$p)")  
  if (length(M$S)!=length(M$off)) stop("M$S and M$off have different lengths")
  if (length(M$S)!=length(M$sp)) stop("M$sp has different length to M$S and M$off")
  
  
  # pack the S array for mgcv call 
  m<-length(M$S)
  Sa<-array(0,0);df<-0
  if (m>0) for (i in 1:m)
  { Sa<-c(Sa,M$S[[i]])
    df[i]<-nrow(M$S[[i]])
    if (M$off[i]+df[i]-1>nar[2]) stop(paste("M$S[",i,"] is too large given M$off[",i,"]",sep=""))
  }
  qra.exist <- FALSE
  if (ncol(M$X)>nrow(M$X)) {
    if (m>0) stop("Penalized model matrix must have no more columns than rows") 
    else { ## absorb M$C constraints
      qra <- qr(t(M$C))
      j <- nrow(M$C);k <- ncol(M$X)
      M$X <- t(qr.qty(qra,t(M$X))[(j+1):k,])
      M$Ain <- t(qr.qty(qra,t(M$Ain))[(j+1):k,])
      M$C <- matrix(0,0,0)
      M$p <- rep(0,ncol(M$X)) 
      nar[2] <- length(M$p)
      nar[4] <- 0
      qra.exist <- TRUE
      if  (ncol(M$X)>nrow(M$X)) stop("Model matrix not full column rank")
    }
  }
  o<-.C(C_RPCLS,as.double(M$X),as.double(M$p),as.double(M$y),as.double(M$w),as.double(M$Ain),as.double(M$bin)
        ,as.double(M$C),as.double(H),as.double(Sa),as.integer(M$off),as.integer(df),as.double(M$sp),
        as.integer(length(M$off)),as.integer(nar))
  p <- array(o[[2]],length(M$p))
  if (qra.exist) p <- qr.qy(qra,c(rep(0,j),p))
  p
}  

mgcv.control<-function(conv.tol=1e-7,max.half=20,target.edf=NULL,min.edf=-1)
# control constants for mgcv
# conv.tol - convergence tolerence for multiple s.p. GCV
# max.half - maximum number of step-length halvings to perform at each newton update
#              of s.p.'s
# min.edf - minimum possible estimated degrees of freedom for model - useful for setting limits
#             on overall smoothing parameter. Set to zero or negative to ignore.
# target.edf - set to negative to ignore. This should only be used if cautious optimization
#                is to be used in mgcv searching. If this is non-negative then the local 
#                minimum closest to the target edf will be returned (which can be the global
#                optimum). Designed for use with non-convergent gams.

{ list(conv.tol=conv.tol,max.half=max.half,target.edf=target.edf,min.edf=min.edf)
}

mgcv <- function(y,X,sp,S,off,C=NULL,w=rep(1,length(y)),H=NULL,scale=1,gcv=TRUE,control=mgcv.control())

# Performs multiple smoothing parameter selection for Generalized ridge regression problems 
# y is the response vector
# X is the model matrix
# sp is an array of smoothing parameters. If control$fixed.sp==TRUE then these are taken as 
#    being the s.p.s. Otherwise any positive are taken as initial estimates and negatives 
#    indicate that autoinitialization should take place.
# S is a list of penalty matrices. S[[i]] is the ith penalty matrix, and is stored as the 
#   smallest sub-matrix of the penalty not excluding any non-zero elements. off[i]
#   indicates which parameter S[[i]][1,1] relates to.
# off off[i] is the first parameter to which S[[i]] applies  
# C is an optional linear constaint matrix.
# w is the weight matrix (often 1/std.dev(y), or an estimate of this).
# H a single fixed penalty matrix to be used instead of the multiple penalties in S
# scale is the scale parameter needed by UBRE
# gcv is true if GCV is to be used (and hence scale ignored), FALSE for UBRE.
# control is a list of control options (see mgcv.control()).
#
# An object is returned 
#
# b estimated parameters
# scale estimated or supplied scale parameter
# score minimising score
# sp estimated smoothing parameters 
# Vb estimated covariance matrix
# hat diagonal of hat matrix
# edf array of edf per parameter
# info - a list of convergence diagnostics
#          g - gradients of gcv/ubre score at termination, h - leading diagonal of Hessian
#          e - eigenvalues of Hessian, iter - iterations taken, init.ok - TRUE if second 
#          autoinitialization guess ok (or intial values supplied), step.fail - TRUE
#          if algorithm terminated on step failure rather than convergence. 
#          edf - array of model edf's from final grid search for overall s.p.
#          score - array of gcv/ubre scores corresponding to edf.
#  
{ if (gcv) scale <- -1
  
  if (!is.null(C)) C.r<-nrow(C)          # number of equality constraints
  else {C.r<-0;C<-0}
  q<-ncol(X)            # number of parameters
  n<-nrow(X)            # number of data
  if (is.null(H))  
  { # pack the S array for mgcv call 
    m<-length(S)
    Sa<-array(0,0);df<-0
    if (m>0) for (i in 1:m)
    { Sa<-c(Sa,S[[i]])
      df[i]<-nrow(S[[i]])
    }
    fixed.sp<-0
  }
  else
  { if (length(S)>0) stop("Can't mix fixed and estimated penalties in mgcv() - use magic()")
    k<-1;while (sum(H[,k]!=0)==0&&sum(H[k,]!=0)==0) k<-k+1
    off <- k
    k<-nrow(H);while (sum(H[,k]!=0)==0&&sum(H[k,]!=0)==0) k<-k-1
    df <- k-off+1
    Sa<-array(H[off:(off+df-1),off:(off+df-1)],df*df)
    m<-1
    fixed.sp<-1
  } 

  # deal with quantities that will be estimated
  p<-matrix(0,q,1)      # set up parameter vector
  Vp<-matrix(0.0,q,q)   # estimated covariance matrix
  edf<-array(0,q)       # estimated degrees of freedom
  hat<-array(0,n)       # elements on leading diagonal of hat matrix
  ddiag<-array(0,3*m)   # array for diagonostics
  idiag<-array(0,3)     # array for diagnostics
  Vp[1,1]<-1.0
  gcv.ubre<-1.0;
  direct.mesh<-100      # number of points for overall s.p. initial direct search
  sdiag<-array(0.0,2*direct.mesh) # array for gcv/ubre vs edf diagnostics
  if (is.null(control$target.edf)) control$target.edf<- -1 # set to signal no target edf

  oo<-.C(C_mgcv,as.double(y),as.double(X),as.double(C),as.double(w^2),as.double(Sa),
         as.double(p),as.double(sp),as.integer(off-1),as.integer(df),as.integer(m),
         as.integer(n),as.integer(q),as.integer(C.r),as.double(scale),as.double(Vp),
		 as.double(edf),as.double(control$conv.tol),as.integer(control$max.half),as.double(ddiag),
                 as.integer(idiag),as.double(sdiag),as.integer(direct.mesh),as.double(control$min.edf),
                 as.double(gcv.ubre),as.double(control$target.edf),as.integer(fixed.sp),
                 as.double(hat))
   
  p<-matrix(oo[[6]],q,1);
  scale<-oo[[14]]
  Vp<-matrix(oo[[15]],q,q)
  sp<-matrix(oo[[7]])
  edf<-oo[[16]]
  ddiag<-oo[[19]]
  idiag<-oo[[20]]
  sdiag<-oo[[21]]
  gcv.ubre<-oo[[24]]
  hat<-oo[[27]]
  conv<-list(edf=sdiag[1:direct.mesh],score=sdiag[direct.mesh+1:direct.mesh],g=ddiag[1:m],h=ddiag[(m+1):(2*m)],
             e=ddiag[(2*m+1):(3*m)],iter=idiag[1],init.ok=as.logical(idiag[2]),step.fail=as.logical(idiag[3]))
  
  list(b=p,scale=scale,score=gcv.ubre,sp=sp,Vb=Vp,hat=hat,edf=edf,info=conv)
 
}




interpret.gam <- function (gf)
# interprets a gam formula of the generic form:
#   y~x0+x1+x3*x4 + s(x5)+ s(x6,x7) ....
# and returns:
# 1. a model formula for the parametric part: pf (and pfok indicating whether it has terms)
# 2. a list of descriptors for the smooths: smooth.spec
{ p.env<-environment(gf) # environment of formula
  tf<-terms.formula(gf,specials=c("s","te")) # specials attribute indicates which terms are smooth
 
  terms<-attr(tf,"term.labels") # labels of the model terms 
  nt<-length(terms) # how many terms?
  
  if (attr(tf,"response")>0)  # start the replacement formulae
  { response<-as.character(attr(tf,"variables")[2])
    pf<-rf<-paste(response,"~",sep="")
  }
  else pf<-rf<-"~"
  sp<-attr(tf,"specials")$s     # array of indices of smooth terms 
  tp<-attr(tf,"specials")$te    # indices of tensor product terms
  off<-attr(tf,"offset") # location of offset in formula
  ## have to translate sp,tp so that they relate to terms,
  ## rather than elements of the formula (26/11/05)
  vtab <- attr(tf,"factors") # cross tabulation of vars to terms
  if (length(sp)>0) for (i in 1:length(sp)) {
    ind <- (1:nt)[as.logical(vtab[sp[i],])]
    sp[i] <- ind # the term that smooth relates to
  }
  if (length(tp)>0) for (i in 1:length(tp)) {
    ind <- (1:nt)[as.logical(vtab[tp[i],])]
    tp[i] <- ind # the term that smooth relates to
  } ## re-referencing is complete

  ns<-length(sp)+length(tp) # number of smooths
  k<-kt<-ks<-kp<-1 # counters for terms in the 2 formulae
  len.sp <- length(sp)
  len.tp <- length(tp)

  smooth.spec<-list()
  if (nt)
  for (i in 1:nt) # work through all terms
  { if (k<=ns&&((ks<=len.sp&&sp[ks]==i)||(kt<=len.tp&&tp[kt]==i))) # it's a smooth
    { st<-eval(parse(text=terms[i]),envir=p.env)
      ##if (k>1||kp>1) rf<-paste(rf,"+",st$full.call,sep="") # add to full formula
      ##else rf<-paste(rf,st$full.call,sep="")
      smooth.spec[[k]]<-st
      if (ks<=len.sp&&sp[ks]==i) ks<-ks+1  # counts s() terms
      else kt<-kt+1              # counts te() terms
      k<-k+1     # counts smooth terms 
    } else          # parametric
    { if (kp>1) pf<-paste(pf,"+",terms[i],sep="") # add to parametric formula
      else pf<-paste(pf,terms[i],sep="")
      ##if (k>1||kp>1) rf<-paste(rf,"+",terms[i],sep="") # add to full formula
      ##else rf<-paste(rf,terms[i],sep="")
      kp<-kp+1    # counts parametric terms
    }
  }    
  if (!is.null(off)) # deal with offset
  { if (kp>1) pf<-paste(pf,"+",sep="")
    if (kp>1||k>1) rf<-paste(rf,"+",sep="")
    pf<-paste(pf,as.character(attr(tf,"variables")[1+off]),sep="")
    ##rf<-paste(rf,as.character(attr(tf,"variables")[1+off]),sep="")
    kp<-kp+1          
  }
  if (attr(tf,"intercept")==0) 
  { pf<-paste(pf,"-1",sep="")##;rf<-paste(rf,"-1",sep="")
    if (kp>1) pfok<-1 else pfok<-0
  } else { 
    pfok<-1;if (kp==1) { 
    pf<-paste(pf,"1"); ##if (k==1) rf<-paste(rf,"1",sep="");
    }
  }
  
  fake.formula<-pf
  if (length(smooth.spec)>0) 
  for (i in 1:length(smooth.spec))
  { nt<-length(smooth.spec[[i]]$term)
    ff1<-paste(smooth.spec[[i]]$term[1:nt],collapse="+")
    fake.formula<-paste(fake.formula,"+",ff1)
    if (smooth.spec[[i]]$by!="NA")
    fake.formula<-paste(fake.formula,"+",smooth.spec[[i]]$by)
  }
  fake.formula<-as.formula(fake.formula,p.env)
  ret<-list(pf=as.formula(pf,p.env),pfok=pfok,smooth.spec=smooth.spec,##full.formula=as.formula(rf,p.env),
            fake.formula=fake.formula,response=response)
  class(ret)<-"split.gam.formula"
  ret
}


fixDependence <- function(X1,X2,tol=.Machine$double.eps^.5)
# model matrix X2 may be linearly dependent on X1. This 
# routine finds which columns of X2 should be zeroed to 
# fix this.
{ qr1 <- qr(X1,LAPACK=TRUE)
  R11 <- abs(qr.R(qr1)[1,1])
  r<-ncol(X1);n<-nrow(X1)
  QtX2 <- qr.qty(qr1,X2)[(r+1):n,] # Q'X2
  qr2 <- qr(QtX2,LAPACK=TRUE)
  R <- qr.R(qr2)
  # now final diagonal block of R may be zero, indicating rank 
  # deficiency. 
  r0<-r<-nrow(R)
  while (mean(abs(R[r0:r,r0:r]))< R11*tol) r0 <- r0 -1
  r0<-r0+1
  if (r0>r) return(NULL) else
  qr2$pivot[r0:r] # the columns of X2 to zero in order to get independence
}


gam.side <- function(sm,Xp,tol=.Machine$double.eps^.5)
# works through a list of smooths, sm, aiming to identify nested or partially
# nested terms, and impose identifiability constraints on them.
# Xp is the parametric model matrix. It is needed in order to check whether
# there is a constant (or equivalent) in the model. If there is then this needs 
# to be included when working out side constraints, otherwise dependencies can be 
# missed. 
{ m <- length(sm)
  if (m==0) return(sm)
  v.names<-array("",0);maxDim<-1
  for (i in 1:m) { ## collect all term names and max smooth `dim'
    vn <- sm[[i]]$term
    ## need to include by variables in names
    if (sm[[i]]$by!="NA") vn <- paste(vn,sm[[i]]$by,sep="")
    ## need to distinguish levels of factor by variables...
    if (!is.null(sm[[i]]$by.level))  vn <- paste(vn,sm[[i]]$by.level,sep="")
    sm[[i]]$vn <- vn ## use this record to identify variables from now
    v.names <- c(v.names,vn)
    if (sm[[i]]$dim > maxDim) maxDim <- sm[[i]]$dim
  } 
  lv <- length(v.names)   
  v.names <- unique(v.names)
  if (lv == length(v.names)) return(sm) ## no repeats => no nesting

  ## Only get this far if there is nesting.
  ## Need to test for intercept or equivalent in Xp
  intercept <- FALSE
  if (ncol(Xp)) {
    f <- rep(1,nrow(Xp))
    ff <- qr.fitted(qr(Xp),f)
    if (max(abs(ff-f))<.Machine$double.eps*100) intercept <- TRUE 
  }

  sm.id <- as.list(v.names)
  names(sm.id) <- v.names
  for (i in 1:length(sm.id)) sm.id[[i]]<-array(0,0)
  sm.dim <- sm.id
  for (d in 1:maxDim) {
    for (i in 1:m) {
      if (sm[[i]]$dim==d) for (j in 1:d) { ## work through terms
        term<-sm[[i]]$vn[j]
        a <- sm.id[[term]]
        la <- length(a)+1
        sm.id[[term]][la] <- i   ## record smooth i.d. for this variable
        sm.dim[[term]][la] <- d  ## ... and smooth dim.
      }
    }
  }
  ## so now each unique variable name has an associated array of 
  ## the smooths of which it is an argument, arranged in ascending 
  ## order of dimension.
  if (maxDim==1) stop("model has repeated 1-d smooths of same variable.")
  for (d in 2:maxDim) { ## work up through dimensions 
    for (i in 1:m) { ## work through smooths
      if (sm[[i]]$dim == d) { ## check for nesting
        X1 <- matrix(1,nrow(sm[[i]]$X),as.integer(intercept))
        for (j in 1:d) { ## work through variables
          b <- sm.id[[sm[[i]]$vn[j]]] # list of smooths dependent on this variable
          k <- (1:length(b))[b==i] ## locate current smooth in list 
          if (k>1) for (l in 1:(k-1)) { ## collect X columns
            X1 <- cbind(X1,sm[[b[l]]]$X)
          }
        } ## Now X1 contains columns for all lower dimensional terms
        if (ncol(X1)==as.integer(intercept)) ind <- NULL else
        ind <- fixDependence(X1,sm[[i]]$X,tol=tol)        
        ## ... the columns to zero to ensure independence
        if (!is.null(ind)) { 
          sm[[i]]$X <- sm[[i]]$X[,-ind]
          for (j in 1:length(sm[[i]]$S)) { 
            sm[[i]]$S[[j]] <- sm[[i]]$S[[j]][-ind,-ind]
            sm[[i]]$rank[j] <- qr(sm[[i]]$S[[j]],tol=tol,LAPACK=FALSE)$rank
          }
          sm[[i]]$df <- ncol(sm[[i]]$X)
          attr(sm[[i]],"del.index") <- ind
        }
        sm[[i]]$vn <- NULL
      } ## end if
    } ## end i in 1:m loop
  }  
  sm
}

clone.smooth.spec <- function(specb,spec) {
## produces a version of base smooth.spec, `specb', but with 
## the variables relating to `spec'. Used by `gam.setup' in handling 
## of linked smooths.
 ## check dimensions same...
 if (specb$dim!=spec$dim) stop("`id' linked smooths must have same number of arguments") 
 ## Now start cloning...
 if (inherits(specb,"tensor.smooth.spec")) { ##`te' generated base smooth.spec
    specb$term <- spec$term
    specb$label <- spec$label 
    specb$by <- spec$by
    k <- 1
    for (i in 1:length(specb$margin)) {
      if (is.null(spec$margin)) { ## sloppy user -- have to construct margin info...
         for (j in 1:length(specb$margin[[i]]$term)) {
           specb$margin[[i]]$term[j] <- spec$term[k]
           k <- k + 1
         }
         specb$margin[[i]]$label <- ""
 
      } else { ## second term was at least `te', so margin cloning is easy
        specb$margin[[i]]$term <- spec$margin[[i]]$term
        specb$margin[[i]]$label <- spec$margin[[i]]$label
        specb$margin[[i]]$xt <- spec$margin[[i]]$xt
      }
    }

  } else { ## `s' generated case
    specb$term <- spec$term
    specb$label <- spec$label 
    specb$by <- spec$by
    specb$xt <- spec$xt ## don't generally know what's in here => don't clone
  }
  specb ## return clone
}


parametricPenalty <- function(pterms,assign,paraPen,sp0) {
## routine to process any penalties on the parametric part of the model.
## paraPen is a list whose items have names corresponding to the 
## term.labels in pterms. Each list item may have named elements 
## L, rank and sp. All other elements should be penalty coefficient matrices.
  S <- list()     ## penalty matrix list
  off <- rep(0,0) ## offset array
  rank <- rep(0,0) ## rank array
  sp <- rep(0,0)    ## smoothing param array
  L <- matrix(0,0,0) 
  k <- 0
  tind <- unique(assign) ## unique term indices
  n.t <- length(tind)
  if (n.t>0) for (j in 1:n.t) if (tind[j]>0) {
    term.label <- attr(pterms[tind[j]],"term.label")
    P <- paraPen[[term.label]] ## get any penalty information for this term
    if (!is.null(P)) { ## then there is information
      ind <- (1:length(assign))[assign==tind[j]] ## index of coefs involved here
      Li <- P$L;P$L <- NULL
      spi <- P$sp;P$sp <- NULL
      ranki <- P$rank;P$rank <- NULL
      ## remaining terms should be penalty matrices...
      np <- length(P)

      if (!is.null(ranki)&&length(ranki)!=np) stop("`rank' has wrong length in `paraPen'") 
      if (np) for (i in 1:np) { ## unpack penalty matrices, offsets and ranks
        k <- k + 1
        S[[k]] <- P[[i]]
        off[k] <- min(ind) ## index of first coef penalized by this term
        if ( ncol(P[[i]])!=nrow(P[[i]])||nrow(P[[i]])!=length(ind)) stop(" a parametric penalty has wrong dimension")
        if (is.null(ranki)) {
          ev <- eigen(S[[k]],symmetric=TRUE,only.values=TRUE)$values
          rank[k] <- sum(ev>max(ev)*.Machine$double.eps*10) ## estimate rank
        } else rank[k] <- ranki[i]
      }
      ## now deal with L matrices
      if (np) { ## only do this stuff if there are any penalties!
        if (is.null(Li)) Li <- diag(np)
        if (nrow(Li)!=np) stop("L has wrong dimension in `paraPen'")
        L <- rbind(cbind(L,matrix(0,nrow(L),ncol(Li))),
                   cbind(matrix(0,nrow(Li),ncol(L)),Li))
        ind <- (length(sp)+1):(length(sp)+ncol(Li))
        if (is.null(spi)) {
          sp[ind] <- -1 ## auto-initialize
        } else {
          if (length(spi)!=ncol(Li)) stop("`sp' dimension wrong in `paraPen'")
          sp[ind] <- spi
        }
        ## add smoothing parameter names....
        if (length(ind)>1) names(sp)[ind] <- paste(term.label,ind-ind[1]+1,sep="") 
        else names(sp)[ind] <- term.label
      }
    } ## end !is.null(P)  
  } ## looped through all terms
  if (k==0) return(NULL)
  if (!is.null(sp0)) {
    if (length(sp0)<length(sp)) stop("`sp' too short")
    sp0 <- sp0[1:length(sp)]
    sp[sp<0] <- sp0[sp<0]
  }
  ## S is list of penalty matrices, off[i] is index of first coefficient penalized by each S[[i]]
  ## sp is array of underlying smoothing parameter (-ve to estimate), L is matrix mapping log
  ## underlying smoothing parameters to log smoothing parameters, rank[i] is the rank of S[[i]].
  list(S=S,off=off,sp=sp,L=L,rank=rank)
}

gam.setup <- function(formula,pterms,data=stop("No data supplied to gam.setup"),knots=NULL,sp=NULL,
                    min.sp=NULL,H=NULL,absorb.cons=TRUE,select=FALSE,idLinksBases=TRUE,
                    scale.penalty=TRUE,paraPen=NULL)
# set up the model matrix, penalty matrices and auxilliary information about the smoothing bases
# needed for a gam fit.
{ # split the formula if the object being passed is a formula, otherwise it's already split

  if (inherits(formula,"formula")) split<-interpret.gam(formula) 
  else if (inherits(formula,"split.gam.formula")) split<-formula
  else stop("First argument is no sort of formula!") 
  if (length(split$smooth.spec)==0)
  { if (split$pfok==0) stop("You've got no model....")
    m<-0
  }  
  else  m<-length(split$smooth.spec) # number of smooth terms
  
  G<-list(m=m,min.sp=min.sp,H=H)
  
  if (is.null(attr(data,"terms"))) # then data is not a model frame
  mf<-model.frame(split$pf,data,drop.unused.levels=FALSE) # must be false or can end up with wrong prediction matrix!
  else mf<-data # data is already a model frame

  G$intercept <-  attr(attr(mf,"terms"),"intercept")>0
  G$offset <- model.offset(mf)   # get the model offset (if any)

  # construct strictly parametric model matrix.... 
  
  X <- model.matrix(pterms,mf)
  G$nsdf <- ncol(X)
  G$contrasts <- attr(X,"contrasts")
  G$xlevels <- .getXlevels(pterms,mf)
  G$assign <- attr(X,"assign") # used to tell which coeffs relate to which pterms

  ## now deal with any user supplied penalties on the parametric part of the model...
  PP <- parametricPenalty(pterms,G$assign,paraPen,sp)
  if (!is.null(PP)) { ## strip out supplied sps already used
    ind <- 1:length(PP$sp)
    if (!is.null(sp)) sp <- sp[-ind]
    if (!is.null(min.sp)) { 
      PP$min.sp <- min.sp[ind]
      min.sp <- min.sp[-ind]
    } 
  }
    
##  if (parametric.only) { G$X<-X;return(G)}
  
  # next work through smooth terms (if any) extending model matrix.....
  
  G$smooth<-list()
  G$S<-list()
 

  if (m>0 && idLinksBases) { ## search smooth.spec[[]] for terms linked by common id's
    id.list <- list() ## id information list
    for (i in 1:m) if (!is.null(split$smooth.spec[[i]]$id)) {
      id <- as.character(split$smooth.spec[[i]]$id)
      if (length(id.list)&&id%in%names(id.list)) { ## it's an existing id
        ni <- length(id.list[[id]]$sm.i) ## number of terms so far with this id
        id.list[[id]]$sm.i[ni+1] <- i    ## adding smooth.spec index to this id's list
        ## clone smooth.spec from base smooth spec....
        base.i <- id.list[[id]]$sm.i[1]
         
        split$smooth.spec[[i]] <- clone.smooth.spec(split$smooth.spec[[base.i]],
                                                      split$smooth.spec[[i]])
        
        ## add data for this term to the data list for basis setup...
        temp.term <- split$smooth.spec[[i]]$term
        for (j in 1:length(temp.term)) id.list[[id]]$data[[j]] <- cbind(id.list[[id]]$data[[j]],
                                                          get.var(temp.term[j],data,vecMat=FALSE))
      } else { ## new id
        id.list[[id]] <- list(sm.i=i) ## start the array of indices of smooths with this id
        id.list[[id]]$data <- list()
        ## need to collect together all data for which this basis will be used,
        ## for basis setup...
        term <- split$smooth.spec[[i]]$term
        for (j in 1:length(term)) id.list[[id]]$data[[j]] <- get.var(term[j],data,vecMat=FALSE)
      } ## new id finished
    }
  } ## id.list complete

  G$off<-array(0,0)
  first.para<-G$nsdf+1
  sm <- list()
  newm <- 0
  if (m>0) for (i in 1:m) 
  { # idea here is that terms are set up in accordance with information given in split$smooth.spec
    # appropriate basis constructor is called depending on the class of the smooth
    # constructor returns penalty matrices model matrix and basis specific information
    ## sm[[i]] <- smoothCon(split$smooth.spec[[i]],data,knots,absorb.cons,scale.penalty=scale.penalty) ## old code
    id <- split$smooth.spec[[i]]$id
    if (is.null(id)||!idLinksBases) { ## regular evaluation
      sml <- smoothCon(split$smooth.spec[[i]],data,knots,absorb.cons,scale.penalty=scale.penalty,
                       null.space.penalty=select) 
    } else { ## it's a smooth with an id, so basis setup data differs from model matrix data
      names(id.list[[id]]$data) <- split$smooth.spec[[i]]$term ## give basis data suitable names
      sml <- smoothCon(split$smooth.spec[[i]],id.list[[id]]$data,knots,
                       absorb.cons,n=nrow(data),dataX=data,scale.penalty=scale.penalty,null.space.penalty=select)
    }
    for (j in 1:length(sml)) {
      newm <- newm + 1
      sm[[newm]] <- sml[[j]]
    }
  }
  
  G$m <- m <- newm ## number of actual smooths

  ## The matrix, L, mapping the underlying log smoothing parameters to the
  ## log of the smoothing parameter multiplying the S[[i]] must be
  ## worked out...
  idx <- list() ## idx[[id]]$c contains index of first col in L relating to id
  L <- matrix(0,0,0)
  sp.names <- rep("",0) ## need a list of names to identify sps in global sp array
  if (m>0) for (i in 1:m) {
    id <- sm[[i]]$id
    ## get the L matrix for this smooth...
    length.S <- length(sm[[i]]$S)
    if (is.null(sm[[i]]$L)) Li <- diag(length.S) else Li <- sm[[i]]$L 
    ## extend the global L matrix...
    if (is.null(id)||is.null(idx[[id]])) { ## new `id'     
      if (!is.null(id)) { ## create record in `idx'
        idx[[id]]$c <- ncol(L)+1   ## starting column in L for this `id'
        idx[[id]]$nc <- ncol(Li)   ## number of columns relating to this `id'
      }
      L <- rbind(cbind(L,matrix(0,nrow(L),ncol(Li))),
                 cbind(matrix(0,nrow(Li),ncol(L)),Li))
      if (length.S > 0) { ## there are smoothing parameters to name
        if (length.S == 1) spn <- sm[[i]]$label else 
          spn <- paste(sm[[i]]$label,1:length.S,sep="")
        sp.names <- c(sp.names,spn) ## extend the sp name vector
      }
    } else { ## it's a repeat id => shares existing sp's
      L0 <- matrix(0,nrow(Li),ncol(L))
      if (ncol(Li)>idx[[id]]$nc) {
        stop("Later terms sharing an `id' can not have more smoothing parameters than the first such term")
      }
      L0[,idx[[id]]$c:(idx[[id]]$c+ncol(Li)-1)] <- Li
      L <- rbind(L,L0)
    }
  }


  ## at this stage, it is neccessary to impose any side conditions required
  ## for identifiability
  if (m>0) sm<-gam.side(sm,X,tol=.Machine$double.eps^.5)

  if (m>0) for (i in 1:m) 
  { n.para<-ncol(sm[[i]]$X)
    # define which elements in the parameter vector this smooth relates to....
    sm[[i]]$first.para<-first.para     
    first.para<-first.para+n.para
    sm[[i]]$last.para<-first.para-1
    ## termwise offset handling ...
    Xoff <- attr(sm[[i]]$X,"offset")
    if (!is.null(Xoff)) { 
      if (is.null(G$offset)) G$offset <- Xoff
      else G$offset <- G$offset + Xoff
    }
    ## model matrix accumulation ...
    X<-cbind(X,sm[[i]]$X);sm[[i]]$X<-NULL
   
    G$smooth[[i]] <- sm[[i]]   
  }
  G$cmX <- colMeans(X) ## useful for componentwise CI construction 
  G$X<-X;rm(X)
  n.p<-ncol(G$X) 
  # deal with penalties


## min.sp must be length nrow(L) to make sense
## sp must be length ncol(L) --- need to partition
## L into columns relating to free log smoothing paramters,
## and columns, L0, corresponding to values supplied in sp.
## lsp0 = L0%*%log(sp[sp>=0]) [need to fudge sp==0 case by
## setting log(0) to, e.g. 10*log(.Machine$double.xmin)]


  ## following deals with supplied and estimated smoothing parameters...
  ## first process the `sp' array supplied to `gam'...
 
  if (!is.null(sp)) # then user has supplied fixed smoothing parameters
  { if (length(sp)!=ncol(L)) { warning("Supplied smoothing parameter vector is too short - ignored.")}
    if (sum(is.na(sp))) { warning("NA's in supplied smoothing parameter vector - ignoring.")}
    G$sp <- sp
  } else { # set up for auto-initialization
    G$sp<-rep(-1,ncol(L))
  }
  
  names(G$sp) <- sp.names

  ## now work through the smooths searching for any `sp' elements
  ## supplied in `s' or `te' terms.... This relies on `idx' created 
  ## above...
  
  k <- 1 ## current location in `sp' array
  if (m>0) for (i in 1:m) {
    id <- sm[[i]]$id
    if (is.null(sm[[i]]$L)) Li <- diag(length(sm[[i]]$S)) else Li <- sm[[i]]$L 
    if (is.null(id)) { ## it's a smooth without an id
      spi <- sm[[i]]$sp
      if (!is.null(spi)) { ## sp supplied in `s' or `te'
        if (length(spi)!=ncol(Li)) stop("incorrect number of smoothing parameters supplied for a smooth term")
        G$sp[k:(k+ncol(Li)-1)] <- spi
      }       
      k <- k + ncol(Li) 
    } else { ## smooth has an id
      spi <- sm[[i]]$sp
      if (is.null(idx[[id]]$sp.done)) { ## not already dealt with these sp's
        if (!is.null(spi)) { ## sp supplied in `s' or `te'
          if (length(spi)!=ncol(Li)) stop("incorrect number of smoothing parameters supplied for a smooth term")
          G$sp[idx[[id]]$c:(idx[[id]]$c+idx[[id]]$nc-1)] <- spi
        }
        idx[[id]]$sp.done <- TRUE ## only makes sense to use supplied `sp' from defining term
        k <- k + idx[[id]]$nc 
      }
    }
  } ## finished processing `sp' vectors supplied in `s' or `te' terms

  ## copy initial sp's back into smooth objects, so there is a record of
  ## fixed and free...
  k <- 1 
  if (length(idx)) for (i in 1:length(idx)) idx[[i]]$sp.done <- FALSE
  if (m>0) for (i in 1:m) { ## work through all smooths
    id <- sm[[i]]$id 
    if (!is.null(id)) { ## smooth with id
      if (idx[[id]]$nc>0) { ## only copy if there are sp's
        G$smooth[[i]]$sp <- G$sp[idx[[id]]$c:(idx[[id]]$c+idx[[id]]$nc-1)]
      }   
      if (!idx[[id]]$sp.done) { ## only update k on first encounter with this smooth
        idx[[id]]$sp.done <- TRUE
        k <- k + idx[[id]]$nc
      }
    
    } else { ## no id, just work through sp 
      if (is.null(sm[[i]]$L)) nc <- length(sm[[i]]$S) else nc <- ncol(sm[[i]]$L)
      if (nc>0) G$smooth[[i]]$sp <- G$sp[k:(k+nc-1)]
      k <- k + nc
    }
  } ## now all elements of G$smooth have a record of initial sp. 


  if (!is.null(min.sp)) # then minimum s.p.'s supplied
  { if (length(min.sp)!=nrow(L)) stop("length of min.sp is wrong.")
    if (sum(is.na(min.sp))) stop("NA's in min.sp.")
    if (sum(min.sp<0)) stop("elements of min.sp must be non negative.")
  }

  k.sp<-0 # count through sp and S
  G$rank<-array(0,0)
  if (m>0) for (i in 1:m)
  { sm<-G$smooth[[i]]
    if (length(sm$S)>0)
    for (j in 1:length(sm$S))  # work through penalty matrices
    { k.sp<-k.sp+1
      G$off[k.sp]<-sm$first.para 
      G$S[[k.sp]]<-sm$S[[j]]
      G$rank[k.sp]<-sm$rank[j]
      if (!is.null(min.sp))
      { if (is.null(H)) H<-matrix(0,n.p,n.p)
        H[sm$first.para:sm$last.para,sm$first.para:sm$last.para]<-
        H[sm$first.para:sm$last.para,sm$first.para:sm$last.para]+min.sp[k.sp]*sm$S[[j]] 
      }           
    } 
  }
 
  ## need to modify L, G$S, G$sp, G$rank and G$off to include any penalties
  ## on parametric stuff, at this point....

  if (!is.null(PP)) { ## deal with penalties on parametric terms
    L <- rbind(cbind(L,matrix(0,nrow(L),ncol(PP$L))),
                 cbind(matrix(0,nrow(PP$L),ncol(L)),PP$L))
    G$off <- c(PP$off,G$off)
    G$S <- c(PP$S,G$S)
    G$rank <- c(PP$rank,G$rank)
    G$sp <- c(PP$sp,G$sp)
    if (!is.null(PP$min.sp)) { ## deal with minimum sps
      if (is.null(H)) H <- matrix(0,n.p,n.p)
      for (i in 1:length(PP$S)) {
        ind <- PP$off[i]:(PP$off[i]+ncol(PP$S[[i]])-1)
        H[ind,ind] <- H[ind,ind] + PP$min.sp[i] * PP$S[[i]]
      }
    } ## min.sp stuff finished
  }


  ## Now remove columns of L and rows of sp relating to fixed 
  ## smoothing parameters, and use removed elements to create lsp0

  fix.ind <- G$sp>=0

  if (sum(fix.ind)) {
    lsp0 <- G$sp[fix.ind]
    ind <- lsp0==0
    lsp0[!ind] <- log(lsp0[!ind])
    lsp0[ind] <- log(.Machine$double.xmin)*1000 ## zero fudge
    lsp0 <- L[,fix.ind,drop=FALSE]%*%lsp0

    L <- L[,!fix.ind,drop=FALSE]  
    G$sp <- G$sp[!fix.ind]
  } else {lsp0 <- rep(0,nrow(L))}

  G$H<-H

  if (ncol(L)==nrow(L)&&!sum(L!=diag(ncol(L)))) L <- NULL ## it's just the identity

  G$L <- L;G$lsp0 <- lsp0

  # deal with constraints 
   
  G$C<-matrix(0,0,n.p)
  if (m>0) 
  { for (i in 1:m)
    { if (is.null(G$smooth[[i]]$C)) n.con<-0 
      else n.con<- nrow(G$smooth[[i]]$C)
      C<-matrix(0,n.con,n.p)
      C[,G$smooth[[i]]$first.para:G$smooth[[i]]$last.para]<-G$smooth[[i]]$C
      G$C<-rbind(G$C,C)
      G$smooth[[i]]$C <- NULL
    }
    rm(C)
  }
  G$y <- data[[split$response]]
         
  ##data[[deparse(split$full.formula[[2]],backtick=TRUE)]]
  
  G$n<-nrow(data)

  if (is.null(data$"(weights)")) G$w<-rep(1,G$n)
  else G$w<-data$"(weights)"  

  # now run some checks on the arguments

  
  ### Should check that there are enough unique covariate combinations to support model dimension

  G
}

formula.gam <- function(x, ...)
# formula.lm and formula.glm reconstruct the formula from x$terms, this is 
# problematic because of the way mgcv handles s() and te() terms 
{ x$formula
}


gam.negbin <- function(lsp,fscale,family,control,method,optimizer,gamma,G,scale,...) {
## negative binomial gam fit, using `negbin' family, when some sort of 
## search for theta parameter is required. If the `theta' parameter to `negbin'
## is a length 2 array with theta[2]>theta[1] then `theta is taken as the 
## search interval over which to optimize theta. Otherwise `theta' is taken
## to be an array giving a discrete set of theta values over which to optimize
## by exhaustive search. Note that AIC is used as the criterion, since the 
## deviance depends on theta, UBRE is not proportional to AIC if theta is varied.
 
  if (method%in%c("ML","REML","P-REML","P-ML")) { use.aic <- FALSE;scoreType=method } 
                                           else { use.aic <- TRUE;scoreType="UBRE"}

  theta <- family$getTheta()
  link <- family$link
  if (length(theta)==2&&(theta[2]>theta[1])) { ## perform interval search
    l.theta <- seq(log(theta[2]),log(theta[1]),length=25) ## grid over which to search
    golden <- TRUE    

  } else { ## perform discrete value search
    l.theta <- log(sort(theta,decreasing=TRUE)) ## the supplied grid
    golden <- FALSE
  }
   
  n.th <- length(l.theta)

  mustart <- list(...)[["mustart"]]

  for (i in 1:n.th) { ## search through theta values
    family <- fix.family.link(negbin(theta=exp(l.theta[i]),link=link))
    if (optimizer[2]=="bfgs") b <- bfgs(
                  lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                  offset=G$offset,U1=G$U1,Mp = G$Mp,
                  family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...) else
    if (optimizer[2]=="newton") b <- newton(lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                offset=G$offset,U1=G$U1,Mp = G$Mp,family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...) else
    b <- simplyFit(lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                offset=G$offset,U1=G$U1,Mp = G$Mp,family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...)
    if (use.aic) score <- b$object$aic + 2*b$object$trA ## AIC
    else score <- b$score ## (P-)(RE)ML
    if (i==1 || score<best.score) {
      best.score <- score
      b.est <- b
      ib <- i
    }
    lsp <- b$lsp 
    mustart <- b$object$fitted.values
  } ## end of discrete search `b.est' contains the best model
   
  if (golden) { ## refine `theta' estimate by golden section search
    tau <- 2/(1+sqrt(5)) ## golden ratio
    ## get bracket ....
    if (ib == 1) { lt0 <- l.theta[2];lt1 <- l.theta[1]} else
    if (ib == n.th) { lt0 <- l.theta[n.th];lt1 <- l.theta[n.th-1]} else
    { lt0 <- l.theta[ib+1];lt1 <- l.theta[ib-1]}
    ## initial evaluations
    lsp <- b.est$lsp 
    mustart <- b.est$object$fitted.values
    lt.tau <- lt0 + tau*(lt1-lt0)
    lt.1tau <- lt0 + (1-tau)*(lt1-lt0)
    for (lt in c(lt.1tau,lt.tau))
    { family <- fix.family.link(negbin(theta=exp(lt),link=link))
      if (optimizer[2]=="bfgs") b <- bfgs(
                  lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                  offset=G$offset,U1=G$U1,Mp = G$Mp,family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...) else 
      if (optimizer[2]=="newton") b <- newton(lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                  offset=G$offset,U1=G$U1,Mp = G$Mp,family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...) else
    b <- simplyFit(lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                offset=G$offset,U1=G$U1,Mp = G$Mp,family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...) 
      if (use.aic) score <- b$object$aic + 2*b$object$trA ## AIC
      else score <- b$score 
      if (lt==lt.tau) f.tau <- score else f.1tau <- score
    }
    lsp <- b$lsp 
    mustart <- b$object$fitted.values
#    tol <- abs(lt1-lt0)*1e-7    
    while (round(exp(lt.tau),digits=3)!=round(exp(lt.1tau),digits=3)) {
      if (f.tau<f.1tau) {
        lt0 <- lt.1tau
        lt.1tau <- lt.tau;f.1tau <- f.tau
        lt.new <- lt.tau <- lt0 + tau*(lt1-lt0)
        f.tau.update <- TRUE
      } else {
        lt1 <- lt.tau
        lt.tau <- lt.1tau;f.tau <- f.1tau
        lt.new <- lt.1tau <- lt0 + (1-tau)*(lt1-lt0)
        f.tau.update <- FALSE 
      }
     
      family <- fix.family.link(negbin(theta=exp(lt.new),link=link))
       if (optimizer[2]=="bfgs") b <- bfgs(
                  lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                  offset=G$offset,U1=G$U1,Mp = G$Mp,family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...) else
      if (optimizer[2]=="newton") b <- newton(lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                  offset=G$offset,U1=G$U1,Mp = G$Mp,family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...)  else
      b <- simplyFit(lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,
                offset=G$offset,U1=G$U1,Mp = G$Mp,family=family,weights=G$w,
                  control=control,gamma=gamma,scale=1,conv.tol=control$newton$conv.tol,
                  maxNstep=control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf,
                  printWarn=FALSE,scoreType=scoreType,mustart=mustart,null.coef=G$null.coef,...)
      if (use.aic) score <- b$object$aic + 2*b$object$trA ## AIC
      else score <- b$score       
      if (f.tau.update) f.tau <- score else f.1tau <- score
    }
    b.est <- b
  }
  object <- b.est$object
  object$GACV <- object$D2 <- object$P2 <- object$UBRE2 <- object$trA2 <- 
  object$GACV1 <- object$GACV2 <- object$GCV2 <- object$D1 <- object$P1 <- NULL
  object$sp <- exp(b.est$lsp)
  b <- list(conv=b$conv,iter=b$iter,grad=b$grad,hess=b$hess,score.hist=b$score.hist) ## return info
  object$outer.info <- b
  object$gcv.ubre <- as.numeric(b.est$score)
  object
}





gam.outer <- function(lsp,fscale,family,control,method,optimizer,criterion,scale,gamma,G,...)
# function for smoothing parameter estimation by outer optimization. i.e.
# P-IRLS scheme iterated to convergence for each trial set of smoothing
# parameters.
# MAJOR STEPS:
#  1. Call appropriate smoothing parameter optimizer, and extract fitted model
#    `object'
#  2. Call `magic.post.proc' to get parameter covariance matrices, edf etc to
#     add to `object' 
{ if (is.null(optimizer[2])) optimizer[2] <- "newton"
  if (!optimizer[2]%in%c("newton","bfgs","nlm","optim","nlm.fd")) stop("unknown outer optimization method.")
  if (length(lsp)==0) { ## no sp estimation to do -- run a fit instead
    optimizer[2] <- "no.sps" ## will cause gam2objective to be called, below
  }
  nbGetTheta <- substr(family$family[1],1,17)=="Negative Binomial" && length(family$getTheta())>1
  if (optimizer[2]=="nlm.fd") {
    if (nbGetTheta) stop("nlm.fd not available with negative binomial Theta estimation")
    if (method%in%c("REML","ML","GACV.Cp","P-ML","P-REML")) stop("nlm.fd only available for GCV/UBRE")
    um<-nlm(full.score,lsp,typsize=lsp,fscale=fscale, stepmax = 
            control$nlm$stepmax, ndigit = control$nlm$ndigit,
	    gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, 
            iterlim = control$nlm$iterlim, G=G,family=family,control=control,
            gamma=gamma,...)
    lsp<-um$estimate
    object<-attr(full.score(lsp,G,family,control,gamma=gamma,...),"full.gam.object")
    object$gcv.ubre <- um$minimum
    object$outer.info <- um
    object$sp <- exp(lsp)
    return(object)
  }
  ## some preparations for the other methods 
 
  family <- fix.family.link(family)
  family <- fix.family.var(family)
  if (method%in%c("REML","ML","P-REML","P-ML")) {reml <- TRUE; family <- fix.family.ls(family)} else reml <- FALSE 
  

  if (nbGetTheta) {
    if (!(optimizer[2]%in%c("newton","bfgs","no.sps"))) {
      warning("only outer methods `newton' & `bfgs' supports `negbin' family and theta selection: reset")
      optimizer[2] <- "newton"
    } 
    object <- gam.negbin(lsp,fscale,family,control,method,optimizer,gamma,G,...)
    ## make sure criterion gets set to UBRE
  } else if (optimizer[2]=="newton"||optimizer[2]=="bfgs"){ ## the gam.fit3 method -- not negbin
    if (optimizer[2]=="bfgs") 
    b <- bfgs(lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,offset=G$offset,U1=G$U1,Mp = G$Mp,
                family=family,weights=G$w,control=control,gamma=gamma,scale=scale,conv.tol=control$newton$conv.tol,
                maxNstep= control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf, 
                printWarn=FALSE,scoreType=criterion,null.coef=G$null.coef,...) else
    b <- newton(lsp=lsp,X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,L=G$L,lsp0=G$lsp0,offset=G$offset,U1=G$U1,Mp=G$Mp,
                family=family,weights=G$w,control=control,gamma=gamma,scale=scale,conv.tol=control$newton$conv.tol,
                maxNstep= control$newton$maxNstep,maxSstep=control$newton$maxSstep,maxHalf=control$newton$maxHalf, 
                printWarn=FALSE,scoreType=criterion,null.coef=G$null.coef,...)                
                
    object <- b$object
    object$REML <- object$REML1 <- object$REML2 <-
    object$GACV <- object$D2 <- object$P2 <- object$UBRE2 <- object$trA2 <- 
    object$GACV1 <- object$GACV2 <- object$GCV2 <- object$D1 <- object$P1 <- NULL
    object$sp <- as.numeric(exp(b$lsp))
    object$gcv.ubre <- as.numeric(b$score)
    b <- list(conv=b$conv,iter=b$iter,grad=b$grad,hess=b$hess,score.hist=b$score.hist) ## return info
    object$outer.info <- b   
  } else { ## methods calling gam.fit3 
    args <- list(X=G$X,y=G$y,Eb=G$Eb,UrS=G$UrS,offset=G$offset,U1=G$U1,Mp=G$Mp,family=family,
             weights=G$w,control=control,scoreType=criterion,gamma=gamma,scale=scale,
             L=G$L,lsp0=G$lsp0,null.coef=G$null.coef)
  
    if (optimizer[2]=="nlm") {
       b <- nlm(gam4objective, lsp, typsize = lsp, fscale = fscale, 
            stepmax = control$nlm$stepmax, ndigit = control$nlm$ndigit,
	    gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, 
            iterlim = control$nlm$iterlim,
	    check.analyticals=control$nlm$check.analyticals,
            args=args,...)
      lsp <- b$estimate
      
    } else if (optimizer[2]=="optim") {
      b<-optim(par=lsp,fn=gam2objective,gr=gam2derivative,method="L-BFGS-B",control=
         list(fnscale=fscale,factr=control$optim$factr,lmm=min(5,length(lsp))),args=args,...)
      lsp <- b$par
    } else b <- NULL
    obj <- gam2objective(lsp,args,printWarn=TRUE,...) # final model fit, with warnings 
    object <- attr(obj,"full.fit")
    object$gcv.ubre <- as.numeric(obj) 
    object$outer.info <- b
    object$sp <- exp(lsp)
  } # end of methods calling gam.fit2
  
  
  if (scale>0) object$scale <- scale else object$scale <- object$scale.est 
  
  mv<-magic.post.proc(G$X,object,w=object$weights)
  object$Vp <- mv$Vb
  object$hat<-mv$hat
  object$Ve <- mv$Ve
  object$edf<-mv$edf
  object$aic <- object$aic + 2*sum(mv$edf)
  object$nsdf <- G$nsdf
  object$GCV<-object$GCV1<-object$UBRE<-object$UBRE1<-object$trA<-
  object$trA1<-object$alpha<-object$alpha1<-object$rV<-object$scale.est<-NULL
  object$sig2 <- object$scale
  
  object
}

get.null.coef <- function(G) {
## Get an estimate of the coefs corresponding to maximum reasonable deviance...
  y <- G$y
  weights <- G$w
  nobs <- G$n
  start <- etastart <- mustart <- NULL
  family <- G$family
  eval(family$initialize) ## have to do this to ensure y numeric
  y <- as.numeric(y)
  null.coef <- qr.coef(qr(G$X),family$linkfun(mean(y)+0*y))
  null.coef[is.na(null.coef)] <- 0;
  null.coef
}

estimate.gam <- function (G,method,optimizer,control,in.out,scale,gamma,...) {
## Do gam estimation and smoothness selection...
  
  if (!optimizer[1]%in%c("perf","outer")) stop("unknown optimizer")
  if (!method%in%c("GCV.Cp","GACV.Cp","REML","P-REML","ML","P-ML")) stop("unknown smoothness selection criterion") 

  G$rS <- mini.roots(G$S,G$off,ncol(G$X))

  if (method%in%c("REML","P-REML","ML","P-ML")) {
    if (optimizer[1]=="perf") {
      warning("Reset optimizer to outer/newton") 
      optimizer <- c("outer","newton")
    }
    reml <- TRUE
  } else reml <- FALSE ## experimental insert
  Ssp <- totalPenaltySpace(G$S,G$H,G$off,ncol(G$X))
  # G$U1 <- Ssp$Y       ## range space of penalty
  G$Eb <- Ssp$E       ## balanced penalty square root for rank determination purrposes 
  G$U1 <- cbind(Ssp$Y,Ssp$Z) ## EXPERIMENTAL: eigen space basis
  G$Mp <- ncol(Ssp$Z) ## null space dimension
  G$UrS <- list()     ## need penalty matrices in overall penalty range space...
  if (length(G$S)>0) for (i in 1:length(G$S)) G$UrS[[i]] <- t(Ssp$Y)%*%G$rS[[i]]
  if (!is.null(G$H)) { ## then the sqrt fixed penalty matrix H is needed for (RE)ML 
      G$UrS[[i+1]] <- t(Ssp$Y)%*%mroot(G$H)
  }


  # is outer looping needed ?
  outer.looping <- ((!G$am && (optimizer[1]=="outer"))||reml||method=="GACV.Cp") ## && length(G$S)>0 && sum(G$sp<0)!=0

  ## sort out exact sp selection criterion to use

  fam.name <- G$family$family[1]

  if (scale==0) ## choose criterion for performance iteration
  { if (fam.name == "binomial"||fam.name == "poisson") G$sig2<-1 #ubre
      else G$sig2 <- -1 #gcv
  } else {G$sig2 <- scale}

  if (reml) { ## then RE(ML) selection, but which variant?
   criterion <- method
   if (fam.name == "binomial"||fam.name == "poisson") scale <- 1
  } else {
    if (scale==0) { 
      if (fam.name=="binomial"||fam.name=="poisson") scale<-1 #ubre
      else scale <- -1 #gcv
    }
    if (scale > 0) criterion <- "UBRE"
    else {
      if (method=="GCV.Cp") criterion <- "GCV" else criterion <- "GACV"
    }
  }

  if (fam.name=="quasi"||fam.name=="quasipoisson"||fam.name=="quasibinomial") {
    ## REML/ML invalid with quasi families
    if (method=="REML") method <- "P-REML"
    if (method=="ML") method <- "P-ML"
  }

  if (substr(fam.name,1,17)=="Negative Binomial") { 
    scale <- 1; ## no choise
    if (method%in%c("GCV.Cp","GACV.Cp")) criterion <- "UBRE"
  }
  ## Reset P-ML or P-REML in known scale parameter cases
  if (scale>0) {
    if (method=="P-ML") criterion <- method <- "ML" else 
    if (method=="P-REML")  criterion <- method <- "REML"
  } 


  # take only a few IRLS steps to get scale estimates for "pure" outer
  # looping...
  family <- G$family  
  if (outer.looping) { 
    fixedSteps <- control$outerPIsteps      ## how many performance iteration steps to use for initialization
    if (substr(G$family$family[1],1,17)=="Negative Binomial") { ## initialize sensibly
      scale <- G$sig2 <- 1
      G$family <- negbin(max(family$getTheta()),link=family$link)
    }
  } else fixedSteps <- control$maxit+2
  
  if (outer.looping && !is.null(in.out)) { # initial s.p.s and scale provided
    ok <- TRUE ## run a few basic checks
    if (is.null(in.out$sp)||is.null(in.out$scale)) ok <- FALSE
    if (length(in.out$sp)!=length(G$sp)) ok <- FALSE
    if (!ok) stop("in.out incorrect: see documentation")
    object<-list() # fake enough of a returned fit object for initialization 
    ##object$sp <- in.out$sp[G$all.sp<0] # only use the values for free s.p.s
    object$gcv.ubre <- in.out$scale
    object$sig2 <- 0 ## just means that in.out$scale acts as total scale
  } else ## do performance iteration.... 
  object <- gam.fit(G,family=G$family,control=control,gamma=gamma,fixedSteps=fixedSteps,...)
  
  G$family <- family ## restore, in case manipulated for negative binomial 
    
  if (outer.looping)
  { # use perf.iter s.p. estimates from gam.fit or supplied initial s.p.s as starting values...
    lsp<-log(object$sp) 
    # don't allow PI initial sp's too far from defaults, otherwise optimizers may
    # get stuck on flat portions of GCV/UBRE score....
    if (is.null(in.out)&&length(lsp)>0) { ## note no checks if supplied 
      lsp2 <- log(initial.sp(G$X,G$S,G$off)) 
      if (!is.null(G$L)) { ## estimate underlying smoothing parameters
        if (is.null(G$lsp0)) G$lsp0 <- rep(0,nrow(G$L))
        lsp2 <- as.numeric(coef(lm(lsp2~G$L-1+offset(G$lsp0))))
      }
      ind <- lsp > lsp2+5;lsp[ind] <- lsp2[ind]+5
      ind <- lsp < lsp2-5;lsp[ind] <- lsp2[ind]-5 
      if (fixedSteps<1) lsp <- lsp2 ## don't use perf iter sp's at all
    }
    mgcv.conv <- object$mgcv.conv  
  
    if (criterion%in%c("REML","ML")&&scale<=0) { ## log(scale) to be estimated as a smoothing parameter
      log.scale <-  log(sum(object$weights*object$residuals^2)/(G$n-sum(object$edf)))
      lsp <- c(lsp,log.scale) ## append log initial scale estimate to lsp
      ## extend G$L, if present...
      if (!is.null(G$L)) G$L <- cbind(rbind(G$L,rep(0,ncol(G$L))),c(rep(0,nrow(G$L)),1))
      if (!is.null(G$lsp0)) G$lsp0 <- c(G$lsp0,0)
    }

    ## Get an estimate of the coefs corresponding to maximum reasonable deviance...

    G$null.coef <- get.null.coef(G)

    object <- gam.outer(lsp,fscale=abs(object$gcv.ubre)+object$sig2/length(G$y),family=G$family,
                        control=control,criterion=criterion,method=method,optimizer=optimizer,scale=scale,gamma=gamma,G=G,...)
    
    if (criterion%in%c("REML","ML")&&scale<=0)  object$sp <- 
                                                object$sp[-length(object$sp)] ## drop scale estimate from sp array
    object$mgcv.conv <- mgcv.conv 

  } else { ## performance iteration already complete, but check for all fixed sp case ...
   # if (!G$am && (optimizer[1]=="outer")) {
   #   ## need to fix up GCV/UBRE score 
   #   if (criterion=="UBRE") object$gcv.ubre <- object$deviance/G$n - scale +
   #                          2 * gamma * scale* sum(object$edf)/G$n else 
   #   if (criterion=="GCV") object$gcv.ubre <- G$n *
   #                     object$deviance/(G$n-sum(object$edf))^2 else 
   #   if (criterion=="GACV") { 
   #     P <- sum(object$weights*object$residuals^2)
   #     tau <- sum(object$edf)
   #     object$gcv.ubre <- object$deviance/G$n + 2 * gamma*tau * P / (G$n*(G$n-tau))
   #   }  
   # }
  }

  ## correct null deviance if there's an offset ....

  if (G$intercept&&any(G$offset!=0)) object$null.deviance <-
                                  glm(G$y~offset(G$offset),family=object$family)$deviance

  object$method <- criterion

  object$smooth<-G$smooth
  # now re-assign variable names to coefficients etc. 
  if (G$nsdf>0) term.names<-colnames(G$X)[1:G$nsdf] else term.names<-array("",0)
  n.smooth<-length(G$smooth)
  if (n.smooth)
  for (i in 1:n.smooth)
  { k<-1
    for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para)
    { term.names[j]<-paste(G$smooth[[i]]$label,".",as.character(k),sep="")
      k<-k+1
    }
  }
  names(object$coefficients) <- term.names  # note - won't work on matrices!!
  names(object$edf) <- term.names
  object
}

gam <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,offset=NULL,
                method="GCV.Cp",optimizer=c("outer","newton"),control=gam.control(),
                scale=0,select=FALSE,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=TRUE,
                paraPen=NULL,G=NULL,in.out=NULL,...)

# Routine to fit a GAM to some data. The model is stated in the formula, which is then 
# interpreted to figure out which bits relate to smooth terms and which to parametric terms.

{  if (is.null(G))
   { # create model frame..... 
    gp<-interpret.gam(formula) # interpret the formula 
    cl<-match.call() # call needed in gam object for update to work
    mf<-match.call(expand.dots=FALSE)
    mf$formula<-gp$fake.formula 
    mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-mf$select <-
               mf$gamma<-mf$method<-mf$fit<-mf$paraPen<-mf$G<-mf$optimizer <- mf$...<-NULL
    mf$drop.unused.levels<-TRUE
    mf[[1]]<-as.name("model.frame")
    pmf <- mf
    mf <- eval(mf, parent.frame()) # the model frame now contains all the data 
    if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf,"terms")
    
    pmf$formula <- gp$pf
    pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part

    pterms <- attr(pmf,"terms") ## pmf only used for this

    if (is.character(family)) family<-eval(parse(text=family))
    if (is.function(family)) family <- family()
    if (is.null(family$family)) stop("family not recognized")
  
    if (family$family[1]=="gaussian" && family$link=="identity") am <- TRUE
    else am <- FALSE
    
    G<-gam.setup(gp,pterms=pterms,data=mf,knots=knots,sp=sp,min.sp=min.sp,
                 H=H,absorb.cons=TRUE,select=select,
                 idLinksBases=control$idLinksBases,scale.penalty=control$scalePenalty,
                 paraPen=paraPen)
    
    G$family <- family
   
    if (ncol(G$X)>nrow(G$X)+nrow(G$C)) stop("Model has more coefficients than data")

    G$terms<-terms;G$pterms<-pterms
    G$mf<-mf;G$cl<-cl;
    G$am <- am

    if (is.null(G$offset)) G$offset<-rep(0,G$n)
     
    G$min.edf<-G$nsdf-dim(G$C)[1]
    if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim

    G$formula<-formula
    environment(G$formula)<-environment(formula)
  }

  if (!fit) return(G)
  
  G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
  G$max.half<-control$mgcv.half # max step halving in Newton update mgcv

  object <- estimate.gam(G,method,optimizer,control,in.out,scale,gamma,...)

  
  if (!is.null(G$L)) object$full.sp <- as.numeric(exp(G$L%*%log(object$sp)+G$lsp0))
  names(object$sp) <- names(G$sp)
  object$formula<-G$formula
  object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
  object$model<-G$mf # store the model frame
  object$na.action <- attr(G$mf,"na.action") # how to deal with NA's
  object$control <- control
  object$terms <- G$terms
  object$pterms <- G$pterms
  object$assign <- G$assign # applies only to pterms
  object$contrasts <- G$contrasts
  object$xlevels <- G$xlevels
  object$offset <- G$offset
  object$data <- data
  object$df.residual <- nrow(G$X) - sum(object$edf)
  object$min.edf<-G$min.edf
  if (G$am&&!(method%in%c("REML","ML","P-ML","P-REML"))) object$optimizer <- "magic" else object$optimizer <- optimizer
  object$call<-G$cl # needed for update() to work
  class(object)<-c("gam","glm","lm")
  object
}


gam.check <- function(b)
# takes a fitted gam object and produces some standard diagnostic plots
{ if (b$method%in%c("GCV","GACV","UBRE","REML","ML","P-ML","P-REML"))
  { old.par<-par(mfrow=c(2,2))
    sc.name<-b$method
    qqnorm(residuals(b))
    plot(b$linear.predictors,residuals(b),main="Resids vs. linear pred.",
         xlab="linear predictor",ylab="residuals");
    hist(residuals(b),xlab="Residuals",main="Histogram of residuals");
    plot(fitted(b),b$y,xlab="Fitted Values",ylab="Response",main="Response vs. Fitted Values")
    
    ## now summarize convergence information 
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
      } else { ## just default print of information...
        cat("\n");print(b$outer.info)
      }
    } else { ## perf iter or AM case
      if (b$mgcv.conv$iter==0) 
      cat("\nModel required no smoothing parameter selection")
      else { 
        cat("\nSmoothing parameter selection converged after",b$mgcv.conv$iter,"iteration")       
        if (b$mgcv.conv$iter>1) cat("s")
         
        if (!b$mgcv.conv$fully.converged)
        cat(" by steepest\ndescent step failure.\n") else cat(".\n")
        cat("The RMS",sc.name,"score gradiant at convergence was",b$mgcv.conv$rms.grad,".\n")
        if (b$mgcv.conv$hess.pos.def)
        cat("The Hessian was positive definite.\n") else cat("The Hessian was not positive definite.\n")
        cat("The estimated model rank was ",b$mgcv.conv$rank,
                   " (maximum possible: ",b$mgcv.conv$full.rank,")\n",sep="")
      }
    }
    cat("\n")
    par(old.par)
  } else ## probably a `gamm' `gam' object
  plot(b$linear.predictor,residuals(b),xlab="linear predictor",ylab="residuals")
}

print.gam<-function (x,...) 
# default print function for gam objects
{ print(x$family)
  cat("Formula:\n")
  print(x$formula)
  n.smooth<-length(x$smooth)
  if (n.smooth==0)
  cat("Total model degrees of freedom",sum(x$edf),"\n")
  else
  { edf<-0
    cat("\nEstimated degrees of freedom:\n")
    for (i in 1:n.smooth)
    edf[i]<-sum(x$edf[x$smooth[[i]]$first.para:x$smooth[[i]]$last.para])
    edf.str <- format(edf,digits=5)
    for (i in 1:n.smooth) {   
    cat(edf.str[i]," ",sep="")
      if (i%%7==0) cat("\n")
    }
    cat(" total =",sum(x$edf),"\n")
  }
  if (!is.null(x$method)&&!(x$method%in%c("PQL","lme.ML","lme.REML")))  
  cat("\n",x$method," score: ",x$gcv.ubre,"\n",sep="")
  invisible(x)
}

gam.control <- function (irls.reg=0.0,epsilon = 1e-06, maxit = 100,
                         mgcv.tol=1e-7,mgcv.half=15,trace =FALSE,
                         rank.tol=.Machine$double.eps^0.5,
                         nlm=list(),optim=list(),newton=list(),outerPIsteps=1,
                         idLinksBases=TRUE,scalePenalty=TRUE) 
# Control structure for a gam. 
# irls.reg is the regularization parameter to use in the GAM fitting IRLS loop.
# epsilon is the tolerance to use in the IRLS MLE loop. maxit is the number 
# of IRLS iterations to use. mgcv.tol is the tolerance to use in the mgcv call within each IRLS. mgcv.half is the 
# number of step halvings to employ in the mgcv search for the optimal GCV score, before giving up 
# on a search direction. trace turns on or off some de-bugging information.
# rank.tol is the tolerance to use for rank determination
# outerPIsteps is the number of performance iteration steps used to intialize
#                         outer iteration
{  
    if (!is.numeric(irls.reg) || irls.reg <0.0) stop("IRLS regularizing parameter must be a non-negative number.")
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    if (rank.tol<0||rank.tol>1) 
    { rank.tol=.Machine$double.eps^0.5
      warning("silly value supplied for rank.tol: reset to square root of machine precision.")
    }
    # work through nlm defaults
    if (is.null(nlm$ndigit)||nlm$ndigit<2) nlm$ndigit <- max(2,ceiling(-log10(epsilon)))
    nlm$ndigit <- round(nlm$ndigit)
    ndigit <- floor(-log10(.Machine$double.eps))
    if (nlm$ndigit>ndigit) nlm$ndigit <- ndigit
    if (is.null(nlm$gradtol)) nlm$gradtol <- epsilon*100
    nlm$gradtol <- abs(nlm$gradtol)
    ## note that nlm will stop after hitting stepmax 5 consecutive times
    ## hence should not be set too small ... 
    if (is.null(nlm$stepmax)||nlm$stepmax==0) nlm$stepmax <- 2
    nlm$stepmax <- abs(nlm$stepmax)
    if (is.null(nlm$steptol)) nlm$steptol <- 1e-4
    nlm$steptol <- abs(nlm$steptol)
    if (is.null(nlm$iterlim)) nlm$iterlim <- 200
    nlm$iterlim <- abs(nlm$iterlim)
    ## Should be reset for a while anytime derivative code altered...
    if (is.null(nlm$check.analyticals)) nlm$check.analyticals <- FALSE
    nlm$check.analyticals <- as.logical(nlm$check.analyticals) 

    # and newton defaults
    if (is.null(newton$conv.tol)) newton$conv.tol <- 1e-6
    if (is.null(newton$maxNstep)) newton$maxNstep <- 5
    if (is.null(newton$maxSstep)) newton$maxSstep <- 2
    if (is.null(newton$maxHalf)) newton$maxHalf <- 30
    if (is.null(newton$use.svd)) newton$use.svd <- FALSE

    # and optim defaults
    if (is.null(optim$factr)) optim$factr <- 1e7
    optim$factr <- abs(optim$factr)

    list(irls.reg=irls.reg,epsilon = epsilon, maxit = maxit,
         trace = trace, mgcv.tol=mgcv.tol,mgcv.half=mgcv.half,
         rank.tol=rank.tol,nlm=nlm,
         optim=optim,newton=newton,outerPIsteps=outerPIsteps,
         idLinksBases=idLinksBases,scalePenalty=scalePenalty)
    
}



mgcv.get.scale<-function(Theta,weights,good,mu,mu.eta.val,G)
# Get scale implied by current fit and trial -ve binom Theta, I've used
# mu and mu.eta.val used in fit rather than implied by it....
{ variance<- negbin(Theta)$variance
  w<-sqrt(weights[good]*mu.eta.val[good]^2/variance(mu)[good])
  wres<-w*(G$y-G$X%*%G$p)
  sum(wres^2)/(G$n-sum(G$edf))
}


mgcv.find.theta<-function(Theta,T.max,T.min,weights,good,mu,mu.eta.val,G,tol)
# searches for -ve binomial theta between given limits to get scale=1 
{ scale<-mgcv.get.scale(Theta,weights,good,mu,mu.eta.val,G)
  T.hi<-T.low<-Theta
  while (scale<1&&T.hi<T.max) 
  { T.hi<-T.hi*2
    T.hi<-min(T.hi,T.max)
    scale<-mgcv.get.scale(T.hi,weights,good,mu,mu.eta.val,G)
  } 
  if (all.equal(T.hi,T.max)==TRUE && scale<1) return(T.hi)
  T.low<-T.hi
  while (scale>=1&&T.low>T.min)
  { T.low<-T.low/2 
    T.low<-max(T.low,T.min)
    scale<-mgcv.get.scale(T.low,weights,good,mu,mu.eta.val,G)
  } 
  if (all.equal(T.low,T.min)==TRUE && scale>1) return(T.low)
  # (T.low,T.hi) now brackets scale=1. 
  while (abs(scale-1)>tol)
  { Theta<-(T.low+T.hi)/2
    scale<-mgcv.get.scale(Theta,weights,good,mu,mu.eta.val,G)
    if (scale<1) T.low<-Theta
    else T.hi<-Theta
  }
  Theta
}


full.score <- function(sp,G,family,control,gamma,...)
# function suitable for calling from nlm in order to polish gam fit
# so that actual minimum of score is found in generalized cases
{ if (is.null(G$L)) {
    G$sp<-exp(sp);
  } else {
    G$sp <- as.numeric(exp(G$L%*%sp + G$lsp0))
  }
  # set up single fixed penalty....
  q<-NCOL(G$X)
  if (is.null(G$H)) G$H<-matrix(0,q,q)
  for (i in 1:length(G$S))
  { j<-ncol(G$S[[i]])
    off1<-G$off[i];off2<-off1+j-1
    G$H[off1:off2,off1:off2]<-G$H[off1:off2,off1:off2]+G$sp[i]*G$S[[i]]
  }
  G$S<-list() # have to reset since length of this is used as number of penalties
  G$L <- NULL
  xx<-gam.fit(G,family=family,control=control,gamma=gamma,...)
  res <- xx$gcv.ubre.dev
  attr(res,"full.gam.object")<-xx
  res
}

gam.fit <- function (G, start = NULL, etastart = NULL, 
    mustart = NULL, family = gaussian(), 
    control = gam.control(),gamma=1,
    fixedSteps=(control$maxit+1),...) 
# fitting function for a gam, modified from glm.fit.
# note that smoothing parameter estimates from one irls iterate are carried over to the next irls iterate
# unless the range of s.p.s is large enough that numerical problems might be encountered (want to avoid 
# completely flat parts of gcv/ubre score). In the latter case autoinitialization is requested.
# fixedSteps < its default causes at most fixedSteps iterations to be taken,
# without warning if convergence has not been achieved. This is useful for
# obtaining starting values for outer iteration.
{   fisher <- TRUE
    if (!fisher) { ## Newton needs extra derivatives...
      family <- fix.family.link(family)
      family <- fix.family.var(family)
      if (family$link==family$canonical) fisher <- TRUE
    }
    intercept<-G$intercept
    conv <- FALSE
    n <- nobs <- NROW(G$y) ## n just there to keep codetools happy
    nvars <- NCOL(G$X) # check this needed
    y<-G$y # original data
    X<-G$X # original design matrix
    if (nvars == 0) stop("Model seems to contain no terms")
    olm <- G$am   # only need 1 iteration as it's a pure additive model. 
    find.theta<-FALSE # any supplied -ve binomial theta treated as known, G$sig2 is scale parameter
    if (substr(family$family[1],1,17)=="Negative Binomial")
    { Theta <- family$getTheta()
      if (length(Theta)==1) { ## Theta fixed
        find.theta <- FALSE
        G$sig2 <- 1
      } else {
        if (length(Theta)>2)
        warning("Discrete Theta search not available with performance iteration")
        Theta <- range(Theta)
        T.max <- Theta[2]          ## upper search limit
        T.min <- Theta[1]          ## lower search limit
        Theta <- sqrt(T.max*T.min) ## initial value
        find.theta <- TRUE
      }
      nb.link<-family$link # negative.binomial family, there's a choise of links
    }

    
    # obtain average element sizes for the penalties
    n.S<-length(G$S)
    if (n.S>0)
    { S.size<-0
      for (i in 1:n.S) S.size[i]<-mean(abs(G$S[[i]])) 
    }
    weights<-G$w # original weights

    n.score <- sum(weights!=0) ## n to use in GCV score (i.e. don't count points with no influence)   

    offset<-G$offset 

    variance <- family$variance;dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv;linkfun <- family$linkfun;mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    if (is.null(mustart))   # new from version 1.5.0 
    { eval(family$initialize)}
    else 
    { mukeep <- mustart
      eval(family$initialize)
      mustart <- mukeep
    }
    if (NCOL(y) > 1) 
        stop("y must be univariate unless binomial")
    
    coefold <- NULL                 # 1.5.0
    eta <- if (!is.null(etastart))  # 1.5.0
        etastart

    else if (!is.null(start)) 
    if (length(start) != nvars) 
    stop(paste("Length of start should equal", nvars,
        "and correspond to initial coefs.")) # 1.5.0
    else 
    { coefold<-start                        #1.5.0
      offset+as.vector(if (NCOL(G$X) == 1)
       G$X * start
       else G$X %*% start)
    }
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) 
        stop("Can't find valid starting values: please specify some")
    devold <- sum(dev.resids(y, mu, weights))
   
    boundary <- FALSE
    scale<-G$sig2

    msp <- G$sp
    magic.control<-list(tol=G$conv.tol,step.half=G$max.half,maxit=control$maxit+control$globit,
                          rank.tol=control$rank.tol)

    for (iter in 1:(control$maxit)) 
    {
        good <- weights > 0
        varmu <- variance(mu)[good]
        if (any(is.na(varmu))) 
            stop("NAs in V(mu)")
        if (any(varmu == 0)) 
            stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good]))) 
            stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0) # note good modified here => must re-calc each iter
        if (all(!good)) {
            conv <- FALSE
            warning(paste("No observations informative at iteration", 
                iter))
            break
        }
   
        mevg <- mu.eta.val[good];mug <- mu[good];yg <- y[good]
        weg <- weights[good];etag <- eta[good]
        var.mug<-variance(mug)

        if (fisher) { ## Conventional Fisher scoring
              G$y <- z <- (eta - offset)[good] + (yg - mug)/mevg
              w <- sqrt((weg * mevg^2)/var.mug)
        } else { ## full Newton (actually this is a problem as w can be negative!!)
              c <- yg - mug
              e <- mevg*(1 + c*(family$dvar(mug)/mevg+var.mug*family$d2link(mug))*mevg/var.mug)
              G$y <- z <- (eta - offset)[good] + c/e ## offset subtracted as eta = X%*%beta + offset
              w <- sqrt(weg*e*mevg/var.mug)
        }
        
        G$w<-w
        G$X<-X[good,,drop=FALSE]  # truncated design matrix       
     
        # must set G$sig2 to scale parameter or -1 here....
        G$sig2<-scale


        if (sum(!is.finite(G$y))+sum(!is.finite(G$w))>0) 
        stop("iterative weights or data non-finite in gam.fit - regularization may help. See ?gam.control.")

        ## solve the working weighted penalized LS problem ...

        mr<-magic(G$y,G$X,msp,G$S,G$off,L=G$L,lsp0=G$lsp0,G$rank,G$H,G$C,G$w,gamma=gamma,G$sig2,G$sig2<0,
                    ridge.parameter=control$irls.reg,control=magic.control,n.score=n.score)
        G$p<-mr$b;msp<-mr$sp;G$sig2<-mr$scale;G$gcv.ubre<-mr$score;

        if (find.theta) {# then family is negative binomial with unknown theta - estimate it here from G$sig2
            ##  need to get edf array
          mv<-magic.post.proc(G$X,mr,w=G$w^2)
          G$edf <- mv$edf

          Theta<-mgcv.find.theta(Theta,T.max,T.min,weights,good,mu,mu.eta.val,G,.Machine$double.eps^0.5)
          family<-do.call("negbin",list(theta=Theta,link=nb.link))
          variance <- family$variance;dev.resids <- family$dev.resids
          aic <- family$aic
          family$Theta <- Theta ## save Theta estimate in family
        }

        if (any(!is.finite(G$p))) {
            conv <- FALSE   
            warning(paste("Non-finite coefficients at iteration",iter))
            break
        }

		
        start <- G$p
        eta <- drop(X %*% start) # 1.5.0
        mu <- linkinv(eta <- eta + offset)
        eta <- linkfun(mu) # force eta/mu consistency even if linkinv truncates
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace) 
            cat("Deviance =", dev, "Iterations -", iter, "\n")
        boundary <- FALSE
        if (!is.finite(dev)) {
            if (is.null(coefold))
            stop("no valid set of coefficients has been found:please supply starting values",
            call. = FALSE)
            warning("Step size truncated due to divergence",call.=FALSE)
            ii <- 1
            while (!is.finite(dev)) {
                if (ii > control$maxit) 
                  stop("inner loop 1; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                eta<-drop(X %*% start)
                mu <- linkinv(eta <- eta + offset)
                eta <- linkfun(mu) 
                dev <- sum(dev.resids(y, mu, weights))
            }
            boundary <- TRUE
            if (control$trace) 
                cat("Step halved: new deviance =", dev, "\n")
        }
        if (!(valideta(eta) && validmu(mu))) {
            warning("Step size truncated: out of bounds.",call.=FALSE)
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit) 
                  stop("inner loop 2; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                eta<-drop(X %*% start)
                mu <- linkinv(eta <- eta + offset)
                eta<-linkfun(mu)
            }
            boundary <- TRUE
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Step halved: new deviance =", dev, "\n")
        }

        ## Test for convergence here ...

        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon || olm ||
            iter >= fixedSteps) {
            conv <- TRUE
            coef <- start #1.5.0
            break
        }
        else {
            devold <- dev
            coefold <- coef<-start
        }
    }
    if (!conv) 
    { warning("Algorithm did not converge") 
    }
    if (boundary) 
        warning("Algorithm stopped at boundary value")
    eps <- 10 * .Machine$double.eps
    if (family$family[1] == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps)) 
            warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family[1] == "poisson") {
        if (any(mu < eps)) 
            warning("fitted rates numerically 0 occurred")
    }
    
    residuals <- rep(NA, nobs)
    residuals[good] <- z - (eta - offset)[good]
    
    wt <- rep(0, nobs)
    wt[good] <- w^2
   
    wtdmu <- if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    
    ## Extract a little more information from the fit....

    mv<-magic.post.proc(G$X,mr,w=G$w^2)
    G$Vp<-mv$Vb;G$hat<-mv$hat;
    G$Ve <- mv$Ve # frequentist cov. matrix
    G$edf<-mv$edf
    G$conv<-mr$gcv.info
    G$sp<-msp
    rank<-G$conv$rank

    aic.model <- aic(y, n, mu, weights, dev) + 2 * sum(G$edf)
    if (scale < 0) { ## deviance based GCV
      gcv.ubre.dev <- n.score*dev/(n.score-gamma*sum(G$edf))^2
    } else { # deviance based UBRE, which is just AIC
      gcv.ubre.dev <- dev/n.score + 2 * gamma * sum(G$edf)/n.score - G$sig2
    }
  
	list(coefficients = as.vector(coef), residuals = residuals, fitted.values = mu, 
        family = family,linear.predictors = eta, deviance = dev,
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,  
        df.null = nulldf, y = y, converged = conv,sig2=G$sig2,edf=G$edf,hat=G$hat,
        boundary = boundary,sp = G$sp,nsdf=G$nsdf,Ve=G$Ve,Vp=G$Vp,mgcv.conv=G$conv,
        gcv.ubre=G$gcv.ubre,aic=aic.model,rank=rank,gcv.ubre.dev=gcv.ubre.dev)
}


model.matrix.gam <- function(object,...)
{ if (!inherits(object,"gam")) stop("`object' is not of class \"gam\"")
  predict.gam(object,type="lpmatrix",...)
}

predict.gam <- function(object,newdata,type="link",se.fit=FALSE,terms=NULL,
                       block.size=1000,newdata.guaranteed=FALSE,na.action=na.pass,...) 
{
# This function is used for predicting from a GAM. object is a gam object, newdata a dataframe to
# be used in prediction......
#
# Type == "link"     - for linear predictor
#      == "response" - for fitted values
#      == "terms"    - for individual terms on scale of linear predictor 
#      == "iterms"   - exactly as "terms" except that se's include uncertainty about mean  
#      == "lpmatrix" - for matrix mapping parameters to l.p.
# Steps are:
#  1. Set newdata to object$model if no newdata supplied
#  2. split up newdata into manageable blocks if too large
#  3. Obtain parametric model matrix (safely!)
#  4. Work through smooths calling prediction.matrix constructors for each term
#  5. Work out required quantities
# 
# The splitting into blocks enables blocks of compiled code to be called efficiently
# using smooth class specific prediciton matrix constructors, without having to 
# build up potentially enormous prediction matrices.
# if newdata.guaranteed == TRUE then the data.frame is assumed complete and
# ready to go, so that only factor levels are checked for sanity.
# 
# if `terms' is non null then it should be a list of terms to be returned 
# when type=="terms" or "iterms". 
# if `object' has an attribute `para.only' then only parametric terms of order
# 1 are returned for type=="terms"/"iterms" : i.e. only what termplot can handle.
#
# if no new data is supplied then na.action does nothing, otherwise 
# if na.action == "na.pass" then NA predictors result in NA predictions (as lm
#                   or glm)
#              == "na.omit" or "na.exclude" then NA predictors result in
#                       dropping

  if (type!="link"&&type!="terms"&&type!="iterms"&&type!="response"&&type!="lpmatrix")  
  { warning("Unknown type, reset to terms.")
    type<-"terms"
  }
  if (!inherits(object,"gam")) stop("predict.gam can only be used to predict from gam objects")
  ## to mimic behaviour of predict.lm, some resetting is required ...
  if (missing(newdata)) na.act <- object$na.action else
  { if (is.null(na.action)) na.act <- NULL 
    else {
      na.txt <- "na.pass"
      if (is.character(na.action))
      na.txt <- substitute(na.action) else
      if (is.function(na.action)) na.txt <- deparse(substitute(na.action))
      if (na.txt=="na.pass") na.act <- "na.exclude" else
      if (na.txt=="na.exclude") na.act <- "na.omit" else
      na.act <- na.action
    }
  } ## ... done
  # get data from which to predict.....  
  nd.is.mf <- FALSE # need to flag if supplied newdata is already a model frame
  if (newdata.guaranteed==FALSE)
  { if (missing(newdata)) # then "fake" an object suitable for prediction 
    { newdata<-object$model
      new.data.ok <- FALSE
      nd.is.mf <- TRUE
    }
    else  # do an R ``standard'' evaluation to pick up data
    { new.data.ok <- TRUE
      if (is.data.frame(newdata)&&!is.null(attr(newdata,"terms"))) # it's a model frame
      { if (sum(!(names(object$model)%in%names(newdata)))) stop(
        "newdata is a model.frame: it should contain all required variables\n")
         nd.is.mf <- TRUE
      } else
      { ## Following is non-standard to allow convenient splitting into blocks
        ## below, and to allow checking that all variables are in newdata ...

        ## get names of required variables, less response, but including offset variable
        Terms <- delete.response(terms(object))
        allNames <- all.vars(Terms)
        ff <- reformulate(allNames)
        if (sum(!(allNames%in%names(newdata)))) { 
        warning("not all required variables have been supplied in  newdata!\n")}
        ## note that `xlev' argument not used here, otherwise `as.factor' in 
        ## formula can cause a problem ... levels reset later.
        newdata <- eval(model.frame(ff,data=newdata,na.action=na.act),parent.frame()) 
        na.act <- attr(newdata,"na.action")
      }
    }
  } else {na.act <- NULL}
  

  if (new.data.ok)
  { ## check factor levels are right ...
    names(newdata)->nn # new data names
    colnames(object$model)->mn # original names
    for (i in 1:length(newdata)) 
    if (nn[i]%in%mn && is.factor(object$model[,nn[i]])) # then so should newdata[[i]] be 
    { newdata[[i]]<-factor(newdata[[i]],levels=levels(object$model[,nn[i]])) # set prediction levels to fit levels
    }

    # split prediction into blocks, to avoid running out of memory
    if (length(newdata)==1) newdata[[2]]<-newdata[[1]] # avoids data frame losing its labels and dimensions below!
    if (is.null(dim(newdata[[1]]))) np<-length(newdata[[1]]) 
    else np<-dim(newdata[[1]])[1] 
    nb<-length(object$coefficients)
    if (block.size<1) block.size <- np
    n.blocks<-np%/%block.size
    b.size<-rep(block.size,n.blocks)
    last.block<-np-sum(b.size)
    if (last.block>0) 
    { n.blocks<-n.blocks+1  
      b.size[n.blocks]<-last.block
    }
  } else # no new data, just use object$model
  { np <- nrow(object$model)
    nb <- length(object$coefficients)
    n.blocks <- 1
    b.size <- array(np,1)
  }
  # setup prediction arrays
  n.smooth<-length(object$smooth)
  if (type=="lpmatrix")
  { H<-matrix(0,np,nb)
  } else
  if (type=="terms"||type=="iterms")
  { term.labels<-attr(object$pterms,"term.labels")
    if (is.null(attr(object,"para.only"))) para.only <-FALSE else
    para.only <- TRUE  # if true then only return information on parametric part
    n.pterms <- length(term.labels)
    fit<-array(0,c(np,n.pterms+as.numeric(!para.only)*n.smooth))
    if (se.fit) se<-fit
    ColNames<-term.labels
  } else
  { fit<-array(0,np)
    if (se.fit) se<-fit
  }
  stop<-0

  Terms <- delete.response(object$pterms)
  s.offset <- NULL # to accumulate any smooth term specific offset
  any.soff <- FALSE # indicator of term specific offset existence
  for (b in 1:n.blocks)  # work through prediction blocks
  { start<-stop+1
    stop<-start+b.size[b]-1
    if (n.blocks==1) data <- newdata else data<-newdata[start:stop,]
    X <- matrix(0,b.size[b],nb)
    Xoff <- matrix(0,b.size[b],n.smooth) ## term specific offsets 
    ## implements safe prediction for parametric part as described in
    ## http://developer.r-project.org/model-fitting-functions.txt
    if (new.data.ok)
    { if (nd.is.mf) mf <- model.frame(data,xlev=object$xlevels) else
      { mf <- model.frame(Terms,data,xlev=object$xlevels)
        if (!is.null(cl <- attr(object$pterms,"dataClasses"))) .checkMFClasses(cl,mf)
      } 
      Xp <- model.matrix(Terms,mf,contrasts=object$contrasts) 
    } else 
    { Xp <- model.matrix(Terms,object$model)
      mf <- newdata # needed in case of offset, below
    }
    
    if (object$nsdf) X[,1:object$nsdf]<-Xp
    if (n.smooth) for (k in 1:n.smooth) 
    { Xfrag <- PredictMat(object$smooth[[k]],data)		 
      X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
      Xfrag.off <- attr(Xfrag,"offset") ## any term specific offsets?
      if (!is.null(Xfrag.off)) { Xoff[,k] <- Xfrag.off; any.soff <- TRUE }
      if (type=="terms"||type=="iterms") ColNames[n.pterms+k]<-object$smooth[[k]]$label
    }
    # have prediction matrix for this block, now do something with it
    if (type=="lpmatrix") { 
      H[start:stop,]<-X
      if (any.soff) s.offset <- rbind(s.offset,Xoff)
    } else 
    if (type=="terms"||type=="iterms")
    {
      ind <- 1:length(object$assign)
      if (n.pterms)  # work through parametric part
      for (i in 1:n.pterms)
      { ii <- ind[object$assign==i]
        fit[start:stop,i] <- as.matrix(X[,ii])%*%object$coefficients[ii]
        if (se.fit) se[start:stop,i]<-
        sqrt(rowSums((as.matrix(X[,ii])%*%object$Vp[ii,ii])*as.matrix(X[,ii])))
      }

      if (n.smooth&&!para.only) 
      { for (k in 1:n.smooth) # work through the smooth terms 
        { first<-object$smooth[[k]]$first.para;last<-object$smooth[[k]]$last.para
          fit[start:stop,n.pterms+k]<-X[,first:last]%*%object$coefficients[first:last] + Xoff[,k]
          if (se.fit) { # diag(Z%*%V%*%t(Z))^0.5; Z=X[,first:last]; V is sub-matrix of Vp
            if (type=="iterms"&&inherits(attr(object$smooth[[k]],"qrc"),"qr")) { ## termwise se to "carry the intercept"
              X1 <- matrix(object$cmX,nrow(X),ncol(X),byrow=TRUE)
              meanL1 <- object$smooth[[k]]$meanL1
              if (!is.null(meanL1)) X1 <- X1 / meanL1              
              X1[,first:last] <- X[,first:last]
              se[start:stop,n.pterms+k] <- sqrt(rowSums((X1%*%object$Vp)*X1))
            } else se[start:stop,n.pterms+k] <- ## terms strictly centred
            sqrt(rowSums((X[,first:last]%*%object$Vp[first:last,first:last])*X[,first:last]))
          } ## end if (se.fit)
        }
        colnames(fit) <- ColNames
        if (se.fit) colnames(se) <- ColNames
      } else { # para.only
        colnames(fit) <- term.labels
        if (se.fit) colnames(se) <- term.labels
        # retain only terms of order 1 - this is to make termplot work
        order <- attr(object$pterms,"order")
        term.labels <- term.labels[order==1]
        fit <- as.matrix(as.matrix(fit)[,order==1])
        colnames(fit) <- term.labels
        if (se.fit) { se <- as.matrix(as.matrix(se)[,order==1])
        colnames(se) <- term.labels } 
      }
      if (!is.null(terms)) # return only terms requested via `terms'
      { if (sum(!(terms %in%colnames(fit)))) 
        warning("non-existent terms requested - ignoring")
        else { names(term.labels) <- term.labels
          term.labels <- term.labels[terms]  # names lost if only one col
          fit <- as.matrix(as.matrix(fit)[,terms])
          colnames(fit) <- term.labels
          if (se.fit) {se <- as.matrix(as.matrix(se)[,terms])
          colnames(se) <- term.labels}
        }
      }
    } else # "link" or "response"
    { k<-attr(attr(object$model,"terms"),"offset")
      fit[start:stop]<-X%*%object$coefficients + rowSums(Xoff)
      if (!is.null(k)) fit[start:stop]<-fit[start:stop]+model.offset(mf) + rowSums(Xoff)
      if (se.fit) se[start:stop]<-sqrt(rowSums((X%*%object$Vp)*X))
      if (type=="response") # transform    
      { fam<-object$family;linkinv<-fam$linkinv;dmu.deta<-fam$mu.eta  
        if (se.fit) se[start:stop]<-se[start:stop]*abs(dmu.deta(fit[start:stop])) 
        fit[start:stop]<-linkinv(fit[start:stop])
      }
    }
  } ## end of prediction block loop
  rn <- rownames(newdata)
  if (type=="lpmatrix") { 
    colnames(H) <- names(object$coefficients);rownames(H)<-rn
    if (!is.null(s.offset)) { 
      s.offset <- napredict(na.act,s.offset)
      attr(H,"offset") <- s.offset
    }
    H <- napredict(na.act,H)
  } else { 
    if (se.fit) { 
      if (is.null(nrow(fit))) {
        names(fit) <- rn
        names(se) <- rn
        fit <- napredict(na.act,fit)
        se <- napredict(na.act,se) 
      } else { 
        rownames(fit)<-rn
        rownames(se)<-rn
        fit <- napredict(na.act,fit)
        se <- napredict(na.act,se)
      }
      H<-list(fit=fit,se.fit=se) 
    } else { 
      H<-fit
      if (is.null(nrow(H))) names(H) <- rn else
      rownames(H)<-rn
      H <- napredict(na.act,H)
    }
  }
  if (type=="terms"||type=="iterms") attr(H,"constant") <- object$coefficients[1]
  H # ... and return
}

plot.gam <- function(x,residuals=FALSE,rug=TRUE,se=TRUE,pages=0,select=NULL,scale=-1,n=100,n2=40,
                     pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,all.terms=FALSE,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,seWithMean=FALSE,by.resids=FALSE,...)

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

{ sub.edf <- function(lab,edf) {
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


  sp.contour <- function(x,y,z,zse,xlab="",ylab="",zlab="",titleOnly=FALSE,
               se.plot=TRUE,se.mult=1,trans=I,shift=0,...)   
  # internal function for contouring 2-d smooths with 1 s.e. limits
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
##  plot(x,y,type="n",xlab="",ylab="",axes=FALSE)
    cs<-(yr/10)/strheight(zlab);if (cs>1) cs<-1 # text scaling based on height  
##  cw<-par()$cxy[1]  
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
##  contour(x,y,z,levels=zlev,lwd=2,labcex=cs*0.65,axes=FALSE,add=TRUE)  
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
#    contour(x,y,z+zse,levels=zlev,add=TRUE,lty=2,col=2,labcex=cs*0.5)  

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
#    contour(x,y,z-zse,levels=zlev,add=TRUE,lty=3,col=3,labcex=cs*0.5)  
    
    if (!titleOnly) {
      xpos<-xrange[2]-xr/5  
      xl<-c(xpos,xpos+xr/10);  
      lines(xl,yl,xpd=TRUE,lty=args$lty,col=args$col)  
      text(xpos+xr/10,ypos,paste("+",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)  
    }
  }  ## end of sp.contour

  # start of main function
  w.resid<-NULL
  if (length(residuals)>1) # residuals supplied 
  { if (length(residuals)==length(x$residuals)) 
    w.resid <- residuals else
    warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  } else partial.resids <- residuals # use working residuals or none
  m<-length(x$smooth) # number of smooth terms
  order <- attr(x$pterms,"order") # array giving order of each parametric term
  if (all.terms) # plot parametric terms as well
  n.para <- sum(order==1) # plotable parametric terms   
  else n.para <- 0 
  if (m+n.para==0) stop("No terms to plot - nothing for plot.gam() to do.")
  if (se)
  { if (is.numeric(se)) se2.mult<-se1.mult<-se else { se1.mult<-2;se2.mult<-1} 
    if (se1.mult<0) se1.mult<-0;if (se2.mult<0) se2.mult<-0
  } else se1.mult<-se2.mult<-1
  
  if (se && x$Vp[1,1]<=0) 
  { se<-FALSE
    warning("No variance estimates available")
  }
  # plot should ignore all "by" variables
  
  # sort out number of pages and plots per page
  n.plots <- m + n.para
  if (pages>n.plots) pages<-n.plots
  if (pages<0) pages<-0
  if (pages!=0)    # figure out how to display things
  { ppp<-n.plots%/%pages
    if (n.plots%%pages!=0) 
    { ppp<-ppp+1
      while (ppp*(pages-1)>=n.plots) pages<-pages-1
 ##     if (n.plots%%pages) last.pages<-0 ##else last.ppp<-n.plots-ppp*pages
    } 
 ## else last.ppp<-0
    # now figure out number of rows and columns
    c<-trunc(sqrt(ppp))
	if (c<1) c<-1
    r<-ppp%/%c
    if (r<1) r<-1
    while (r*c<ppp) r<-r+1
    while (r*c-ppp >c && r>1) r<-r-1
    while (r*c-ppp >r && c>1) c<-c-1 
    oldpar<-par(mfrow=c(r,c))
  
  } else
  { ppp<-1;oldpar<-par()}
  
  if ((pages==0&&prod(par("mfcol"))<n.plots&&dev.interactive())||
       pages>1&&dev.interactive()) ask <- TRUE else ask <- FALSE 

  ##if (pages==0&&is.null(select)) par(mfrow=par("mfrow")) ## new display

  if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

  # work through all smooth terms assembling the plot data list pd with elements
  # dim, x, fit, se, ylab, xlab for 1-d terms;
  # dim, xm, ym, fit, se, ylab, xlab, title for 2-d terms;
  # and dim otherwise
  if (partial.resids) 
  { fv.terms <- predict(x,type="terms")
    if (is.null(w.resid)) w.resid<-x$residuals*sqrt(x$weights) # weighted working residuals
  }
  pd<-list();
  i<-1 # needs a value if no smooths, but parametric terms ...
  if (m>0) for (i in 1:m) # work through smooth terms
  { if (x$smooth[[i]]$dim==1)
    { raw<-x$model[x$smooth[[i]]$term]
      xx<-seq(min(raw),max(raw),length=n)   # generate x sequence for prediction
      if (x$smooth[[i]]$by!="NA")         # deal with any by variables
      { by<-rep(1,n);dat<-data.frame(x=xx,by=by)
        names(dat)<-c(x$smooth[[i]]$term,x$smooth[[i]]$by)
      } else
      { dat<-data.frame(x=xx);names(dat)<-x$smooth[[i]]$term}  # prediction data.frame
      X <- PredictMat(x$smooth[[i]],dat)   # prediction matrix from this term
      first<-x$smooth[[i]]$first.para;last<-x$smooth[[i]]$last.para
      p<-x$coefficients[first:last]       # relevent coefficients 
      offset <- attr(X,"offset")
      if (is.null(offset)) 
      fit <- X%*%p else fit<-X%*%p + offset       # fitted values
      if (se) {
        ## test whether mean variability to be added to variability (only for centred terms)
        if (seWithMean&&inherits(attr(x$smooth[[i]],"qrc"),"qr")) {
          X1 <- matrix(x$cmX,nrow(X),ncol(x$Vp),byrow=TRUE)
          meanL1 <- x$smooth[[i]]$meanL1
          if (!is.null(meanL1)) X1 <- X1 / meanL1
          X1[,first:last] <- X
          se.fit <- sqrt(rowSums((X1%*%x$Vp)*X1))
        } else se.fit <- ## se in centred (or anyway unconstained) space only
        sqrt(rowSums((X%*%x$Vp[first:last,first:last])*X))
      }
      edf<-sum(x$edf[first:last])
      xterm <- x$smooth[[i]]$term
      if (is.null(xlab)) xlabel<- xterm else xlabel <- xlab
      if (is.null(ylab)) ylabel <- sub.edf(x$smooth[[i]]$label,edf) else
                         ylabel <- ylab
      pd.item<-list(fit=fit,dim=1,x=xx,ylab=ylabel,xlab=xlabel,raw=raw[[1]])
      if (partial.resids) {pd.item$p.resid <- fv.terms[,length(order)+i]+w.resid}
      if (se) pd.item$se=se.fit*se1.mult  # Note multiplier
      pd[[i]]<-pd.item;rm(pd.item)
    } else 
    if (x$smooth[[i]]$dim==2)
    { xterm <- x$smooth[[i]]$term[1]
      if (is.null(xlab)) xlabel <- xterm else xlabel <- xlab
      yterm <- x$smooth[[i]]$term[2]
      if (is.null(ylab)) ylabel <- yterm else ylabel <- ylab
      raw<-data.frame(x=as.numeric(x$model[xterm][[1]]),
                      y=as.numeric(x$model[yterm][[1]]))
      n2<-max(10,n2)
      xm<-seq(min(raw$x),max(raw$x),length=n2)
      ym<-seq(min(raw$y),max(raw$y),length=n2)  
      xx<-rep(xm,n2)
      yy<-rep(ym,rep(n2,n2))
      if (too.far>0)
      exclude <- exclude.too.far(xx,yy,raw$x,raw$y,dist=too.far) else
      exclude <- rep(FALSE,n2*n2)
      if (x$smooth[[i]]$by!="NA")         # deal with any by variables
      { by<-rep(1,n2^2);dat<-data.frame(x=xx,y=yy,by=by)
        names(dat)<-c(xterm,yterm,x$smooth[[i]]$by)
      } else
      { dat<-data.frame(x=xx,y=yy);names(dat)<-c(xterm,yterm)}  # prediction data.frame
      X <- PredictMat(x$smooth[[i]],dat)   # prediction matrix for this term
      first<-x$smooth[[i]]$first.para;last<-x$smooth[[i]]$last.para
      p<-x$coefficients[first:last]      # relevent coefficients 
      offset <- attr(X,"offset")
      if (is.null(offset)) 
      fit <- X%*%p else fit<-X%*%p + offset       # fitted values
      fit[exclude] <- NA                 # exclude grid points too far from data
      if (se) {  
        if (seWithMean&&inherits(attr(x$smooth[[i]],"qrc"),"qr")) { ## then se to include uncertainty in overall mean
          X1 <- matrix(x$cmX,nrow(X),ncol(x$Vp),byrow=TRUE)
          meanL1 <- x$smooth[[i]]$meanL1
          if (!is.null(meanL1)) X1 <- X1 / meanL1
          X1[,first:last] <- X
          se.fit <- sqrt(rowSums((X1%*%x$Vp)*X1))
        } else se.fit <- ## se in centred space only
        sqrt(rowSums((X%*%x$Vp[first:last,first:last])*X))

        se.fit[exclude] <- NA # exclude grid points too distant from data
      }
      edf<-sum(x$edf[first:last])
      if (is.null(main)) 
      { title <- sub.edf(x$smooth[[i]]$label,edf)
      }
      else title <- main
      pd.item<-list(fit=fit,dim=2,xm=xm,ym=ym,ylab=ylabel,xlab=xlabel,title=title,raw=raw)
      if (is.null(ylim)) pd.item$ylim <- range(ym) else pd.item$ylim <- ylim
      if (is.null(xlim)) pd.item$xlim <- range(xm) else pd.item$xlim <- xlim
      if (se) pd.item$se=se.fit*se2.mult  # Note multiplier
      pd[[i]]<-pd.item;rm(pd.item)
    } else
    { pd[[i]]<-list(dim=x$smooth[[i]]$dim)}
  }

  
  # now plot .....
  if (se)   # pd$fit and pd$se
  { k<-0
    if (scale==-1&&is.null(ylim)) # getting common scale for 1-d terms
    if (m>0) for (i in 1:m)
    { if (pd[[i]]$dim==1)
      { ul<-pd[[i]]$fit+pd[[i]]$se
        ll<-pd[[i]]$fit-pd[[i]]$se
        if (k==0) 
        { ylim<-c(min(ll),max(ul));k<-1;
        } else
        { if (min(ll)<ylim[1]) ylim[1]<-min(ll)
	  if (max(ul)>ylim[2]) ylim[2]<-max(ul)
        }
        if (partial.resids)
        { ul <- max(pd[[i]]$p.resid,na.rm=TRUE)
          if (ul > ylim[2]) ylim[2] <- ul
          ll <-  min(pd[[i]]$p.resid,na.rm=TRUE)
          if (ll < ylim[1]) ylim[1] <- ll
        }
      }
    }
    j<-1
    if (m>0) for (i in 1:m)
    { if (is.null(select)||i==select)
      { ##if (interactive()&& is.null(select) && pd[[i]]$dim<3 && i>1&&(i-1)%%ppp==0) 
        ##readline("Press return for next page....")
        if (pd[[i]]$dim==1)
        { ul<-pd[[i]]$fit+pd[[i]]$se
          ll<-pd[[i]]$fit-pd[[i]]$se
          if (scale==0&&is.null(ylim)) 
          { ylimit<-c(min(ll),max(ul))
            if (partial.resids)
            { max.r <- max(pd[[i]]$p.resid,na.rm=TRUE)
              if ( max.r> ylimit[2]) ylimit[2] <- max.r
              min.r <-  min(pd[[i]]$p.resid,na.rm=TRUE)
              if (min.r < ylimit[1]) ylimit[1] <- min.r
            }
          }
          if (!is.null(ylim)) ylimit <- ylim
          if (shade)
          { plot(pd[[i]]$x,trans(pd[[i]]$fit+shift),type="n",xlab=pd[[i]]$xlab,ylim=trans(ylimit+shift),
                 xlim=xlim,ylab=pd[[i]]$ylab,main=main,...)
            polygon(c(pd[[i]]$x,pd[[i]]$x[n:1],pd[[i]]$x[1]),
                     trans(c(ul,ll[n:1],ul[1])+shift),col = shade.col,border = NA)
            lines(pd[[i]]$x,trans(pd[[i]]$fit+shift))
          } else
          { plot(pd[[i]]$x,trans(pd[[i]]$fit+shift),type="l",xlab=pd[[i]]$xlab,ylim=trans(ylimit+shift),xlim=xlim,
                 ylab=pd[[i]]$ylab,main=main,...)
	    if (is.null(list(...)[["lty"]]))
            { lines(pd[[i]]$x,trans(ul+shift),lty=2,...)
              lines(pd[[i]]$x,trans(ll+shift),lty=2,...)
            } else
            { lines(pd[[i]]$x,trans(ul+shift),...)
              lines(pd[[i]]$x,trans(ll+shift),...)
            }
          } 
          if (partial.resids&&(by.resids||x$smooth[[i]]$by=="NA"))
          { if (length(pd[[i]]$raw)==length(pd[[i]]$p.resid)) {
              if (is.null(list(...)[["pch"]]))
              points(pd[[i]]$raw,trans(pd[[i]]$p.resid+shift),pch=".",...) else
              points(pd[[i]]$raw,trans(pd[[i]]$p.resid+shift),...) 
            } else {
              warning("Partial residuals do not have a natural x-axis location for linear functional terms")
            }
          }
	  if (rug) 
          { if (jit) rug(jitter(as.numeric(pd[[i]]$raw)),...)
             else rug(as.numeric(pd[[i]]$raw),...)
	  }
        } else if (pd[[i]]$dim==2)
        { 
          if (pers) 
          { if (!is.null(main)) pd[[i]]$title <- main
            persp(pd[[i]]$xm,pd[[i]]$ym,matrix(trans(pd[[i]]$fit+shift),n2,n2),xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,
                  zlab=pd[[i]]$title,ylim=pd[[i]]$ylim,xlim=pd[[i]]$xlim,theta=theta,phi=phi,...)
          } else
          { sp.contour(pd[[i]]$xm,pd[[i]]$ym,matrix(pd[[i]]$fit,n2,n2),matrix(pd[[i]]$se,n2,n2),
                     xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,zlab=pd[[i]]$title,titleOnly=!is.null(main),
                     se.mult=se2.mult,trans=trans,shift=shift,...)
            if (rug) { 
              if (is.null(list(...)[["pch"]]))
              points(pd[[i]]$raw$x,pd[[i]]$raw$y,pch=".",...) else
              points(pd[[i]]$raw$x,pd[[i]]$raw$y,...) 
            }
          } 
        } else
        { warning("no automatic plotting for smooths of more than two variables")
        }
      }  
      j<-j+pd[[i]]$dim
    }
  } else # don't plot confidence limits
  { k<-0
    if (scale==-1&&is.null(ylim))
    if (m>0) for (i in 1:m)
    { if (pd[[i]]$dim==1)
      { if (k==0) { 
          if (partial.resids) ylim <- range(pd[[i]]$p.resid,na.rm=TRUE) else 
          ylim<-range(pd[[i]]$fit);k<-1 
        } else
        { if (partial.resids)
          { if (min(pd[[i]]$p.resid)<ylim[1]) ylim[1]<-min(pd[[i]]$p.resid,na.rm=TRUE)
	    if (max(pd[[i]]$p.resid)>ylim[2]) ylim[2]<-max(pd[[i]]$p.resid,na.rm=TRUE)
          } else
          { if (min(pd[[i]]$fit)<ylim[1]) ylim[1]<-min(pd[[i]]$fit)
	    if (max(pd[[i]]$fit)>ylim[2]) ylim[2]<-max(pd[[i]]$fit)
          }
	}
      }
    }
    j<-1
    if (m>0) for (i in 1:m)
    { if (is.null(select)||i==select)
      {### if (interactive() && is.null(select) && pd[[i]]$dim<3 && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
        if (pd[[i]]$dim==1)
        { if (scale==0&&is.null(ylim)) 
          { if (partial.resids) ylimit <- range(pd[[i]]$p.resid,na.rm=TRUE) else ylimit <-range(pd[[i]]$fit)}
          if (!is.null(ylim)) ylimit <- ylim
          plot(pd[[i]]$x,trans(pd[[i]]$fit+shift),type="l",,xlab=pd[[i]]$xlab,
               ylab=pd[[i]]$ylab,ylim=trans(ylimit+shift),xlim=xlim,main=main,...)
          if (rug) 
	  { if (jit) rug(jitter(as.numeric(pd[[i]]$raw)),...)
            else rug(as.numeric(pd[[i]]$raw),...) 
          }
          if (partial.resids&&(by.resids||x$smooth[[i]]$by=="NA"))
          { if (is.null(list(...)[["pch"]]))
            points(pd[[i]]$raw,trans(pd[[i]]$p.resid+shift),pch=".",...) else
            points(pd[[i]]$raw,trans(pd[[i]]$p.resid+shift),...)
          }
        } else if (pd[[i]]$dim==2)
        { if (!is.null(main)) pd[[i]]$title <- main
          if (pers) 
          { persp(pd[[i]]$xm,pd[[i]]$ym,matrix(trans(pd[[i]]$fit+shift),n2,n2),xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,
                          zlab=pd[[i]]$title,theta=theta,phi=phi,xlim=pd[[i]]$xlim,ylim=pd[[i]]$ylim,...)
          }
          else
          { contour(pd[[i]]$xm,pd[[i]]$ym,matrix(trans(pd[[i]]$fit+shift),n2,n2),xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,
                    main=pd[[i]]$title,xlim=pd[[i]]$xlim,ylim=pd[[i]]$ylim,...)
            if (rug) 
            {  if (is.null(list(...)[["pch"]])) points(pd[[i]]$raw$x,pd[[i]]$raw$y,pch=".",...) else
               points(pd[[i]]$raw$x,pd[[i]]$raw$y,...)
            }
          }  

        } else
        { warning("no automatic plotting for smooths of more than one variable")}
      }
      j<-j+pd[[i]]$dim
    } 
  }
  if (n.para>0) # plot parameteric terms
  { class(x) <- c("gam","glm","lm") # needed to get termplot to call model.frame.glm 
    if (is.null(select)) {
      attr(x,"para.only") <- TRUE
    #  if (interactive() && m && i%%ppp==0) 
    #  readline("Press return for next page....")
      termplot(x,se=se,rug=rug,col.se=1,col.term=1)
    } else { # figure out which plot is required
      if (select > m) { 
        select <- select - m # i.e. which parametric term
        term.labels <- attr(x$pterms,"term.labels")
        term.labels <- term.labels[order==1]
        if (select <= length(term.labels)) {
        if (interactive() && m &&i%%ppp==0) 
##        readline("Press return for next page....")
        termplot(x,terms=term.labels[select],se=se,rug=rug,col.se=1,col.term=1)
        }  
      }
    }
  }
  if (pages>0) par(oldpar)
}


residuals.gam <-function(object, type = c("deviance", "pearson","scaled.pearson", "working", "response"),...)
# calculates residuals for gam object - default for glm (from which this is adapted) seems to be buggy
{ type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
##  family <- object$family
  wts <- object$prior.weights
  res<- switch(type,working = object$residuals,
         scaled.pearson = (y-mu)*sqrt(wts)/sqrt(object$sig2*object$family$variance(mu)),
              pearson = (y-mu)*sqrt(wts)/sqrt(object$family$variance(mu)),
              deviance = { d.res<-sqrt(pmax(object$family$dev.resids(y,mu,wts),0))
                           ifelse(y>mu , d.res, -d.res)             
                         },
              response = y - mu)
  res <- naresid(object$na.action,res)
  res
}


## Start of anova and summary (with contributions from Henric Nilsson) ....

eigXVX <- function(X,V,rank=NULL,tol=.Machine$double.eps^.5) {
## forms truncated eigen-decomposition of XVX', efficiently,
## where V is symmetric, and X has more rows than columns
## first `rank' eigen values/vectors are returned, where `rank'
## is the smaller of any non-NULL supplied value, and the rank
## estimated using `tol'
  qrx <- qr(X)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  ind <- abs(ed$values) > max(abs(ed$values))*tol
  erank <- sum(ind) ## empirical rank
  if (is.null(rank)) {
    rank <- erank
  } else { if (rank<erank) ind <- 1:rank else rank <- erank }
  vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
  list(values=ed$values[ind],vectors=vec[,ind],rank=rank)
}

pinvXVX <- function(X,V,rank=NULL) {
## Routine for forming fractionally trunctated
## pseudoinverse of XVX'. Returns as D where
## DD' gives the pseudoinverse itself.
  k <- floor(rank)
  nu <- rank - k
#  if (k < 1) { k <- 1; nu <- 0}
  if (nu>0) k1 <- k+1 else k1 <- k

  qrx <- qr(X)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)

  ## Get the eigenvectors...
  vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
  if (k1<ncol(vec)) vec <- vec[,1:k1,drop=FALSE]
  if (k==0) {
     vec <- t(t(vec)*sqrt(nu/ed$val[1]))
     return(vec)
  }
 
  ## deal with the fractional part of the pinv...
  if (nu>0) {
     if (k>1) vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
     b12 <- .5*nu*(1-nu)
     if (b12<0) b12 <- 0
     b12 <- sqrt(b12)
     B <- matrix(c(1,b12,b12,nu),2,2)
     ev <- diag(ed$values[k:k1]^-.5)
     B <- ev%*%B%*%ev
     eb <- eigen(B,symmetric=TRUE)
     rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
     vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
  } else {
    vec <- t(t(vec)/sqrt(ed$val[1:k]))
  }
  vec ## vec%*%t(vec) is the pseudoinverse
}




summary.gam <- function (object, dispersion = NULL, freq = FALSE,alpha=0, ...) 
# summary method for gam object - provides approximate p values for terms + other diagnostics
# Improved by Henric Nilsson
{ pinv<-function(V,M,rank.tol=1e-6)
  { D<-La.svd(V)
    M1<-length(D$d[D$d>rank.tol*D$d[1]])
    if (M>M1) M<-M1 # avoid problems with zero eigen-values
    if (M+1<=length(D$d)) D$d[(M+1):length(D$d)]<-1
    D$d<- 1/D$d
    if (M+1<=length(D$d)) D$d[(M+1):length(D$d)]<-0
    res <- D$u%*%diag(D$d)%*%D$v
    attr(res,"rank") <- M
    res
  }
  p.table <- pTerms.table <- s.table <- NULL
  if (freq) covmat <- object$Ve else covmat <- object$Vp
  name <- names(object$edf)
  dimnames(covmat) <- list(name, name)
  covmat.unscaled <- covmat/object$sig2
  est.disp <- TRUE
  if(object$method == "UBRE") est.disp <- FALSE
  if (!is.null(dispersion)) { 
    covmat <- dispersion * covmat.unscaled
    est.disp <- FALSE
  } else dispersion <- object$sig2
  se<-0;for (i in 1:length(object$coefficients)) se[i] <- covmat[i,i]^0.5
  residual.df<-length(object$y)-sum(object$edf)
  if (object$nsdf>0) # individual parameters
  { p.coeff<-object$coefficients[1:object$nsdf]
    p.se <- se[1:object$nsdf]
    p.t<-p.coeff/p.se
    if (!est.disp) {
      p.pv<-2*pnorm(abs(p.t),lower.tail=FALSE)
      p.table<-cbind(p.coeff, p.se, p.t, p.pv)   
      dimnames(p.table) <- list(names(p.coeff), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    } else {
      p.pv<-2*pt(abs(p.t),df=residual.df,lower.tail=FALSE)
      p.table<-cbind(p.coeff, p.se, p.t, p.pv)
      dimnames(p.table) <- list(names(p.coeff), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    }    
  } else {p.coeff<-p.t<-p.pv<-array(0,0)}
  
  term.labels<-attr(object$pterms,"term.labels")
  nt<-length(term.labels)
  if (nt>0) # individual parametric terms
  { np<-length(object$assign)
    Vb<-matrix(covmat[1:np,1:np],np,np)
    bp<-array(object$coefficients[1:np],np)
    pTerms.pv <- array(0,nt)
    attr(pTerms.pv,"names") <- term.labels
    pTerms.df <- pTerms.chi.sq <- pTerms.pv
    for (i in 1:nt)
    { ind <- object$assign==i
      b <- bp[ind];V <- Vb[ind,ind]
      ## pseudo-inverse needed in case of truncation of parametric space 
      if (length(b)==1) { 
        V <- 1/V 
        pTerms.df[i] <- nb <- 1      
        pTerms.chi.sq[i] <- V*b*b
      } else {
        V <- pinv(V,length(b),rank.tol=.Machine$double.eps^.5)
        pTerms.df[i] <- nb <- attr(V,"rank")      
        pTerms.chi.sq[i] <- t(b)%*%V%*%b
      }
      if (!est.disp)
      pTerms.pv[i]<-pchisq(pTerms.chi.sq[i],df=nb,lower.tail=FALSE)
      else
      pTerms.pv[i]<-pf(pTerms.chi.sq[i]/nb,df1=nb,df2=residual.df,lower.tail=FALSE)      
    }
    if (!est.disp) {      
      pTerms.table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)   
      dimnames(pTerms.table) <- list(term.labels, c("df", "Chi.sq", "p-value"))
    } else {
      pTerms.table <- cbind(pTerms.df, pTerms.chi.sq/pTerms.df, pTerms.pv)   
      dimnames(pTerms.table) <- list(term.labels, c("df", "F", "p-value"))
    }
  } else { pTerms.df<-pTerms.chi.sq<-pTerms.pv<-array(0,0)}

  m<-length(object$smooth) # number of smooth terms
  df <- edf <- s.pv <- chi.sq <- array(0, m)
  if (m>0) # form test statistics for each smooth
  { if (!freq) { 
      X <- model.matrix(object)
      ## get corrected edf
      ##  edf1 <- 2*object$edf - rowSums(object$Ve*(t(X)%*%X))/object$sig2
    }
    for (i in 1:m)
    { start<-object$smooth[[i]]$first.para;stop<-object$smooth[[i]]$last.para
      V <- covmat[start:stop,start:stop] # cov matrix for smooth
      p<-object$coefficients[start:stop]  # params for smooth
      edf[i]<-sum(object$edf[start:stop]) # edf for this smooth
      if (freq) { ## old style frequentist
        M1<-object$smooth[[i]]$df
        M<-min(M1,ceiling(2*sum(object$edf[start:stop]))) ## upper limit of 2*edf on rank
        V<-pinv(V,M) # get rank M pseudoinverse of V
        chi.sq[i]<-t(p)%*%V%*%p
        df[i] <- attr(V, "rank")
      } else { ## Nychka statistics
        Xt <- X[,start:stop] 
        ft <- Xt%*%p
        if (FALSE) { ## current release version
          trial.rank <- ceiling(edf[i]) ## R 2.7.0 ceiling is not as advertised!
          if (edf[i]-trial.rank>0) trial.rank <- trial.rank+1
     
          ed <- eigXVX(Xt,V,trial.rank)
          
          iv <- 1/ed$values
       
          ## t(ft)%*%ginv(Ats)%*%ft where Ats = Xt%*%Vt%*%t(Xt), efficiently calculated...
          chi.sq[i] <- sum(((t(ed$vectors)%*%ft)*sqrt(iv))^2)
          df[i] <- edf[i] + .5
        } else { ## experimental version --- problematic
          df[i] <- min(ncol(Xt),edf[i])
          D <- pinvXVX(Xt,V,df[i])
          df[i] <- df[i]+alpha*sum(object$smooth[[i]]$sp<0) ## i.e. alpha * (number free sp's)
          chi.sq[i] <- sum((t(D)%*%ft)^2)   
          if (FALSE) { ## full interval inversion  
            if (inherits(attr(object$smooth[[i]],"qrc"),"qr")) {
              X1 <- matrix(object$cmX,nrow(Xt),ncol(object$Vp),byrow=TRUE)
              meanL1 <- object$smooth[[i]]$meanL1
              if (!is.null(meanL1)) X1 <- X1 / meanL1
              X1[,start:stop] <- Xt
              df[i] <- edf[i]
              D <- pinvXVX(X1,object$Vp,df[i]+1)
              #se.fit <- sqrt(rowSums((X1%*%x$Vp)*X1))
            } else { ## se in centred (or anyway unconstained) space only
              df[i] <- edf[i]
              D <- pinvXVX(Xt,V,df[i])
            }
            fm <- sum(D%*%(t(D)%*%ft))/sum(colSums(D)^2) ## re-centering
            chi.sq[i] <- sum((t(D)%*%(ft-fm))^2)  
          }
        } ## end experimental version
      }
      names(chi.sq)[i]<- object$smooth[[i]]$label
      if (!est.disp)
      s.pv[i]<-pchisq(chi.sq[i], df = df[i], lower.tail = FALSE)
      else
      s.pv[i] <- pf(chi.sq[i]/df[i], df1 = df[i], df2 = residual.df, lower.tail = FALSE)
      ## p-values are meaningless for very small edf. Need to set to NA
      if (df[i] < 0.5) s.pv[i] <- NA
    }
    if (!est.disp) {
      if (freq) {
        s.table <- cbind(edf, df, chi.sq, s.pv)      
        dimnames(s.table) <- list(names(chi.sq), c("edf", "Est.rank", "Chi.sq", "p-value"))
      } else {
        s.table <- cbind(edf, df, chi.sq, s.pv)      
        dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "Chi.sq", "p-value"))
      }
    } else {
      if (freq) {
        s.table <- cbind(edf, df, chi.sq/df, s.pv)      
        dimnames(s.table) <- list(names(chi.sq), c("edf", "Est.rank", "F", "p-value"))
      } else {
        s.table <- cbind(edf, df, chi.sq/df, s.pv)      
        dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "F", "p-value"))
      }
    }
  }
  w <- object$prior.weights
  nobs <- nrow(object$model)
  r.sq<- 1 - var(w*(object$y-object$fitted.values))*(nobs-1)/(var(w*object$y)*residual.df) 
  dev.expl<-(object$null.deviance-object$deviance)/object$null.deviance
  ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,
       s.pv=s.pv,scale=dispersion,r.sq=r.sq,family=object$family,formula=object$formula,n=nobs,
       dev.expl=dev.expl,edf=edf,dispersion=dispersion,pTerms.pv=pTerms.pv,pTerms.chi.sq=pTerms.chi.sq,
       pTerms.df = pTerms.df, cov.unscaled = covmat.unscaled, cov.scaled = covmat, p.table = p.table,
       pTerms.table = pTerms.table, s.table = s.table,method=object$method,sp.criterion=object$gcv.ubre)
 
  class(ret)<-"summary.gam"
  ret
}

print.summary.gam <- function(x, digits = max(3, getOption("digits") - 3), 
                              signif.stars = getOption("show.signif.stars"), ...)
# print method for gam summary method. Improved by Henric Nilsson
{ print(x$family)
  cat("Formula:\n")
  print(x$formula)
  if (length(x$p.coeff)>0)
  { cat("\nParametric coefficients:\n")
    printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n")
  if(x$m>0)
  { cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
  }
  cat("\nR-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5))
  if (length(x$dev.expl)>0) cat("   Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%\n",sep="")
  
  if (!is.null(x$method)&&!(x$method%in%c("PQL","lme.ML","lme.REML")))  
    cat(x$method," score = ",formatC(x$sp.criterion,digits=5),sep="")
 
  cat("  Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
  invisible(x)
}


anova.gam <- function (object, ..., dispersion = NULL, test = NULL, alpha=0, freq=FALSE)
# improved by Henric Nilsson
{   # adapted from anova.glm: R stats package
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named))
        warning("The following arguments to anova.glm(..) are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.glm <- unlist(lapply(dotargs, function(x) inherits(x,
        "glm")))
    dotargs <- dotargs[is.glm]
    if (length(dotargs) > 0)
        return(anova.glmlist(c(list(object), dotargs), dispersion = dispersion,
            test = test))
    if (!is.null(test)) warning("test argument ignored")
    if (!inherits(object,"gam")) stop("anova.gam called with non gam object")
    sg <- summary(object, dispersion = dispersion, freq = freq,alpha = alpha)
    class(sg) <- "anova.gam"
    sg
}


print.anova.gam <- function(x, digits = max(3, getOption("digits") - 3), ...)
{ # print method for class anova.gam resulting from single
  # gam model calls to anova. Improved by Henric Nilsson.
  print(x$family)
  cat("Formula:\n")
  print(x$formula)
  if (length(x$pTerms.pv)>0)
  { cat("\nParametric Terms:\n")
    printCoefmat(x$pTerms.table, digits = digits, signif.stars = FALSE, has.Pvalue = TRUE, na.print = "NA", ...)
  }
  cat("\n")
  if(x$m>0)
  { cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$s.table, digits = digits, signif.stars = FALSE, has.Pvalue = TRUE, na.print = "NA", ...)
  }
  invisible(x)
}

## End of improved anova and summary code. 



cooks.distance.gam <- function(model,...)
{ res <- residuals(model,type="pearson")
  dispersion <- model$sig2
  hat <- model$hat
  p <- sum(model$edf)
  (res/(1 - hat))^2 * hat/(dispersion * p)
}


sp.vcov <- function(x) {
## get cov matrix of smoothing parameters, if available
  if (!inherits(x,"gam")) stop("argument is not a gam object")
  if (x$method%in%c("ML","P-ML","REML","P-REML")&&!is.null(x$outer.info$hess)) {
    return(solve(x$outer.info$hess))
  } else return(NULL)
}

vcov.gam <- function(object, freq = FALSE, dispersion = NULL, ...)
## supplied by Henric Nilsson <henric.nilsson@statisticon.se> 
{ if (freq)
    vc <- object$Ve
  else vc <- object$Vp
  if (!is.null(dispersion))
    vc <- dispersion * vc / object$sig2
  name <- names(object$edf)
  dimnames(vc) <- list(name, name)
  vc
}




influence.gam <- function(model,...) { model$hat }




logLik.gam <- function (object, ...)
{  # based on logLik.glm - is ordering of p correction right???
    if (length(list(...)))
        warning("extra arguments discarded")
    fam <- family(object)$family
    p <- sum(object$edf)
    if (fam %in% c("gaussian", "Gamma", "inverse.gaussian"))
        p <- p + 1
    val <- p - object$aic/2
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
}




exclude.too.far<-function(g1,g2,d1,d2,dist)
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
  o<-.C(C_MinimumSeparation,as.double(g1),as.double(g2),as.integer(n),as.double(d1),as.double(d2),
         as.integer(m),distance=as.double(distance))  
  res<-rep(FALSE,n)
  res[o$distance > dist] <-TRUE
  res
}

strip.offset <- function(x)
# sole purpose is to take a model frame and rename any "offset(a.name)"
# columns "a.name"
{ na <- names(x)
  for (i in 1:length(na)) {
    if (substr(na[i],1,7)=="offset(") 
      na[i] <- substr(na[i],8,nchar(na[i])-1)
  }
  names(x) <- na
  x
}


vis.gam <- function(x,view=NULL,cond=list(),n.grid=30,too.far=0,col=NA,color="heat",
           contour.col=NULL,se=-1,type="link",plot.type="persp",zlim=NULL,nCol=50,...)
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

  #x$model <- strip.offset(x$model) 
  ## ... remove "offset(" and ")" from offset column name

  v.names <- row.names(attr(delete.response(x$terms),"factors"))

  if (is.null(view)) # get default view if none supplied
  { # v.names<-attr(attr(x$model,"terms"),"term.labels") # BUG... too many of these!!
   
    if (length(v.names)<2) stop("Model doesn't seem to have enough terms to do anything useful")
    view<-v.names[1:2]
  }
  if (!sum(view%in%names(x$model))) stop(
  paste(c("view variables must be one of",v.names),collapse=", "))
  if (length(unique(x$model[,view[1]]))<=1||length(unique(x$model[,view[2]]))<=1) 
  stop(paste("View variables must contain more than one value. view = c(",view[1],",",view[2],").",sep=""))

  # now get the values of the variables which are not the arguments of the plotted surface
  marg<-x$model[1,]
  m.name<-names(x$model)
  for (i in 1:length(marg))
  { ma<-cond[[m.name[i]]][1]
    if (is.null(ma)) 
    { if (is.factor(x$model[[i]]))
      marg[i]<-factor(levels(x$model[[i]])[1],levels(x$model[[i]]))
      else marg[i]<-mean(x$model[[i]]) 
    } else
    { if (is.factor(x$model[[i]]))
      marg[i]<-factor(ma,levels(x$model[[i]]))
      else marg[i]<-ma
    }
  }
  # marg includes conditioning values for view variables, but these will be ignored
  
  # Make dataframe....
  if (is.factor(x$model[,view[1]]))
  m1<-fac.seq(x$model[,view[1]],n.grid)
  else { r1<-range(x$model[,view[1]]);m1<-seq(r1[1],r1[2],length=n.grid)}
  if (is.factor(x$model[,view[2]]))
  m2<-fac.seq(x$model[,view[2]],n.grid)
  else {r2<-range(x$model[,view[2]]);m2<-seq(r2[1],r2[2],length=n.grid)}
  v1<-rep(m1,n.grid);v2<-rep(m2,rep(n.grid,n.grid))
  newd<-data.frame(v1=rep(marg[[1]],n.grid*n.grid))
  for (i in 2:dim(x$model)[2]) newd[[i]]<-rep(marg[[i]],n.grid*n.grid)
  row.names <- attr(newd,"row.names")
  attributes(newd) <- attributes(x$model) # done so that handling of offsets etc. works
  attr(newd,"row.names") <- row.names
  newd[[view[1]]]<-v1
  newd[[view[2]]]<-v2
  # call predict.gam to get predictions.....
  if (type=="link") zlab<-paste("linear predictor")
  else if (type=="response") zlab<-type
  else stop("type must be \"link\" or \"response\"")
  ## turn newd into a model frame, so that names and averages are valid
  attributes(newd)<-attributes(x$model)
  attr(newd,"row.names")<-as.character(1:(n.grid*n.grid))
  fv<-predict.gam(x,newdata=newd,se=TRUE,type=type)
  z<-fv$fit # store NA free copy now
  if (too.far>0) # exclude predictions too far from data
  { ex.tf<-exclude.too.far(v1,v2,x$model[,view[1]],x$model[,view[2]],dist=too.far)
    fv$se.fit[ex.tf]<-fv$fit[ex.tf]<-NA
  }
  # produce a continuous scale in place of any factors
  if (is.factor(m1)) 
  { m1<-as.numeric(m1);m1<-seq(min(m1)-0.5,max(m1)+0.5,length=n.grid) }
  if (is.factor(m2)) 
  { m2<-as.numeric(m2);m2<-seq(min(m1)-0.5,max(m2)+0.5,length=n.grid) }
  if (se<=0)
  { old.warn<-options(warn=-1)
    av<-matrix(c(0.5,0.5,rep(0,n.grid-1)),n.grid,n.grid-1)
    options(old.warn)
    # z is without any exclusion of gridpoints, so that averaging works nicely
    max.z <- max(z,na.rm=TRUE)
    z[is.na(z)] <- max.z*10000 # make sure NA's don't mess it up
    z<-matrix(z,n.grid,n.grid) # convert to matrix
    surf.col<-t(av)%*%z%*%av   # average over tiles  
    surf.col[surf.col>max.z*2] <- NA # restore NA's
    # use only non-NA data to set colour limits
    if (!is.null(zlim))
    { if (length(zlim)!=2||zlim[1]>=zlim[2]) stop("Something wrong with zlim")
      min.z<-zlim[1]
      max.z<-zlim[2]
    } else
    { min.z<-min(fv$fit,na.rm=TRUE)
      max.z<-max(fv$fit,na.rm=TRUE)
    }
    surf.col<-surf.col-min.z
    surf.col<-surf.col/(max.z-min.z)  
    surf.col<-round(surf.col*nCol)
    con.col <-1
    if (color=="heat") { pal<-heat.colors(nCol);con.col<-3;}
    else if (color=="topo") { pal<-topo.colors(nCol);con.col<-2;}
    else if (color=="cm") { pal<-cm.colors(nCol);con.col<-1;}
    else if (color=="terrain") { pal<-terrain.colors(nCol);con.col<-2;}
    else if (color=="gray"||color=="bw") {pal <- gray(seq(0.1,0.9,length=nCol));con.col<-1}
    else stop("color scheme not recognised")
    if (is.null(contour.col)) contour.col<-con.col   # default colour scheme
    surf.col[surf.col<1]<-1;surf.col[surf.col>nCol]<-nCol # otherwise NA tiles can get e.g. -ve index
    if (is.na(col)) col<-pal[as.array(surf.col)]
    z<-matrix(fv$fit,n.grid,n.grid)
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
                    ifelse("main" %in% dnm, "" , ",zlab=zlab"),",...)",sep="")
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
    { subs <- paste("grey are +/-",se,"s.e.") 
      lo.col <- "gray" ## ignore codetools claims about this
      hi.col <- "gray" ## ignore codetools 
    } else
    { subs<-paste("red/green are +/-",se,"s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim))
    { if (length(zlim)!=2||zlim[1]>=zlim[2]) stop("Something wrong with zlim")
      min.z<-zlim[1]
      max.z<-zlim[2]
    } else
    { z.max<-max(fv$fit+fv$se.fit*se,na.rm=TRUE)
      z.min<-min(fv$fit-fv$se.fit*se,na.rm=TRUE)
    }
    zlim<-c(z.min,z.max)
    z<-fv$fit-fv$se.fit*se;z<-matrix(z,n.grid,n.grid)
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
    z<-fv$fit;z<-matrix(z,n.grid,n.grid)

    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim",
                 ifelse("border" %in% dnm, "" ,",border=\"black\""),
                 stub,sep="")
    eval(parse(text=txt))

    par(new=TRUE) # don't clean device
    z<-fv$fit+se*fv$se.fit;z<-matrix(z,n.grid,n.grid)
    
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim",
                 ifelse("border" %in% dnm, "" ,",border=hi.col"),
                 stub,sep="")
    eval(parse(text=txt))

  }
}

# From here on is the code for magic.....

mroot <- function(A,rank=NULL,method="chol")
# finds the smallest square root of A, or the best approximate square root of 
# given rank. B is returned where BB'=A. A assumed non-negative definite. 
# Current methods "chol", "svd". "svd" is much slower, but much better at getting the 
# correct rank if it isn't known in advance. 
{ if (!all.equal(A,t(A))) stop("Supplied matrix not symmetric")
  if (method=="svd")
  { um<-La.svd(A)
    if (sum(um$d!=sort(um$d,decreasing=TRUE))>0) 
    stop("singular values not returned in order")
    if (is.null(rank)) # have to work out rank
    { rank<-dim(A)[1]
      if (um$d[1]<=0) rank <- 1 else
      while (rank>0&&(um$d[rank]/um$d[1]<.Machine$double.eps||
                           all.equal(um$u[,rank],um$vt[rank,])!=TRUE)) rank<-rank-1 
      if (rank==0) stop("Something wrong - matrix probably not +ve semi definite")    
    }
    d<-um$d[1:rank]^0.5
    return(t(t(um$u[,1:rank])*as.vector(d))) # note recycling rule used for efficiency
  } else
  if (method=="chol")
  { op<-options(warn=-1) ## don't want to be warned it's not +ve def
    L<-chol(A,pivot=TRUE)
    options(op) ## reset default warnings
    piv<-order(attr(L,"pivot"))
    if (is.null(rank)) rank<-attr(L,"rank")
    L<-L[,piv];L<-t(L[1:rank,])
    if (rank <= 1) dim(L) <- c(nrow(A),1)
    return(L)
  } else
  stop("method not recognised.")
}



magic.post.proc <- function(X,object,w=NULL)
# routine to take list returned by magic and extract:
# Vb the estimated bayesian parameter covariance matrix. rV%*%t(rV)*scale
# Ve the frequentist parameter estimator covariance matrix.
# edf the array of estimated degrees of freedom per parameter Vb%*%t(X)%*%W%*%X /scale
# hat the leading diagonal of the hat/influence matrix 
# NOTE: W=diag(w) if w non-matrix, otherwise w is a matrix square root. 
# flop count is O(nq^2) if X is n by q... this is why routine not part of magic
{ V<-object$rV%*%t(object$rV)
  if (!is.null(w)) 
  { if (is.matrix(w)) WX <- X <- w%*%X else 
    WX <- as.vector(w)*X # use recycling rule to form diag(w)%*%X cheaply 
    
  } else {WX <- X}
  M <- WX%*%V
  Ve <- (V%*%t(X))%*%M*object$scale # frequentist cov. matrix
  B <- X*M
  rm(M)
  hat <- apply(B,1,sum) # diag(X%*%V%*%t(WX))
  edf <- apply(B,2,sum) # diag(V%*%t(X)%*%WX)
  Vb <- V*object$scale;rm(V)
  list(Ve=Ve,Vb=Vb,hat=hat,edf=edf)
}

single.sp <- function(X,S,target=.5,tol=.Machine$double.eps*100)
## function to find smoothing parameter corresponding to particular 
## target e.d.f. for a single smoothing parameter problem. 
## X is model matrix; S is penalty matrix; target is target 
## average e.d.f. per penalized term.
{ R <- qr.R(qr(X)) ### BUG? pivoting?
  te <- try(RS <- backsolve(R,S,transpose=TRUE),silent=TRUE)
  if (inherits(te,"try-error")) return(-1)
  te <- try(RSR <- backsolve(R,t(RS),transpose=TRUE),silent=TRUE)
  if (inherits(te,"try-error")) return(-1)
  RSR <- (RSR+t(RSR))/2
  d <- eigen(RSR,symmetric=TRUE)$values
  d <- d[d>max(d)*tol]
  ff <- function(lambda,d,target) { 
    mean(1/(1+exp(lambda)*d))-target
  }
  lower <- 0
  while (ff(lower,d,target) <= 0) lower <- lower - 1
  upper <- lower
  while (ff(upper,d,target) > 0) upper <- upper + 1
  exp(uniroot(ff,c(lower,upper),d=d,target=target)$root)
}


initial.sp <- function(X,S,off,expensive=FALSE)
# Find initial smoothing parameter guesstimates based on model matrix X 
# and penalty list S. off[i] is the index of the first parameter to
# which S[[i]] applies, since S[[i]]'s only store non-zero submatrix of 
# penalty coefficient matrix.
{ n.p <- length(S) 
  def.sp <- array(0,n.p)
  if (n.p) { 
    ldxx <- colSums(X*X) # yields diag(t(X)%*%X)
    ldss <- ldxx*0       # storage for combined penalty l.d. 
    if (expensive) St <- matrix(0,ncol(X),ncol(X)) 
    pen <- rep(FALSE,length(ldxx)) # index of what actually gets penalized
    for (i in 1:n.p) { # loop over penalties
      maS <- max(abs(S[[i]])) 
      rsS <- rowMeans(abs(S[[i]]))
      csS <- colMeans(abs(S[[i]]))
      thresh <- .Machine$double.eps*maS*10
      ind <- rsS > thresh & csS > thresh # only these columns really penalize
      ss <- diag(S[[i]])[ind] # non-zero elements of l.d. S[[i]]
      start <- off[i];finish <- start+ncol(S[[i]])-1
      xx <- ldxx[start:finish]
      xx <- xx[ind]
      pen[start:finish] <- pen[start:finish]|ind
      sizeXX <- mean(xx)
      sizeS <- mean(ss)
      if (sizeS <= 0) stop(paste("S[[",i,"]] matrix is not +ve definite.",sep=""))
      def.sp[i] <- sizeXX/ sizeS # relative s.p. estimate
      ## accumulate leading diagonal of \sum sp[i]*S[[i]]
      ldss[start:finish] <- ldss[start:finish] + def.sp[i]*diag(S[[i]]) 
      
      if (expensive) St[start:finish,start:finish] <- 
                     St[start:finish,start:finish] + def.sp[i]*S[[i]]
    }
    if (expensive) { ## does full search for overall s.p.
      msp <- single.sp(X,St)           
      if (msp>0) def.sp <- def.sp*msp  
    } else {
      ind <- ldss>0&pen # base following only on penalized terms
      ldxx<-ldxx[ind];ldss<-ldss[ind]
      while (mean(ldxx/(ldxx+ldss))>.4) { def.sp <- def.sp*10;ldss <- ldss*10 }
      while (mean(ldxx/(ldxx+ldss))<.4) { def.sp <- def.sp/10;ldss <- ldss/10 }
    }
  } 
  as.numeric(def.sp)
}




magic <- function(y,X,sp,S,off,L=NULL,lsp0=NULL,rank=NULL,H=NULL,C=NULL,w=NULL,gamma=1,scale=1,gcv=TRUE,
                ridge.parameter=NULL,control=list(maxit=50,tol=1e-6,step.half=25,
                rank.tol=.Machine$double.eps^0.5),extra.rss=0,n.score=length(y))
# Wrapper for C routine magic. Deals with constraints weights and square roots of 
# penalties. 
# y is data vector, X is model matrix, sp is array of smoothing parameters,
# S is list of penalty matrices stored as smallest square submatrix excluding no 
# non-zero entries, off[i] is the location on the leading diagonal of the
# total penalty matrix of element (1,1) of S[[i]], rank is an array of penalty 
# ranks, L is a matrix mapping the log underlying smoothing parameters to the 
# smoothing parameters that actually multiply the penalties. i.e. the 
# log smoothing parameters are L%*%sp + lsp0
# H is any fixed penalty, C is a linear constraint matrix and w is the 
# weight vector. gamma is the dof inflation factor, scale is the scale parameter, only 
# used with UBRE, gcv TRUE means use GCV, if false, use UBRE.  
# Return list includes rV such that cov(b)=rV%*%t(rV)*scale and the leading diagonal
# of rV%*%t(rV)%*%t(X)%*%X gives the edf for each parameter.
# NOTE: W is assumed to be square root of inverse of covariance matrix. i.e. if
# W=diag(w) RSS is ||W(y-Xb||^2  
# If `ridge.parameter' is a positive number then then it is assumed to be the multiplier
# for a ridge penalty to be applied during fitting. 
# `extra.rss' is an additive constant by which the RSS is modified in the
#  GCV/UBRE or scale calculations, n.score is the `n' to use in the GCV/UBRE
#  score calcualtions (Useful for dealing with huge datasets).
{ n.p<-length(S)
  n.b<-dim(X)[2] # number of parameters
  # get initial estimates of smoothing parameters, using better method than is
  # built in to C code. This must be done before application of general 
  # constraints.
  if (n.p) def.sp <- initial.sp(X,S,off) else def.sp <- sp

  if (!is.null(L)) { ## have to estimate appropriate starting coefs
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=n.p||ncol(L)!=length(sp)) stop("L has inconsistent dimensions.")
    if (is.null(lsp0)) lsp0 <- rep(0,nrow(L))
    if (ncol(L)) def.sp <- exp(as.numeric(coef(lm(log(def.sp)~L-1+offset(lsp0)))))
  }

  # get square roots of penalties using supplied ranks or estimated 
  if (n.p>0)
  { for (i in 1:n.p) 
    { if (is.null(rank)) B<-mroot(S[[i]],method="svd") 
      else B<-mroot(S[[i]],rank=rank[i],method="chol")
      m<-dim(B)[2]
      R<-matrix(0,n.b,m)
      R[off[i]:(off[i]+dim(B)[1]-1),]<-B
      S[[i]]<-R
    }
    rm(B);rm(R)
  }
  # if there are constraints then need to form null space of constraints Z 
  # (from final columns of Q, from QR=C'). Then form XZ and Z'S_i^0.5 for all i 
  # and Z'HZ.
  # On return from mgcv2 set parameters to Zb (apply Q to [0,b']').   
  ##Xo<-X
  if (!is.null(C)) # then impose constraints 
   { n.con<-dim(C)[1]
    ns.qr<-qr(t(C)) # last n.b-n.con columns of Q are the null space of C
    X<-t(qr.qty(ns.qr,t(X)))[,(n.con+1):n.b] # last n.b-n.con cols of XQ (=(Q'X')')
    # need to work through penalties forming Z'S_i^0.5 's
    if (n.p>0) for (i in 1:n.p) { 
      S[[i]]<-qr.qty(ns.qr,S[[i]])[(n.con+1):n.b,,drop=FALSE]
      ## following essential given assumptions of the C code...
      if (ncol(S[[i]])>nrow(S[[i]])) { ## no longer have a min col square root.
        S[[i]] <- t(qr.R(qr(t(S[[i]])))) ## better!
      }
    }
    # and Z'HZ too
    if (!is.null(H))
    { H<-qr.qty(ns.qr,H)[(n.con+1):n.b,] # Z'H
      H<-t(qr.qty(ns.qr,t(H))[(n.con+1):n.b,]) # Z'HZ = (Z'[Z'H]')' 
    }
    full.rank=n.b-n.con
  } else full.rank=n.b
  # now deal with weights....
  if (!is.null(w))
  { if (is.matrix(w))
    { if (dim(w)[1]!=dim(w)[2]||dim(w)[2]!=dim(X)[1]) stop("dimensions of supplied w wrong.")
      y<-w%*%y
      X<-w%*%X
    } else
    { if (length(y)!=length(w)) stop("w different length from y!")
      y<-y*w
      X<-as.vector(w)*X # use recycling rule to form diag(w)%*%X cheaply
    }
  }
  if (is.null(dim(X))) # lost dimensions as result of being single columned! 
  { n<-length(y)
    if (n!=length(X)) stop("X lost dimensions in magic!!")
    dim(X)<-c(n,1)
  }
  # call real mgcv engine...
  Si<-array(0,0);cS<-0
  if (n.p>0) for (i in 1:n.p) 
  { Si<-c(Si,S[[i]]);
    cS[i]<-dim(S[[i]])[2]
  }
  icontrol<-as.integer(gcv);icontrol[2]<-length(y);q<-icontrol[3]<-dim(X)[2];
  if (!is.null(ridge.parameter)&&ridge.parameter>0)
  { if(is.null(H)) H<-diag(ridge.parameter,q) else H<-H+diag(ridge.parameter,q)}
  icontrol[4]<-as.integer(!is.null(H));icontrol[5]<- n.p;icontrol[6]<-control$step.half
  if (is.null(L)) { icontrol[7] <- -1;L <- diag(n.p) } else icontrol[7]<-ncol(L)
  if (is.null(lsp0)) lsp0 <- rep(0,nrow(L))
 
  b<-array(0,icontrol[3])
  # argument names in call refer to returned values.
  um<-.C(C_magic,as.double(y),as.double(X),sp=as.double(sp),as.double(def.sp),as.double(Si),as.double(H),as.double(L),
          lsp0=as.double(lsp0),score=as.double(gamma),scale=as.double(scale),info=as.integer(icontrol),as.integer(cS),
          as.double(control$rank.tol),rms.grad=as.double(control$tol),b=as.double(b),rV=double(q*q),
          as.double(extra.rss),as.integer(n.score))
  res<-list(b=um$b,scale=um$scale,score=um$score,sp=um$sp,sp.full=as.numeric(exp(L%*%log(um$sp))))
  res$rV<-matrix(um$rV[1:(um$info[1]*q)],q,um$info[1])
  gcv.info<-list(full.rank=full.rank,rank=um$info[1],fully.converged=as.logical(um$info[2]),
      hess.pos.def=as.logical(um$info[3]),iter=um$info[4],score.calls=um$info[5],rms.grad=um$rms.grad)
  res$gcv.info<-gcv.info
  if (!is.null(C)) # need image of constrained parameter vector in full space
  { b<-c(rep(0,n.con),res$b)
    res$b<-qr.qy(ns.qr,b) # Zb 
    b<-matrix(0,n.b,dim(res$rV)[2])
    b[(n.con+1):n.b,]<-res$rV 
    res$rV<-qr.qy(ns.qr,b)# ZrV
  } 
  res
}


print.mgcv.version <- function()
{ library(help=mgcv)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  cat(paste("This is mgcv ",version,". For overview type `help(\"mgcv-package\")'.\n"),sep="")
}

set.mgcv.options <- function()
## function used to set optional value used in notLog
## and notExp...
{ runif(1) ## ensure there is a seed 
  options(mgcv.vc.logrange=25)
}

.onAttach <- function(...) { 
  print.mgcv.version()
  set.mgcv.options()
}

.onUnload <- function(libpath) library.dynam.unload("mgcv", libpath)

.First.lib <- function(lib, pkg) {
  library.dynam("mgcv", pkg, lib)
  print.mgcv.version()
  set.mgcv.options()
}


###############################################################################
### ISSUES.....






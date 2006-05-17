##  R routines for the package mgcv (c) Simon Wood 2000-2006
##  With contributions from Henric Nilsson

mono.con<-function(x,up=TRUE,lower=NA,upper=NA)
# Takes the knot sequence x for a cubic regression spline and returns a list with 
# 2 elements matrix A and array b, such that if p is the vector of coeffs of the
# spline, then Ap>b ensures monotonicity of the spline.
# up=TRUE gives monotonic increase, up=FALSE gives decrease.
# lower and upper are the optional lower and upper bounds on the spline.
{
  if (is.na(lower)) {lo<-0;lower<-0;} else lo<-1
  if (is.na(upper)) {hi<-0;upper<-0;} else hi<-1
  if (up) inc<-1 else inc<-0
  control<-4*inc+2*lo+hi
  n<-length(x)
  if (n<4) stop("At least three knots required in call to mono.con.")
  A<-matrix(0,4*(n-1)+lo+hi,n)
  b<-array(0,4*(n-1)+lo+hi)
  if (lo*hi==1&&lower>=upper) stop("lower bound >= upper bound in call to mono.con()")
  oo<-.C(C_RMonoCon,as.double(A),as.double(b),as.double(x),as.integer(control),as.double(lower),
         as.double(upper),as.integer(n))
  A<-matrix(oo[[1]],dim(A)[1],dim(A)[2])
  b<-array(oo[[2]],dim(A)[1])
  list(A=A,b=b)
}  


uniquecombs<-function(x) {
# takes matrix x and counts up unique rows
if (is.null(x)) stop("x is null")
if (is.null(nrow(x))) stop("x has no row attribute")
if (is.null(ncol(x))) stop("x has no col attribute")
res<-.C(C_RuniqueCombs,as.double(x),as.integer(nrow(x)),as.integer(ncol(x)))
n<-res[[2]]*res[[3]]
x<-matrix(res[[1]][1:n],res[[2]],res[[3]])
x
}


null.space.dimension<-function(d,m)
# vectorized function for calculating null space dimension for penalties of order m
# for dimension d data M=(m+d+1)!/(d!(m-d)!). Any m not satisfying 2m>d is reset so 
# that 2m>d+1 (assuring "visual" smoothness) 
{ if (sum(d<0)) stop("d can not be negative in call to null.space.dimension().")
  ind<-2*m<d+1
  if (sum(ind)) # then default m required for some elements
  { m[ind]<-1;ind<-2*m<d+2
    while (sum(ind)) { m[ind]<-m[ind]+1;ind<-2+m<d+2;}
  }
  M<-m*0+1;ind<-M==1;i<-0
  while(sum(ind))
  { M[ind]<-M[ind]*(d[ind]+m[ind]-1-i);i<-i+1;ind<-i<d
  }
  ind<-d>1;i<-2
  while(sum(ind))
  { M[ind]<-M[ind]/i;ind<-d>i;i<-i+1   
  }
  M
}







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

  o<-.C(C_RPCLS,as.double(M$X),as.double(M$p),as.double(M$y),as.double(M$w),as.double(M$Ain),as.double(M$bin)
        ,as.double(M$C),as.double(H),as.double(Sa),as.integer(M$off),as.integer(df),as.double(M$sp),
        as.integer(length(M$off)),as.integer(nar))
  array(o[[2]],length(M$p))
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

mgcv<-function(y,X,sp,S,off,C=NULL,w=rep(1,length(y)),H=NULL,scale=1,gcv=TRUE,control=mgcv.control())

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
  
  ret<-list(b=p,scale=scale,score=gcv.ubre,sp=sp,Vb=Vp,hat=hat,edf=edf,info=conv)
 
}

tensor.prod.model.matrix<-function(X)
# X is a list of model matrices, from which a tensor product model matrix is to be produced.
# e.g. ith row is basically X[[1]][i,]%x%X[[2]][i,]%x%X[[3]][i,], but this routine works 
# column-wise, for efficiency
{ m<-length(X)
  X1<-X[[m]]
  n<-nrow(X1)
  if (m>1) for (i in (m-1):1)
  { X0<-X1;X1<-matrix(0,n,0)
    for (j in 1:ncol(X[[i]]))
    X1<-cbind(X1,X[[i]][,j]*X0)
  }
  X1
}

tensor.prod.penalties <- function(S)
# Given a list S of penalty matrices for the marginal bases of a tensor product smoother
# this routine produces the resulting penalties for the tensor product basis. 
# e.g. if S_1, S_2 and S_3 are marginal penalties and I_1, I_2, I_3 are identity matrices 
# of the same dimensions then the tensor product penalties are:
#   S_1 %x% I_2 %x% I_3, I_1 %x% S_2 %x% I_3 and I_1 %*% I_2 %*% S_3
# Note that the penalty list must be in the same order as the model matrix list supplied
# to tensor.prod.model() when using these together.
{ m<-length(S)
  I<-list(); for (i in 1:m) { 
    n<-ncol(S[[i]])
    I[[i]]<-diag(n)
  #  I[[i]][1,1] <- I[[i]][n,n]<-.5 
  }
  TS<-list()
  if (m==1) TS[[1]]<-S[[1]] else
  for (i in 1:m)
  { if (i==1) M0<-S[[1]] else M0<-I[[1]]
    for (j in 2:m)
    { if (i==j) M1<-S[[i]] else M1<-I[[j]] 
      M0<-M0%x%M1
    }
    TS[[i]]<- (M0+t(M0))/2 # ensure exactly symmetric 
  }
  TS
}


get.var<-function(txt,data)
# txt contains text that may be a variable name and may be an expression 
# for creating a variable. get.var first tries data[[txt]] and if that 
# fails tries evaluating txt within data (only). Routine returns NULL
# on failure, or if result is not numeric or a factor.
{ x <- data[[txt]]
  if (is.null(x)) 
  { x <- try(eval(parse(text=txt),data,enclos=NULL),silent=TRUE)
    if (inherits(x,"try-error")) x <- NULL
  }
  if (!is.numeric(x)&&!is.factor(x)) x <- NULL  
  x
}



te <- function(..., k=NA,bs="cr",m=0,d=NA,by=NA,fx=FALSE,mp=TRUE,np=TRUE)
# function for use in gam formulae to specify a tensor product smooth term.
# e.g. te(x0,x1,x2,k=c(5,4,4),bs=c("tp","cr","cr"),m=c(1,1,2),by=x3) specifies a rank 80 tensor  
# product spline. The first basis is rank 5, t.p.r.s. basis penalty order 1, and the next 2 bases
# are rank 4 cubic regression splines with m ignored.  
# k, bs,m,d and fx can be supplied as single numbers or arrays with an element for each basis.
# Returns a list consisting of:
# * margin - a list of smooth.spec objects specifying the marginal bases
# * term   - array of covariate names
# * by     - the by variable name
# * fx     - array indicating which margins should be treated as fixed (i.e unpenalized).
# * full.call - fully exapanded call for this term
# * label  - label for this term
# * mp - TRUE to use a penalty per dimension, FALSE to use a single penalty
{ vars<-as.list(substitute(list(...)))[-1] # gets terms to be smoothed without evaluation
  dim<-length(vars) # dimension of smoother
  by.var<-deparse(substitute(by),backtick=TRUE) #getting the name of the by variable
  term<-deparse(vars[[1]],backtick=TRUE) # first covariate
  if (dim>1) # then deal with further covariates
  for (i in 2:dim)
  { term[i]<-deparse(vars[[i]],backtick=TRUE)
  }
  for (i in 1:dim) term[i] <- attr(terms(reformulate(term[i])),"term.labels")
  # term now contains the names of the covariates for this model term
  
  # check d - the number of covariates per basis
  if (sum(is.na(d))||is.null(d)) { n.bases<-dim;d<-rep(1,dim)} # one basis for each dimension
  else  # array d supplied, the dimension of each term in the tensor product 
  { d<-round(d)
    ok<-TRUE
    if (sum(d<=0)) ok<-FALSE 
    if (sum(d)!=dim) ok<-FALSE
    if (ok)
    n.bases<-length(d)
    else 
    { warning("something wrong with argument d.")
      n.bases<-dim;d<-rep(1,dim)
    }     
  }
  
  # now evaluate k 
  if (sum(is.na(k))||is.null(k)) k<-5^d 
  else 
  { k<-round(k);ok<-TRUE
    if (sum(k<3)) { ok<-FALSE;warning("one or more supplied k too small - reset to default")}
    if (length(k)==1&&ok) k<-rep(k,n.bases)
    else if (length(k)!=n.bases) ok<-FALSE
    if (!ok) k<-5^d 
  }
  # evaluate fx
  if (sum(is.na(fx))||is.null(fx)) fx<-rep(FALSE,n.bases)
  else if (length(fx)==1) fx<-rep(fx,n.bases)
  else if (length(fx)!=n.bases)
  { warning("dimension of fx is wrong") 
    fx<-rep(FALSE,n.bases)
  }

  # now check the basis types
  if (length(bs)==1) bs<-rep(bs,n.bases)
  if (length(bs)!=n.bases) {warning("bs wrong length and ignored.");bs<-rep("cr",n.bases)}
  bs[d>1&bs!="tp"&bs!="ts"]<-"tp"
  # finally the penalty orders
  if (length(m)==1) m<-rep(m,n.bases)
  if (length(m)!=n.bases) 
  { warning("m wrong length and ignored.");m<-rep(0,n.bases)}
  m[m<0]<-0
  # check for repeated variables in function argument list
  if (length(unique(term))!=dim) stop("Repeated variables as arguments of a smooth are not permitted")
  # Now construct smooth.spec objects for the margins
  j<-1 # counter for terms
  margin<-list()
  for (i in 1:n.bases)
  { j1<-j+d[i]-1
    stxt<-"s("
    for (l in j:j1) stxt<-paste(stxt,term[l],",",sep="")
    stxt<-paste(stxt,"k=",deparse(k[i],backtick=TRUE),",bs=",deparse(bs[i],backtick=TRUE),
                ",m=",deparse(m[i],backtick=TRUE),")")
    margin[[i]]<- eval(parse(text=stxt))  # NOTE: fx and by not dealt with here!
    j<-j1+1
  }
  # assemble version of call with all options expanded as text
  if (mp) mp <- TRUE else mp <- FALSE
  if (np) np <- TRUE else np <- FALSE
  full.call<-paste("te(",term[1],sep="")
  if (dim>1) for (i in 2:dim) full.call<-paste(full.call,",",term[i],sep="")
  label<-paste(full.call,")",sep="")   # label for parameters of this term
  full.call<-paste(full.call,",k=",deparse(k,backtick=TRUE),",bs=",deparse(bs,backtick=TRUE),
                   ",m=",deparse(m,backtick=TRUE),",d=",deparse(d,backtick=TRUE),
                   ",by=",by.var,",fx=",deparse(fx,backtick=TRUE),",mp=",deparse(mp,backtick=TRUE),
                   ",np=",deparse(np,backtick=TRUE),")",sep="")
  ret<-list(margin=margin,term=term,by=by.var,fx=fx,full.call=full.call,label=label,dim=dim,mp=mp,np=np)
  class(ret)<-"tensor.smooth.spec"
  ret
}





s <- function (..., k=-1,fx=FALSE,bs="tp",m=0,by=NA)
# function for use in gam formulae to specify smooth term, e.g. s(x0,x1,x2,k=40,m=3,by=x3) specifies 
# a rank 40 thin plate regression spline of x0,x1 and x2 with a third order penalty, to be multiplied by
# covariate x3, when it enters the model.
# Returns a list of consisting of the names of the covariates, and the name of any by variable
# of a model formula term representing the smooth, the basis dimension, the type of basis
# (1 -t.p.r.s.; 0 - cubic), whether it is fixed or penalized and the order of the penalty (0 for auto).
# backwards compatibility with versions < 0.6 is maintained so that terms like: s(x,6|f)
# still work correctly.
# The returned full.call is a string with the call made fully explicit, so that when pasted
# into a formula, that formula can be parsed without reference to any variables that may
# have been used to define k,fx,bs or m.
{ vars<-as.list(substitute(list(...)))[-1] # gets terms to be smoothed without evaluation
 # call<-match.call() # get function call
  d<-length(vars) # dimension of smoother
  term<-deparse(vars[[d]],backtick=TRUE) # last term in the ... arguments
  by.var<-deparse(substitute(by),backtick=TRUE) #getting the name of the by variable
  if (by.var==".") stop("by=. not allowed")
  term<-deparse(vars[[1]],backtick=TRUE) # first covariate
  if (term[1]==".") stop("s(.) not yet supported.")
  if (d>1) # then deal with further covariates
  for (i in 2:d)
  { term[i]<-deparse(vars[[i]],backtick=TRUE)
    if (term[i]==".") stop("s(.) not yet supported.")
  }
  for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])),"term.labels")
  # term now contains the names of the covariates for this model term
  # now evaluate all the other 
  k.new<-round(k) # in case user has supplied non-integer basis dimension
  if (!all.equal(k.new,k)) {warning("argument k of s() should be integer and has been rounded")}
  k<-k.new
  if (k==-1) k<-10*3^(d-1) # auto-initialize basis dimension
  if (k < 2) 
  { k <- 2
    warning("meaninglessly low k; reset to 2\n")
  }
  if (bs=="cr") # set basis types
  { if (d>1) { warning("cr basis only works with 1-d smooths!\n");bs<-"tp";}
  } 
  m[m<0]<-0
  # check for repeated variables in function argument list
  if (length(unique(term))!=d) stop("Repeated variables as arguments of a smooth are not permitted")
  # assemble version of call with all options expanded as text
  full.call<-paste("s(",term[1],sep="")
  if (d>1) for (i in 2:d) full.call<-paste(full.call,",",term[i],sep="")
  label<-paste(full.call,")",sep="") # used for labelling parameters
  full.call<-paste(full.call,",k=",deparse(k,backtick=TRUE),",fx=",deparse(fx,backtick=TRUE),",bs=",
             deparse(bs,backtick=TRUE),",m=",deparse(m,backtick=TRUE),
                   ",by=",by.var,")",sep="")
  ret<-list(term=term,bs.dim=k,fixed=fx,dim=d,p.order=m,by=by.var,full.call=full.call,label=label)
  class(ret)<-paste(bs,".smooth.spec",sep="")
  ret
}

smooth.construct.tensor.smooth.spec<-function(object,data,knots)
# the constructor for a tensor product basis object
{ m<-length(object$margin)  # number of marginal bases
  Xm<-list();Sm<-list();nr<-r<-d<-array(0,m)
  for (i in 1:m)
  { object$margin[[i]]<-smooth.construct(object$margin[[i]],data,knots)
    Xm[[i]]<-object$margin[[i]]$X
    Sm[[i]]<-object$margin[[i]]$S[[1]]
    d[i]<-nrow(Sm[[i]])
    r[i]<-object$margin[[i]]$rank
    nr[i]<-object$margin[[i]]$null.space.dim
  }
  XP <- list()
  if (object$np) # reparameterize 
  for (i in 1:m)
  { if (object$margin[[i]]$dim==1) {
      if (!inherits(object$margin[[i]],c("cs.smooth","cr.smooth"))) { # these classes already optimal
        x <- get.var(object$margin[[i]]$term,data)
        np <- ncol(object$margin[[i]]$X) ## number of params
        ## note: to avoid extrapolating wiggliness measure
        ## must include extremes as eval points
        knt <- quantile(unique(x),(0:(np-1))/(np-1)) ## evaluation points
#        knt <- seq(min(x),max(x),length=np) 
        pd <- data.frame(knt)
        names(pd) <- object$margin[[i]]$term
        XP[[i]] <- solve(Predict.matrix(object$margin[[i]],pd),tol=0)
        Xm[[i]] <- Xm[[i]]%*%XP[[i]]
        Sm[[i]] <- t(XP[[i]])%*%Sm[[i]]%*%XP[[i]]
      } else XP[[i]]<-NULL
    } else XP[[i]]<-NULL
  }
  # scale `nicely' - mostly to avoid problems with lme ...
  for (i in 1:m)  Sm[[i]] <- Sm[[i]]/svd(Sm[[i]])$d[1] 
  max.rank<-prod(d)
  r<-max.rank*r/d # penalty ranks
  X<-tensor.prod.model.matrix(Xm)
  if (object$mp) # multiple penalties
  { S<-tensor.prod.penalties(Sm)
    for (i in m:1) if (object$fx[i]) S[[i]]<-NULL # remove penalties for un-penalized margins
  } else # single penalty
  { S<-Sm[[1]];r<-object$margin[[i]]$rank
    if (m>1) for (i in 2:m) 
    { S<-S%x%Sm[[i]]
      r<-r*object$margin[[i]]$rank
    } 
    if (sum(object$fx)==m) 
    { S <- list();object$fixed=TRUE } else
    { S<-list(S);object$fixed=FALSE }
    nr <- max.rank-r
    object$bs.dim<-max.rank
  }
  C<-matrix(colSums(X),1,ncol(X))
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    X<-as.numeric(by)*X
  }
  object$X<-X;object$S<-S;object$C<-C
  object$df<-ncol(X)-1
  object$null.space.dim <- prod(nr) # penalty null space rank 
  object$rank<-r
  object$XP <- XP
  class(object)<-"tensor.smooth"
  object
}

 


smooth.construct.tp.smooth.spec<-function(object,data,knots)
# The constructor for a t.p.r.s. basis object.
{ shrink <- attr(object,"shrink")
  x<-array(0,0)
  shift<-array(0,object$dim)
  for (i in 1:object$dim) 
  { xx <- get.var(object$term[[i]],data)
    shift[i]<-mean(xx)  # centre covariates
    xx <- xx - shift[i]
    x<-c(x,xx)
  }
  if (is.null(knots)) {knt<-0;nk<-0}
  else 
  { knt<-array(0,0)
    for (i in 1:object$dim) 
    { dum <- get.var(object$term[[i]],knots)-shift[i]
      if (is.null(dum)) {knt<-0;nk<-0;break} # no valid knots for this term
      knt <- c(knt,dum)
      nk0 <- length(dum)
      if (i > 1 && nk != nk0) 
      stop("components of knots relating to a single smooth must be of same length")
      nk <- nk0
    }
  }
  n<-nrow(data)
  if (nk>n) { nk <- 0
  warning("more knots than data in a tp term: knots ignored.")}
  mtk <- attr(object,"max.tprs.knots")
  if (is.null(mtk)) mtk <- n 
  k<-object$bs.dim 
  M<-null.space.dimension(object$dim,object$p.order) 
  if (k<M+1) # essential or construct_tprs will segfault, as tprs_setup does this
  { k<-M+1
    object$bs.dim<-k
    warning("basis dimension, k, increased to minimum possible\n")
  }
  X<-array(0,n*k)
  S<-array(0,k*k)
 
  UZ<-array(0,(n+M)*k)
  Xu<-x
  C<-array(0,k)
  nXu<-0  
  oo<-.C(C_construct_tprs,as.double(x),as.integer(object$dim),as.integer(n),as.double(knt),as.integer(nk),
               as.integer(object$p.order),as.integer(object$bs.dim),X=as.double(X),S=as.double(S),
               UZ=as.double(UZ),Xu=as.double(Xu),n.Xu=as.integer(nXu),C=as.double(C),
               max.knots=as.integer(mtk))
  object$X<-matrix(oo$X,n,k)                   # model matrix
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    object$X<-as.numeric(by)*object$X
  }
  object$S<-list()
  if (!object$fixed) 
  { object$S[[1]]<-matrix(oo$S,k,k)         # penalty matrix
    object$S[[1]]<-(object$S[[1]]+t(object$S[[1]]))/2 # ensure exact symmetry
    if (!is.null(shrink)) # then add shrinkage term to penalty 
    { norm <- mean(object$S[[1]]^2)^0.5
      object$S[[1]] <- object$S[[1]] + diag(k)*norm*abs(shrink)
    }
  }
  UZ.len <- (oo$n.Xu+M)*k
  object$UZ<-matrix(oo$UZ[1:UZ.len],oo$n.Xu+M,k)         # truncated basis matrix
  Xu.len <- oo$n.Xu*object$dim
  object$Xu<-matrix(oo$Xu[1:Xu.len],oo$n.Xu,object$dim)  # unique covariate combinations
  object$C<-matrix(oo$C,1,k)                   # constraints
  object$df<-object$bs.dim-1                   # DoF given constraint
  object$shift<-shift                          # covariate shifts
  if (is.null(shrink)) { 
    object$rank <- k-M 
  } else object$rank <- k                             # penalty rank
  object$null.space.dim<-M

  class(object)<-"tprs.smooth"
  object
}

smooth.construct.ts.smooth.spec<-function(object,data,knots)
# implements a class of tprs like smooths with an additional shrinkage
# term in the penalty... this allows for fully integrated GCV model selection
{ attr(object,"shrink") <- 1e-4
  object <- smooth.construct.tp.smooth.spec(object,data,knots)
  class(object) <- "ts.smooth"
  object
}

smooth.construct.cr.smooth.spec<-function(object,data,knots)
# this routine is the constructor for cubic regression spline basis objects
# It takes a cubic regression spline specification object and returns the 
# corresponding basis object.
{ shrink <- attr(object,"shrink")
  x <- get.var(object$term,data)
  nx<-length(x)
  if (is.null(knots)) ok <- FALSE
  else 
  { k <- get.var(object$term,knots)
    if (is.null(k)) ok <- FALSE
    else ok<-TRUE
  }
  nk <- object$bs.dim
  if (!ok) { k <- rep(0,nk);k[2]<- -1}
  
  if (length(k)!=nk) stop("number of supplied knots != k for a cr smooth")

  X <- rep(0,nx*nk);S<-rep(0,nk*nk);C<-rep(0,nk);control<-0
  
  if (length(unique(x))<nk) 
  { msg <- paste(object$term," has insufficient unique values to support ",
                 nk," knots: reduce k.",sep="")
    stop(msg)
  }

  oo <- .C(C_construct_cr,as.double(x),as.integer(nx),as.double(k),
           as.integer(nk),as.double(X),as.double(S),
           as.double(C),as.integer(control))

  object$X <- matrix(oo[[5]],nx,nk)
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    object$X <- as.numeric(by)*object$X
  }
  object$S<-list()     # only return penalty if term not fixed
  if (!object$fixed) 
  { object$S[[1]] <- matrix(oo[[6]],nk,nk)
    object$S[[1]]<-(object$S[[1]]+t(object$S[[1]]))/2 # ensure exact symmetry
    if (!is.null(shrink)) # then add shrinkage term to penalty 
    { norm <- mean(object$S[[1]]^2)^0.5
      object$S[[1]] <- object$S[[1]] + diag(nk)*norm*abs(shrink)
    }
  }
  if (is.null(shrink)) { 
  object$rank<-nk-2 
  } else object$rank <- nk   # penalty rank
  object$C <- matrix(oo[[7]],1,nk)  # constraint
  object$df<-object$bs.dim-1 # degrees of freedom, given constraint
  object$null.space.dim <- 2
  object$xp <- oo[[3]]  # knot positions 
  class(object) <- "cr.smooth"
  object
}

smooth.construct.cs.smooth.spec<-function(object,data,knots)
# implements a class of cr like smooths with an additional shrinkage
# term in the penalty... this allows for fully integrated GCV model selection
{ attr(object,"shrink") <- 1e-4
  object <- smooth.construct.cr.smooth.spec(object,data,knots)
  class(object) <- "cs.smooth"
  object
}

place.knots<-function(x,nk)
# knot placement code. x is a covariate array, nk is the number of knots,
# and this routine spaces nk knots evenly throughout the x values, with the 
# endpoints at the extremes of the data.
{ x<-sort(unique(x));n<-length(x)
  if (nk>n) stop("more knots than unique data values is not allowed")
  if (nk<2) stop("too few knots")
  if (nk==2) return(range(x))
  delta<-(n-1)/(nk-1) # how many data steps per knot
  lbi<-floor(delta*1:(nk-2))+1 # lower interval bound index
  frac<-delta*1:(nk-2)+1-lbi # left over proportion of interval  
  x.shift<-x[-1]
  knot<-array(0,nk)
  knot[nk]<-x[n];knot[1]<-x[1]
  knot[2:(nk-1)]<-x[lbi]*(1-frac)+x.shift[lbi]*frac
  knot
}

smooth.construct.cc.smooth.spec<-function(object,data,knots)
# constructor function for cyclic cubic splines
{ getBD<-function(x)
  # matrices B and D in expression Bm=Dp where m are s"(x_i) and 
  # p are s(x_i) and the x_i are knots of periodic spline s(x)
  # B and D slightly modified (for periodicity) from Lancaster 
  # and Salkauskas (1986) Curve and Surface Fitting section 4.7.
  { n<-length(x)
    h<-x[2:n]-x[1:(n-1)]
    n<-n-1
    D<-B<-matrix(0,n,n)
    B[1,1]<-(h[n]+h[1])/3;B[1,2]<-h[1]/6;B[1,n]<-h[n]/6
    D[1,1]<- -(1/h[1]+1/h[n]);D[1,2]<-1/h[1];D[1,n]<-1/h[n]
    for (i in 2:(n-1))
    { B[i,i-1]<-h[i-1]/6
      B[i,i]<-(h[i-1]+h[i])/3
      B[i,i+1]<-h[i]/6
      D[i,i-1]<-1/h[i-1]
      D[i,i]<- -(1/h[i-1]+1/h[i])
      D[i,i+1]<- 1/h[i]
    }
    B[n,n-1]<-h[n-1]/6;B[n,n]<-(h[n-1]+h[n])/3;B[n,1]<-h[n]/6
    D[n,n-1]<-1/h[n-1];D[n,n]<- -(1/h[n-1]+1/h[n]);D[n,1]<-1/h[n]
    list(B=B,D=D)
  } # end of getBD local function
  # evaluate covariate, x, and knots, k.
  x <- get.var(object$term,data)
  nx<-length(x)
  nk <- object$bs.dim
  if (!is.null(knots))  k <- get.var(object$term,knots)
  else k<-NULL
  if (is.null(k)) k<-place.knots(x,nk)   

  if (length(k)!=nk) stop("number of supplied knots != k for a cc smooth")

  um<-getBD(k)
  BD<-solve(um$B,um$D) # s"(k)=BD%*%s(k) where k are knots minus last knot
  if (!object$fixed)
  { object$S<-list(t(um$D)%*%BD)      # the penalty
    object$S[[1]]<-(object$S[[1]]+t(object$S[[1]]))/2 # ensure exact symmetry
  }
  object$BD<-BD # needed for prediction
  object$xp<-k  # needed for prediction   
  X<-Predict.matrix.cyclic.smooth(object,data) 
  C<-matrix(colSums(X),1,ncol(X))
  object$X<-X
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    object$X<-as.numeric(by)*object$X
  }
  object$C<-C
  object$rank<-ncol(X)-1  # rank of smoother matrix
  object$df<-object$bs.dim-2 # degrees of freedom, accounting for centring and cycling
  object$null.space.dim <- 1  
  class(object)<-"cyclic.smooth"
  object
}

Predict.matrix.tensor.smooth<-function(object,data)
# the prediction method for a tensor product smooth
{ m<-length(object$margin)
  X<-list()
  for (i in 1:m) X[[i]]<-Predict.matrix(object$margin[[i]],data)
  mxp <- length(object$XP)
  if (mxp>0) 
  for (i in 1:mxp) if (!is.null(object$XP[[i]])) X[[i]] <- X[[i]]%*%object$XP[[i]]
  T <- tensor.prod.model.matrix(X)
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    T <- as.numeric(by)*T
  }
  T
}

Predict.matrix.cyclic.smooth<-function(object,data)
# this is the prediction method for a cyclic cubic regression spline
{ pred.mat<-function(x,knots,BD)
  # BD is B^{-1}D. Basis as given in Lancaster and Salkauskas (1986)
  # Curve and Surface fitting, but wrapped to give periodic smooth.
  { j<-x
    n<-length(knots)
    h<-knots[2:n]-knots[1:(n-1)]
    if (max(x)>max(knots)||min(x)<min(knots)) 
    stop("can't predict outside range of knots with periodic smoother")
    for (i in n:2) j[x<=knots[i]]<-i
    j1<-hj<-j-1
    j[j==n]<-1
    I<-diag(n-1)
    X<-BD[j1,]*as.numeric(knots[j1+1]-x)^3/as.numeric(6*h[hj])+
       BD[j,]*as.numeric(x-knots[j1])^3/as.numeric(6*h[hj])-
       BD[j1,]*as.numeric(h[hj]*(knots[j1+1]-x)/6)-
       BD[j,]*as.numeric(h[hj]*(x-knots[j1])/6) +
       I[j1,]*as.numeric((knots[j1+1]-x)/h[hj]) +
       I[j,]*as.numeric((x-knots[j1])/h[hj])
    X
  }
  x <- get.var(object$term,data)
  if (length(x)<1) stop("no data to predict at")
  X <- pred.mat(x,object$xp,object$BD)
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    X<-as.numeric(by)*X
  }
  X
}


Predict.matrix.cr.smooth<-function(object,data)
# this is the prediction method for a cubic regression spline
{ x <- get.var(object$term,data)
  if (length(x)<1) stop("no data to predict at")
  nx<-length(x)
  nk<-object$bs.dim
  X <- rep(0,nx*nk);S<-rep(0,nk*nk);C<-rep(0,nk);control<-0

  oo <- .C(C_construct_cr,as.double(x),as.integer(nx),as.double(object$xp),
            as.integer(object$bs.dim),as.double(X),as.double(S),
                   as.double(C),as.integer(control))
  X<-matrix(oo[[5]],nx,nk) # the prediction matrix
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    X<-as.numeric(by)*X
  }
  X
}

Predict.matrix.cs.smooth<-function(object,data)
# this is the prediction method for a cubic regression spline 
# with shrinkage
{ Predict.matrix.cr.smooth(object,data)
}

Predict.matrix.tprs.smooth<-function(object,data)
# prediction matrix method for a t.p.r.s. term 
{ x<-array(0,0)
  for (i in 1:object$dim) 
  { xx <- get.var(object$term[[i]],data)
    xx <- xx - object$shift[i]
    if (length(xx)<1) stop("no data to predict at")
    x<-c(x,xx)
  }
  n<-nrow(data)
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    by.exists<-TRUE
  } else
  { by<-0;by.exists<-FALSE}
  X<-matrix(0,n,object$bs.dim)
  oo<-.C(C_predict_tprs,as.double(x),as.integer(object$dim),as.integer(n),as.integer(object$p.order),
      as.integer(object$bs.dim),as.integer(object$null.space.dim),as.double(object$Xu),
      as.integer(nrow(object$Xu)),as.double(object$UZ),as.double(by),as.integer(by.exists),X=as.double(X))
  X<-matrix(oo$X,n,object$bs.dim)
}

Predict.matrix.ts.smooth<-function(object,data)
# this is the prediction method for a t.p.r.s
# with shrinkage
{ Predict.matrix.tprs.smooth(object,data)
}

smooth.construct <- function(object,data,knots) UseMethod("smooth.construct")


Predict.matrix <- function(object,data) UseMethod("Predict.matrix")


smoothCon <- function(object,data,knots,absorb.cons=FALSE,scale.penalty=TRUE)
## wrapper function which calls smooth.construct methods, but can then modify
## the parameterizaion used. If absorb.cons==TRUE then a constraint free
## parameterization is used. 
{ sm <- smooth.construct(object,data,knots)
  if (!is.null(attr(sm,"qrc"))) warning("smooth objects should not have a qrc attribute.")
 
  ## following is intended to make scaling `nice' for better gamm performance
  if (scale.penalty && length(sm$S)>0) # then the penalty coefficient matrix is rescaled
  { maXX <- mean(abs(t(sm$X)%*%sm$X)) # `size' of X'X
    for (i in 1:length(sm$S)) {
      maS <- mean(abs(sm$S[[i]]))
      sm$S[[i]] <- sm$S[[i]] * maXX / maS
    }
  } 

  if (absorb.cons)
  { k<-ncol(sm$X)
    j<-nrow(sm$C)
    if (j>0) # there are constraints
    { qrc<-qr(t(sm$C))
      if (length(sm$S)>0)
      for (l in 1:length(sm$S)) # tensor product terms have > 1 penalty 
      { ZSZ<-qr.qty(qrc,sm$S[[l]])[(j+1):k,]
        sm$S[[l]]<-t(qr.qty(qrc,t(ZSZ))[(j+1):k,])
      }
      sm$X<-t(qr.qy(qrc,t(sm$X))[(j+1):k,])
      #sm$qrc<-qrc
      attr(sm,"qrc") <- qrc
      attr(sm,"nCons") <- j;
      sm$C <- NULL
      sm$rank <- pmin(sm$rank,k-j)
      ## ... so qr.qy(sm$qrc,c(rep(0,nrow(sm$C)),b)) gives original para.'s
    }  
  } else attr(sm,"qrc") <-NULL

  sm 
}

PredictMat <- function(object,data)
## wrapper function which calls Predict.matrix and imposes same constraints as 
## smoothCon on resulting Prediction Matrix
{ X <- Predict.matrix(object,data)
  qrc <- attr(object,"qrc")
  if (!is.null(qrc)) { ## then a transform has been applied by smoothCon
    j<-attr(object,"nCons");k<-ncol(X)
    X<-t(qr.qy(qrc,t(X))[(j+1):k,])
  }
  del.index <- attr(object,"del.index")
  if (!is.null(del.index)) X <- X[,-del.index]
  X
}







interpret.gam<-function (gf)
# interprets a gam formula of the generic form:
#   y~x0+x1+x3*x4 + s(x5)+ s(x6,x7) ....
# and returns:
# 1. a model formula for the parametric part: pf (and pfok indicating whether it has terms)
# 2. a list of descriptors for the smooths: smooth.spec
# 3. a full version of the formulae with all terms expanded in full
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
#  if (!is.null(off)) 
#  { sp[sp>off]<-sp[sp>off]-1 # have to remove the offset from this index list 
#    tp[tp>off]<-tp[tp>off]-1
#  } # commented out 26/11/05 after reworking of term handling above
  ns<-length(sp)+length(tp) # number of smooths
  k<-kt<-ks<-kp<-1 # counters for terms in the 2 formulae
  len.sp <- length(sp)
  len.tp <- length(tp)

  smooth.spec<-list()
  if (nt)
  for (i in 1:nt) # work through all terms
  { if (k<=ns&&((ks<=len.sp&&sp[ks]==i)||(kt<=len.tp&&tp[kt]==i))) # it's a smooth
    { st<-eval(parse(text=terms[i]),envir=p.env)
      if (k>1||kp>1) rf<-paste(rf,"+",st$full.call,sep="") # add to full formula
      else rf<-paste(rf,st$full.call,sep="")
      smooth.spec[[k]]<-st
      if (ks<=len.sp&&sp[ks]==i) ks<-ks+1  # counts s() terms
      else kt<-kt+1              # counts te() terms
      k<-k+1     # counts smooth terms 
    } else          # parametric
    { if (kp>1) pf<-paste(pf,"+",terms[i],sep="") # add to parametric formula
      else pf<-paste(pf,terms[i],sep="")
      if (k>1||kp>1) rf<-paste(rf,"+",terms[i],sep="") # add to full formula
      else rf<-paste(rf,terms[i],sep="")
      kp<-kp+1    # counts parametric terms
    }
  }    
  if (!is.null(off)) # deal with offset
  { if (kp>1) pf<-paste(pf,"+",sep="")
    if (kp>1||k>1) rf<-paste(rf,"+",sep="")
    pf<-paste(pf,as.character(attr(tf,"variables")[1+off]),sep="")
    rf<-paste(rf,as.character(attr(tf,"variables")[1+off]),sep="")
    kp<-kp+1          
  }
  if (attr(tf,"intercept")==0) 
  {pf<-paste(pf,"-1",sep="");rf<-paste(rf,"-1",sep="");if (kp>1) pfok<-1 else pfok<-0}
  else { pfok<-1;if (kp==1) { pf<-paste(pf,"1"); if (k==1) rf<-paste(rf,"1",sep="");}}
  
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
  ret<-list(pf=as.formula(pf,p.env),pfok=pfok,smooth.spec=smooth.spec,full.formula=as.formula(rf,p.env),
            fake.formula=fake.formula,response=response)
  class(ret)<-"split.gam.formula"
  ret
}

fixDependence <- function(X1,X2,tol=.Machine$double.eps^.5)
# model matrix X2 may be linearly dependent on X1. This 
# routine finds which columns of X2 should be zeroed to 
# fix this.
{ qr1 <- qr(X1,LAPACK=TRUE)
  r<-ncol(X1);n<-nrow(X1)
  QtX2 <- qr.qty(qr1,X2)[(r+1):n,] # Q'X2
  qr2 <- qr(QtX2,LAPACK=TRUE)
  R <- qr.R(qr2)
  # now final diagonal block of R may be zero, indicating rank 
  # deficiency. 
  r0<-r<-nrow(R)
  while (mean(abs(R[r0:r,r0:r]))<abs(R[1,1])*tol) r0 <- r0 -1
  r0<-r0+1
  if (r0>r) return(NULL) else
  qr2$pivot[r0:r] # the columns of X2 to zero in order to get independence
}


gam.side <- function(sm,tol=.Machine$double.eps^.5)
# works through a list of smooths, aiming to identify nested or partially
# nested terms, and impose identifiability constraints on them
{ m <- length(sm)
  if (m==0) return(sm)
  v.names<-array("",0);maxDim<-1
  for (i in 1:m) { ## collect all term names and max smooth `dim'
    vn <- sm[[i]]$term
    ## need to include by variables in names
    if (sm[[i]]$by!="NA") vn <- paste(vn,sm[[i]]$by,sep="")
    sm[[i]]$vn <- vn ## use this record to identify variables from now
    v.names <- c(v.names,vn)
    if (sm[[i]]$dim > maxDim) maxDim <- sm[[i]]$dim
  } 
  lv <- length(v.names)   
  v.names <- unique(v.names)
  if (lv == length(v.names)) return(sm) ## no repeats => no nesting
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
        X1 <- matrix(0,nrow(sm[[i]]$X),0)
        for (j in 1:d) { ## work through variables
          b <- sm.id[[sm[[i]]$vn[j]]] # list of smooths dependent on this variable
          k <- (1:length(b))[b==i] ## locate current smooth in list 
          if (k>1) for (l in 1:(k-1)) { ## collect X columns
            X1 <- cbind(X1,sm[[b[l]]]$X)
          }
        } ## Now X1 contains columns for all lower dimensional terms
        if (ncol(X1)==0) ind <- NULL else
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



gam.setup <- function(formula,pterms,data=stop("No data supplied to gam.setup"),knots=NULL,sp=NULL,
                    min.sp=NULL,H=NULL,fit.method="magic",parametric.only=FALSE,absorb.cons=FALSE,
                    max.tprs.knots = NULL)
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
  
  G<-list(m=m,full.formula=split$full.formula,min.sp=min.sp,H=H)

  

  if (fit.method=="fastest") 
  { if (G$m==1) G$fit.method<-"mgcv" else G$fit.method<-"magic"
  } else G$fit.method<-fit.method

  if (is.null(attr(data,"terms"))) # then data is not a model frame
  mf<-model.frame(split$pf,data,drop.unused.levels=FALSE) # must be false or can end up with wrong prediction matrix!
  else mf<-data # data is already a model frame

  G$intercept <-  attr(attr(mf,"terms"),"intercept")>0
  G$offset <- model.offset(mf)   # get the model offset (if any)

  # construct model matrix.... 
  
  X <- model.matrix(pterms,mf)
  G$nsdf <- ncol(X)
  G$contrasts <- attr(X,"contrasts")
  G$xlevels <- .getXlevels(pterms,mf)
  G$assign <- attr(X,"assign") # used to tell which coeffs relate to which pterms

  if (parametric.only) { G$X<-X;return(G)}
  
  # next work through smooth terms (if any) extending model matrix.....
  
  G$smooth<-list()
  G$S<-list()
 
  G$off<-array(0,0)
  first.para<-G$nsdf+1
  sm <- list()
  if (m>0) for (i in 1:m) 
  { # idea here is that terms are set up in accordance with information given in split$smooth.spec
    # appropriate basis constructor is called depending on the class of the smooth
    # constructor returns penalty matrices model matrix and basis specific information
    attr(split$smooth.spec[[i]],"max.tprs.knots") <- max.tprs.knots
    sm[[i]] <- smoothCon(split$smooth.spec[[i]],data,knots,absorb.cons)
  }
  
  ## at this stage, it is neccessary to impose any side conditions required
  ## for identifiability
  if (m>0) sm<-gam.side(sm,tol=.Machine$double.eps^.5)

  if (m>0) for (i in 1:m) 
  { n.para<-ncol(sm[[i]]$X)
    # define which elements in the parameter vector this smooth relates to....
    sm[[i]]$first.para<-first.para     
    first.para<-first.para+n.para
    sm[[i]]$last.para<-first.para-1
    X<-cbind(X,sm[[i]]$X);sm[[i]]$X<-NULL
   
    G$smooth[[i]] <- sm[[i]]   
  }
  G$X<-X;rm(X)
  n.p<-ncol(G$X) 
  # deal with penalties
  n.smooths<-0 # number of `free' penalties
  k.sp<-0 # count through sp and min.sp
  G$rank<-array(0,0)
  G$sp<-array(0,0)
  if (m>0) for (i in 1:m)
  { sm<-G$smooth[[i]]
    if (length(sm$S)>0)
    for (j in 1:length(sm$S))  # work through penalty matrices
    { k.sp<-k.sp+1
      if (is.null(sp)||sp[k.sp]<0) # s.p.s to be estimated
      { n.smooths<-n.smooths+1
        G$off[n.smooths]<-sm$first.para 
        G$S[[n.smooths]]<-sm$S[[j]]
        G$rank[n.smooths]<-sm$rank[j]
        G$sp[n.smooths] <- -1 # confirms that it's to be estimated
      } else
      if (sp[k.sp]>0) # add penalty to fixed penalty H
      { if (is.null(H)) H<-matrix(0,n.p,n.p)
        H[sm$first.para:sm$last.para,sm$first.para:sm$last.para]<-
        H[sm$first.para:sm$last.para,sm$first.para:sm$last.para]+sp[k.sp]*sm$S[[j]]  
      } # if sp is zero then ignore term
      if (!is.null(min.sp))
      { if (is.null(H)) H<-matrix(0,n.p,n.p)
        H[sm$first.para:sm$last.para,sm$first.para:sm$last.para]<-
        H[sm$first.para:sm$last.para,sm$first.para:sm$last.para]+min.sp[k.sp]*sm$S[[j]] 
      }           
    } 
  }
  if (!is.null(sp)&&length(sp)!=k.sp) stop("supplied sp has wrong length")  
  if (!is.null(min.sp)&&length(min.sp)!=k.sp) stop("supplied min.sp has wrong length")
  G$H<-H

  if (!is.null(sp)) # then user has supplied fixed smoothing parameters
  { ok<-TRUE
    if (length(sp)!=k.sp) { ok<-FALSE;warning("Supplied smoothing parameter vector is too short - ignored.")}
    if (sum(is.na(sp))) { ok<-FALSE;warning("NA's in supplied smoothing parameter vector - ignoring.")}
  } else # set up for auto-initialization
  G$sp<-rep(-1,k.sp) # is this really needed?
  if (is.null(sp)) G$all.sp<-G$sp else G$all.sp <- sp
  if (!is.null(min.sp)) # then minimum s.p.'s supplied
  { if (length(min.sp)!=n.smooths) stop("length of min.sp is wrong.")
    if (sum(is.na(min.sp))) stop("NA's in min.sp.")
    if (sum(min.sp<0)) stop("elements of min.sp must be non negative.")
  }


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
  G$y <- data[[deparse(split$full.formula[[2]],backtick=TRUE)]]
  
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

gam.method.description <- function(method,am=TRUE)
## produces short fitting method description string
{ if (am) return(method$am)
  if (method$gam=="perf.magic") return("performance iteration - magic")
  if (method$gam=="perf.mgcv") return("performance iteration - mgcv")
  if (method$gam=="perf.outer") return(paste("perf. iter. magic + outer",method$outer))
  if (method$pearson==FALSE) {
    if (method$outer=="nlm") return("deviance based outer iter. - nlm exact derivs.")
    if (method$outer=="optim")  return("deviance based outer iter. - Quasi-Newton exact derivs.")
    if (method$outer=="nlm.fd") return("deviance based outer iter. - nlm with finite differences.")
  } else {
    if (method$outer=="nlm") return("Pearson based outer iter. - nlm exact derivs.")
    if (method$outer=="optim")  return("Pearson based outer iter. - Quasi-Newton exact derivs.")
    if (method$outer=="nlm.fd") return("Pearson based outer iter. - nlm with finite differences.")
  }
}

gam.method <- function(am="magic",gam="outer",outer="nlm",pearson=FALSE,family=NULL)
# Function for returning fit method control list for gam.
# am controls the fitting method to use for pure additive models.
# gam controls the type of iteration to use for Gams.
# outer controls the optimization method to use when using outer
# looping with gams.
# pearson determined whether to base scores on Pearson statistic or deviance.
{ if (sum(am==c("mgcv","magic"))==0) stop("Unknown additive model fit method.") 
  if (sum(gam==c("perf.magic","perf.mgcv","perf.outer","outer","outer"))==0) 
  stop("Unknown *generalized* additive model fit method.") 
  if (sum(outer==c("optim","nlm","nlm.fd"))==0) 
  stop("Unknown GAM outer optimizing method.") 
  if (!is.logical(pearson)) {pearson <- FALSE;
    warning("pearson should be TRUE or FALSE - set to FALSE.")
  } 
  if (!is.null(family)&&substr(family$family,1,17)=="Negative Binomial" 
       &&gam!="perf.magic"&&gam!="perf.mgcv") gam <- "perf.magic"  
  list(am=am,gam=gam,outer=outer,pearson=pearson)
}

gam.outer <- function(lsp,fscale,family,control,method,gamma,G)
# function for smoothing parameter estimation by outer optimization. i.e.
# P-IRLS scheme iterated to convergence for each trial set of smoothing
# parameters.
{ if (method$outer=="nlm.fd") {
    um<-nlm(full.score,lsp,typsize=lsp,fscale=fscale, stepmax = 
            control$nlm$stepmax, ndigit = control$nlm$ndigit,
	    gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, 
            iterlim = control$nlm$iterlim, G=G,family=family,control=control,
            gamma=gamma,pearson=method$pearson)
    lsp<-um$estimate
    object<-attr(full.score(lsp,G,family,control,gamma=gamma,pearson=method$pearson),"full.gam.object")
    object$gcv.ubre <- um$minimum
    object$outer.info <- um
  } else { ## methods calling gam.fit2
    if (substr(family$family,1,17)=="Negative Binomial") 
    stop("Negative binomial family not (yet) usable with type 2 iteration methods.")
    if (is.null(attr(G$smooth[[1]],"qrc"))) 
    stop("Must use gam.control(absorb.cons=TRUE), for type 2 iteration
    methods.")
    family <- fix.family.link(family)
    family <- fix.family.var(family)
    G$rS <- mini.roots(G$S,G$off,ncol(G$X))
    if (G$sig2>0) {criterion <- "UBRE";scale <- G$sig2} else { criterion <- "GCV";scale<-1}
    args <- list(X=G$X,y=G$y,S=G$S,rS=G$rS,off=G$off,H=G$H,offset=G$offset,family=family,
             weights=G$w,control=control,scoreType=criterion,gamma=gamma,scale=scale,pearson=method$pearson)
    if (method$outer=="nlm") {
      b <- nlm(gam3objective, lsp, typsize = lsp, fscale = fscale, 
            stepmax = control$nlm$stepmax, ndigit = control$nlm$ndigit,
	    gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, 
            iterlim = control$nlm$iterlim,
	    check.analyticals=control$nlm$check.analyticals,
            args=args)
      lsp <- b$estimate

    } else if (method$outer=="optim") {
      b<-optim(par=lsp,fn=gam2objective,gr=gam2derivative,method="L-BFGS-B",control=
         list(fnscale=fscale,factr=control$optim$factr,lmm=min(5,length(lsp))),args=args)
      lsp <- b$par
    }
    obj <- gam2objective(lsp,args,printWarn=TRUE) # final model fit, with warnings 
    object <- attr(obj,"full.fit")
    object$outer.info <- b
    object$gcv.ubre <- as.numeric(obj)
    if (args$scoreType=="GCV") object$scale <- object$scale.est else
    object$scale <- G$sig2
    mv<-magic.post.proc(G$X,object,w=sqrt(object$weights))
    object$Vp <- mv$Vb
    object$hat<-mv$hat
    object$Ve <- mv$Ve
    object$edf<-mv$edf
    object$aic <- object$aic + 2*sum(mv$edf)
    object$nsdf <- G$nsdf
    object$GCV<-object$GCV1<-object$UBRE<-object$UBRE1<-object$trA<-
    object$trA1<-object$alpha<-object$alpha1<-object$rV<-object$scale.est<-NULL
    object$sig2 <- object$scale
  }
  object$sp <- exp(lsp)
  object
}


gam <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,offset=NULL,
                control=gam.control(),method=gam.method(),
                scale=0,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=TRUE,G=NULL,...)

# Routine to fit a GAM to some data. The model is stated in the formula, which is then 
# interpreted to figure out which bits relate to smooth terms and which to parametric terms.

{  if (is.null(G))
   { # create model frame..... 
    gp<-interpret.gam(formula) # interpret the formula 
    cl<-match.call() # call needed in gam object for update to work
    mf<-match.call(expand.dots=FALSE)
    mf$formula<-gp$fake.formula 
    mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-
               mf$gamma<-mf$method<-mf$fit<-mf$G<-mf$...<-NULL
    mf$drop.unused.levels<-TRUE
    mf[[1]]<-as.name("model.frame")
    pmf <- mf
    mf <- eval(mf, parent.frame()) # the model frame now contains all the data 
    if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf,"terms")
    
    pmf$formula <- gp$pf
    pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part

    pterms <- attr(pmf,"terms")

    if (is.character(family)) family<-eval(parse(text=family))
    if (is.function(family)) family <- family()
    if (is.null(family$family)) stop("family not recognized")
  
    if (family$family=="gaussian" && family$link=="identity") am <- TRUE
    else am <- FALSE
    
    if (am) fit.method <- method$am else { 
      if (method$gam=="perf.mgcv") fit.method <- "mgcv" else fit.method <- "magic"}

    G<-gam.setup(gp,pterms=pterms,data=mf,knots=knots,sp=sp,min.sp=min.sp,
                 H=H,fit.method=fit.method,parametric.only=FALSE,absorb.cons=control$absorb.cons,
                 max.tprs.knots=control$max.tprs.knots)
    
    if (ncol(G$X)>nrow(G$X)+nrow(G$C)) stop("Model has more coefficients than data")

    G$terms<-terms;G$pterms<-pterms
    G$mf<-mf;G$cl<-cl;
    G$am <- am

    if (is.null(G$offset)) G$offset<-rep(0,G$n)
     
    method <- gam.method(method$am,method$gam,method$outer,method$pearson,family) # checking it's ok

    if (scale==0) 
    { if (family$family=="binomial"||family$family=="poisson") scale<-1 #ubre
      else scale <- -1 #gcv
    }
  
    G$sig2<-scale

    G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
    G$max.half<-control$mgcv.half # max step halving in Newton update mgcv
    G$min.edf<-G$nsdf-dim(G$C)[1]
    if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim

    environment(G$full.formula)<-environment(formula) 
    G$formula<-formula
    environment(G$formula)<-environment(formula)
  }

  if (!fit) return(G)

  # is outer looping needed ?
  outer.looping <- !G$am && (method$gam=="perf.outer"||method$gam=="outer") &&
                    length(G$S)>0 && sum(G$sp<0)!=0

  # take only a few IRLS step to get scale estimates for "pure" outer
  # looping...
    
  if (outer.looping && method$gam=="outer") fixedSteps <- control$outerPIsteps else 
      fixedSteps <- control$globit+control$maxit+2
  
  object<-gam.fit(G,family=family,control=control,gamma=gamma,fixedSteps=fixedSteps,...)
  
  if (!is.null(sp)) # fill returned s.p. array with estimated and supplied terms
  { temp.sp<-object$sp
    object$sp<-sp
    object$sp[sp<0]<-temp.sp
  }
   
  # mgcv.conv<-object$mgcv.conv
   
  if (outer.looping)
  { # use perf.iter s.p. estimates from gam.fit as starting values...
    lsp<-log(object$sp[G$all.sp<0]) # make sure only free s.p.s are optimized!
    # don't allow initial sp's too far from defaults, otherwise optimizers may
    # get stuck on flat portions of GCV/UBRE score....
    if (method$gam!="perf.outer") { 
      lsp2 <- log(initial.sp(G$X,G$S,G$off)) 
      ind <- lsp > lsp2+5;lsp[ind] <- lsp2[ind]+5
      ind <- lsp < lsp2-5;lsp[ind] <- lsp2[ind]-5 
      if (fixedSteps<1) lsp <- lsp2 ## don't use perf iter sp's at all
    }
   # temp.sp <-object$sp # keep copy off all sp's
    mgcv.conv <- object$mgcv.conv  
  
    object <- gam.outer(lsp,fscale=abs(object$gcv.ubre)+object$sig2/length(object$y),family=family,
                        control=control,method=method,gamma=gamma,G=G)
    object$mgcv.conv <- mgcv.conv 
    temp.sp <- G$all.sp
    temp.sp[G$all.sp<0] <- object$sp # copy estimated sp's into whole vector
    object$sp <- temp.sp   # correct object sp vector
  }

  ## correct null deviance if there's an offset ....

  if (G$intercept&&any(G$offset)) object$null.deviance <-
                                  glm(G$y~offset(G$offset),family=family)$deviance

  if (G$sig2<0) object$method <- "GCV" else object$method <- "UBRE"

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
  object$full.formula<-as.formula(G$full.formula)
  environment(object$full.formula)<-environment(G$formula) 
  object$formula<-G$formula
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
  object$fit.method <- gam.method.description(method,G$am)
  object$call<-G$cl # needed for update() to work
  class(object)<-c("gam","glm","lm")
  object
}


gam.check <- function(b)
# takes a fitted gam object and produces some standard diagnostic plots
{ if (b$fit.method=="mgcv"||b$fit.method=="performance iteration - mgcv")
  fit.method <- "mgcv" else fit.method <- "other"
  if (b$method=="GCV"||b$method=="UBRE")
  { old.par<-par(mfrow=c(2,2))
    sc.name<-b$method
    if (fit.method=="mgcv")
    { if (b$mgcv.conv$iter>0)
      { plot(b$mgcv.conv$edf,b$mgcv.conv$score,xlab="Estimated Degrees of Freedom",
         ylab=paste(sc.name,"Score"),main=paste(sc.name,"w.r.t. model EDF"),type="l")
        points(b$nsdf+sum(b$edf),b$gcv.ubre,col=2,pch=20)
      }
    } else qqnorm(residuals(b))
    plot(b$linear.predictors,residuals(b),main="Resids vs. linear pred.",
         xlab="linear predictor",ylab="residuals");
    hist(residuals(b),xlab="Residuals",main="Histogram of residuals");
    plot(fitted(b),b$y,xlab="Fitted Values",ylab="Response",main="Response vs. Fitted Values")
    if (b$mgcv.conv$iter>0)   
    cat("\nSmoothing parameter selection converged after",b$mgcv.conv$iter,"iteration")
    else 
    cat("\nModel required no smoothing parameter selection")
    if (b$mgcv.conv$iter>1) cat("s")
    
    if ((fit.method=="mgcv"&&b$mgcv.conv$step.fail)||(b$fit.method=="magic"&&!b$mgcv.conv$fully.converged)) 
    cat(" by steepest\ndescent step failure.\n") else cat(".\n")
    if (fit.method=="mgcv")
    { if (length(b$smooth)>1&&b$mgcv.conv$iter>0)
      { cat("The mean absolute",sc.name,"score gradient at convergence was ",mean(abs(b$mgcv.conv$g)),".\n")
        if (sum(b$mgcv.conv$e<0)) cat("The Hessian of the",sc.name ,"score at convergence was not positive definite.\n")
        else cat("The Hessian of the",sc.name,"score at convergence was positive definite.\n")
      }
      if (!b$mgcv.conv$init.ok&&(b$mgcv.conv$iter>0)) cat("Note: the default second smoothing parameter guess failed.\n")
    } else
    { cat("The RMS",sc.name,"score gradiant at convergence was",b$mgcv.conv$rms.grad,".\n")
      if (b$mgcv.conv$hess.pos.def)
      cat("The Hessian was positive definite.\n") else cat("The Hessian was not positive definite.\n")
      cat("The estimated model rank was ",b$mgcv.conv$rank," (maximum possible: ",b$mgcv.conv$full.rank,")\n",sep="")
    }
    cat("\n")
    par(old.par)
  }
  else
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
    for (i in 1:n.smooth)
    edf[i]<-sum(x$edf[x$smooth[[i]]$first.para:x$smooth[[i]]$last.para])
    cat("\nEstimated degrees of freedom:\n",edf,"  total = ",sum(x$edf),"\n")
  }
  if (x$method=="GCV")  
  cat("\nGCV score: ",x$gcv.ubre,"\n")
  else if (x$method=="UBRE")
  cat("\nUBRE score: ",x$gcv.ubre,"\n")
}

gam.control <- function (irls.reg=0.0,epsilon = 1e-06, maxit = 100,globit = 20,
                         mgcv.tol=1e-7,mgcv.half=15,nb.theta.mult=10000,trace =FALSE,
                         rank.tol=.Machine$double.eps^0.5,absorb.cons=TRUE,
                         max.tprs.knots=5000,nlm=list(),optim=list(),outerPIsteps=4) 
# Control structure for a gam. 
# irls.reg is the regularization parameter to use in the GAM fitting IRLS loop.
# epsilon is the tolerance to use in the IRLS MLE loop. maxit is the number 
# of IRLS iterations to use with local search for optimal s.p. after globit iterations have used global 
# searches. mgcv.tol is the tolerance to use in the mgcv call within each IRLS. mgcv.half is the 
# number of step halvings to employ in the mgcv search for the optimal GCV score, before giving up 
# on a search direction. trace turns on or off some de-bugging information.
# nb.theta.mult controls the upper and lower limits on theta estimates - for use with negative binomial  
# for single s.p. case and "magic" otherwise. 
# rank.tol is the tolerance to use for rank determination
# outerPIsteps is the number of performance iteration steps used to intialize
#                         outer iteration
{
    if (!is.numeric(irls.reg) || irls.reg <0.0) stop("IRLS regularizing parameter must be a non-negative number.")
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(globit) || globit <= 0) 
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(nb.theta.mult)||nb.theta.mult<2) 
        stop("nb.theta.mult must be >= 2")
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

    # and optim defaults
    if (is.null(optim$factr)) optim$factr <- 1e7
    optim$factr <- abs(optim$factr)

    list(irls.reg=irls.reg,epsilon = epsilon, maxit = maxit,globit = globit,
         trace = trace, mgcv.tol=mgcv.tol,mgcv.half=mgcv.half,nb.theta.mult=nb.theta.mult,
         rank.tol=rank.tol,absorb.cons=absorb.cons,max.tprs.knots=max.tprs.knots,nlm=nlm,
         optim=optim,outerPIsteps=outerPIsteps)
    
}



mgcv.get.scale<-function(Theta,weights,good,mu,mu.eta.val,G)
# Get scale implied by current fit and trial -ve binom Theta, I've used
# mu and mu.eta.val used in fit rather than implied by it....
{ variance<-neg.bin(Theta)$variance
  w<-sqrt(weights[good]*mu.eta.val[good]^2/variance(mu)[good])
  wres<-w*(G$y-G$X%*%G$p)
  scale<-sum(wres^2)/(G$n-sum(G$edf)-G$nsdf)
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


full.score <- function(sp,G,family,control,gamma,pearson)
# function suitable for calling from nlm in order to polish gam fit
# so that actual minimum of score is found in generalized cases
{ G$sp<-exp(sp);
  # set up single fixed penalty....
  q<-NCOL(G$X)
  if (is.null(G$H)) G$H<-matrix(0,q,q)
  for (i in 1:length(G$S))
  { j<-ncol(G$S[[i]])
    off1<-G$off[i];off2<-off1+j-1
    G$H[off1:off2,off1:off2]<-G$H[off1:off2,off1:off2]+G$sp[i]*G$S[[i]]
  }
  G$S<-list() # have to reset since length of this is used as number of penalties
  xx<-gam.fit(G,family=family,control=control,gamma=gamma)
  if (pearson) res <- xx$gcv.ubre else res <- xx$gcv.ubre.dev
  attr(res,"full.gam.object")<-xx
  res
}

gam.fit <- function (G, start = NULL, etastart = NULL, 
    mustart = NULL, family = gaussian(), 
    control = gam.control(),gamma=1,
    fixedSteps=(control$maxit+control$globit+1)) 
# fitting function for a gam, modified from glm.fit.
# note that smoothing parameter estimates from one irls iterate are carried over to the next irls iterate
# unless the range of s.p.s is large enough that numerical problems might be encountered (want to avoid 
# completely flat parts of gcv/ubre score). In the latter case autoinitialization is requested.
# fixedSteps < its default causes at most fixedSteps iterations to be taken,
# without warning if convergence has not been achieved. This is useful for
# obtaining starting values for outer iteration.
{
    intercept<-G$intercept
    conv <- FALSE
    nobs <- NROW(G$y)
    nvars <- NCOL(G$X) # check this needed
    y<-G$y # original data
    X<-G$X # original design matrix
    if (nvars == 0) stop("Model seems to contain no terms")
    olm <- G$am   # only need 1 iteration as it's a pure additive model. 
    find.theta<-FALSE # any supplied -ve binomial theta treated as known, G$sig2 is scale parameter
    if (substr(family$family,1,17)=="Negative Binomial")
    { if (G$sig2<=0) find.theta<-TRUE # find theta by GCV
      # now get theta/initial theta
      V<-mu<-0.5
      while(all.equal(V,mu)==TRUE)
      { mu<-mu*2;V<-family$variance(mu)
        if (all.equal(V,mu)!=TRUE) Theta<-mu^2/(V-mu)
      }
      T.max<-Theta*control$nb.theta.mult;T.min<-Theta/control$nb.theta.mult
      if (family$family=="Negative Binomial") nb.link<-NULL # neg.bin family, no link choises
      else nb.link<-family$link # negative.binomial family, there's a choise of links
    }

    
    # obtain average element sizes for the penalties
    n.free<-length(G$S)
    if (n.free>0)
    { S.size<-0
      for (i in 1:n.free) S.size[i]<-mean(abs(G$S[[i]])) 
    }
    weights<-G$w # original weights
   
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
    if (G$fit.method=="magic") 
    { msp<-rep(-1,n.free) # free smoothing parameter vector for magic
      magic.control<-list(tol=G$conv.tol,step.half=G$max.half,maxit=control$maxit+control$globit,
                          rank.tol=control$rank.tol)
    } else
    { mgcv.control<-list(conv.tol=G$conv.tol,max.half=G$max.half,min.edf=G$min.edf,target.edf=-1)
    }
    for (iter in 1:(control$maxit+control$globit)) 
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
   
        z<-G$y <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w<- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        
        G$w<-w
        G$X<-X[good,]  # truncated design matrix       
		if (dim(X)[2]==1) dim(G$X)<-c(length(X[good,]),1) # otherwise dim(G$X)==NULL !!
        ngoodobs <- as.integer(nobs - sum(!good))
        ncols <- as.integer(1)
        # must set G$sig2 to scale parameter or -1 here....
        G$sig2<-scale

        if (G$fit.method=="mgcv"&&n.free>0) # check that s.p.'s haven't drifted too far apart
        { temp.sp<-G$sp;temp.S.size<-S.size*temp.sp
          # check if there is a danger of getting stuck on a flat section of gcv/ubre score...
          if (min(temp.sp)>0 && min(temp.S.size)<.Machine$double.eps^0.5*max(temp.S.size)) 
          G$sp<- rep(-1.0,n.free) # .... if so use use auto-initialization in mgcv
          if (control$trace) cat("Re-initializing smoothing parameters\n") 
          if (iter>control$globit) # solution could be cycling - use more cautious optimization approach
          { mgcv.control$target.edf<-G$nsdf+sum(G$edf)
          } else
          mgcv.control$target.edf<- -1 # want less cautious optimization - better at local minimum avoidance
        }
        if (sum(!is.finite(G$y))+sum(!is.finite(G$w))>0) 
        stop("iterative weights or data non-finite in gam.fit - regularization may help. See ?gam.control.")

        if (G$fit.method=="mgcv") 
        { mr<-mgcv(G$y,G$X,G$sp,G$S,G$off,G$C,G$w,H=G$H,scale=G$sig2,gcv=(G$sig2<0),control=mgcv.control)
          G$p<-mr$b;G$sp<-mr$sp;G$sig2<-mr$scale;G$gcv.ubre<-mr$score
          G$Vp<-mr$Vb;G$hat<-mr$hat;G$edf<-mr$edf;G$conv<-mr$info
        }
        else
        { mr<-magic(G$y,G$X,msp,G$S,G$off,G$rank,G$H,G$C,G$w,gamma=gamma,G$sig2,G$sig2<0,
                    ridge.parameter=control$irls.reg,control=magic.control)
          G$p<-mr$b;msp<-mr$sp;G$sig2<-mr$scale;G$gcv.ubre<-mr$score;
         
        }

        if (find.theta) # then family is negative binomial with unknown theta - estimate it here from G$sig2
        { Theta<-mgcv.find.theta(Theta,T.max,T.min,weights,good,mu,mu.eta.val,G,.Machine$double.eps^0.5)
          if (is.null(nb.link)) family<-neg.bin(Theta)
          else family<-do.call("negative.binomial",list(theta=Theta,link=nb.link))
          variance <- family$variance;dev.resids <- family$dev.resids
          aic <- family$aic
        }

        if (control$trace&&G$fit.method=="mgcv")
        { cat("sp: ",G$sp,"\n")
          plot(G$conv$edf,G$conv$score,xlab="EDF",ylab="GCV/UBRE score",type="l");
          points(G$nsdf+sum(G$edf),G$gcv.ubre,pch=20,col=2)
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
            iter > fixedSteps) {
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
    if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps)) 
            warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
        if (any(mu < eps)) 
            warning("fitted rates numerically 0 occurred")
    }
    
    residuals <- rep(NA, nobs)
    residuals[good] <- z - (eta - offset)[good]
    
    nr <- min(sum(good), nvars)

    wt <- rep(0, nobs)
    wt[good] <- w^2
   
    wtdmu <- if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    if (G$fit.method=="magic") # then some post processing is needed to extract covariance matrix etc...
    { mv<-magic.post.proc(G$X,mr,w=G$w)
      G$Vp<-mv$Vb;G$hat<-mv$hat;
      G$Ve <- mv$Ve # frequentist cov. matrix
      G$edf<-mv$edf
      G$conv<-mr$gcv.info
      G$sp<-msp
      rank<-G$conv$rank
    } else { 
      X <- as.vector(G$w) * G$X 
      X <- G$Vp %*% t(X)
      G$Ve <- X%*%t(X)/G$sig2 # frequentist cov. matrix 
      rank <- ncol(G$X)-ncol(G$C)
    }
    aic.model <- aic(y, n, mu, weights, dev) + 2 * sum(G$edf)
    if (scale < 0) { ## deviance based GCV
      gcv.ubre.dev <- length(y)*dev/(length(y)-gamma*sum(G$edf))^2
    } else { # deviance based UBRE, which is just AIC
      gcv.ubre.dev <- dev/length(y) + 2 * gamma * sum(G$edf)/length(y) - G$sig2
    }
  
	list(coefficients = as.vector(coef), residuals = residuals, fitted.values = mu, 
        family = family,linear.predictors = eta, deviance = dev,
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,  
        df.null = nulldf, y = y, converged = conv,sig2=G$sig2,edf=G$edf,hat=G$hat,
        boundary = boundary,sp = G$sp,nsdf=G$nsdf,Ve=G$Ve,Vp=G$Vp,mgcv.conv=G$conv,
        gcv.ubre=G$gcv.ubre,aic=aic.model,rank=rank,gcv.ubre.dev=gcv.ubre.dev)
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
#      == "lpmatrix" - for matrix mapping parameters to l.p.
# Steps are:
#  1. Set newdata to object$model if no newdata supplied
#  2. split up newdata into manageable blocks if too large
#  3. Obtain parametric model matrix by call to gam.setup (NOTE: take care to use number of 
#     levels in original data!)
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
# when type=="terms". 
# if `object' has an attribute `para.only' then only parametric terms of order
# 1 are returned for type=="terms" : i.e. only what termplot can handle.
#
# if no new data is supplied then na.action does nothing, otherwise 
# if na.action == "na.pass" then NA predictors result in NA predictions (as lm
#                   or glm)
#              == "na.omit" or "na.exclude" then NA predictors result in
#                       dropping

  if (type!="link"&&type!="terms"&&type!="response"&&type!="lpmatrix")  
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
  if (type=="terms")
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

  for (b in 1:n.blocks)  # work through prediction blocks
  { start<-stop+1
    stop<-start+b.size[b]-1
    if (n.blocks==1) data <- newdata else data<-newdata[start:stop,]
    X<-matrix(0,b.size[b],nb)

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
    { X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para]<-
                              PredictMat(object$smooth[[k]],data)
      if (type=="terms") ColNames[n.pterms+k]<-object$smooth[[k]]$label
    }
    # have prediction matrix for this block, now do something with it
    if (type=="lpmatrix") H[start:stop,]<-X else 
    if (type=="terms")
    { ##
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
          fit[start:stop,n.pterms+k]<-X[,first:last]%*%object$coefficients[first:last]
          if (se.fit) # diag(Z%*%V%*%t(Z))^0.5; Z=X[,first:last]; V is sub-matrix of Vp
          se[start:stop,n.pterms+k]<-
          sqrt(rowSums((X[,first:last]%*%object$Vp[first:last,first:last])*X[,first:last]))
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
      fit[start:stop]<-X%*%object$coefficients
      if (!is.null(k)) fit[start:stop]<-fit[start:stop]+model.offset(mf)
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
  if (type=="terms") attr(H,"constant") <- object$coefficients[1]
  H # ... and return
}

plot.gam<-function(x,residuals=FALSE,rug=TRUE,se=TRUE,pages=0,select=NULL,scale=-1,n=100,n2=40,
                   pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                   ylim=NULL,xlim=NULL,too.far=0.1,all.terms=FALSE,shade=FALSE,shade.col="gray80",
                   shift=0,trans=I,...)

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

{ sp.contour <- function(x,y,z,zse,xlab="",ylab="",zlab="",titleOnly=FALSE,
               se.plot=TRUE,se.mult=1,trans=I,shift=0,...)   
  # internal function for contouring 2-d smooths with 1 s.e. limits
  { gap<-median(zse,na.rm=TRUE)  
    zr<-max(trans(z+zse+shift),na.rm=TRUE)-min(trans(z-zse+shift),na.rm=TRUE) # plotting range  
    n<-10  
    while (n>1 && zr/n<2.5*gap) n<-n-1    
    zrange<-c(min(trans(z-zse+shift),na.rm=TRUE),max(trans(z+zse+shift),na.rm=TRUE))  
    zlev<-pretty(zrange,n)  
    yrange<-range(y);yr<-yrange[2]-yrange[1]  
    xrange<-range(x);xr<-xrange[2]-xrange[1]  
    ypos<-yrange[2]+yr/10
    args <- as.list(substitute(list(...)))[-1]
    args$x <- substitute(x);args$y <- substitute(y)
    args$type="n";args$xlab<-args$ylab<-"";args$axes<-FALSE
    do.call("plot",args)
#    plot(x,y,type="n",xlab="",ylab="",axes=FALSE)
    cs<-(yr/10)/strheight(zlab);if (cs>1) cs<-1 # text scaling based on height  
    cw<-par()$cxy[1]  
    tl<-strwidth(zlab);  
    if (tl*cs>3*xr/10) cs<-(3*xr/10)/tl  
    args <- as.list(substitute(list(...)))[-1]
    n.args <- names(args)
    zz <- trans(z+shift)
    args$x<-substitute(x);args$y<-substitute(y);args$z<-substitute(zz)
    if (!"levels"%in%n.args) args$levels<-substitute(zlev)
    if (!"lwd"%in%n.args) args$lwd<-2
    if (!"labcex"%in%n.args) args$labcex<-cs*.65
    if (!"axes"%in%n.args) args$axes <- FALSE
    if (!"add"%in%n.args) args$add <- TRUE
    do.call("contour",args)
    #contour(x,y,z,levels=zlev,lwd=2,labcex=cs*0.65,axes=FALSE,add=TRUE)  
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
      if (n.plots%%pages) last.pages<-0 else last.ppp<-n.plots-ppp*pages
    } 
    else last.ppp<-0
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
      p<-x$coefficients[first:last]      # relevent coefficients 
      fit<-X%*%p                         # fitted values
      if (se) se.fit<-sqrt(rowSums((X%*%x$Vp[first:last,first:last])*X))
      edf<-sum(x$edf[first:last])
      xterm <- x$smooth[[i]]$term
      if (is.null(xlab)) xlabel<- xterm else xlabel <- xlab
      if (is.null(ylab)) 
      ylabel<-paste("s(",xterm,",",as.character(round(edf,2)),")",sep="") else
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
      raw<-data.frame(x=x$model[xterm][[1]],y=x$model[yterm][[1]])
      n2<-max(10,n2)
      xm<-seq(min(raw$x),max(raw$x),length=n2)
      ym<-seq(min(raw$y),max(raw$y),length=n2)  
      xx<-rep(xm,n2)
      yy<-rep(ym,rep(n2,n2))
      if (too.far>0)
      exclude <- exclude.too.far(xx,yy,raw$x,raw$y,dist=too.far) else
      exclude <- rep(FALSE,n2*n2)
      if (x$smooth[[i]]$by!="NA")         # deal with any by variables
      { by<-rep(1,n);dat<-data.frame(x=xx,y=yy,by=by)
        names(dat)<-c(xterm,yterm,x$smooth[[i]]$by)
      } else
      { dat<-data.frame(x=xx,y=yy);names(dat)<-c(xterm,yterm)}  # prediction data.frame
      X <- PredictMat(x$smooth[[i]],dat)   # prediction matrix for this term
      first<-x$smooth[[i]]$first.para;last<-x$smooth[[i]]$last.para
      p<-x$coefficients[first:last]      # relevent coefficients 
      fit<-X%*%p                         # fitted values
      fit[exclude] <- NA                 # exclude grid points too far from data
      if (se) 
      { se.fit<-sqrt(rowSums((X%*%x$Vp[first:last,first:last])*X))
        se.fit[exclude] <- NA # exclude grid points too distant from data
      }
      edf<-sum(x$edf[first:last])
      if (is.null(main)) 
      { if (is.null(x$smooth[[i]]$margin))
        title<-paste("s(",xterm,",",yterm,",",as.character(round(edf,2)),")",sep="") else
        title<-paste("te(",xterm,",",yterm,",",as.character(round(edf,2)),")",sep="")
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
      { if (interactive()&& is.null(select) && pd[[i]]$dim<3 && i>1&&(i-1)%%ppp==0) 
        readline("Press return for next page....")
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
          if (partial.resids)
          { if (is.null(list(...)[["pch"]]))
            points(pd[[i]]$raw,trans(pd[[i]]$p.resid+shift),pch=".",...) else
            points(pd[[i]]$raw,trans(pd[[i]]$p.resid+shift),...)
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
      { if (interactive() && pd[[i]]$dim<3 && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
        if (pd[[i]]$dim==1)
        { if (scale==0&&is.null(ylim)) 
          { if (partial.resids) ylimit <- range(pd[[i]]$p.resid,na.rm=TRUE) else ylimit <-range(pd[[i]]$fit)}
          if (!is.null(ylim)) ylimit <- ylim
          plot(pd[[i]]$x,trans(pd[[i]]$fit+shift),type="l",,xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,ylim=trans(ylimit+shift),xlim=xlim,main=main,...)
          if (rug) 
	  { if (jit) rug(jitter(as.numeric(pd[[i]]$raw)),...)
            else rug(as.numeric(pd[[i]]$raw),...) 
          }
          if (partial.resids)
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
      if (interactive() && m && i%%ppp==0) 
      readline("Press return for next page....")
      termplot(x,se=se,rug=rug,col.se=1,col.term=1)
    } else { # figure out which plot is required
      if (select > m) { 
        select <- select - m # i.e. which parametric term
        term.labels <- attr(x$pterms,"term.labels")
        term.labels <- term.labels[order==1]
        if (select <= length(term.labels)) {
        if (interactive() && m &&i%%ppp==0) 
        readline("Press return for next page....")
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
  family <- object$family
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

## old summary and anova code starts here .....
## functions and classes renamed summary2 and anova2.


summary2.gam <- function (object,freq=TRUE,...) 
# summary method for gam object - provides approximate p values for terms + other diagnostics
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
  if (freq) Vp <- object$Ve else Vp <- object$Vp
  se<-0;for (i in 1:length(object$coefficients)) se[i] <- Vp[i,i]^0.5
  residual.df<-length(object$y)-sum(object$edf)
  if (object$nsdf>0) # individual parameters
  { p.coeff<-object$coefficients[1:object$nsdf]
    p.t<-p.coeff/se[1:object$nsdf]
    if (object$method=="UBRE") 
    p.pv<-2*pnorm(abs(p.t),lower.tail=FALSE) else
    p.pv<-2*pt(abs(p.t),df=residual.df,lower.tail=FALSE)
  } else {p.coeff<-p.t<-p.pv<-array(0,0)}
  
  term.labels<-attr(object$pterms,"term.labels")
  nt<-length(term.labels)
  if (nt>0) # individual parametric terms
  { np<-length(object$assign)
    Vb<-matrix(Vp[1:np,1:np],np,np)
    bp<-array(object$coefficients[1:np],np)
    pTerms.pv <- array(0,nt)
    attr(pTerms.pv,"names") <- term.labels
    pTerms.df <- pTerms.chi.sq <- pTerms.pv
    for (i in 1:nt)
    { ind <- object$assign==i
      b <- bp[ind];V <- Vb[ind,ind]
      pTerms.df[i] <- nb <- length(b)
      pTerms.chi.sq[i] <- b%*%solve(V,b)
      if (object$method=="UBRE")
      pTerms.pv[i]<-pchisq(pTerms.chi.sq[i],df=nb,lower.tail=FALSE)
      else     
      pTerms.pv[i]<-pf(pTerms.chi.sq[i]/nb,df1=nb,df2=residual.df,lower.tail=FALSE)
    }
  } else { pTerms.df<-pTerms.chi.sq<-pTerms.pv<-array(0,0)}

  m<-length(object$smooth) # number of smooth terms
  edf<-s.pv<-chi.sq<-array(0,m)
  if (m>0) # form test statistics for each smooth
  { for (i in 1:m)
    { start<-object$smooth[[i]]$first.para;stop<-object$smooth[[i]]$last.para
      V <- Vp[start:stop,start:stop] # cov matrix for smooth
      p<-object$coefficients[start:stop]  # params for smooth
      M1<-object$smooth[[i]]$df
      M<-round(sum(object$edf[start:stop]))
      V<-pinv(V,M1) # get rank M pseudoinverse of V
      chi.sq[i]<-t(p)%*%V%*%p
      er<-names(object$coefficients)[start]
      er<-substring(er,1,nchar(er)-2)
      if (object$smooth[[i]]$by!="NA") 
      { er<-paste(er,":",object$smooth[[i]]$by,sep="")} 
      names(chi.sq)[i]<-er
      edf[i]<-sum(object$edf[start:stop])
      if (freq) df <- attr(V,"rank") else df <- edf[i]
      if (object$method=="UBRE")
      s.pv[i]<-pchisq(chi.sq[i],df=df,lower.tail=FALSE)
      else     
      s.pv[i]<-pf(chi.sq[i]/df,df1=df,df2=residual.df,lower.tail=FALSE) 
   }
  }
  w <- object$prior.weights
  nobs <- nrow(object$model)
  r.sq<- 1 - var(w*(object$y-object$fitted.values))*(nobs-1)/(var(w*object$y)*residual.df) 
  dev.expl<-(object$null.deviance-object$deviance)/object$null.deviance
  ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,
       s.pv=s.pv,scale=object$sig2,r.sq=r.sq,family=object$family,formula=object$formula,n=nobs,
       dev.expl=dev.expl,edf=edf,dispersion=object$sig2,pTerms.pv=pTerms.pv,pTerms.chi.sq=pTerms.chi.sq,
       pTerms.df=pTerms.df)
  if (object$method=="GCV") ret$gcv<-object$gcv.ubre else if (object$method=="UBRE") ret$ubre<-object$gcv.ubre
  class(ret)<-"summary.gam2"
  ret
}

print.summary.gam2<-function(x,...)
# print method for gam summary method.
{ print(x$family)
  cat("Formula:\n")
  print(x$formula)
  if (length(x$p.coeff)>0)
  { cat("\nParametric coefficients:\n")
    width<-max(nchar(names(x$p.coeff)))
    cat(rep(" ",width),"   Estimate  std. err.    t ratio    Pr(>|t|)\n",sep="")
    for (i in 1:length(x$p.coeff))
    cat(formatC(names(x$p.coeff)[i],width=width)," ",formatC(x$p.coeff[i],width=10,digits=5)," ",
    formatC(x$se[i],width=10,digits=4)," ",formatC(x$p.t[i],width=10,digits=4),"    ",format.pval(x$p.pv[i]),"\n",sep="")
  }
  cat("\n")
  if(x$m>0)
  { cat("Approximate significance of smooth terms:\n")
    width<-max(nchar(names(x$chi.sq)))
    cat(rep(" ",width),"        edf       chi.sq     p-value\n",sep="")
    for (i in 1:x$m)
    cat(formatC(names(x$chi.sq)[i],width=width)," ",formatC(x$edf[i],width=10,digits=4),"   ",
    formatC(x$chi.sq[i],width=10,digits=5),"     ",format.pval(x$s.pv[i]),"\n",sep="")
  }
  cat("\nR-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5))
  if (length(x$dev.expl)>0) cat("   Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%\n",sep="")
  if (!is.null(x$ubre)) cat("UBRE score = ",formatC(x$ubre,digits=5),sep="")
  if (!is.null(x$gcv)) cat("GCV score = ",formatC(x$gcv,digits=5)," ",sep="")
  cat("  Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
}


anova.gam2 <- function (object, ..., dispersion = NULL, test = NULL)
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
    if (!is.null(dispersion)) warning("dispersion argument ignored")
    if (!is.null(test)) warning("test argument ignored")
    if (!inherits(object,"gam")) stop("anova.gam called with non gam object")
    sg <- summary(object,freq=TRUE) 
    class(sg) <- "anova.gam2"
    sg
}


print.anova.gam2 <- function(x,...)
{ # print method for class anova.gam resulting from single
  # gam model calls to anova.
  print(x$family)
  cat("Formula:\n")
  print(x$formula)
  if (length(x$pTerms.pv)>0)
  { cat("\nParametric Terms:\n")
    term.names <- names(x$pTerms.pv)
    width<-max(nchar(term.names))
    cat(rep(" ",width),"         df       chi.sq     p-value\n",sep="")
    for (i in 1:length(term.names))
    cat(formatC(term.names[i],width=width)," ",formatC(x$pTerms.df[i],width=10,digits=4),"   ",
    formatC(x$pTerms.chi.sq[i],width=10,digits=5),"     ",format.pval(x$pTerms.pv[i]),"\n",sep="")
  }
  cat("\n")
  if(x$m>0)
  { cat("Approximate significance of smooth terms:\n")
    width<-max(nchar(names(x$chi.sq)))
    cat(rep(" ",width),"        edf       chi.sq     p-value\n",sep="")
    for (i in 1:x$m)
    cat(formatC(names(x$chi.sq)[i],width=width)," ",formatC(x$edf[i],width=10,digits=4),"   ",
    formatC(x$chi.sq[i],width=10,digits=5),"     ",format.pval(x$s.pv[i]),"\n",sep="")
  }
}

## end of old summary and anova code

## Start of anova and summary code as improved by Henric Nilsson ....
## Added 10/8/05...

summary.gam <- function (object, dispersion = NULL, freq = TRUE, ...) 
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
      pTerms.df[i] <- nb <- length(b)      
      pTerms.chi.sq[i] <- b%*%solve(V,b)
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
  { for (i in 1:m)
    { start<-object$smooth[[i]]$first.para;stop<-object$smooth[[i]]$last.para
      V <- covmat[start:stop,start:stop] # cov matrix for smooth
      p<-object$coefficients[start:stop]  # params for smooth
      M1<-object$smooth[[i]]$df
      M<-min(M1,ceiling(2*sum(object$edf[start:stop]))) ## upper limit of 2*edf on rank
      V<-pinv(V,M) # get rank M pseudoinverse of V
      chi.sq[i]<-t(p)%*%V%*%p
      er<-names(object$coefficients)[start]
      er<-substring(er,1,nchar(er)-2)
      if (object$smooth[[i]]$by!="NA") 
      { er<-paste(er,":",object$smooth[[i]]$by,sep="")} 
      names(chi.sq)[i]<-er
      edf[i]<-sum(object$edf[start:stop])
      if (freq) df[i] <- attr(V, "rank") else df[i] <- edf[i]
      if (!est.disp)
      s.pv[i]<-pchisq(chi.sq[i], df = df[i], lower.tail = FALSE)
      else
      s.pv[i] <- pf(chi.sq[i]/df[i], df1 = df[i], df2 = residual.df, lower.tail = FALSE)
    }
    if (!est.disp) {
      if (freq) {
        s.table <- cbind(edf, df, chi.sq, s.pv)      
        dimnames(s.table) <- list(names(chi.sq), c("edf", "Est.rank", "Chi.sq", "p-value"))
      } else {
        s.table <- cbind(edf, chi.sq, s.pv)      
        dimnames(s.table) <- list(names(chi.sq), c("edf", "Chi.sq", "p-value"))
      }
    } else {
      if (freq) {
        s.table <- cbind(edf, df, chi.sq/df, s.pv)      
        dimnames(s.table) <- list(names(chi.sq), c("edf", "Est.rank", "F", "p-value"))
      } else {
        s.table <- cbind(edf, chi.sq/df, s.pv)      
        dimnames(s.table) <- list(names(chi.sq), c("edf", "F", "p-value"))
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
       pTerms.table = pTerms.table, s.table = s.table)
  if (object$method=="GCV") ret$gcv<-object$gcv.ubre else if (object$method=="UBRE") ret$ubre<-object$gcv.ubre
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
  if (!is.null(x$ubre)) cat("UBRE score = ",formatC(x$ubre,digits=5),sep="")
  if (!is.null(x$gcv)) cat("GCV score = ",formatC(x$gcv,digits=5)," ",sep="")
  cat("  Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
}


anova.gam <- function (object, ..., dispersion = NULL, test = NULL)
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
    sg <- summary(object, dispersion = dispersion, freq = TRUE)
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
}

## End of improved anova and summary code. 



cooks.distance.gam <- function(model,...)
{ res <- residuals(model,type="pearson")
  dispersion <- model$sig2
  hat <- model$hat
  p <- sum(model$edf)
  (res/(1 - hat))^2 * hat/(dispersion * p)
}

vcov.gam <- function(object, freq = TRUE, dispersion = NULL, ...)
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
    z<-matrix(z,n.grid,n.grid) # convert to matrix
    surf.col<-t(av)%*%z%*%av   # average over tiles  
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
      lo.col <- "gray"
      hi.col <- "gray"
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
# NOTE: W=diag(w^2). 
# flop count is O(nq^2) if X is n by q... this is why routine not part of magic
{ V<-object$rV%*%t(object$rV)
  if (!is.null(w)) 
  { if (is.matrix(w))
    X<-w%*%X else 
    X<-as.vector(w)*X # use recycling rule to form diag(w)%*%X cheaply 
    
  }
  M <- V%*%t(X);B <- X*t(M)
  Ve <- M%*%t(M)*object$scale # frequentist cov. matrix
  rm(M)
  hat <- apply(B,1,sum) # diag(X%*%V%*%t(X))
  edf <- apply(B,2,sum) # diag(V%*%t(X)%*%X)
  Vb <- V*object$scale;rm(V)
  list(Ve=Ve,Vb=Vb,hat=hat,edf=edf)
}

single.sp <- function(X,S,target=.5,tol=.Machine$double.eps*100)
## function to find smoothing parameter corresponding to particular 
## target e.d.f. for a single smoothing parameter problem. 
## X is model matrix; S is penalty matrix; target is target 
## average e.d.f. per penalized term.
{ R <- qr.R(qr(X))
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
  def.sp
}




magic <- function(y,X,sp,S,off,rank=NULL,H=NULL,C=NULL,w=NULL,gamma=1,scale=1,gcv=TRUE,
                ridge.parameter=NULL,control=list(maxit=50,tol=1e-6,step.half=25,
                rank.tol=.Machine$double.eps^0.5),extra.rss=0,n.score=length(y))
# Wrapper for C routine magic. Deals with constraints weights and square roots of 
# penalties. Currently only a diagonal weight matrix is allowed, but this 
# is easy to change.
# y is data vector, X is model matrix, sp is array of smoothing parameters,
# S is list of penalty matrices stored as smallest square submatrix excluding no 
# non-zero entries, off[i] is the location on the leading diagonal of the
# total penalty matrix of element (1,1) of S[[i]], rank is an array of penalty 
# ranks, H is any fixed penalty, C is a linear constraint matrix and w is the 
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
  Xo<-X
  if (!is.null(C)) # then impose constraints 
   { n.con<-dim(C)[1]
    ns.qr<-qr(t(C)) # last n.b-n.con columns of Q are the null space of C
    X<-t(qr.qty(ns.qr,t(X)))[,(n.con+1):n.b] # last n.b-n.con cols of XQ (=(Q'X')')
    # need to work through penalties forming Z'S_i^0.5 's
    if (n.p>0) for (i in 1:n.p) S[[i]]<-qr.qty(ns.qr,S[[i]])[(n.con+1):n.b,,drop=FALSE]
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
  icontrol[4]<-as.integer(!is.null(H));icontrol[5]<-length(sp);icontrol[6]<-control$step.half
  icontrol[7]<-control$maxit
  b<-array(0,icontrol[3])
  # argument names in call refer to returned values.
  um<-.C(C_magic,as.double(y),as.double(X),sp=as.double(sp),as.double(def.sp),as.double(Si),as.double(H),
          score=as.double(gamma),scale=as.double(scale),info=as.integer(icontrol),as.integer(cS),
          as.double(control$rank.tol),rms.grad=as.double(control$tol),b=as.double(b),rV=double(q*q),
          as.double(extra.rss),as.integer(n.score))
  res<-list(b=um$b,scale=um$scale,score=um$score,sp=um$sp)
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
  cat(paste("This is mgcv",version,"\n"))
}

set.mgcv.options <- function()
## function used to set optional value used in notLog
## and notExp...
{ options(mgcv.vc.logrange=25)
}

.onAttach <- function(...) { 
  print.mgcv.version()
  set.mgcv.options()
}


.First.lib <- function(lib, pkg) {
  library.dynam("mgcv", pkg, lib)
  print.mgcv.version()
  set.mgcv.options()
}


###############################################################################
### ISSUES.....






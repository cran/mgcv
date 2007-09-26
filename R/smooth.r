##  R routines for the package mgcv (c) Simon Wood 2000-2007

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
# for dimension d data M=(m+d-1)!/(d!(m-1)!). Any m not satisfying 2m>d is reset so 
# that 2m>d+1 (assuring "visual" smoothness) 
{ if (sum(d<0)) stop("d can not be negative in call to null.space.dimension().")
  ind<-2*m<d+1
  if (sum(ind)) # then default m required for some elements
  { m[ind]<-1;ind<-2*m<d+2
    while (sum(ind)) { m[ind]<-m[ind]+1;ind<-2*m<d+2;}
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



te <- function(..., k=NA,bs="cr",m=0,d=NA,by=NA,fx=FALSE,mp=TRUE,np=TRUE,xt=NA)
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

  # deal with `xt' extras list
  xtra <- list()
  if (length(xt)==1) for (i in 1:n.bases) xtra[[i]] <- xt else
  if (length(xt)==n.bases) xtra <- xt else
  stop("xt argument is faulty.")

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
                ",m=",deparse(m[i],backtick=TRUE),",xt=xtra[[i]]", ")")
    margin[[i]]<- eval(parse(text=stxt))  # NOTE: fx and by not dealt with here!
    j<-j1+1
  }
  # assemble term.label 
  if (mp) mp <- TRUE else mp <- FALSE
  if (np) np <- TRUE else np <- FALSE
  full.call<-paste("te(",term[1],sep="")
  if (dim>1) for (i in 2:dim) full.call<-paste(full.call,",",term[i],sep="")
  label<-paste(full.call,")",sep="")   # label for parameters of this term
 ## full.call<-paste(full.call,",k=",deparse(k,backtick=TRUE),",bs=",deparse(bs,backtick=TRUE),
 ##                  ",m=",deparse(m,backtick=TRUE),",d=",deparse(d,backtick=TRUE),
 ##                  ",by=",by.var,",fx=",deparse(fx,backtick=TRUE),",mp=",deparse(mp,backtick=TRUE),
 ##                  ",np=",deparse(np,backtick=TRUE),")",sep="")
  ret<-list(margin=margin,term=term,by=by.var,fx=fx,label=label,dim=dim,mp=mp,np=np)##full.call=full.call)
  class(ret)<-"tensor.smooth.spec"
  ret
}





s <- function (..., k=-1,fx=FALSE,bs="tp",m=0,by=NA,xt=NA)
# function for use in gam formulae to specify smooth term, e.g. s(x0,x1,x2,k=40,m=3,by=x3) specifies 
# a rank 40 thin plate regression spline of x0,x1 and x2 with a third order penalty, to be multiplied by
# covariate x3, when it enters the model.
# Returns a list consisting of the names of the covariates, and the name of any by variable,
# a model formula term representing the smooth, the basis dimension, the type of basis
# , whether it is fixed or penalized and the order of the penalty (0 for auto).
# xt contains information to be passed straight on to the basis constructor
{ vars<-as.list(substitute(list(...)))[-1] # gets terms to be smoothed without evaluation
 # call<-match.call() # get function call
  d<-length(vars) # dimension of smoother
  term<-deparse(vars[[d]],backtick=TRUE,width.cutoff=500) # last term in the ... arguments
  by.var<-deparse(substitute(by),backtick=TRUE,width.cutoff=500) #getting the name of the by variable
  if (by.var==".") stop("by=. not allowed")
  term<-deparse(vars[[1]],backtick=TRUE,width.cutoff=500) # first covariate
  if (term[1]==".") stop("s(.) not yet supported.")
  if (d>1) # then deal with further covariates
  for (i in 2:d)
  { term[i]<-deparse(vars[[i]],backtick=TRUE,width.cutoff=500)
    if (term[i]==".") stop("s(.) not yet supported.")
  }
  for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])),"term.labels")
  # term now contains the names of the covariates for this model term
  # now evaluate all the other 
  k.new<-round(k) # in case user has supplied non-integer basis dimension
  if (!all.equal(k.new,k)) {warning("argument k of s() should be integer and has been rounded")}
  k<-k.new
  if (length(k)==1 && k==-1) k<-10*3^(d-1) # auto-initialize basis dimension
  ind <- k<2
  if (sum(ind)) 
  { k[ind] <- 2
    warning("meaninglessly low k; reset to 2\n")
  }
  if (bs=="cr"||bs=="cc") # a check
  { if (d>1) { warning("cr/cc basis only works with 1-d smooths!\n");bs<-"tp";}
  } 
  m[m<0]<-0
  # check for repeated variables in function argument list
  if (length(unique(term))!=d) stop("Repeated variables as arguments of a smooth are not permitted")
  # assemble label for term
  full.call<-paste("s(",term[1],sep="")
  if (d>1) for (i in 2:d) full.call<-paste(full.call,",",term[i],sep="")
  label<-paste(full.call,")",sep="") # used for labelling parameters
##  full.call<-paste(full.call,",k=",deparse(k,backtick=TRUE,width.cutoff=500),",fx=",
##                   deparse(fx,backtick=TRUE,width.cutoff=500),",bs=",
##                   deparse(bs,backtick=TRUE,width.cutoff=500),",m=",deparse(m,backtick=TRUE,width.cutoff=500),
##                   ",by=",by.var,")",sep="")
  ret<-list(term=term,bs.dim=k,fixed=fx,dim=d,p.order=m,by=by.var,label=label,xt=xt)
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
      if (!inherits(object$margin[[i]],c("cs.smooth","cr.smooth","cyclic.smooth"))) { # these classes already optimal
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
  ## deal with possible extra arguments of "tp" type smooth
  xtra <- list()
  object$xt <- as.list(object$xt)
  if (is.null(object$xt$max.knots)) xtra$max.knots <- 3000 
  else xtra$max.knots <- object$xt$max.knots 
  if (is.null(object$xt$seed)) xtra$seed <- 1 
  else xtra$seed <- object$xt$seed 

  ## now collect predictors
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
  ## deal with possibility of large data set
  if (nk==0 && n>xtra$max.knots) { ## then there *may* be too many data  
    xu <- uniquecombs(matrix(x,n,object$dim)) ## find the unique `locations'
    nu <- nrow(xu)  ## number of unique locations
    if (nu>xtra$max.knots) { ## then there is really a problem 
      seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed
      kind <- RNGkind(NULL)
      RNGkind("default","default")
      set.seed(xtra$seed) ## ensure repeatability
      nk <- xtra$max.knots ## going to create nk knots
      ind <- sample(1:nu,nk,replace=FALSE)  ## by sampling these rows from xu
      knt <- as.numeric(xu[ind,])  ## ... like this
      RNGkind(kind[1],kind[2])
      assign(".Random.seed",seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
    }
  } ## end of large data set handling
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
               UZ=as.double(UZ),Xu=as.double(Xu),n.Xu=as.integer(nXu),C=as.double(C))
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
  
  if (object$bs.dim <3) { object$bs.dim <- 3
    warning("basis dimension, k, increased to minimum possible\n")
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
  object$C <- matrix(colSums(object$X),1,ncol(object$X))
##  object$C <- matrix(oo[[7]],1,nk)  # constraint
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
##  nx<-length(x)

  if (object$bs.dim <4) { object$bs.dim <- 4
    warning("basis dimension, k, increased to minimum possible\n")
  }

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
    } else {
      attr(sm,"qrc") <- "no constraints"
      attr(sm,"nCons") <- 0;
    } 
  } else attr(sm,"qrc") <-NULL

  sm 
}

PredictMat <- function(object,data)
## wrapper function which calls Predict.matrix and imposes same constraints as 
## smoothCon on resulting Prediction Matrix
{ X <- Predict.matrix(object,data)
  offset <- attr(X,"offset")
  qrc <- attr(object,"qrc")
  if (!is.null(qrc)) { ## then smoothCon absorbed constraints
    j <- attr(object,"nCons")
    if (j>0) { ## there were constraints to absorb - need to untransform
      k<-ncol(X)
      if (sum(is.na(X))) {
        ind <- !is.na(rowSums(X))
        X1 <- t(qr.qy(qrc,t(X[ind,]))[(j+1):k,])
        X <- matrix(NA,nrow(X),ncol(X1))
        X[ind,] <- X1
      } else {
        X <- t(qr.qy(qrc,t(X))[(j+1):k,])
      }
    }
  }
  ## drop columns eliminated by side-conditions...
  del.index <- attr(object,"del.index") 
  if (!is.null(del.index)) X <- X[,-del.index]
  attr(X,"offset") <- offset
  X
}


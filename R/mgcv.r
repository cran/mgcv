# These are the R routines for the package mgcv (c) Simon Wood 2000-2004



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
  oo<-.C("RMonoCon",as.double(A),as.double(b),as.double(x),as.integer(control),as.double(lower),
         as.double(upper),as.integer(n),PACKAGE="mgcv")
  A<-matrix(oo[[1]],dim(A)[1],dim(A)[2])
  b<-array(oo[[2]],dim(A)[1])
  list(A=A,b=b)
}  


uniquecombs<-function(x) {
# takes matrix x and counts up unique rows
if (is.null(x)) stop("x is null")
if (is.null(nrow(x))) stop("x has no row attribute")
if (is.null(ncol(x))) stop("x has no col attribute")
res<-.C("RuniqueCombs",as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),PACKAGE="mgcv")
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


null.space.basis.powers<-function(m,d)
# generates the sequence of powers required to specify the M polynomials spanning the 
# null space of the penalty of a d-dimensional tps with wiggliness penalty order d 
# So, if x_i are the co-ordinates the kth polynomial is x_1^pi[k][1]*x_2^pi[k][2] ....
{ if (d<0) stop("d can not be negative in call to null.space.basis.powers()")
  M<-null.space.dimension(d=d,m=m)
  if (2*m<=d) {m<-1;while (2*m<d+2) m<-m+1;}
  index<-array(0,d)
  pi<-array(0,c(M,d))
  for (i in 1:M)
  { # copy index to pi 
    pi[i,]<-index
    # update index.... 
    sumi<-sum(index)
    if (sumi<m-1) # then increase lowest index 
    index[1]<-index[1]+1
    else          # pass the problem up 
    { sumi<-sumi-index[1];
      index[1]<-0;
      for (j in 2:d)
      { index[j]<-index[j]+1;sumi<-sumi+1;
        if (sumi==m) { sumi<-sumi-index[j];index[j]<-0;}
        else break; # problem resolved! 
      } 
    }
  }
  rm(index)
  pi
}

null.space.basis.labels<-function(names,m,by="NA")
# given the list of covariates of one spline term, and the order of the penalty for the
# term this function returns a list of names for the polynomial terms spanning the null
# space of the penalty, in the order in which the terms appear in the model design matrix.
# The routine ensures that the same term always gets the same label, irrespective of the
# ordering of the input names list.... this lets the labels be used for testing for
# redundancy in a GAM design matrix. If a names is supplied in `by' then it  multiplies
# all terms in the basis, except the constant. 
{ d<-length(names[names!=""])
  if (length(unique(names))<d) stop("Repeated variables in a smooth term (null space basis)")
  pi<-null.space.basis.powers(m=m,d=d)
  M<-dim(pi)[1]
  term<-array("",M)
  for (i in 1:M)
  { label<-""
    for (j in 1:d)
    if (pi[i,j]==0) label[j]<-""
    else if (pi[i,j]==1) label[j]<-names[j]
    else label[j]<-paste(names[j],"^",pi[i,j],collapse="",sep="") 
    label<-label[label!=""]
    if (length(label))
    { if (by!="NA") label[length(label)+1]<-by
      label<-sort(label)
      term[i]<-paste(label,collapse="*",sep="*")
    }
    else term[i]<-"1"
  }  
  term
}  



gam.side.conditions<-function(G)
# This routine is designed to:
#  i) check the smooth part of a gam for identifiability
#  ii) impose side conditions to impose identifiability, where possible
# some by-variable lack of identifiability may be missed. Smooths remain centred even with by-vars.
{ gam.repeat.check<-function(v.names,by.names)
  # takes a 2-d array of variable names and checks for repeated rows (order unimportant)
  { n<-dim(v.names)[1];m<-dim(v.names)[2]
    if (n<=1) return(FALSE) # nothing to do
    for (i in 1:n) 
    { vn<-v.names[i,]<-sort(v.names[i,])
      vn<-vn[vn!=""];
      if (length(unique(vn))!=length(vn)) return(TRUE) # there is a within term repeat
    }
    for (i in 1:(n-1))
    { vn<-c(v.names[i,],by.names[i])
      for (j in (i+1):n) 
      if (sum(vn==c(v.names[j,],by.names[j]))==(m+1)) return(TRUE) # there is a repeat in this list of terms
    }
    return(FALSE) # no repeats
  }
  # end of local function declaration, start of main code.....
  if (G$m==0) return(G$C) # nothing to do
  v.names<-array("",0)
  ok.classes<-c("cr.smooth","tprs.smooth") # classes of smooths that can be handled here
  max.dim<-0
  for (i in 1:G$m) # collect all smoothed terms except from tensor product smooths
  { if (sum(class(G$smooth[[i]])==ok.classes)) # then this term is one of the class of smooths that can be handled
    { v.names<-c(v.names,G$smooth[[i]]$term)
      max.dim<-max(G$smooth[[i]]$dim,max.dim)
    }
  }
  n.vars<-length(v.names)
  if (length(unique(v.names))==length(v.names))
  return(G$C) # nothing to do - there are no repeated variable names in the spline part
  # first assemble the array of term labels
  term.names<-array("",c(G$m,max.dim))
  by.names<-rep("NA",G$m) 
  off<-p.order<-bs.dim<-dim<-array(-1,G$m)
 
  for (i in 1:G$m)   # term[i,j] holds jth variable name for ith term
  { if (sum(class(G$smooth[[i]])==ok.classes))
    { dim[i]<-G$smooth[[i]]$dim
      term.names[i,1:dim[i]]<-G$smooth[[i]]$term
      by.names[i]<-G$smooth[[i]]$by
      off[i]<-G$smooth[[i]]$first.para
      p.order[i]<-G$smooth[[i]]$p.order
      bs.dim[i]<-G$smooth[[i]]$bs.dim
    }
  }  
  # now work through the term dimensions
  all.ind<-1:G$m # index vector
  var.stack<-array("",0)  # stores names for null space basis vectors used so far - checked to avoid repeats
  stack.size<-0  # size of var.stack

  for (d in 1:max.dim) # note if d=1 were not first then would have to deal with "cr" terms specially 
  { # assemble all terms of dimension d, need variable names, p.order and parameter vector offsets
    ind<-all.ind[dim==d]
    n.d<-length(ind)
    if (n.d>0) # then there is something to do
    { d.vnames<-term.names[ind,]
      dim(d.vnames)<-c(n.d,dim(term.names)[2]) # otherwise vectors lose dimension attribute!
      d.by.names<-by.names[ind]
      d.off<-off[ind]
      d.dim<-dim[ind]
      d.p.order<-p.order[ind] # order of penalty 
      d.k<-bs.dim[ind] # basis dimension
      # check these terms for repeats within and between terms
      if(gam.repeat.check(v.names=d.vnames,by.names=d.by.names)) 
      stop("The smooth part of this model is not identifiable")
      # work through the terms imposing constraints on redundant null basis  
      # parameters and updating the null basis stack.....
      for (j in 1:n.d) # work through all n.d d-dimensional terms 
      { #cat(d.vnames[j,]," ",d.p.order[j]," ")
        by<-d.by.names[j]
        bs.label<-null.space.basis.labels(d.vnames[j,],d.p.order[j],by=by) # get unique null basis vector names
        bs.label<-bs.label[bs.label!="1"] # constants already constrained by sum to zero constraints
        #print(bs.label)
        for (i in 1:length(bs.label)) 
        if (length(var.stack[var.stack==bs.label[i]])>0) # already on the stack
        { # locate the parameter for this term 
          k<-d.off[j]+d.k[j]-null.space.dimension(d.dim[j],d.p.order[j])+i
          rows<-dim(G$C)[1];cols<-dim(G$C)[2]
          CC<-matrix(0,rows+1,cols)
          CC[1:rows,1:cols]<-G$C
          CC[rows+1,k]<-1
          G$C<-CC;rm(CC)  # parameter constrained to be zero
        } else # pop it on the stack
        { stack.size<-stack.size+1
          var.stack[stack.size]<-bs.label[i]
        } 
      }
    }
  }   
  G$C # the final constraint matrix
}  



pcls<-function(M)
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
   
  # pack the S array for mgcv call 
  m<-length(M$S)
  Sa<-array(0,0);df<-0
  if (m>0) for (i in 1:m)
  { Sa<-c(Sa,M$S[[i]])
    df[i]<-nrow(M$S[[i]])
  }

  o<-.C("RPCLS",as.double(M$X),as.double(M$p),as.double(M$y),as.double(M$w),as.double(M$Ain),as.double(M$bin)
        ,as.double(M$C),as.double(H),as.double(Sa),as.integer(M$off),as.integer(df),as.double(M$sp),
        as.integer(length(M$off)),as.integer(nar),PACKAGE="mgcv")
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
#          autonitialization guess ok (or intial values supplied), step.fail - TRUE
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

  oo<-.C("mgcv",as.double(y),as.double(X),as.double(C),as.double(w^2),as.double(Sa),
         as.double(p),as.double(sp),as.integer(off-1),as.integer(df),as.integer(m),
         as.integer(n),as.integer(q),as.integer(C.r),as.double(scale),as.double(Vp),
		 as.double(edf),as.double(control$conv.tol),as.integer(control$max.half),as.double(ddiag),
                 as.integer(idiag),as.double(sdiag),as.integer(direct.mesh),as.double(control$min.edf),
                 as.double(gcv.ubre),as.double(control$target.edf),as.integer(fixed.sp),
                 as.double(hat),PACKAGE="mgcv")
   
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

tensor.prod.penalties<-function(S)
# Given a list S of penalty matrices for the marginal bases of a tensor product smoother
# this routine produces the resulting penalties for the tensor product basis. 
# e.g. if S_1, S_2 and S_3 are marginal penalties and I_1, I_2, I_3 are identity matrices 
# of the same dimensions then the tensor product penalties are:
#   S_1 %x% I_2 %x% I_3, I_1 %x% S_2 %x% I_3 and I_1 %*% I_2 %*% S_3
# Note that the penalty list must be in the same order as the model matrix list supplied
# to tensor.prod.model() when using these together.
{ m<-length(S)
  I<-list(); for (i in 1:m) I[[i]]<-diag(ncol(S[[i]]))
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
# fails tries evaluating txt within data (only)
{ x <- data[[txt]]
  if (is.null(x)) 
  { x <- try(eval(parse(text=txt),data,enclos=NULL),silent=TRUE)
    if (inherits(x,"try-error")) x <- NULL
  }
  x
}



te<-function(..., k=NA,bs="cr",m=0,d=NA,by=NA,fx=FALSE,mp=TRUE)
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
  bs[d>1]<-"tp"
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
  full.call<-paste("te(",term[1],sep="")
  if (dim>1) for (i in 2:dim) full.call<-paste(full.call,",",term[i],sep="")
  label<-paste(full.call,")",sep="")   # label for parameters of this term
  full.call<-paste(full.call,",k=",deparse(k,backtick=TRUE),",bs=",deparse(bs,backtick=TRUE),
                   ",m=",deparse(m,backtick=TRUE),",d=",deparse(d,backtick=TRUE),
                   ",by=",by.var,",fx=",deparse(fx,backtick=TRUE),",mp=",deparse(mp,backtick=TRUE),")",sep="")
  ret<-list(margin=margin,term=term,by=by.var,fx=fx,full.call=full.call,label=label,dim=dim,mp=mp)
  class(ret)<-"tensor.smooth.spec"
  ret
}





s<-function (..., k=-1,fx=FALSE,bs="tp",m=0,by=NA)
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
  term<-deparse(vars[[1]],backtick=TRUE) # first covariate
  if (d>1) # then deal with further covariates
  for (i in 2:d)
  { term[i]<-deparse(vars[[i]],backtick=TRUE)
  }
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
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    X<-by*X
  }
  C<-matrix(colSums(X),1,ncol(X))
  object$X<-X;object$S<-S;object$C<-C
  object$df<-ncol(X)-1
  object$null.space.dim <- prod(nr) # penalty null space rank 
  object$rank<-r
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
  oo<-.C("construct_tprs",as.double(x),as.integer(object$dim),as.integer(n),as.double(knt),as.integer(nk),
               as.integer(object$p.order),as.integer(object$bs.dim),X=as.double(X),S=as.double(S),
               UZ=as.double(UZ),Xu=as.double(Xu),n.Xu=as.integer(nXu),C=as.double(C),PACKAGE="mgcv")
  object$X<-matrix(oo$X,n,k)                   # model matrix
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    object$X<-by*object$X
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
  X <- rep(0,nx*nk);S<-rep(0,nk*nk);C<-rep(0,nk);control<-0
  
  if (length(unique(x))<nk) 
  { msg <- paste(object$term," has insufficient unique values to support ",
                 nk," knots: reduce k.",sep="")
    stop(msg)
  }

  oo <- .C("construct_cr",as.double(x),as.integer(nx),as.double(k),
           as.integer(nk),as.double(X),as.double(S),
           as.double(C),as.integer(control),PACKAGE="mgcv")

  object$X <- matrix(oo[[5]],nx,nk)
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    object$X<-by*object$X
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
    object$X<-by*object$X
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
  tensor.prod.model.matrix(X)
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
  X <- pred.mat(x,object$xp,object$BD)
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    X<-by*X
  }
  X
}


Predict.matrix.cr.smooth<-function(object,data)
# this is the prediction method for a cubic regression spline
{ x <- get.var(object$term,data)
  nx<-length(x)
  nk<-object$bs.dim
  X <- rep(0,nx*nk);S<-rep(0,nk*nk);C<-rep(0,nk);control<-0

  oo <- .C("construct_cr",as.double(x),as.integer(nx),as.double(object$xp),
            as.integer(object$bs.dim),as.double(X),as.double(S),
                   as.double(C),as.integer(control),PACKAGE="mgcv")
  X<-matrix(oo[[5]],nx,nk) # the prediction matrix
  if (object$by!="NA")  # deal with "by" variable 
  { by <- get.var(object$by,data)
    if (is.null(by)) stop("Can't find by variable")
    X<-by*X
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
  oo<-.C("predict_tprs",as.double(x),as.integer(object$dim),as.integer(n),as.integer(object$p.order),
      as.integer(object$bs.dim),as.integer(object$null.space.dim),as.double(object$Xu),
      as.integer(nrow(object$Xu)),as.double(object$UZ),as.double(by),as.integer(by.exists),X=as.double(X)
      ,PACKAGE="mgcv")
  X<-matrix(oo$X,n,object$bs.dim)
}

Predict.matrix.ts.smooth<-function(object,data)
# this is the prediction method for a t.p.r.s
# with shrinkage
{ Predict.matrix.tprs.smooth(object,data)
}

smooth.construct<-function(object,data,knots) UseMethod("smooth.construct")


Predict.matrix<-function(object,data) UseMethod("Predict.matrix")


### the following two functions are for use in place of log and exp
### in positivity ensuring re-parameterization.... they have much better 
### over/underflow characteristics, but are still continuous to second
### derivative. 

notExp <- function(x)
# overflow avoiding C2 function for ensuring positivity
{ f <- x
  ind <- x > 1
  f[ind] <- exp(1)*(x[ind]^2+1)/2
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  x[ind] <- -x[ind] ;f[ind] <-  exp(1)*(x[ind]^2+1)/2; f[ind]<-1/f[ind]
  f
}


notLog <- function(x)
# inverse function of notExp
{ f <- x
  ind <- x> exp(1)
  f[ind] <- sqrt(2*x[ind]/exp(1)-1)
  ind <- !ind & x > exp(-1)
  f[ind] <- log(x[ind])
  ind <- x <= exp(-1)
  x[ind]<- 1/x[ind]; f[ind] <- sqrt(2*x[ind]/exp(1)-1);f[ind] <- -f[ind]
 f
}




#### pdMat class definitions, to enable tensor product smooths to be employed with gamm()
#### Based on various Pinheiro and Bates pdMat classes.

pdTens <- function(value = numeric(0), form = NULL, nam = NULL, data = sys.frame(sys.parent()))
## Constructor for the pdTens pdMat class: 
# the inverse of the scaled random effects covariance matrix for this class
# is given by a weighted sum of the matrices in the list that is the "S" attribute of 
# a pdTens formula. The weights are the exponentials of the class parameters.
# i.e. the inverse of the r.e. covariance matrix is 
#   \sum_i \exp(\theta_i) S_i / \sigma^2 
# The class name relates to the fact that these objects are used with tensor product smooths.  
{
  object <- numeric(0)
  class(object) <- c("pdTens", "pdMat")
  pdConstruct(object, value, form, nam, data)
}

## Methods for local generics


pdConstruct.pdTens <-
  function(object, value = numeric(0), form = formula(object),
	   nam = Names(object), data = sys.frame(sys.parent()), ...)
# used to initialize pdTens objects. Note that the initialization matrices supplied
# are (factors of) trial random effects covariance matrices and not their inverses
{
  val <- NextMethod()
  if (length(val) == 0) {               # uninitiliazed object
     if ((ncol <- length(Names(val))) > 0) {
      attr(val, "ncol") <- ncol
    }
    return(val)
  }
  if (is.matrix(val)) {			# initialize from a positive definite
    S <- attr(form,"S")
    m <- length(S)
    y <- as.numeric(solve(crossprod(val)))     # it's a factor that gets passed in
    lform <- "y ~ as.numeric(S[[1]])"
    for (i in 2:m) lform <- paste(lform," + as.numeric(S[[",i,"]])",sep="")
    lform <- formula(paste(lform,"-1"))
    value <- coef(lm(lform))  
    value[value <=0] <- mean(abs(value))/1000
    value <- notLog(value)
    attributes(value) <- attributes(val)[names(attributes(val)) != "dim"]
    class(value) <- c("pdTens", "pdMat")
    return(value)
  }
  m <- length(attr(form,"S"))
  if ((aux <- length(val)) > 0) {
    if (aux && (aux != m)) {
      stop(paste("An object of length", aux,
		 "does not match the required parameter size"))
    }
  }
  val
}


pdFactor.pdTens <-
  function(object)
# The factor of the inverse of the scaled r.e. covariance matrix is returned here
{ sp <- as.vector(object)
  m <- length(sp)
  S <- attr(formula(object),"S")
  value <- S[[1]]*notExp(sp[1])
  if (m>1) for (i in 2:m) value <- value + notExp(sp[i])*S[[i]] 
  if (sum(is.na(value))>0) warning("NA's in pdTens factor")
  value <- (value+t(value))/2
  mroot(value,rank=nrow(value))
}


pdMatrix.pdTens <-
  function(object, factor = FALSE) 
# the inverse of the scaled random effect covariance matrix is returned here, or
# its factor if factor==TRUE
{
  if (!isInitialized(object)) {
    stop("Cannot extract the matrix from an uninitialized object")
  }
  sp <- as.vector(object)
  m <- length(sp)
  S <- attr(formula(object),"S")
  value <- S[[1]]*notExp(sp[1])   
  if (m>1) for (i in 2:m) value <- value + notExp(sp[i])*S[[i]]  
  value <- (value + t(value))/2 # ensure symmetry
  if (sum(is.na(value))>0) warning("NA's in pdTens matrix")
  if (factor) {
    value <- mroot(value,rank=nrow(value))
    attr(value, "logDet") <- sum(log(diag(value)))
  } 
  dimnames(value) <- attr(object, "Dimnames")
  value
}

#### Methods for standard generics

coef.pdTens <-
  function(object, unconstrained = TRUE, ...)
{
  if (unconstrained) NextMethod()
  else {
    val <- notExp(as.vector(object))
    names(val) <- paste("sp.",1:length(val), sep ="")
    val
  }
}

summary.pdTens <-
  function(object, structName = "Tensor product smooth term", ...)
{
  # summary.pdMat(object, structName, noCorrelation = TRUE)

  ## ... summary.pdMat is not exported in the nlme NAMESPACE file, so....

  NextMethod(object, structName, noCorrelation=TRUE)
}


# .... end of pdMat definitions for tensor product smooths


### pdIdnot: multiple of the identity matrix - the parameter is
### the notLog of the multiple. This is directly modified form 
### Pinheiro and Bates pdIdent class.

####* Constructor

pdIdnot <-
  ## Constructor for the pdIdnot class
  function(value = numeric(0), form = NULL, nam = NULL, data = sys.frame(sys.parent()))
{
  object <- numeric(0)
  class(object) <- c("pdIdnot", "pdMat")
  pdConstruct(object, value, form, nam, data)
}

####* Methods for local generics

corMatrix.pdIdnot <-
  function(object, ...)
{
  if (!isInitialized(object)) {
    stop("Cannot extract the matrix from an uninitialized pdMat object")
  }
  if (is.null(Ncol <- attr(object, "ncol"))) {
    stop(paste("Cannot extract the matrix with uninitialized dimensions"))
  }
  val <- diag(Ncol)
  attr(val, "stdDev") <- rep(notExp(as.vector(object)), Ncol)
  if (length(nm <- Names(object)) == 0) {
    nm <- paste("V", 1:len, sep = "")
    dimnames(val) <- list(nm, nm)
  }
  names(attr(val, "stdDev")) <- nm
  val
}

pdConstruct.pdIdnot <-
  function(object, value = numeric(0), form = formula(object),
	   nam = Names(object), data = sys.frame(sys.parent()), ...)
{
  val <- NextMethod()
  if (length(val) == 0) {			# uninitialized object
    if ((ncol <- length(Names(val))) > 0) {
      attr(val, "ncol") <- ncol
    }
    return(val)
  }
  if (is.matrix(val)) {
    value <- notLog(sqrt(mean(diag(crossprod(val)))))
    attributes(value) <- attributes(val)[names(attributes(val)) != "dim"]
    attr(value, "ncol") <- dim(val)[2]
    class(value) <- c("pdIdnot", "pdMat")
    return(value)
  }
  if (length(val) > 1) {
    stop(paste("An object of length", length(val),
	       "does not match the required parameter size"))
  }
  if (((aux <- length(Names(val))) == 0) && is.null(formula(val))) {
    stop(paste("Must give names when initializing pdIdnot from parameter.",
	       "without a formula"))
  } else {
    attr(val, "ncol") <- aux
  }
  val
}

pdFactor.pdIdnot <-
  function(object)
{
  notExp(as.vector(object)) * diag(attr(object, "ncol"))
}


pdMatrix.pdIdnot <-
  function(object, factor = FALSE)
{
  if (!isInitialized(object)) {
    stop("Cannot extract the matrix from an uninitialized pdMat object")
  }
  if (is.null(Ncol <- attr(object, "ncol"))) {
    stop(paste("Cannot extract the matrix with uninitialized dimensions"))
  }
  value <- diag(Ncol)
  
  if (factor) {
   
    value <- notExp(as.vector(object)) * value
     attr(value, "logDet") <- Ncol * log(notExp(as.vector(object)))
   
  } else {
    value <- notExp(as.vector(object))^2 * value
  }
  dimnames(value) <- attr(object, "Dimnames")
  value
}

####* Methods for standard generics

coef.pdIdnot <-
  function(object, unconstrained = TRUE, ...)
{
  if (unconstrained) NextMethod()
  else structure(notExp(as.vector(object)),
           names = c(paste("sd(", deparse(formula(object)[[2]],backtick=TRUE),")",sep = "")))
}

Dim.pdIdnot <-
  function(object, ...)
{
  if (!is.null(val <- attr(object, "ncol"))) {
    c(val, val)
  } else {
    stop("Cannot extract the dimensions")
  }
}

logDet.pdIdnot <-
  function(object, ...)
{
  attr(object, "ncol") * log(notExp(as.vector(object)))
}

solve.pdIdnot <-
  function(a, b, ...)
{
  if (!isInitialized(a)) {
    stop("Cannot extract the inverse from an uninitialized object")
  }
  coef(a) <- -coef(a, TRUE)
  a
}

summary.pdIdnot <-
  function(object, structName = "Multiple of an Identity", ...)
{
  # summary.pdMat(object, structName, noCorrelation = TRUE)

  ## ... summary.pdMat is not exported in the nlme NAMESPACE file, so....

  NextMethod(object, structName, noCorrelation=TRUE)
}

### end of pdIdnot class


gamm.setup<-function(formula,pterms,data=stop("No data supplied to gam.setup"),knots=NULL,parametric.only=FALSE)
# set up the model matrix, penalty matrices and auxilliary information about the smoothing bases
# needed for a gamm fit.

# NOTE: lme can't deal with offset terms.
# There is an implicit assumption that any rank deficient penalty does not penalize 
# the constant term in a basis. 
{ # split the formula if the object being passed is a formula, otherwise it's already split
  if (inherits(formula,"formula")) split<-interpret.gam(formula) 
  else if (inherits(formula,"split.gam.formula")) split<-formula
  else stop("First argument is no sort of formula!") 
  if (length(split$smooth.spec)==0)
  { if (split$pfok==0) stop("You've got no model....")
    m<-0
  }  
  else  m<-length(split$smooth.spec) # number of smooth terms
  
  G<-list(m=m,full.formula=split$full.formula)

  if (is.null(attr(data,"terms"))) # then data is not a model frame
  mf<-model.frame(split$pf,data,drop.unused.levels=FALSE) # FALSE or can end up with wrong prediction matrix!
  else mf<-data # data is already a model frame
  
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
 
  G$off<-array(0,0)
  first.f.para<-G$nsdf+1
  first.r.para<-1
 
  G$Xf <- X # matrix in which to accumulate full GAM model matrix, treating
            # smooths as fixed effects
  random<-list()
  random.i<-0
  k.sp <- 0  # counter for penalties
  if (m)
  for (i in 1:m) 
  { # idea here is that terms are set up in accordance with information given in split$smooth.spec
    # appropriate basis constructor is called depending on the class of the smooth
    # constructor returns penalty matrices model matrix and basis specific information
    sm<-smooth.construct(split$smooth.spec[[i]],data,knots)
    # incorporate any constraints, and store null space information in smooth object
    if (length(sm$S)>1)
    { if (sum(sm$fx)==length(sm$fx)) sm$fixed <- TRUE
      else 
      { sm$fixed <- FALSE
        if (sum(sm$fx)!=0) warning("gamm can not fix only some margins of tensor product.")
      }
    } 
    if (!sm$fixed) random.i <- random.i+1
    k<-ncol(sm$X)
    j<-nrow(sm$C)
    ZSZ <- list()
    if (nrow(sm$C)) # there are constraints
    { qrc<-qr(t(sm$C))
      if (!sm$fixed)
      for (l in 1:length(sm$S)) # tensor product terms have > 1 penalty 
      { ZSZ[[l]]<-qr.qty(qrc,sm$S[[l]])[(j+1):k,]
        ZSZ[[l]]<-t(qr.qty(qrc,t(ZSZ[[l]]))[(j+1):k,])
        k.sp <- k.sp+1
      }
      XZ<-t(qr.qy(qrc,t(sm$X))[(j+1):k,])
      sm$qrc<-qrc
    } else 
    { if (!sm$fixed) 
      for (l in 1:length(sm$S)) ZSZ[[l]]<-sm$S[[l]]
      XZ<-sm$X
    }
    G$Xf<-cbind(G$Xf,XZ)   # accumulate model matrix that treats all smooths as fixed
    if (!sm$fixed) 
    { sm$ZSZ <- ZSZ            # store these too - for construction of Vp matrix
      if (length(ZSZ)>1) # tensor product term - need to find null space from sum of penalties
      { sum.ZSZ <- ZSZ[[1]]/mean(abs(ZSZ[[1]]))
        null.rank <- sm$margin[[1]]$bs.dim-sm$margin[[1]]$rank
        for (l in 2:length(ZSZ)) 
        { sum.ZSZ <- sum.ZSZ + ZSZ[[l]]/mean(abs(ZSZ[[l]]))
          null.rank <- # the rank of the null space of the penalty 
                       null.rank * (sm$margin[[l]]$bs.dim-sm$margin[[l]]$rank)
        }
        sum.ZSZ <- (sum.ZSZ+t(sum.ZSZ))/2 # ensure symmetry
        ev <- eigen(sum.ZSZ,symmetric=TRUE)
        mult.pen <- TRUE
      } else            # regular s() term
      { ZSZ[[1]] <- (ZSZ[[1]]+t(ZSZ[[1]]))/2
        ev<-eigen(ZSZ[[1]],symmetric=TRUE)
        null.rank <- sm$bs.dim - sm$rank
        mult.pen <- FALSE
      }
      p.rank <- ncol(sm$X) - null.rank
      if (p.rank>ncol(XZ)) p.rank <- ncol(XZ)
      U<-ev$vectors
      D<-1/sqrt(ev$values[1:p.rank])
      if (sum(is.na(U))||sum(is.na(D))) stop(
      "NA returned from eigen: please email simon@stats.gla.ac.uk with as much detail as possible ")
      XZU<-XZ%*%U
      if (p.rank<k-j) Xf<-XZU[,(p.rank+1):(k-j)]
      else Xf<-matrix(0,nrow(sm$X),0) # no fixed terms left
      if (mult.pen) 
      { Xr <- XZU[,1:p.rank] # tensor product case
        for (l in 1:length(ZSZ))   # transform penalty explicitly
        ZSZ[[l]] <- (t(U)%*%ZSZ[[l]]%*%U)[1:p.rank,1:p.rank]
      }
      else Xr<-t(t(XZU[,1:p.rank])*D)
      n.para<-k-j-p.rank # indices for fixed parameters
      sm$first.f.para<-first.f.para
      first.f.para<-first.f.para+n.para
      sm$last.f.para<-first.f.para-1
      n.para<-ncol(Xr) # indices for random parameters
      sm$first.r.para<-first.r.para
      first.r.para<-first.r.para+n.para
      sm$last.r.para<-first.r.para-1
    
      sm$D<-D;sm$U<-U # information (with qrc) for backtransforming to original space 

      term.name <- paste("Xr.",random.i,sep="")
      term.name <- new.name(term.name,names(data))
      form <- as.formula(paste("~",term.name,"-1",sep=""))
      if (mult.pen)  # tensor product case
      { attr(form,"S") <- ZSZ
        random[[random.i]] <- pdTens(form)
      } else  # single penalty smooth
      random[[random.i]] <- pdIdnot(form)
      names(random)[random.i] <- term.name
      eval(parse(text=paste("G$",term.name,"<-Xr",sep="")))
    } else # term is fixed, so model matrix appended to fixed matrix
    { Xf <- XZ # whole term goes to fixed 
      n.para <- ncol(Xf)       # now define where the parameters of this term live 
      sm$first.f.para <- first.f.para
      first.f.para <- first.f.para+n.para
      sm$last.f.para <- first.f.para-1
    }
    X<-cbind(X,Xf) # add fixed model matrix to overall X
  
    sm$X <- NULL
  
    G$smooth[[i]] <- sm   
  }
  G$m<-m
  G$random<-random
  G$X<-X  

  G$y <- data[[deparse(split$full.formula[[2]],backtick=TRUE)]]
  
  G$n<-nrow(data)

  if (is.null(data$"(weights)")) G$w<-rep(1,G$n)
  else G$w<-data$"(weights)"  

  # now run some checks on the arguments 

  ### Should check that there are enough unique covariate combinations to support model dimension

  G
}




extract.lme.cov2<-function(b,data,start.level=1)
# function to extract the response data covariance matrix from an lme fitted
# model object b, fitted to the data in data. "inner" == "finest" grouping 
# start.level is the r.e. grouping level at which to start the construction, 
# levels outer to this will not be included in the calculation - this is useful
# for gamm calculations
# 
# This version aims to be efficient, by not forming the complete matrix if it
# is diagonal or block diagonal. To this end the matrix is returned in a form
# that relates to the data re-ordered according to the coarsest applicable 
# grouping factor. ind[i] gives the row in the original data frame
# corresponding to the ith row/column of V. 
# V is either returned as an array, if it's diagonal, a matrix if it is
# a full matrix or a list of matrices if it is block diagonal.
{ if (!inherits(b,"lme")) stop("object does not appear to be of class lme")
  grps<-getGroups(b) # labels of the innermost groupings - in data frame order
  n<-length(grps)    # number of data
  n.levels <- length(b$groups) # number of levels of grouping
  if (n.levels >= start.level)
  { Cgrps <- getGroups(b,level=start.level) # labels for outer grouping (df order) 
    Cind <- sort(as.numeric(Cgrps),index.return=TRUE)$ix
    # Cind[i] is where row i of sorted Cgrps is in original data frame order 
    rCind <- 1:n; rCind[Cind] <- 1:n
    # rCind[i] is location of ith original datum in the coarse ordering
    CFgrps <- grps[Cind] # fine group levels in coarse group order
    Clevel <- levels(Cgrps) # levels of coarse grouping factor
    n.cg <- length(Clevel)  # number of outer groups
    size.cg <- array(0,n.cg)  
    for (i in 1:n.cg) size.cg[i] <- sum(Cgrps==Clevel[i]) # size of each coarse group
    ## Cgrps[Cind] is sorted by coarsest grouping factor level
    ## so e.g. y[Cind] would be data in c.g.f. order
  } else {n.cg <- 1;Cind<-1:n}
  if (is.null(b$modelStruct$varStruct)) w<-rep(b$sigma,n) ### 
  else 
  { w<-1/varWeights(b$modelStruct$varStruct) 
    # w is not in data.frame order - it's in inner grouping level order
    group.name<-names(b$groups) # b$groups[[i]] doesn't always retain factor ordering
    order.txt <- paste("ind<-order(data[[\"",group.name[1],"\"]]",sep="")
    if (length(b$groups)>1) for (i in 2:length(b$groups)) 
    order.txt <- paste(order.txt,",data[[\"",group.name[i],"\"]]",sep="")
    order.txt <- paste(order.txt,")")
    eval(parse(text=order.txt))
    w[ind] <- w # into data frame order
    w<-w*b$sigma
  }
  w <- w[Cind] # re-order in coarse group order
  if (is.null(b$modelStruct$corStruct)) V<-array(1,n) 
  else
  { c.m<-corMatrix(b$modelStruct$corStruct) # correlation matrices for each innermost group
    if (!is.list(c.m)) { # copy and re-order into coarse group order
      V <- c.m;V[Cind,] -> V;V[,Cind] -> V 
    } else { 
      V<-list()   # V[[i]] is cor matrix for ith coarse group
      ind <- list() # ind[[i]] is order index for V[[i]] 
      for (i in 1:n.cg) { 
        V[[i]] <- matrix(0,size.cg[i],size.cg[i]) 
        ind[[i]] <- 1:size.cg[i]
      }
      # Voff[i] is where, in coarse order data, first element of V[[i]]
      # relates to ... 
      Voff <- cumsum(c(1,size.cg)) 
      gr.name <- names(c.m) # the names of the innermost groups
      n.g<-length(c.m)   # number of innermost groups
      j0<-rep(1,n.cg) # place holders in V[[i]]'s
      ii <- 1:n
      for (i in 1:n.g) # work through innermost groups
      { # first identify coarse grouping
        Clev <- unique(Cgrps[grps==gr.name[i]])  # level for coarse grouping factor
        if (length(Clev)>1) stop("inner groupings not nested in outer!!")
        k <- (1:n.cg)[Clevel==Clev] # index of coarse group - i.e. update V[[k]] 
        # now need to get c.m into right place within V[[k]]
        j1<-j0[k]+nrow(c.m[[i]])-1
        V[[k]][j0[k]:j1,j0[k]:j1]<-c.m[[i]]
        ind1 <- ii[grps==gr.name[i]] 
        # ind1 is the rows of original data.frame to which c.m[[i]] applies 
        # assuming that data frame order is preserved at the inner grouping
        ind2 <- rCind[ind1] 
        # ind2 contains the rows of the coarse ordering to which c.m[[i]] applies
        ind[[k]][j0[k]:j1] <- ind2 - Voff[k] + 1
        # ind[k] accumulates rows within coarse group k to which V[[k]] applies
        j0[k]<-j1+1  
      }
      for (k in 1:n.cg) { # pasting correlations into right place in each matrix
        V[[k]][ind[[k]],]<-V[[k]];V[[k]][,ind[[k]]]<-V[[k]] 
      }
    }
  } 
  # now form diag(w)%*%V%*%diag(w), depending on class of V
  if (is.list(V)) # it's  a block diagonal structure
  { for (i in 1:n.cg)
    { wi <- w[Voff[i]:(Voff[i]+size.cg[i]-1)] 
      V[[i]] <- wi*t(wi*V[[i]])
    }
  } else
  if (is.matrix(V))
  { V <- w*t(w*V)
  } else # it's a diagonal matrix
  { V <- w^2*V
  }
  # ... covariance matrix according to fitted correlation structure in coarse
  # group order
  
  ## Now work on the random effects ..... 
  X<-list()
  grp.dims<-b$dims$ncol # number of Zt columns for each grouping level (inner levels first)
  # inner levels are first in Zt
  Zt<-model.matrix(b$modelStruct$reStruct,data)  # a sort of proto - Z matrix
  # b$groups and cov (defined below have the inner levels last)
  cov<-as.matrix(b$modelStruct$reStruct) # list of estimated covariance matrices (inner level last)
  i.col<-1
  Z<-matrix(0,n,0) # Z matrix
  if (start.level<=n.levels)
  { for (i in 1:(n.levels-start.level+1)) # work through the r.e. groupings inner to outer
    { # get matrix with columns that are indicator variables for ith set of groups...
      # groups has outer levels first 
      X[[1]] <- model.matrix(~b$groups[[n.levels-i+1]]-1,
                contrasts.arg=c("contr.treatment","contr.treatment"))
      # Get `model matrix' columns relevant to current grouping level...
      X[[2]] <- as.matrix(Zt[,i.col:(i.col+grp.dims[i]-1)])
      i.col <- i.col+grp.dims[i]
      # tensor product the X[[1]] and X[[2]] rows...
      Z <- cbind(Z,tensor.prod.model.matrix(X))
    } # so Z assembled from inner to outer levels
    # Now construct overall ranef covariance matrix
    Vr <- matrix(0,ncol(Z),ncol(Z))
    start <- 1
    for (i in 1:(n.levels-start.level+1))
    { k <- n.levels-i+1
      for (j in 1:b$dims$ngrps[i]) 
      { stop <- start+ncol(cov[[k]])-1
        Vr[start:stop,start:stop]<-cov[[k]]
        start <- stop+1
      }
    }
    Vr <- Vr*b$sigma^2
    ## Now re-order Z into coarse group order
    Z <- Z[Cind,]
    ## Now Z %*% Vr %*% t(Z) is block diagonal: if Z' = [Z1':Z2':Z3': ... ]
    ## where Zi contains th rows of Z for the ith level of the coarsest
    ## grouping factor, then the ith block of (Z Vr Z') is (Zi Vr Zi')
    if (n.cg == 1) { 
      if (is.matrix(V)) { 
        V <- V+Z%*%Vr%*%t(Z)
      } else V <- diag(V) + Z%*%Vr%*%t(Z) 
    } else { # V has a block - diagonal structure
      j0<-1
      Vz <- list()
      for (i in 1:n.cg) {
        j1 <- size.cg[i] + j0 -1
        Zi <- Z[j0:j1,]
        Vz[[i]] <- Zi %*% Vr %*% t(Zi) 
        j0 <- j1+1
      }
      if (is.list(V)) {
        for (i in 1:n.cg) V[[i]] <- V[[i]]+Vz[[i]] 
      } else { 
        j0 <-1
        for (i in 1:n.cg) {
          j1 <- size.cg[i] + j0 -1
          Vz[[i]] <- Vz[[i]] + diag(V[j0:j1])
          j0 <- j1+1
        }
        V <- Vz
      }
    }
  }
  list(V=V,ind=Cind)
}


extract.lme.cov<-function(b,data,start.level=1)
# function to extract the response data covariance matrix from an lme fitted
# model object b, fitted to the data in data. "inner" == "finest" grouping 
# start.level is the r.e. grouping level at which to start the construction, 
# levels outer to this will not be included in the calculation - this is useful
# for gamm calculations
{ if (!inherits(b,"lme")) stop("object does not appear to be of class lme")
  grps<-getGroups(b) # labels of the innermost groupings - in data frame order
  n<-length(grps)    # number of data
  if (is.null(b$modelStruct$varStruct)) w<-rep(b$sigma,n) ### 
  else 
  { w<-1/varWeights(b$modelStruct$varStruct) 
    # w is not in data.frame order - it's in inner grouping level order
    group.name<-names(b$groups) # b$groups[[i]] doesn't always retain factor ordering
    order.txt <- paste("ind<-order(data[[\"",group.name[1],"\"]]",sep="")
    if (length(b$groups)>1) for (i in 2:length(b$groups)) 
    order.txt <- paste(order.txt,",data[[\"",group.name[i],"\"]]",sep="")
    order.txt <- paste(order.txt,")")
    eval(parse(text=order.txt))
    w[ind] <- w # into data frame order
    w<-w*b$sigma
  }
  if (is.null(b$modelStruct$corStruct)) V<-diag(n) #*b$sigma^2
  else
  { c.m<-corMatrix(b$modelStruct$corStruct) # correlation matrices for each group
    if (!is.list(c.m)) V<-c.m
    else
    { V<-matrix(0,n,n)   # data cor matrix
      gr.name <- names(c.m) # the names of the groups
      n.g<-length(c.m)   # number of innermost groups
      j0<-1
      ind<-ii<-1:n
      for (i in 1:n.g) 
      { j1<-j0+nrow(c.m[[i]])-1
        V[j0:j1,j0:j1]<-c.m[[i]]
        ind[j0:j1]<-ii[grps==gr.name[i]]
        j0<-j1+1  
      }
      V[ind,]<-V;V[,ind]<-V # pasting correlations into right place in overall matrix
      # V<-V*b$sigma^2
    }
  }  
  V <- w*t(w*V) # diag(w)%*%V%*%diag(w)
  # ... covariance matrix according to fitted correlation structure
  X<-list()
  grp.dims<-b$dims$ncol # number of Zt columns for each grouping level (inner levels first)
  # inner levels are first in Zt
  Zt<-model.matrix(b$modelStruct$reStruct,data)  # a sort of proto - Z matrix
  # b$groups and cov (defined below have the inner levels last)
  cov<-as.matrix(b$modelStruct$reStruct) # list of estimated covariance matrices (inner level last)
  i.col<-1
  n.levels<-length(b$groups)
  Z<-matrix(0,n,0) # Z matrix
  if (start.level<=n.levels)
  { for (i in 1:(n.levels-start.level+1)) # work through the r.e. groupings inner to outer
    { # get matrix with columns that are indicator variables for ith set of groups...
      # groups has outer levels first 
      X[[1]] <- model.matrix(~b$groups[[n.levels-i+1]]-1,
                contrasts.arg=c("contr.treatment","contr.treatment"))
      # Get `model matrix' columns relevant to current grouping level...
      X[[2]] <- as.matrix(Zt[,i.col:(i.col+grp.dims[i]-1)])
      i.col <- i.col+grp.dims[i]
      # tensor product the X[[1]] and X[[2]] rows...
      Z <- cbind(Z,tensor.prod.model.matrix(X))
    } # so Z assembled from inner to outer levels
    # Now construct overall ranef covariance matrix
    Vr <- matrix(0,ncol(Z),ncol(Z))
    start <- 1
    for (i in 1:(n.levels-start.level+1))
    { k <- n.levels-i+1
      for (j in 1:b$dims$ngrps[i]) 
      { stop <- start+ncol(cov[[k]])-1
        Vr[start:stop,start:stop]<-cov[[k]]
        start <- stop+1
      }
    }
    Vr <- Vr*b$sigma^2
    V <- V+Z%*%Vr%*%t(Z)
  }
  V
}

formXtViX <- function(V,X)
## forms X'V^{-1}X as efficiently as possible given the structure of
## V (diagonal, block-diagonal, full)
{ X <- X[V$ind,] # have to re-order X according to V ordering
  if (is.list(V$V)) {     ### block diagonal case
    Z <- X
    j0<-1
    for (i in 1:length(V$V))
    { Cv <- chol(V$V[[i]])
      j1 <- j0+nrow(V$V[[i]])-1
      Z[j0:j1,]<-backsolve(Cv,X[j0:j1,],transpose=TRUE)
      j0 <- j1 + 1
    }
    res <- t(Z)%*%Z
  } else if (is.matrix(V$V)) { ### full matrix case
    Cv<-chol(V$V)
    Z<-backsolve(Cv,X,transpose=TRUE)
    res <- t(Z)%*%Z
  } else {                ### diagonal matrix case
    res <- t(X)%*%(X/as.numeric(V$V))
  }
  res
}


new.name <- function(proposed,old.names)
# finds a name based on proposed, that is not in old.names
# if the proposed name is in old.names then ".xx" is added to it 
# where xx is a number chosen to ensure the a unique name
{ prop <- proposed
  k <- 0
  while (sum(old.names==prop))
  { prop<-paste(proposed,".",k,sep="")
    k <- k + 1
  }
  prop
}



gamm <- function(formula,random=NULL,correlation=NULL,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,
               knots=NULL,control=lmeControl(niterEM=3),niterPQL=20,verbosePQL=TRUE,...)
# Routine to fit a GAMM to some data. Fixed and smooth terms are defined in the formula, but the wiggly 
# parts of the smooth terms are treated as random effects. The onesided formula random defines additional 
# random terms. correlation describes the correlation structure. This routine is basically an interface
# between the bases constructors provided in mgcv and the glmmPQL routine used to estimate the model.

# NOTE: need to fill out the gam object properly

{   if (!require(nlme)) stop("gamm() requires package nlme to be installed")
    if (!require(MASS)) stop("gamm() requires package MASS to be installed")
    # check that random is a named list
    if (!is.null(random))
    { if (is.list(random)) 
      { r.names<-names(random)
        if (is.null(r.names)) stop("random argument must be a *named* list.")
        else if (sum(r.names=="")) stop("all elements of random list must be named")
      }
      else stop("gamm() can only handle random effects defined as named lists")
      random.vars<-c(unlist(lapply(random, function(x) all.vars(formula(x)))),r.names)
    } else random.vars<-NULL
    if (!is.null(correlation))
    { cor.for<-attr(correlation,"formula")
      if (!is.null(cor.for))
      cor.vars<-all.vars(cor.for)
    } else cor.vars<-NULL
    # create model frame.....
    gp<-interpret.gam(formula) # interpret the formula 
    cl<-match.call() # call needed in gamm object for update to work
    mf<-match.call(expand.dots=FALSE)
  
    allvars <- c(cor.vars,random.vars)
    if (length(allvars))
    mf$formula<-as.formula(paste(paste(deparse(gp$fake.formula,backtick=TRUE),collapse=""),
                           "+",paste(allvars,collapse="+")))
    else mf$formula<-gp$fake.formula
    mf$correlation<-mf$random<-mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-
    mf$min.sp<-mf$H<-mf$gamma<-mf$fit<-mf$niterPQL<-mf$verbosePQL<-mf$G<-mf$...<-NULL
    mf$drop.unused.levels<-TRUE
    mf[[1]]<-as.name("model.frame")
    pmf <- mf
    mf <- eval(mf, parent.frame()) # the model frame now contains all the data 

    Terms <- attr(mf,"terms")    
  
    pmf$formula <- gp$pf
    pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part 

    pTerms <- attr(pmf,"terms")

    if (is.character(family)) family<-eval(parse(text=family))
    if (is.function(family)) family <- family()
    if (is.null(family$family)) stop("family not recognized")
  
    # now call gamm.setup 

    G<-gamm.setup(gp,pterms=pTerms,data=mf,knots=knots,parametric.only=FALSE)

    if (is.null(random)&&G$m==0) 
    stop("gamm models must have at least 1 smooth with unknown smoothing paraemter or at least one other random effect")

    g<-as.factor(G$y*0+1)

    offset.name <- attr(mf,"names")[attr(attr(mf,"terms"),"offset")]

    yname <- new.name("y",names(mf))
    eval(parse(text=paste("mf$",yname,"<-G$y",sep="")))
    Xname <- new.name("X",names(mf))
    eval(parse(text=paste("mf$",Xname,"<-G$X",sep="")))
    
    fixed.formula <- paste(yname,"~",Xname,"-1")
    if (length(offset.name)) 
    { fixed.formula <- paste(fixed.formula,"+",offset.name) 
    }
    fixed.formula <- as.formula(fixed.formula)
    
    n.sr <- length(G$random) # number of random smooths (i.e. s(...,fx=FALSE,...) terms)
    group.name<-rep("",n.sr)
    r.name <- names(G$random) 
    if (n.sr) for (i in 1:n.sr) # adding the constructed variables to the model frame avoiding name duplication
    { mf[[r.name[i]]] <- G[[r.name[i]]]
      group.name[i] <- new.name(paste("g.",i,sep=""),names(mf))
      eval(parse(text=paste("mf$",group.name[i]," <- g",sep="")))
    }
    ret<-list()
    rand<-G$random;names(rand)<-group.name  
    if (!is.null(random)) # add explicit random effects
    { r.m<-length(random)
      r.names<-c(names(rand),names(random))
      for (i in 1:r.m) rand[[G$m+i]]<-random[[i]]    
      names(rand)<-r.names
    }
 
    ### Actually do fitting ....
    if (family$family=="gaussian"&&family$link=="identity"&&
    length(offset.name)==0) reml.used <- TRUE else reml.used <- FALSE
    
    if (reml.used)
    { ret$lme<-lme(fixed.formula,random=rand,data=mf,correlation=correlation,control=control)
    } else
    { ret$lme<-glmmPQL(fixed.formula,random=rand,data=mf,family=family,correlation=correlation,
                       control=control,niter=niterPQL,verbose=verbosePQL)
    }

    ### .... fitting finished

    # now fake a gam object 
    
    object<-list(model=mf,formula=formula,smooth=G$smooth,nsdf=G$nsdf,family=family,
                 df.null=nrow(G$X),y=G$y,terms=Terms,pterms=pTerms,xlevels=G$xlevels,
                 contrasts=G$contrasts,assign=G$assign)
    # Transform  parameters back to the original space....
    bf<-as.numeric(ret$lme$coefficients$fixed)
    br<-as.numeric(unlist(ret$lme$coefficients$random))
    if (G$nsdf) p<-bf[1:G$nsdf] else p<-array(0,0)
    n.pen <- 0 # count up the penalties
    if (G$m>0) for (i in 1:G$m)
    { fx <- G$smooth[[i]]$fixed 
      first<-G$smooth[[i]]$first.f.para;last<-G$smooth[[i]]$last.f.para
      if (first <=last) beta<-bf[first:last] else beta<-array(0,0)
      if (fx) b <- beta 
      else # not fixed so need to undo transform of random effects etc. 
      { dum <- length(G$smooth[[i]]$S) # number of penalties for this term
        n.pen <- n.pen + dum
        if (dum > 1) mult.pen <- TRUE else mult.pen <- FALSE
        b<-br[G$smooth[[i]]$first.r.para:G$smooth[[i]]$last.r.para]     
        if (mult.pen) b <- c(b,beta) # multiple penalties not reduced to identity
        else b <- c(G$smooth[[i]]$D*b,beta) # single penalty case
        b<-G$smooth[[i]]$U%*%b 
      }
      nc <- nrow(G$smooth[[i]]$C) 
      if (nc) b <- qr.qy(G$smooth[[i]]$qrc,c(rep(0,nc),b))
      object$smooth[[i]]$first.para<-length(p)+1
      p<-c(p,b)
      object$smooth[[i]]$last.para<-length(p)
    }
 
    cov<-as.matrix(ret$lme$modelStruct$reStruct)
    var.param <- coef(ret$lme$modelStruct$reStruct)
    n.v <- length(var.param) 
    k <- 1
    if (G$m>0) for (i in 1:G$m) # var.param in reverse term order, but forward order within terms!!
    { n.sp <- length(object$smooth[[i]]$S) # number of s.p.s for this term 
      if (inherits(object$smooth[[i]],"tensor.smooth")&&n.sp>1) 
      object$sp[k:(k+n.sp-1)] <- notExp(var.param[(n.v-n.sp+1):n.v])*2
      else object$sp[k:(k+n.sp-1)] <- 1/notExp(var.param[(n.v-n.sp+1):n.v])^2
      k <- k + n.sp
      n.v <- n.v - n.sp

    }
   
    object$coefficients<-p
    
    V<-extract.lme.cov2(ret$lme,mf,n.sr+1) # the data covariance matrix, excluding smooths
    XVX <- formXtViX(V,G$Xf)
#    Cv<-chol(V)
#    X<-G$Xf
#    Z<-backsolve(Cv,X,transpose=TRUE)
    S<-matrix(0,ncol(G$Xf),ncol(G$Xf)) # penalty matrix
    first <- G$nsdf+1
    k <- 1
    if (G$m>0) for (i in 1:G$m) # Accumulate the total penalty matrix
    { n.para <- object$smooth[[i]]$last.para - object$smooth[[i]]$first.para + 1 - 
                nrow(G$smooth[[i]]$C)
      last <- first + n.para - 1 
      if (!object$smooth[[i]]$fixed)
      { for (l in 1:length(object$smooth[[i]]$ZSZ))
         { S[first:last,first:last] <- S[first:last,first:last] + 
                  object$smooth[[i]]$ZSZ[[l]]*object$sp[k]
          k <- k+1
        }
      }
      first <- last + 1 
    }
    S<-S/ret$lme$sigma^2 # X'V^{-1}X divided by \sigma^2, so should S be
#    Z<-t(Z)%*%Z # X'V^{-1}X # this was XVX
    Vb <- chol2inv(chol(XVX+S)) # covariance matrix - in constraint space
    # need to project out of constraint space
    Vp <- matrix(Vb[1:G$nsdf,],G$nsdf,ncol(Vb))
    X <- matrix(XVX[1:G$nsdf,],G$nsdf,ncol(XVX))
    first <- G$nsdf+1
    if (G$m >0) for (i in 1:G$m)
    { nc <- nrow(G$smooth[[i]]$C)
      last <- first+object$smooth[[i]]$df-1 # old code: nrow(object$smooth[[i]]$ZSZ)-1
      if (nc)
      { V <- rbind(matrix(0,nc,ncol(Vp)),Vb[first:last,])
        V <- qr.qy(G$smooth[[i]]$qrc,V)
        Vp <- rbind(Vp,V) # cov matrix
        V <- rbind(matrix(0,nc,ncol(Vp)),XVX[first:last,])
        V <- qr.qy(G$smooth[[i]]$qrc,V)
        X <- rbind(X,V)  # X'V^{-1}X
      }
      first <- last+1
    }
    Vb<-matrix(Vp[,1:G$nsdf],nrow(Vp),G$nsdf)
    Z <- matrix(X[,1:G$nsdf],nrow(X),G$nsdf)
    first<-G$nsdf+1
    if (G$m>0) for (i in 1:G$m)
    { last <- first + object$smooth[[i]]$df-1   # old code: nrow(object$smooth[[i]]$ZSZ)-1
      nc <- nrow(G$smooth[[i]]$C)
      if (nc)
      { V <- cbind(matrix(0,nrow(Vb),nc),Vp[,first:last])
        V <- t(qr.qy(G$smooth[[i]]$qrc,t(V)))
        Vb<-cbind(Vb,V) 
        V <- cbind(matrix(0,nrow(Vb),nc),X[,first:last])
        V <- t(qr.qy(G$smooth[[i]]$qrc,t(V)))
        Z<-cbind(Z,V) 
      }
      first <- last+1
    }
    
    # then get edf diag(Vb%*%Z)
    
    object$edf<-rowSums(Vb*t(Z))
    
    
    object$sig2 <- ret$lme$sigma^2
    if (reml.used) { object$fit.method <- "lme";object$method <- "REML"} 
    else { object$fit.method <- "glmmPQL";object$method <- "PQL"}

    if (!reml.used) Vb<-Vb*length(G$y)/(length(G$y)-G$nsdf)
    object$Vp<- Vb
    
    object$prior.weights <- weights
    class(object)<-"gam"
    object$full.formula <- G$full.formula

    object$fitted.values <- predict.gam(object,type="response")
    object$residuals <- residuals(ret$lme) #as.numeric(G$y) - object$fitted.values

    if (G$nsdf>0) term.names<-colnames(G$X)[1:G$nsdf] else term.names<-array("",0)
    n.smooth<-length(G$smooth) 
    if (n.smooth)
    for (i in 1:n.smooth)
    { k<-1
      for (j in object$smooth[[i]]$first.para:object$smooth[[i]]$last.para)
      { term.names[j]<-paste(object$smooth[[i]]$label,".",as.character(k),sep="")
        k<-k+1
      }
    }
    names(object$coefficients)<-term.names  # note - won't work on matrices!!
    if (is.null(weights))
    object$prior.weights <- object$y*0+1
    else object$prior.weights <- weights 
    object$weights<-object$prior.weights   

    ret$gam<-object
    ret

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
  sp<-attr(tf,"specials")$s     # array of indeces of smooth terms 
  tp<-attr(tf,"specials")$te    # indeces of tensor product terms
  off<-attr(tf,"offset") # location of offset in formula
  if (!is.null(off)) 
  { sp[sp>off]<-sp[sp>off]-1 # have to remove the offset from this index list 
    tp[tp>off]<-tp[tp>off]-1
  }
  ns<-length(sp)+length(tp) # number of smooths
  k<-kt<-ks<-kp<-1 # counters for terms in the 2 formulae
  len.sp <- length(sp)
  len.tp <- length(tp)

  smooth.spec<-list()
  if (nt)
  for (i in 1:nt) # work through all terms
  { if (k<=ns&&((ks<=len.sp&&sp[ks]==i+1)||(kt<=len.tp&&tp[kt]==i+1))) # it's a smooth
    { st<-eval(parse(text=terms[i]),envir=p.env)
      if (k>1||kp>1) rf<-paste(rf,"+",st$full.call,sep="") # add to full formula
      else rf<-paste(rf,st$full.call,sep="")
      smooth.spec[[k]]<-st
      if (ks<=len.sp&&sp[ks]==i+1) ks<-ks+1  # counts s() terms
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




gam.setup<-function(formula,pterms,data=stop("No data supplied to gam.setup"),knots=NULL,sp=NULL,
                    min.sp=NULL,H=NULL,fit.method="magic",parametric.only=FALSE)
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
  if (m>0) for (i in 1:m) 
  { # idea here is that terms are set up in accordance with information given in split$smooth.spec
    # appropriate basis constructor is called depending on the class of the smooth
    # constructor returns penalty matrices model matrix and basis specific information
    sm<-smooth.construct(split$smooth.spec[[i]],data,knots)
    n.para<-ncol(sm$X)
    # define which elements in the parameter vector this smooth relates to....
    sm$first.para<-first.para     
    first.para<-first.para+n.para
    sm$last.para<-first.para-1
    X<-cbind(X,sm$X);sm$X<-NULL
   
    G$smooth[[i]] <- sm   
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
        G$sp[n.smooths] <- -1 # confirms that it's to be estiamted
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
    if (length(sp)!=k.sp) { ok<-FALSE;warning("Fixed smoothing parameter vector is too short - ignored.")}
    if (sum(is.na(sp))) { ok<-FALSE;warning("NA's in fixed smoothing parameter vector - ignoring.")}
  } else # set up for auto-initialization
  G$sp<-rep(-1,k.sp) # is this really needed?
  
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



gam <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,
                control=gam.control(),
                scale=0,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=TRUE,G=NULL,...)

# Routine to fit a GAM to some data. The model is stated in the formula, which is then 
# parsed and interpreted to figure out which bits relate to smooth terms and which to parametric terms.

{  if (is.null(G))
   { # create model frame..... 
    gp<-interpret.gam(formula) # interpret the formula 
    cl<-match.call() # call needed in gam object for update to work
    mf<-match.call(expand.dots=FALSE)
    mf$formula<-gp$fake.formula 
    mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-mf$gamma<-mf$fit<-mf$G<-mf$...<-NULL
    mf$drop.unused.levels<-TRUE
    mf[[1]]<-as.name("model.frame")
    pmf <- mf
    mf <- eval(mf, parent.frame()) # the model frame now contains all the data 

    terms <- attr(mf,"terms")
    
    pmf$formula <- gp$pf
    pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part

    pterms <- attr(pmf,"terms")

    if (is.character(family)) family<-eval(parse(text=family))
    if (is.function(family)) family <- family()
    if (is.null(family$family)) stop("family not recognized")
  
    if (sum(control$fit.method==c("mgcv","magic","fastest"))==0) stop("Unknown fit method.") 
    G<-gam.setup(gp,pterms=pterms,data=mf,knots=knots,sp=sp,min.sp=min.sp,
                 H=H,fit.method=control$fit.method,parametric.only=FALSE)
    G$terms<-terms;G$pterms<-pterms
    G$mf<-mf;G$cl<-cl;

    G$C<-gam.side.conditions(G) # check for identifiability of smooth part, constrain if poss.
 
    if (is.null(G$offset)) G$offset<-rep(0,G$n)
  
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

  # take only one IRLS step to get initial sp estimates for "pure" outer
  # looping...
    
  if (control$spIterType=="outer") oneStep <- TRUE else oneStep <- FALSE
  
  object<-gam.fit(G,family=family,control=control,gamma=gamma,oneStep=oneStep,...)

  if (!is.null(sp)) # fill returned s.p. array with estimated and supplied terms
  { temp.sp<-object$sp
    object$sp<-sp
    object$sp[sp<0]<-temp.sp
  }
 
  mgcv.conv<-object$mgcv.conv
  # need to check that i) following is wanted; (ii) there are free s.p.s to estimate...
  if ((control$spIterType!="perf")&&(length(G$S)>0)&&(G$fit.method=="magic"||sum(G$sp<0)!=0) )
  { lsp<-log(object$sp[G$sp<0]) # make sure only free s.p.s are optimized!
    um<-nlm(full.score,lsp,typsize=lsp,fscale=abs(object$gcv.ubre),stepmax=1,
            ndigit=12,gradtol=1e-4,steptol=0.01,G=G,family=family,control=control,gamma=gamma)
    lsp<-um$estimate
    
    object<-attr(full.score(lsp,G,family,control,gamma=gamma),"full.gam.object")
    object$mgcv.conv<-mgcv.conv # want info on power iteration, not single evaluation calls used by nlm
    object$mgcv.conv$nlm.iterations<-um$iterations
  }
  object$fit.method<-G$fit.method

  if (object$fit.method=="magic") object$rank <- object$mgcv.conv$rank else {
  object$rank <- ncol(G$X)-nrow(G$C) }  # model rank

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
  names(object$coefficients)<-term.names  # note - won't work on matrices!!
  object$full.formula<-as.formula(G$full.formula)
  environment(object$full.formula)<-environment(G$formula) 
  object$formula<-G$formula
  object$model<-G$mf # store the model frame
  object$control <- control
  object$terms <- G$terms
  object$pterms <- G$pterms
  object$assign <- G$assign # applies only to pterms
  object$contrasts <- G$contrasts
  object$xlevels <- G$xlevels
  object$offset <- G$offset
  object$data <- data
  object$df.residual <- object$df.null - sum(object$edf)
  object$min.edf<-G$min.edf
  object$call<-G$cl # needed for update() to work
  class(object)<-c("gam","glm","lm")
  object
}

gam.check<-function(b)
# takes a fitted gam object and produces some standard diagnostic plots
{ if (b$method=="GCV"||b$method=="UBRE")
  { old.par<-par(mfrow=c(2,2))
    sc.name<-b$method
    if (b$fit.method=="mgcv")
    { if (b$mgcv.conv$iter>0)
      { plot(b$mgcv.conv$edf,b$mgcv.conv$score,xlab="Estimated Degrees of Freedom",
         ylab=paste(sc.name,"Score"),main=paste(sc.name,"w.r.t. model EDF"),type="l")
        points(b$nsdf+sum(b$edf),b$gcv.ubre,col=2,pch=20)
      }
    } else qqnorm(residuals(b))
    plot(fitted(b),residuals(b),main="Residuals vs. Fitted",xlab="Fitted Values",ylab="Residuals");
    hist(residuals(b),xlab="Residuals",main="Histogram of residuals");
    plot(fitted(b),b$y,xlab="Fitted Values",ylab="Response",main="Response vs. Fitted Values")
    if (b$mgcv.conv$iter>0)   
    cat("\nSmoothing parameter selection converged after",b$mgcv.conv$iter,"iteration")
    else 
    cat("\nModel required no smoothing parameter selection")
    if (b$mgcv.conv$iter>1) cat("s")
    
    if ((b$fit.method=="mgcv"&&b$mgcv.conv$step.fail)||(b$fit.method=="magic"&&!b$mgcv.conv$fully.converged)) 
    cat(" by steepest\ndescent step failure.\n") else cat(".\n")
    if (b$fit.method=="mgcv")
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
  plot(fitted(b),residuals(b),xlab="fitted values",ylab="residuals")
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

gam.control<-function (irls.reg=0.0,epsilon = 1e-04, maxit = 20,globit = 20,mgcv.tol=1e-6,mgcv.half=15, 
                       nb.theta.mult=10000,trace =FALSE,fit.method="magic",perf.iter=NULL,spIterType="perf",
                       rank.tol=.Machine$double.eps^0.5) 
# Control structure for a gam. 
# irls.reg is the regularization parameter to use in the GAM fitting IRLS loop.
# epsilon is the tolerance to use in the IRLS MLE loop. maxit is the number 
# of IRLS iterations to use with local search for optimal s.p. after globit iterations have used global 
# searches. mgcv.tol is the tolerance to use in the mgcv call within each IRLS. mgcv.half is the 
# number of step halvings to employ in the mgcv search for the optimal GCV score, before giving up 
# on a search direction. trace turns on or off some de-bugging information.
# nb.theta.mult controls the upper and lower limits on theta estimates - for use with negative binomial  
# fit.method can be "magic" for QR/SVD method or "mgcv" for Wood (2000) method, or "fastest" to use "mgcv" 
# for single s.p. case and "magic" otherwise. 
# perf.iter (deprecated) TRUE to use Gu's performance iteration, FALSE to follow it up with
#           O'Sullivan's slower approach.
# spIterType is one of "perf" for Gu's performance iteration, "outer" for the
#             slower O'Sullivan type approach, "perf+outer" for one after the other 
# rank.tol is the tolerance to use for rank determination
{   if (!is.null(perf.iter)) {
      warning("perf.iter is deprecated: use spIterType")
      if (!is.logical(perf.iter)) stop("power.iter must be one of TRUE or FALSE.")
      if (perf.iter) spIterType <- "perf" else
      spIterType <- "perf+outer"
    }
    if (!(spIterType%in%c("perf","perf+outer","outer"))) stop(
    "spIterType must be one of \"perf\", \"perf+outer\" or \"outer\".")
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
    list(irls.reg=irls.reg,epsilon = epsilon, maxit = maxit,globit = globit, trace = trace, mgcv.tol=mgcv.tol,
         mgcv.half=mgcv.half,nb.theta.mult=nb.theta.mult,fit.method=fit.method,perf.iter=perf.iter,
         rank.tol=rank.tol,spIterType=spIterType)
    
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


full.score<-function(sp,G,family,control,gamma)
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
  res<-xx$gcv.ubre
  attr(res,"full.gam.object")<-xx
  res
}

gam.fit<-function (G, start = NULL, etastart = NULL, 
    mustart = NULL, family = gaussian(), 
    control = gam.control(),gamma=1,oneStep=FALSE) 
# fitting function for a gam, modified from glm.fit.
# note that smoothing parameter estimates from one irls iterate are carried over to the next irls iterate
# unless the range of s.p.s is large enough that numerical problems might be encountered (want to avoid 
# completely flat parts of gcv/ubre score). In the latter case autoinitialization is requested.
# oneStep == TRUE causes only a single IRLS step to be taken
{
    intercept<-G$intercept
    conv <- FALSE
    nobs <- NROW(G$y)
    nvars <- NCOL(G$X) # check this needed
    y<-G$y # original data
    X<-G$X # original design matrix
    if (nvars == 0) stop("Model seems to contain no terms")
    if (family$family=="gaussian" && family$link=="identity") olm<-TRUE # i.e. only one iteration needed
    else olm<-FALSE 
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

        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon || olm || oneStep) {
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
    nulldf <- n.ok 
    if (G$fit.method=="magic") # then some post processing is needed to extract covariance matrix etc...
    { mv<-magic.post.proc(G$X,mr,w=G$w)
      G$Vp<-mv$Vb;G$hat<-mv$hat;
      #G$edf<-array(0,0) # Now some edf's for each term....
      #if (G$m) for (i in 1:G$m) G$edf[i]<-sum(mv$edf[G$off[i]:(G$off[i]+G$df[i]-1)])
      G$edf<-mv$edf
      G$conv<-mr$gcv.info
      G$sp<-msp
    }

    aic.model <- aic(y, n, mu, weights, dev) + 2 * sum(G$edf)

	list(coefficients = as.vector(coef), residuals = residuals, fitted.values = mu, 
        family = family,linear.predictors = eta, deviance = dev,
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,  
        df.null = nulldf, y = y, converged = conv,sig2=G$sig2,edf=G$edf,hat=G$hat,
        boundary = boundary,sp = G$sp,nsdf=G$nsdf,Vp=G$Vp,mgcv.conv=G$conv,
        gcv.ubre=G$gcv.ubre,aic=aic.model)
}


predict.gam<-function(object,newdata,type="link",se.fit=FALSE,terms=NULL,
                       block.size=1000,newdata.guaranteed=FALSE,...) 
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

  if (type!="link"&&type!="terms"&&type!="response"&&type!="lpmatrix")  
  { warning("Unknown type, reset to terms.")
    type<-"terms"
  }
  if (!inherits(object,"gam")) stop("predict.gam can only be used to predict from gam objects")

  # get data from which to predict.....  
  if (newdata.guaranteed==FALSE)
  { if (missing(newdata)) # then "fake" an object suitable for prediction 
    { newdata<-object$model
      new.data.ok <- FALSE
    }
    else  # do an R ``standard'' evaluation to pick up data
    { new.data.ok <- TRUE
      if (is.data.frame(newdata)&&!is.null(attr(newdata,"terms"))) # it's a model frame
      { if (sum(!(names(object$model)%in%names(newdata)))) stop(
        "newdata is a model.frame it should contain all required variables\n")
      } else
      { ff <- interpret.gam(object$full.formula)$fake.formula[-2]
        if (sum(!(all.vars(ff)%in%names(newdata)))) { 
        warning("not all required variables have been supplied in newdata: bad practice!\n")}
        newdata <-
        eval(model.frame(ff,data=newdata),parent.frame()) 
      }
    }
  }
  # check that factor levels match for prediction and original fit 
  if (new.data.ok)
  { names(newdata)->nn # new data names
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
    { mf <- model.frame(Terms,data,xlev=object$xlevels)
      if (!is.null(cl <- attr(object$pterms,"dataClasses"))) .checkMFClasses(cl,mf)
      Xp <- model.matrix(Terms,mf,contrasts=object$contrasts)
    } else Xp <- model.matrix(Terms,object$model)
    
    if (object$nsdf) X[,1:object$nsdf]<-Xp
    if (n.smooth) for (k in 1:n.smooth) 
    { X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para]<-
                              Predict.matrix(object$smooth[[k]],data)
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
  }
  if (type=="lpmatrix") { colnames(H) <- names(object$coefficients)} else { 
    if (se.fit) H<-list(fit=fit,se.fit=se) else H<-fit
  }
  H # ... and return
}

plot.gam<-function(x,residuals=FALSE,rug=TRUE,se=TRUE,pages=0,select=NULL,scale=-1,n=100,n2=40,
                   pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                   ylim=NULL,xlim=NULL,too.far=0.1,all.terms=FALSE,shade=FALSE,col.shade="gray80",...)

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

{ sp.contour<-function (x,y,z,zse,xlab="",ylab="",zlab="main",se.plot=TRUE,se.mult=1,...)   
  # internal function for contouring 2-d smooths with 1 s.e. limits
  { gap<-median(zse,na.rm=TRUE)  
    zr<-max(z+zse,na.rm=TRUE)-min(z-zse,na.rm=TRUE) # plotting range  
    n<-10  
    while (n>1 && zr/n<2.5*gap) n<-n-1    
    zrange<-c(min(z-zse,na.rm=TRUE),max(z+zse,na.rm=TRUE))  
    zlev<-pretty(zrange,n)  
    yrange<-range(y);yr<-yrange[2]-yrange[1]  
    xrange<-range(x);xr<-xrange[2]-xrange[1]  
    ypos<-yrange[2]+yr/10
    plot(x,y,type="n",xlab="",ylab="",axes=FALSE)
    cs<-(yr/10)/strheight(zlab);if (cs>1) cs<-1 # text scaling based on height  
    cw<-par()$cxy[1]  
    tl<-strwidth(zlab);  
    if (tl*cs>3*xr/10) cs<-(3*xr/10)/tl  
    contour(x,y,z,levels=zlev,lwd=2,labcex=cs*0.65,axes=FALSE,add=TRUE)  
    axis(1,cex.axis=cs);axis(2,cex.axis=cs);box();  
    mtext(xlab,1,2.5,cex=cs);mtext(ylab,2,2.5,cex=cs)  
    contour(x,y,z+zse,levels=zlev,add=TRUE,lty=2,col=2,labcex=cs*0.5)  
    contour(x,y,z-zse,levels=zlev,add=TRUE,lty=3,col=3,labcex=cs*0.5)  
    
    xpos<-xrange[1]  
    xl<-c(xpos,xpos+xr/10);yl<-c(ypos,ypos)  
    lines(xl,yl,xpd=TRUE,lty=3,col=3)  
    text(xpos+xr/10,ypos,paste("-",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs,off=0.5*cs)  
  
    xpos<-xpos+3*xr/10  
    xl<-c(xpos,xpos+xr/10);  
    lines(xl,yl,xpd=TRUE,lwd=2)  
    text(xpos+xr/10,ypos,zlab,xpd=TRUE,pos=4,cex=cs,off=0.5*cs)  
    
    xpos<-xrange[2]-xr/5  
    xl<-c(xpos,xpos+xr/10);  
    lines(xl,yl,xpd=TRUE,lty=2,col=2)  
    text(xpos+xr/10,ypos,paste("+",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs,off=0.5*cs)  
  }   


  # start of main function
  w.resid<-NULL
  if (length(residuals)>1) # residuals supplied 
  { if (length(residuals)==length(x$residuals)) 
    w.resid <- residuals else
    warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  } else partial.resids <- residuals # use working residuals or none
  m<-length(x$smooth) # number of smooth terms
  if (all.terms) # plot parametric terms as well
  { order <- attr(x$pterms,"order")
    n.para <- sum(order==1) # plotable parametric terms   
  } else n.para <- 0
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
      X <- Predict.matrix(x$smooth[[i]],dat)   # prediction matrix from this term
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
      X <- Predict.matrix(x$smooth[[i]],dat)   # prediction matrix for this term
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
        { ul <- max(pd[[i]]$p.resid)
          if (ul > ylim[2]) ylim[2] <- ul
          ll <-  min(pd[[i]]$p.resid)
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
            { max.r <- max(pd[[i]]$p.resid)
              if ( max.r> ylimit[2]) ylimit[2] <- max.r
              min.r <-  min(pd[[i]]$p.resid)
              if (min.r < ylimit[1]) ylimit[1] <- min.r
            }
          }
          if (!is.null(ylim)) ylimit <- ylim
          if (shade)
          { plot(pd[[i]]$x,pd[[i]]$fit,type="n",xlab=pd[[i]]$xlab,ylim=ylimit,
                 xlim=xlim,ylab=pd[[i]]$ylab,main=main,...)
            polygon(c(pd[[i]]$x,pd[[i]]$x[n:1],pd[[i]]$x[1]),
                     c(ul,ll[n:1],ul[1]),col = col.shade,border = NA)
            lines(pd[[i]]$x,pd[[i]]$fit)
          } else
          { plot(pd[[i]]$x,pd[[i]]$fit,type="l",xlab=pd[[i]]$xlab,ylim=ylimit,xlim=xlim,
                 ylab=pd[[i]]$ylab,main=main,...)
	    if (is.null(list(...)[["lty"]]))
            { lines(pd[[i]]$x,ul,lty=2,...)
              lines(pd[[i]]$x,ll,lty=2,...)
            } else
            { lines(pd[[i]]$x,ul,...)
              lines(pd[[i]]$x,ll,...)
            }
          } 
          if (partial.resids)
          { if (is.null(list(...)[["pch"]]))
            points(pd[[i]]$raw,pd[[i]]$p.resid,pch=".",...) else
            points(pd[[i]]$raw,pd[[i]]$p.resid,...)
          }
	  if (rug) 
          { if (jit) rug(jitter(as.numeric(pd[[i]]$raw)),...)
             else rug(as.numeric(pd[[i]]$raw),...)
	  }
        } else if (pd[[i]]$dim==2)
        { 
          if (pers) 
          { if (!is.null(main)) pd[[i]]$title <- main
            persp(pd[[i]]$xm,pd[[i]]$ym,matrix(pd[[i]]$fit,n2,n2),xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,
                  zlab=pd[[i]]$title,ylim=pd[[i]]$ylim,xlim=pd[[i]]$xlim,theta=theta,phi=phi,...)
          } else
          { sp.contour(pd[[i]]$xm,pd[[i]]$ym,matrix(pd[[i]]$fit,n2,n2),matrix(pd[[i]]$se,n2,n2),
                     xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,zlab=pd[[i]]$title,se.mult=se2.mult,...)
            if (rug) points(pd[[i]]$raw$x,pd[[i]]$raw$y,pch=".",...)
          } 
        } else
        { warning("no automatic plotting for smooths of more than one variable")
        }
      }  
      j<-j+pd[[i]]$dim
    }
  } else # don't plot confidence limits
  { k<-0
    if (scale==-1&&is.null(ylim))
    if (m>0) for (i in 1:m)
    { if (pd[[i]]$dim==1)
      { if (k==0) 
        { if (partial.resids) ylim <- range(pd[[i]]$p.resid) else ylim<-range(pd[[i]]$fit);k<-1 }
	else
        { if (partial.resids)
          { if (min(pd[[i]]$p.resid)<ylim[1]) ylim[1]<-min(pd[[i]]$p.resid)
	    if (max(pd[[i]]$p.resid)>ylim[2]) ylim[2]<-max(pd[[i]]$p.resid)
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
          { if (partial.resids) ylimit <- range(pd[[i]]$p.resid) else ylimit <-range(pd[[i]]$fit)}
          if (!is.null(ylim)) ylimit <- ylim
          plot(pd[[i]]$x,pd[[i]]$fit,type="l",,xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,ylim=ylimit,xlim=xlim,main=main,...)
          if (rug) 
	  { if (jit) rug(jitter(as.numeric(pd[[i]]$raw)),...)
            else rug(as.numeric(pd[[i]]$raw),...) 
          }
          if (partial.resids)
          { if (is.null(list(...)[["pch"]]))
            points(pd[[i]]$raw,pd[[i]]$p.resid,pch=".",...) else
            points(pd[[i]]$raw,pd[[i]]$p.resid,...)
          }
        } else if (pd[[i]]$dim==2)
        { if (!is.null(main)) pd[[i]]$title <- main
          if (pers) 
          { persp(pd[[i]]$xm,pd[[i]]$ym,matrix(pd[[i]]$fit,n2,n2),xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,
                          zlab=pd[[i]]$title,theta=theta,phi=phi,xlim=pd[[i]]$xlim,ylim=pd[[i]]$ylim,...)
          }
          else
          { contour(pd[[i]]$xm,pd[[i]]$ym,matrix(pd[[i]]$fit,n2,n2),xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,
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
      if (interactive() && m && pd[[i]]$dim<3 && i>1&&(i-1)%%ppp==0) 
      readline("Press return for next page....")
      termplot(x,se=se,rug=rug,col.se=1,col.term=1)
    } else { # figure out which plot is required
      if (select > m) { 
        select <- select - m # i.e. which parametric term
        term.labels <- attr(x$pterms,"term.labels")
        term.labels <- term.labels[order==1]
        if (select <= length(term.labels)) {
        if (interactive() && m && pd[[i]]$dim<3 && i>1&&(i-1)%%ppp==0) 
        readline("Press return for next page....")
        termplot(x,terms=term.labels[select],se=se,rug=rug,col.se=1,col.term=1)
        }  
      }
    }
  }
  if (pages>0) par(oldpar)
}


residuals.gam <-function(object, type = c("deviance", "pearson","scaled.pearson", "working", "response"),...)
# calculates residuals for gam object - defualt for glm (from which this is adapted) seems to be buggy
{ type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  family <- object$family
  wts <- object$prior.weights
  switch(type,working = object$residuals,
       scaled.pearson = (y-mu)*sqrt(wts)/sqrt(object$sig2*object$family$variance(mu)),
              pearson = (y-mu)*sqrt(wts)/sqrt(object$family$variance(mu)),
              deviance = { d.res<-sqrt(pmax(object$family$dev.resids(y,mu,wts),0))
                           ifelse(y>mu , d.res, -d.res)             
                         },
              response = y - mu)
}

summary.gam<-function (object,...) 
# summary method for gam object - provides approximate p values for terms + other diagnostics
{ pinv<-function(V,M)
  { D<-La.svd(V)
    M1<-length(D$d[D$d>1e-12*D$d[1]]);if (M>M1) M<-M1 # avoid problems with zero eigen-values
    D$d[(M+1):length(D$d)]<-1
    D$d<- 1/D$d
    D$d[(M+1):length(D$d)]<-0
    D$u%*%diag(D$d)%*%D$v
  } 
  se<-0;for (i in 1:length(object$coefficients)) se[i]<-object$Vp[i,i]^0.5
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
    Vp<-matrix(object$Vp[1:np,1:np],np,np)
    bp<-array(object$coefficients[1:np],np)
    pTerms.pv <- array(0,nt)
    attr(pTerms.pv,"names") <- term.labels
    pTerms.df <- pTerms.chi.sq <- pTerms.pv
    for (i in 1:nt)
    { ind <- object$assign==i
      b <- bp[ind];V <- Vp[ind,ind]
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
      V<-object$Vp[start:stop,start:stop] # cov matrix for smooth
      p<-object$coefficients[start:stop]  # params for smooth
      # now get null space dimension for this term
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
      if (object$method=="UBRE")
      s.pv[i]<-pchisq(chi.sq[i],df=max(1,edf[i]),lower.tail=FALSE)
      else     
      s.pv[i]<-pf(chi.sq[i]/edf[i],df1=max(1,edf[i]),df2=residual.df,lower.tail=FALSE) 
    }
  }
  r.sq<- 1 - var(object$y-object$fitted.values)*(object$df.null-1)/(var(object$y)*residual.df) 
  dev.expl<-(object$null.deviance-object$deviance)/object$null.deviance
  ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,
       s.pv=s.pv,scale=object$sig2,r.sq=r.sq,family=object$family,formula=object$formula,n=object$df.null,
       dev.expl=dev.expl,edf=edf,dispersion=object$sig2,pTerms.pv=pTerms.pv,pTerms.chi.sq=pTerms.chi.sq,
       pTerms.df=pTerms.df)
  if (object$method=="GCV") ret$gcv<-object$gcv.ubre else if (object$method=="UBRE") ret$ubre<-object$gcv.ubre
  class(ret)<-"summary.gam"
  ret
}

anova.gam <- function (object, ..., dispersion = NULL, test = NULL)
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
    sg <- summary(object) 
    class(sg) <- "anova.gam"
    sg
}

influence.gam <- function(model,...) { model$hat }

print.anova.gam <- function(x,...)
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


print.summary.gam<-function(x,...)
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
  o<-.C("MinimumSeparation",as.double(g1),as.double(g2),as.integer(n),as.double(d1),as.double(d2),as.integer(m),
         distance=as.double(distance),PACKAGE="mgcv")  
  res<-rep(FALSE,n)
  res[o$distance > dist] <-TRUE
  res
}

vis.gam<-function(x,view=NULL,cond=list(),n.grid=30,too.far=0,col=NA,color="heat",
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

  if (is.null(view)) # get default view if none supplied
  { v.names<-attr(attr(x$model,"terms"),"term.labels")
    if (length(v.names)<2) stop("Model doesn't seem to have enough terms to do anything useful")
    view<-v.names[1:2]
  }
  if (!sum(view%in%names(x$model))) stop(
  paste("view variables must be one of",names(x$model)))
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


mroot<-function(A,rank=NULL,method="chol")
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
      while (rank>0&&(um$d[rank]/um$d[1]<.Machine$double.eps||
                           all.equal(um$u[,rank],um$vt[rank,])!=TRUE)) rank<-rank-1 
      if (rank==0) stop("Something wrong - matrix probably not +ve semi definite")    
    }
    d<-um$d[1:rank]^0.5
    return(t(t(um$u[,1:rank])*d)) # note recycling rule used for efficiency
  } else
  if (method=="chol")
  { L<-chol(A,pivot=TRUE)
    piv<-order(attr(L,"pivot"))
    if (is.null(rank)) rank<-attr(L,"rank")
    L<-L[,piv];L<-t(L[1:rank,])
    return(L)
  } else
  stop("method not recognised.")
}



magic.post.proc<-function(X,object,w=NULL)
# routine to take list returned by magic and extract:
# Vb the estimated parameter covariance matrix. rV%*%t(rV)*scale
# edf the array of estimated degrees of freedom per parameter Vb%*%t(X)%*%W%*%X /scale
# hat the leading diagonal of the hat/influence matrix 
# NOTE: W=diag(w^2). 
# flop count is O(nq^2) if X is n by q... this is why routine not part of magic
{ V<-object$rV%*%t(object$rV)
  if (!is.null(w)) 
  { if (is.matrix(w))
    X<-w%*%X else 
    X<-w*X # use recycling rule to form diag(w)%*%X cheaply 
    
  }
  M<-V%*%t(X);B<-X*t(M);rm(M)
  hat<-apply(B,1,sum) # diag(X%*%V%*%t(X))
  edf<-apply(B,2,sum) # diag(V%*%t(X)%*%X)
  Vb<-V*object$scale;rm(V)
  list(Vb=Vb,hat=hat,edf=edf)
}


magic<-function(y,X,sp,S,off,rank=NULL,H=NULL,C=NULL,w=NULL,gamma=1,scale=1,gcv=TRUE,
                ridge.parameter=NULL,control=list(maxit=50,tol=1e-6,step.half=25,
                rank.tol=.Machine$double.eps^0.5))
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
# If ridge.parameter is a positive number then then it is assumed to be the multiplier
# for a ridge penalty to be applied during fitting. 
{ n.p<-length(S)
  n.b<-dim(X)[2] # number of parameters
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
    if (n.p>0) for (i in 1:n.p) S[[i]]<-qr.qty(ns.qr,S[[i]])[(n.con+1):n.b,]
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
      X<-w*X # use recycling rule to form diag(w)%*%X cheaply
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
  um<-.C("magic",as.double(y),as.double(X),sp=as.double(sp),as.double(Si),as.double(H),
          score=as.double(gamma),scale=as.double(scale),info=as.integer(icontrol),as.integer(cS),
          as.double(control$rank.tol),rms.grad=as.double(control$tol),b=as.double(b),rV=double(q*q),
          PACKAGE="mgcv")
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



.onAttach <- function(...) cat("This is mgcv 1.1-2 \n")


.First.lib <- function(lib, pkg) {
    library.dynam("mgcv", pkg, lib)
    cat("This is mgcv 1.1-2 \n")
}


###############################################################################
### ISSUES.....






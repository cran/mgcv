# These are the R routines for the package mgcv (c) Simon Wood 2000-2001

QT<-function(A) {

# This routine performs a QT factorization of matrix A. A QT factorization
# is of the form AQ=[0,T] where T is reverse lower triangular (top left zero)
# nrow(A)<ncol(A) and the first (ncol(A)-nrow(A)) columns of Q are the null
# space of A. Q is actually stored as Householder rotations in over A. 
# Specifically if u_i is row i of A then let H_i=(I-u_i u_i') then 
# Q = H_1 H_2 ... there are nrow(A) u_i's and H_i's. The last i-1 elements
# of u_i are always zero (i>=1). The main use of this routine is to find
# the null space of a constraint matrix C, as used in mgcv().
# The return value is the over-written A. The factorization is performed
# using compiled code.

oo<-.C("RQT",as.double(A),as.integer(nrow(A)),as.integer(ncol(A)))
A<-matrix(oo[[1]],nrow(A),ncol(A))
A
}


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
  A<-matrix(0,4*(n-1)+lo+hi,n)
  b<-array(0,4*(n-1)+lo+hi)
  if (lo*hi==1&&lower>=upper) stop("lower bound >= upper bound in call to mono.con()")
  oo<-.C("RMonoCon",as.double(A),as.double(b),as.double(x),as.integer(control),as.double(lower),as.double(upper),as.integer(n))
  A<-matrix(oo[[1]],dim(A)[1],dim(A)[2])
  b<-array(oo[[2]],dim(A)[1])
  list(A=A,b=b)
}  


uniquecombs<-function(x) {
# takes matrix x and counts up unique rows
res<-.C("RuniqueCombs",as.double(x),as.integer(nrow(x)),as.integer(ncol(x)))
n<-res[[2]]*res[[3]]
x<-matrix(res[[1]][1:n],res[[2]],res[[3]])
x
}



null.space.dimension<-function(d,m)

{ if (2*m<=d) {m<-1;while (2*m<d+2) m<-m+1;}
  M<-1;     # dimension of penalty null space
  for (i in 0:(d-1)) M<-M*(d+m-1-i);
  if (d>1) for (i in 2:d) M<-M/i;     # M = (m+d+1)!/(d!(m-d!)
  M
}


GAMsetup<-function(G) {
# This function sets up the design matrix and penalty matrices for a GAM formulated with
# one dimensional cubic penalized regression splines, or as thin plate regression splines.
# It takes a list G as argument which  should have the following named elements:
#
# G$m the number of smooth terms in the model
# G$n the number of data to be modelled
# G$nsdf the number of user supplied columns of the design matrix for any parametric model parts
#        including the constant (if any)
# G$df an array of G$m integers specifying the maximum d.f. for each spline term
# G$dim an array containing the dimensions of each smooth
# G$s.type an array containing a code for the type of each smooth 0 - cubic regression spline
#          1 - tprs.
# G$p.order is an array giving the orders of the spline penalties to use 
# G$x an array of G$n element arrays of data and (optionally) design matrix columns. The first
#     G$nsdf elements of G$x should contain the elements of the columns of the design matrix
#     corresponding to the parametric part of the model. The remaining G$m elements of G$x 
#     are the values of the covariates that are arguments of the spline terms.
#
# The function returns a list H containing all the above named elements plus the following:
#
# H$X the full design matrix.
# H$S contains the elements of the G$m penalty matrices. let start.k=sum(G$df[1:(k-1)]^2) and start_1=0
#     Element i,j of the kth penalty matrix is S[start.k+i+G$df[k]*(j-1)]. Note however that this
#     matrix is infact the smallest matrix containing all the non-zero elements of the full
#     penalty. Paste it into location starting at M$off[k], M$off[k] of an appropriate
#     matrix of zeroes to get the full penalty.
# H$off is an array of offsets, used to facilitate efficient storage of the penalty matrices
#       and to indicate where in the overall parameter vector the parameters of the ith 
#       spline reside (e.g. first parameter of ith spline is at p[off[i]]).
# H$C - the matrix of linear equality constraints on the parameters .  
# H$covariate.shift - the shift applied to the covariates in order to centre them
# H$UZ - Array containing matrices giving the truncated basis for the ith tprs term.
#        Let start[k]=sum_{i=1}^{k-1} (n+M[i])*df[i],  where n is number of data
#        M[i] is null space dimension and df[i] the basis dimension for term i. start[0]=0. 
#        Then row i col j of basis matrix k is at H$UZ[start+i+(j-1)*(n+M[k])] 
# H$Xu - This is an array used for storing matrices containing all the unique covariate combinations for
#        each model term (when it's a t.p.r.s.). The packing convention is a s follows:
#        let start[1]=0 and start[k+1]=start[k]+xu[k]*dim[k] where xu[k] and dim[k] give the
#        number of unique covariate combinations and the number of coivariates for the kth term.
#        then the jth covariate for the ith unique combination for the kth model term is at:
#        h$Xu[start[k]+i+(j-1)*xu[k]].
# H$xu.length - H$xu.length[i] is the number of unique covariate values for the ith tprs term
# H$xp - if the spline term is a cubic regression spline then this is  the matrix whose 
#        rows contain the covariate values corresponding to the parameters
#        of each spline - the splines are parameterized using their y- values at a series
#        of x values - these vectors contain those x values!

  q<-G$nsdf
  if (G$m) for (i in 1:G$m) { q<-q+G$df[i] }  # q stores total number of parameters
  X<-matrix(0,G$n,q)             # design matrix
  if (G$m>0)
  { mdf<-max(G$df)
    # centre the covariates to improve numerical stability
    nc<-0;for (i in 1:G$m) nc<-nc+G$dim[i] # count up covariates
    G$covariate.shift<-array(0,nc)
    for (i in 1:nc) G$covariate.shift[i]<-mean(G$x[i+G$nsdf,])
    for (i in 1:nc) G$x[i+G$nsdf,]<-G$x[i+G$nsdf,]-G$covariate.shift[i]
    xp<-matrix(0,G$m,mdf)          # covariate matrix
    S<-array(0,sum(G$df^2))   # array of wiggliness penalty matrices
    C<-matrix(0,G$m,q)             # fixed constraint matrix (one per spline term)
    m.type1<-sum(G$s.type)
    M<-array(0,G$m);for (i in 1:G$m) M[i]<-null.space.dimension(G$dim[i],G$p.order[i])
    if (m.type1>0) UZ.length<-sum((G$n+M[G$s.type==1])*G$df[G$s.type==1])
    else UZ.length<-0
    UZ<-array(0,UZ.length) # storage for truncated basis
    if (m.type1>0) Xu.length<-sum(G$n*G$dim[G$s.type==1])
    else Xu.length<-0
    Xu<-array(0,Xu.length) # storage for set of unique covariate combinations
    xu.length<-array(0,m.type1)  # storage for number of unique covariate combinations for each tprs term
    off<-array(0,c(G$m))           # array for storing offsets for Wiggliness penalties
  } else mdf<-G$covariate.shift<-xp<-S<-C<-m.type1<-M<-UZ<-Xu<-xu.length<-off<-UZ.length<-Xu.length<-0;
  o<-.C("RGAMsetup",as.double(X),as.double(C),as.double(S),as.double(UZ),as.double(Xu),as.integer(xu.length),as.double(xp),
        as.integer(off),as.double(G$x),as.integer(G$m),as.integer(G$n),as.integer(G$df),
        as.integer(G$nsdf),as.integer(G$dim),as.integer(G$s.type),as.integer(G$p.order)) # compiled code to set up matrices for model
  G$X<-matrix(o[[1]],G$n,q);
  G$C<-matrix(o[[2]],G$m,q);
  G$S<-array(o[[3]],sum(G$df^2));             #dim=c(G$m,mdf,mdf));
  G$UZ<-array(o[[4]],UZ.length)  #dim=c(m.type1,G$n+maxM,mdf.type1))
  G$xu.length<-array(o[[6]],m.type1)
  if (m.type1>0) Xu.length<-sum(G$xu.length[G$s.type==1]*G$dim[G$s.type==1])
  else Xu.length<-0
  G$Xu<-array(o[[5]],Xu.length)    #dim=c(m.type1,G$n,mdim))
  G$xp<-matrix(o[[7]],G$m,mdf);
  G$off<-array(o[[8]],c(G$m));
  G # ... and return
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
# M$S  - array packed with non-zero coefficients of S_i's see mgcv or GAMsetup for details of packing.
# M$off - used for unpacking M$S
# M$df  - used for unpacking M$S
# M$sp - array of theta_i's 
# Ain, bin and p are not in the object needed to call mgcv....
#
{ nar<-c(length(M$y),length(M$p),dim(M$Ain)[1],dim(M$C)[1],0)
  H<-0
  o<-.C("RPCLS",as.double(M$X),as.double(M$p),as.double(M$y),as.double(M$w),as.double(M$Ain),as.double(M$bin)
        ,as.double(M$C),as.double(H),as.double(M$S),as.integer(M$off),as.integer(M$df),as.double(M$sp),
        as.integer(length(M$off)),as.integer(nar))
  array(o[[2]],length(M$p))
}  

mgcv<-function(M) {

# Performs multiple smoothing parameter selection for Generalized ridge regression problems 
# M is a list which must containing the following named elements:
#
# M$y  - the response variable vector
# M$X  - the design matrix (declared as a matrix, with correct numbers of rows and columns,
#        as number of data and number of parameters are read from this)
# M$C  - matrix defining linear equality constraints on parameters (Cp=0). 
#        Number of rows is number of constraints.
# M$w  - weight vector (often proportional to inverse of variance)
# M$S  -  contains the elements of the G$m penalty matrices. let start.k=sum(M$df[1:(k-1)]^2) and start_1=0
#         Element i,j of the kth penalty matrix is S[start.k+i+M$df[k]*(j-1)]. Note however that this
#         matrix is infact the smallest matrix containing all the non-zero elements of the full
#         penalty. Paste it into location starting at M$off[k], M$off[k] of an appropriate
#         matrix of zeroes to get the full penalty.
# M$off - array of offsets used as described for M$S and also indicates where in the parameter
#         vector (p, say) the parameters penalized by the ith penalty start.
# M$df  - array containing actual dimensions of non-zero part of S[i], i.e. S[i] contains only
#         an M$df square block of non-zero elements.
# M$fix - array containing TRUE if smooth i is to have fixed d.f. and FALSE otherwise
# M$sig2 - the error variance if it is known and to be used for UBRE - must be <=0 to use GCV.
#          Note the convention used here is that var(y_i) = sig2/w_i.
# M$sp   - the initial smoothing parameter estimates: if any of these are <=0 then automatic
#           initialisation is used
# M$conv.tol - convergence tolerence for multiple s.p. GCV
# M$max.half - maximum number of step-length halvings to perform at each newton update
#              of s.p.'s
#
# The routine returns M with the following elements added (or reset):
#
# M$p    - the best fit parameter vector given the selected smoothing parameters
# M$sp   - the vector of smoothing parameters (\theta_i/\rho) estimated by the method
# M$sig2 - the estimate of the error variance (if GCV used)
# M$Vp   - the estimated covariance matrix of the parameters set Vp[1,1] <0 to not calculate
# M$edf  - the estimated degrees of freedom for the ith smooth term if Vp calculated
  
  if (is.null(M$sig2)) M$sig2<-0
  C.r<-nrow(M$C)          # number of equality constraints
  if (is.null(C.r)) C.r<-0
  q<-ncol(M$X)            # number of parameters
  
  n<-nrow(M$X)            # number of data
  # need to strip out penalties for fixed df smooths..... 
  #k<-dim(M$S)[1]             # number of penalties (before stripping)
  k<-length(M$df)
  if (is.null(k)) k<-0
  if (k!=0)  m<-k-sum(M$fix) # count penalty terms to include
  else m<-0
  #S<-array(0,c(m,dim(M$S)[2],dim(M$S)[3])) # reduced smoothness array
  if (sum(!M$fix)) S<-array(0,sum(M$df[!M$fix]^2))
  else S<-array(0,0)
  off<-array(0,m)  # reduced offset array
  df<-array(0,m)   # reduced penalty size array
  j<-1             # element index for off[], df and S[]
  stop<-0          # where to stop in copying from M$S 
  stopj<-0         # where to stop within S
  if (k!=0)
  for (i in 1:k)   # do actual stripping
  { start<-stop;stop<-start+M$df[i]^2;start<-start+1
    if (M$fix[i]==FALSE) 
    { #S[j,,]<-M$S[i,,]
      startj<-stopj;stopj<-startj+M$df[i]^2;startj<-startj+1
      S[startj:stopj]<-M$S[start:stop]
      off[j]<-M$off[i]
      df[j]<-M$df[i]
      j<-j+1
    }  
  }
  # deal with quantities that will be estimated
  p<-matrix(0,q,1)      # set up parameter vector
  Vp<-matrix(0.0,q,q)   # estimated covariance matrix
  edf<-array(0,m)       # estimated degrees of freedom
  Vp[1,1]<-1.0

  oo<-.C("mgcv",as.double(M$y),as.double(M$X),as.double(M$C),as.double(M$w),as.double(S),
         as.double(p),as.double(M$sp),as.integer(off),as.integer(df),as.integer(m),
         as.integer(n),as.integer(q),as.integer(C.r),as.double(M$sig2),as.double(Vp),
		 as.double(edf),as.double(M$conv.tol),as.integer(M$max.half))
   
  p<-matrix(oo[[6]],q,1);
  sig2<-oo[[14]]
  Vp<-matrix(oo[[15]],q,q)
  sp<-matrix(oo[[7]])
  edf<-oo[[16]]
  # unpack results back to correct place in output (which includes fixed d.f. and free d.f. terms)
  M$sp<-array(0,k)
  M$edf<-array(0,k)
  j<-1
  if (k!=0)
  for (i in 1:k)
  { if (M$fix[i]==TRUE)
    { M$sp[i]<-0
	  M$edf[i]<- -1 # must be filled in outside this routine, with reference to constraints
    }
	else
	{ M$sp[i]<-sp[j]
	  M$edf[i]<-edf[j]
      j<-j+1
    }
  }
  M$Vp<-Vp
  M$p<-p;
  M$sig2<-sig2;
  M    # return list
}


SANtest<-function(n=100,sig2=-1) { 

# Example function that simulates data from a simple additive normal truth and then attempts
# to reconstruct it using a GAM constructed from penalized regression splines. 
# n is number of data. Magnitude of sig2 is error variance used for simulation.
# sig2<=0 => use GCV, otherwise use UBRE for smoothing parameter estimation.

  if (sig2<=0) noquote("GCV used") else noquote("UBRE used")
  if (n<60) 
  { n<-60
    warning("n reset to 60")
  }
# simulate some data to which GAM can be fitted......

  x<-runif(5*n,0,1)      # simulate covariates
  x<-array(x,dim=c(5,n)) # load into array
  pi<-asin(1)*2          
  y<-2*sin(pi*x[2,])  
  y<-y+exp(2*x[3,])-3.75887
  y<-y+0.2*x[4,]^11*(10*(1-x[4,]))^6+10*(10*x[4,])^3*(1-x[4,])^10-1.396 # simulate truth
  e<-rnorm(n,0,sqrt(abs(sig2)))  # noise to add to truth
  y<-y+e              # create data from truth + noise   
  w<-matrix(1,n,1)    # weight matrix
  

# display scatter plots of data against covariates.....

  par(mfrow=c(2,2))   # set display environment to 2 by 2
  plot(x[2,],y)      
  plot(x[3,],y)
  plot(x[4,],y)
  plot(x[5,],y)

# set up input list for GAMsetup and call GAMsetup.......

  x[1,]<-1 # set up dummy variable for intercept term
  G<-list(m=4,n=n,nsdf=1,df=c(15,15,15,15),dim=c(1,1,1,1),s.type=c(0,0,0,0),p.order=c(0,0,0,0),x=x) # input list for GAMsetup
  H<-GAMsetup(G)  
  H$fix<-array(FALSE,H$m)
  H$y<-y;H$sig2<-sig2;H$w<-w # add data, variance and weights to mgcv input list
  H$sp<-array(-1,H$m)
  H$conv.tol<-1e-6;H$max.half<-15
  H<-mgcv(H)   # fit model
  xp<-H$xp
  p<-H$p;
  sig2<-H$sig2

  readline("Hit return to see fitted model and truth")
# plot reconstruction and truth....
  i<-c((H$off[1]+1):H$off[2]);yy<-p[i];yy<-yy-mean(yy);i<-i-H$off[1];xx<-xp[1,i]+H$covariate.shift[1];
  plot(xx,yy,type="l");t1<-2*sin(pi*xx);t1<-t1-mean(t1);lines(xx,t1); # t1 is truth
  i<-c((H$off[2]+1):H$off[3]);yy<-p[i];yy<-yy-mean(yy);i<-i-H$off[2];xx<-xp[2,i]+H$covariate.shift[2];
  plot(xx,yy,type="l");t2<-exp(2*xx)-3.75887;t2<-t2-mean(t2);lines(xx,t2);
  i<-c((H$off[3]+1):H$off[4]);yy<-p[i];yy<-yy-mean(yy);i<-i-H$off[3];xx<-xp[3,i]+H$covariate.shift[3];
  plot(xx,yy,type="l"); 
  t3<-0.2*xx^11*(10*(1-xx))^6+10*(10*xx)^3*(1-xx)^10-1.396;t3<-t3-mean(t3); lines(xx,t3);
  i<-c((H$off[4]+1):ncol(H$X));yy<-p[i];yy<-yy-mean(yy);i<-i-H$off[4];xx<-xp[4,i]+H$covariate.shift[4];
  plot(xx,yy,type="l");t4<-xx*0.0;lines(xx,t4);
  sig2 # the return value
}



########################################################################################
## The following provide a gam() function and associated gubbins for R......
########################################################################################



s<-function (..., k=-1,fx=FALSE,bs="tp",m=0)
# function for use in gam formulae to specify smooth term, e.g. s(x0,x1,x2,k=40,m=3) specifies 
# a rank 40 thin plate regression spline of x0,x1 and x2 witha third order penalty.
# Returns a list of consisting
# of a model formula term representing the smooth, the basis dimension, the type of basis
# (1 -t.p.r.s.; 0 - cubic), whether it is fixed or penalized and the order of the penalty (0 for auto).
# backwards compatibility with versions < 0.6 is maintained so that terms like: s(x,6|f)
# still work correctly.
# The returned full.call is a string with the call made fully explicit, so that when pasted
# into a formula, that formula can be parsed without reference to any variables that may
# have been used to define k,fx,bs or m.
{ vars<-as.list(substitute(list(...)))[-1] # gets terms to be smoothed without evaluation
  call<-match.call() # get function call
  d<-length(vars) # dimension of smoother
  term<-as.character(call[d+1]) # last term in the ... arguments
  if (term[1]=="$") term<-paste(term[2],term[1],term[3],sep="") # allow "a$b" type terms
  # deal with old style formula
  c<-substring(term,1,1) # first character of last ... argument
  old.warn<-options(warn=-1)       # suppress warnings when testing for numerals
  if (is.na(as.numeric(c))==FALSE) # then it's an old style formula smooth term (version's before 0.6)
  { options(old.warn) # turn on warnings again
    ss<-c;i<-2
    ns<-nchar(term)   # length of dimension specifying term
    if (ns>1)
    while (c!="|"&&i<=ns)
    { c<-" ";while(c==" "&&i<=ns) { c<-substring(term,i,i);i<-i+1}
      if (c!="|") ss<-paste(ss,c,sep="")  # then add digit to basis dimension string
    }
    k<-as.numeric(ss)
    if (c=="|")  # check whether term specified as fised d.f.
    { if (i>ns) error("Syntax error in s() term: missing f")
      c<-" ";while(c==" "&&i<=ns) { c<-substring(term,i,i);i<-i+1}
      if (c=="f"||c=="F") fx<-TRUE 
      else warning("Strange argument in s() call ignored")
    }
    d<-d-1 # reduce m since last ... argument specified basis dimension not variable
    if (d==1) bs<-"cr" # maintain strict backwards compatibility of code
  } else options(old.warn) # turn on regular warnings again
  term<-as.character(vars[[1]]) # first covariate
  if (term[1]=="$") term<-paste(term[2],term[1],term[3],sep="")  
  if (d>1) # then deal with further covariates
  for (i in 2:d)
  { t<-as.character(call[i+1])
    if (t[1]=="$") t<-paste(t[2],t[1],t[3],sep="")
    #term <- paste(term, ":", t, sep = "")
    term[i]<-t
  }
  # term now contains the names of the covariates for this model term
  # now evaluate all the other 
  if (k==-1) k<-10*3^{d-1} # auto-initialize basis dimension
  if (bs=="cr") # set basis types
  { bs.type<-0
    if (d>1) { warning("cr basis only works with 1-d smooths!");bs.type<-1;}
  } else 
  { if (bs!="tp") warning("unrecognised basis specifier - set to t.p.r.s. basis.");
    bs.type<-1
  }
  if (m<0) m<-0
  # assemble version of call with all options expanded as text
  full.call<-paste("s(",term[1],sep="")
  if (d>1) for (i in 2:d) full.call<-paste(full.call,",",term[i],sep="")
  full.call<-paste(full.call,",k=",deparse(k),",fx=",deparse(fx),",bs=",deparse(bs),",m=",deparse(m),")",sep="")
  ret<-list(term=term,bs.dim=k,bs.type=bs.type,fixed=fx,dim=d,p.order=m,full.call=full.call)
  ret
}


gam.parser<-function (gf,parent.level=1)
# parses a gam formula of the generic form:
#   y~x0+x1+x3*x4 + s(x5)+ s(x6,x7) ....
# and returns 2 model formulae: 1 containing only the non-smooth stuff
# and the other relating to the smooths in a conventional model notation 
# without s(,,,)'s e.g:
#   y~x5+x6:x7
# On exit a list is returned (see ret below). The parametric part may also 
# include any model offsets.
# See s() above for syntax of smooth terms - s() is called for each smooth term and
# returns information on the arguments of the smooth and it's properties. 
# Note that this routine does not evaluate model variables.
# parent.level specifies how far back to go from frame of this function to find
# arguments of s() in GAM formula
{ tf<-terms.formula(gf,specials=c("s")) # specials attribute indicates which terms are smooth
  terms<-attr(tf,"term.labels") # labels of the model terms
  nt<-length(terms) # how many terms?
  pf<-"~";          # start of parametric model formula
  if (attr(tf,"response")>0)  # start the replacement formula
  rf<-paste(as.character(attr(tf,"variables")[2]),"~",sep="")
  else rf<-"~"
  sp<-attr(tf,"specials")$s     # array of indeces of smooth terms 
  ns<-length(sp) # number of smooths
  ks<-1;kp<-1 # counters for terms in the 2 formulae
  df<- -1     # max d.f. array 
  fix<-FALSE  # fixed d.f. array 
  bs.type<-1  # basis type array
  dim<-0      # dimension of smoother array
  p.order<-0  # order of the penalties
  v.names<-as.character(attr(tf,"variables")[2])  # names of covariates for smooths starting with response
  n.cov<-1     # total number of covariates for smooths
  for (i in 1:nt) # work through all terms
  { if (ks<=ns&&sp[ks]==i+1) # it's a smooth
    { #stxt<-paste(substring(terms[i],1,nchar(terms[i])-1),",parent.level=",deparse(parent.level+1),")",sep="")
      st<-eval(parse(text=terms[i]),envir=sys.frame(sys.parent(n=parent.level))) # get smooth term information
      if (ks>1||kp>1) rf<-paste(rf,"+",st$full.call,sep="") # add to smooth formula
      else rf<-paste(rf,st$full.call,sep="")
      for (i in 1:st$dim) v.names[n.cov+i]<-st$term[i]
      n.cov<-n.cov+st$dim
      df[ks]<-st$bs.dim;          # getting other info on smooth
      fix[ks]<-st$fixed
      bs.type[ks]<-st$bs.type
      p.order[ks]<-st$p.order
      dim[ks]<-st$dim
      ks<-ks+1
    } else          # parametric
    { if (kp>1) pf<-paste(pf,"+",terms[i],sep="") # add to parametric formula
      else pf<-paste(pf,terms[i],sep="")
      if (ks>1||kp>1) rf<-paste(rf,"+",terms[i],sep="") # add to smooth formula
      else rf<-paste(rf,terms[i],sep="")
      kp<-kp+1
    }
  }
  off<-attr(tf,"offset")  
  if (length(off)>0) # deal with offset
  { if (kp>1) pf<-paste(pf,"+",sep="")
    if (kp>1||ks>1) rf<-paste(rf,"+",sep="")
    pf<-paste(pf,as.character(attr(tf,"variables")[1+off]),sep="")
    rf<-paste(rf,as.character(attr(tf,"variables")[1+off]),sep="")
    kp<-kp+1          
  }
  if (attr(tf,"intercept")==0) {pf<-paste(pf,"-1",sep="");rf<-paste(rf,"-1",sep="");pfok<-0}
  else { pfok<-1;if (kp==1) pf<-"~1"}
  sfok<-0;if (ks>1) sfok<-1;
  ret<-list(pftext=pf,pf=as.formula(pf),pfok=pfok,v.names=v.names,fix=fix,df=df,bs.type=bs.type,s.dim=dim,p.order=p.order,full.formula=as.formula(rf))
#  if (sfok!=0) ret$sf<-as.formula(sf)
  ret
}



gam.setup<-function(formula,data=list(),gam.call,predict=TRUE,parent.level=1,nsdf=-1)

# This gets the data referred to in the model formula, either from data frame "data"
# or from the level parent.level parent (e.g. when called from gam() this will be
# the parent that called gam(). 
# G$names[i] contains the names associated with each column of the design matrix,
# followed by the name of each smooth term - the constant is un-named! 
# Note that it is assumed that the default parent.level assumes that this function
# is called from the function that was called by the user! i.e. the user's calling 
# environment is the grandparent of the environment of gam.setup() 
# if nsdf==-1 then it assumed that the design matrix for the non-spline part of the 
# model is to be found here. Otherwise it is assumed that nsdf is the known 
# number of degrees of freedom for the non-spline part (including intercept), but
# the columns corresponding to the non-spline part of the model and the offset 
#  are set to zero.
# 
{ # now split the formula
  split<-gam.parser(formula,parent.level=parent.level+1) 
  dmiss<-missing(data)  
  if (dmiss) data<-sys.frame(sys.parent(n=parent.level))
  if (split$df[1]==-1)
  { if (split$pfok==0) stop("You've got no model....")
    m<-0
  }  
  else  m<-length(split$df) # number of smooth terms
  G<-list(m=m,df=split$df,full.formula=split$full.formula)
  G$fix<-split$fix
  add.constant<-FALSE
  if (split$pfok)   # deal with the strictly parametric part of the model 
  { if (length(attr(terms(split$pf),"variables"))==1) # then the formula is ~ 1 and model.frame(), etc can't cope 
    { G$nsdf<-1
      add.constant<-TRUE  # must construct model matrix later when length of response is known
    } else
    { if (nsdf>=0)  # set G$nsdf to supplied value, but deal with X later
      { G$nsdf<-nsdf
      } else 
      { if (predict) # then user MUST have supplied all data in data frame and searching is easy
        { mf<-model.frame(split$pf,data)
        } else   # all sorts of stupid practices are allowed and must be catered for.... 
        { mf<-gam.call # assemble a call to model frame which will allow searching in the calling env. 
          mf[[1]]<-as.name("model.frame")      # replace "gam" with "model.frame" in call
          mf$family<-mf$weights<-mf$control<-mf$scale<-NULL
          mf$formula<-split$pf                 # set formula to parametric only part     
          # evaluate in the grandparent environment to get model frame......
          mf<-eval(mf,sys.frame(sys.parent(n = parent.level))) 
        }
        G$offset <- model.offset(mf)   # get the model offset (if any)
        if (!is.null(G$offset) && length(attr(terms(split$pf),"variables"))<=2) # then there is an offest but no terms other than "+1" or "-1"
        { if (length(grep("-",as.character(split$pf[2])))>=1) # then there is a "-1" term in formula
          X<-matrix(0,length(G$offset),0)
          else X <- model.matrix(split$pf,mf)  # then there is a constant in the model and model.matrix() can cope 
        }
        else X <- model.matrix(split$pf,mf)       # getting the model matrix
        G$nsdf <- dim(X)[2]                  # extract the total non- smooth df
      }        
    }	
  } else
  { if (m==0) Stop("Model appears to have no terms") 
    G$nsdf<-0
  }
  #xvars <- as.character(attr(mt, "variables"))[-1]
  if (!predict) # obtain the "y variable"  
  { # evaluates rhs of model using data and failing that the calling environment...
    #G$y<-eval(terms(split$sf)[[2]],data,sys.frame(sys.parent(n = parent.level)))
    G$y<-eval(parse(text=split$v.names[1]),data,sys.frame(sys.parent(n = parent.level)))
    # checks that this has worked.....
    if (is.null(G$y)) stop(paste("Failed to find variable for",split.v.names))
  }
  if (m>0) 
  { # now find number of data without looking at y (needed for predicting)
    v.names<-split$v.names          #all.vars(split$sf)  # getting the list of variable names for the smooths
    # need to know number of rows in design matrix .....
    z<-eval(parse(text=v.names[2]),data,sys.frame(sys.parent(n = parent.level))) 
    G$n<-NROW(z)
  } else
  { if (add.constant) # model only has ~1
    { if (predict) stop("Model only has a constant - prediction doesn't need this call!!")
      G$n<-length(G$y) # get n from response
    } else G$n<-dim(X)[1] # get n from covariates
  }  
  if (add.constant)          # then X must be created by hand for case where parametric model is ~1
  { X<-matrix(1,G$n,G$nsdf)
    colnames(X)<-"constant"
  }
  if (nsdf>=0) # then X still needs to be filled
  { X<-matrix(0,G$n,G$nsdf)
    colnames(X)<-rep("",G$nsdf)         # don't need them, but will fail below without!  
    G$offset<-array(0,G$n)
  }
  xcols<-G$nsdf
  if (m>0) xcols<-xcols+length(v.names)-1 # -1 is to remove ~y   
  G$x<-array(0,dim=c(xcols,G$n)) # explanatory data array
  G$names<-array("",G$nsdf+m)
  G$s.type<-array(0,m)  # array to store smooth types
  G$dim<-split$s.dim     # gam parser produces this array of smooth dimensions
  G$p.order<-split$p.order # gam.parser's array of penalty orders
  if (G$nsdf>0) for (k in 1:G$nsdf) 
  { G$x[k,]<-X[,k] # copy columns of parametric design matrix
    G$names[k]<-colnames(X)[k] # copy column names for later use as labels
  }
  kk<-k<-G$nsdf+1
  jj<-1
  v.off<-2  # must be 2 to avoid response!!
  if (m>0)
  for (jj in 1:m)  # loop through smooths
  { vlist<-v.names[v.off:(v.off+G$dim[jj]-1)];v.off<-v.off+G$dim[jj];
    if (jj==1) 
    {  if (G$nsdf>0) G$vnames<-c(colnames(X),vlist) else
       G$vnames<-vlist 
    } else G$vnames<-c(G$vnames,vlist)
    G$names[kk]<-"s("
    for (i in 1:G$dim[jj]) # produces a column for each variable in this smooth
    { z<-eval(parse(text=vlist[i]),data,sys.frame(sys.parent(n = parent.level)))
     
      if (is.null(z)) stop(paste("Failed to find variable",vlist[i]))  
      G$x[k,]<-z
      if (i<length(vlist)) G$names[kk]<-paste(G$names[kk],vlist[i],",",sep="")
      else G$names[kk]<-paste(G$names[kk],vlist[i],sep="")
      k<-k+1
    }
    G$s.type[jj]<-split$bs.type[jj] 
    G$names[kk]<-paste(G$names[kk],")",sep="")
    kk<-kk+1
  }
  # now sort out max dfs per term (only if get.y - otherwise call from predict)
  if (!predict&&m>0)
  { k<-0;off<-G$nsdf
    for (i in 1:G$m)
    { xx<-G$x[(off+1):(off+G$dim[i]),]
      M<-null.space.dimension(G$dim[i],G$p.order[i])
      if (G$dim[i]==1) n<-length(unique(xx))
      else n<-dim(uniquecombs(xx))[2]
      if (M>=n) stop("There are too few unique covariate combinations to support one of your model terms")
      if (G$df[i]>n) 
      { warning("Max d.f. reset to number of unique covariate combinations for one of your model terms")
        G$df[i]<-n
      } else
      if (G$df[i]<=M)
      { warning("Max d.f. for a term must be greater than dimension of null space of penalty - max. d.f. has been raised for a term");
        G$df[i]<-M+1
      }
      k<-k+G$df[i]
      off<-off+G$dim[i] 
    }
    xx<-G$x[(G$nsdf+1):off,]
    if (off-G$nsdf==1) n<-length(unique(xx))
    else n<-dim(uniquecombs(xx))[2]
    if (k>n) stop("You have specified more degrees of freedom than you have unique covariate combinations - use uniquecombs() to investigate further") 
  }
  G
}


gam<-function(formula,family=gaussian(),data=list(),weights=NULL,control=gam.control(),scale=0)

# Routine to fit a GAM to some data. The model is stated in the formula, which is then 
# parsed to figure out which bits relate to smooth terms and which to parametric terms.
# -ve binomial stuff by Mike Lonergan, rest by Simon Wood.

{ gam.call<-match.call()  # store the call to facilitate searching in gam.setup()

  family<-get.family(family) # deals with -ve binomial handling as well as any coercion

  if (missing(data)) G<-gam.setup(formula,gam.call=gam.call,parent.level=2,predict=FALSE)
  else G<-gam.setup(formula,data,gam.call=gam.call,parent.level=2,predict=FALSE)
 
  G<-GAMsetup(G) 
 
  if (is.null(G$offset)) G$offset<-rep(0,G$n)
  if (is.null(weights)) G$w<-rep(1,G$n) else G$w<-weights
  
  if (scale==0) 
  { if (family$family=="binomial"||family$family=="poisson"
         || substr(family$family,1,17)=="Negative Binomial") scale<-1 #ubre
    else scale <- -1 #gcv
  }

  G$sig2<-scale
  G$sp<-array(-1,G$m) # set up smoothing parameters for autoinitialization at first run
  G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
  G$max.half<-control$mgcv.max.half # max step halving in Newton update mgcv
  if(family$family=="Negative Binomial(NA)") object<-gam.nbut(G, family$link, control, scale)
  else object<-gam.fit(G,family=family,control=control)
  object$covariate.shift<-G$covariate.shift # basis stabilizing linear translations  
  object$xp<-G$xp  # position of knots for cubic regression splines
  object$UZ<-G$UZ  # truncated basis for tprs terms
  object$Xu<-G$Xu # unique covariate combinations for tprs terms
  object$xu.length<-G$xu.length # array of numbers of unique covariate combinations for tprs terms  
  # return correct d.f. for fixed d.f. terms
  if (G$m>0) for (i in 1:G$m) if (G$fix[i]) object$edf[i]<-G$df[i]-1
  # now re-assign variable names to coefficients etc. 
  term.names<-array("",length(object$coefficients))
  if (G$nsdf>0) for (i in 1:G$nsdf) term.names[i]<-G$names[i]
  i<-G$nsdf+1
  if (G$m>0)
  for (k in 1:G$m)
  for (j in 1:G$df[k])
  { term.names[i]<-paste(G$names[G$nsdf+k],".",as.character(j),sep="")
    i<-i+1
  }
  names(object$coefficients)<-term.names  # note - won't work on matrices!!
  object$full.formula<-G$full.formula
  object$formula<-formula
  object$x<-G$x
  object$s.type<-G$s.type
  object$p.order<-G$p.order
  object$dim<-G$dim
  object$call<-gam.call
  rownames(object$x)<-G$vnames
  names(object$x)<-NULL
  class(object)<-"gam"

  object
}

print.gam<-function (x,...) 
# default print function for gam objects
{ print(x$family)
  cat("Formula:\n")
  print(x$formula)
  if (x$dim==0)
  cat("Total model degrees of freedom",x$nsdf,"\n")
  else
  cat("\nEstimated degrees of freedom:\n",x$edf,"  total = ",sum(x$edf)+x$nsdf,"\n")
  gcv<-x$df.null*x$sig2/(x$df.null-sum(x$edf)-x$nsdf)
  cat("\nGCV score: ",gcv,"\n")
}

gam.control<-function (epsilon = 1e-04, maxit = 30,mgcv.tol=1e-6,mgcv.max.half=15, trace = FALSE) 
# control structure for a gam
{   if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace, mgcv.tol=mgcv.tol,mgcv.max.half=mgcv.max.half)
}

gam.fit<-function (G, start = NULL, etastart = NULL, 
    mustart = NULL, family = gaussian(), 
    control = gam.control()) 
# fitting function for a gam, modified from glm.fit
{
    #x <- as.matrix(x)
    #xnames <- dimnames(x)[[2]]
    #ynames <- names(y)
    conv <- FALSE
    nobs <- NROW(G$y)
    nvars <- NCOL(G$X) # check this needed
    y<-G$y # original data
    X<-G$X # original design matrix
    if (nvars == 0) stop("Model seems to contain no terms")
    if (family$family=="gaussian" && family$link=="identity") olm<-TRUE # i.e. only one iteration needed
    else olm<-FALSE 

    weights<-G$w # original weights
   
    offset<-G$offset 

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
    if (is.null(mustart)) 
        eval(family$initialize, sys.frame(sys.nframe()))
    
    if (NCOL(y) > 1) 
        stop("y must be univariate unless binomial")
    eta <- if (!is.null(etastart) && valideta(etastart)) 
        etastart
    else if (!is.null(start)) 
        if (length(start) != nvars) 
            stop(paste("Length of start should equal", nvars ))
        else as.vector(if (NCOL(G$X) == 1) 
            G$X * start
        else G$X %*% start)
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) 
        stop("Can't find valid starting values: please specify some")
    devold <- sum(dev.resids(y, mu, weights))
    coefold <- start
    boundary <- FALSE
    scale<-G$sig2
    
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
        good <- (weights > 0) & (mu.eta.val != 0) # note good modified here => must re-calc each iter
        if (all(!good)) {
            conv <- FALSE
            warning(paste("No observations informative at iteration", 
                iter))
            break
        }
   
        z<-G$y <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w<-G$w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        G$w<-G$w^2 # this line is somewhat important
        G$X<-X[good,]  # truncated design matrix       
        ngoodobs <- as.integer(nobs - sum(!good))
        ncols <- as.integer(1)
        # must set G$sig2 to scale parameter or -1 here....
        G$sig2<-scale
        G<-mgcv(G) 
		
        start <- coef <- G$p
       
        eta[good] <- drop(X[good, , drop = FALSE] %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace) 
            cat("Deviance =", dev, "Iterations -", iter, "\n")
        boundary <- FALSE
        if (is.na(dev) || any(is.na(coef))) {
            warning("Step size truncated due to divergence")
            ii <- 1
            while ((is.na(dev) || any(is.na(start)))) {
                if (ii > control$maxit) 
                  stop("inner loop 1; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                eta[good] <- drop(X[good, , drop = FALSE] %*% 
                  start)
                mu <- linkinv(eta <- eta + offset)
                dev <- sum(dev.resids(y, mu, weights))
            }
            boundary <- TRUE
            coef <- start
            if (control$trace) 
                cat("New Deviance =", dev, "\n")
        }
        if (!(valideta(eta) && validmu(mu))) {
            warning("Step size truncated: out of bounds.")
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit) 
                  stop("inner loop 2; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                eta[good] <- drop(X[good, , drop = FALSE] %*% 
                  start)
                mu <- linkinv(eta <- eta + offset)
            }
            boundary <- TRUE
            coef <- start
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("New Deviance =", dev, "\n")
        }
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon || olm) {
            conv <- TRUE
            break
        }
        else {
            devold <- dev
            coefold <- coef
        }
    }
    if (!conv) 
        warning("Algorithm did not converge")
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
   
    wtdmu <- sum(weights * y)/sum(weights)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok 

	list(coefficients = as.vector(coef), residuals = residuals, fitted.values = mu, 
        family = family, 
        linear.predictor = eta, deviance = dev,
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        #df.residual = resdf, 
        df.null = nulldf, y = y, converged = conv,sig2=G$sig2,edf=G$edf,
        boundary = boundary,sp = G$sp,df=G$df,nsdf=G$nsdf,Vp=G$Vp)
}


predict.gam<-function(object,newdata,type="link",se.fit=FALSE,plot.call=FALSE,...) {

# This function is used for predicting from a GAM. object is a gam object, newdata a dataframe to
# be used in prediction..... it's all done via a call to the compiled C routine RGAMpredict().
#
# Type == "link"     - for linear predictor
#      == "response" - for fitted values
#      == "terms"    - for individual terms on scale of linear predictor 
#      == "lpmatrix" - for matrix mapping parameters to l.p.
#  Steps are:
#  1. Construct x[i,j] - the explanatory data array for which predictions are required, this means
#     using the model formulae from the gam object, "object" to extract the model variables from 
#     the data-frame "newdata". Note np - the number of elements to predict for.
#  2. Convert the knot position data to a form usable by the C routine RGAMpredict(). 
#  3. Initialize storage space for eta and se.   
#  4. Call RGAMpredict()
#  5. Use eta and se to construct the returned vector, matrix or list.
#  6. Tidy up and return.  

  if (type!="link"&&type!="terms"&&type!="response"&&type!="lpmatrix")  
  { warning("Unknown type, reset to terms.")
    type<-"terms"
  }
  if (class(object)!="gam") stop("predict.gam can only be used to predict from gam objects")

  # get data from which to predict.....  
  
  if (missing(newdata)) # then "fake" an object suitable for prediction 
  { if (object$dim==0) m<-0
    else m<-length(object$sp)
    n<-length(object$y)
    x<-object$x
    for (i in 1:m) x[object$nsdf+i,]<-x[object$nsdf+i,]+object$covariate.shift[i]
    G<-list(x=x,nsdf=object$nsdf,m=m,n=n,dim=object$dim)
    no.data<-FALSE
  }
  else 
  { if (plot.call) # then it's a call from plot, and only spline parts are to be evaluated
    G<-gam.setup(object$full.formula,newdata,gam.call=object$call,parent.level=3,predict=TRUE,nsdf=object$nsdf)
    else G<-gam.setup(object$full.formula,newdata,gam.call=object$call,parent.level=2,predict=TRUE)    
    no.data<-FALSE
  }
  np<-G$n
  
  # now set up the other information for prediction.....
  control<-0
  if (type=="lpmatrix")
  { X<-matrix(0,np,length(object$coefficients))
    eta<-se<-0;control<-4
  } else 
  if (type=="terms") 
  { control<-2
    eta<-array(0,dim=c(G$nsdf+G$m,np))
    se<-array(0,dim=c(G$nsdf+G$m,np))
    X<-0
  } else
  { eta<-array(0,dim=c(np))
    se<-array(0,dim=c(np))
    X<-0
  }
  if (se.fit&&object$Vp[1,1]<=0)
  { se.fit<-FALSE
    warning("No variance estimates available")
  }
  if (se.fit&&control<4) control<-control+1
  if (G$m)
  { nc<-0;for (i in 1:G$m) nc<-nc+G$dim[i]
    for (i in 1:nc) G$x[i+G$nsdf,] <- G$x[i+G$nsdf,] - object$covariate.shift[i]  # stabilizing shift
  }
  # call compiled code for prediction ....... NOTE: change for new basis....
  o<-.C("RGAMpredict", as.integer(object$xu.length),as.double(object$Xu),as.double(object$UZ),
             as.double(object$xp),as.integer(G$nsdf),as.integer(object$dim),as.integer(object$s.type),
            as.integer(object$df),as.integer(object$p.order), as.integer(G$m),as.integer(length(object$y)),as.double(G$x),
            as.integer(np),as.double(object$coefficients),as.double(object$Vp),as.double(eta),as.double(se),as.double(X),
            as.integer(control))
  # now eta contains linear predictor (terms) and se may contain corresponding standard errors, or
  # X contains the matrix mapping params to l.p.    
  if (type=="lpmatrix")
  { H<-matrix(o[[18]],np,length(object$coefficients)) 
  } else  
  { if (type=="terms")
    { eta<-array(o[[16]],c(G$nsdf+G$m,np));
      se<-array(o[[17]],c(G$nsdf+G$m,np));
    } else
    { eta<-array(o[[16]],c(np));
      se<-array(o[[17]],c(np));
      if (type=="response") # transform onto scale of data
      { fam<-object$family;linkinv<-fam$linkinv;dmu.deta<-fam$mu.eta
        se<-se*abs(dmu.deta(eta)) 
        eta<-linkinv(eta)
      }
    }
    if (se.fit)
    { H<-list(fit=eta,se.fit=se) }
    else
    { H<-eta}
  }
  # tidy up? 
   
  H # ... and return
}

plot.gam<-function(x,rug=TRUE,se=TRUE,pages=0,scale=-1,n=100,n2=40,pers=FALSE,theta=30,phi=30,...)

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
{ sp.contour<-function (x,y,z,zse,xlab="",ylab="",zlab="main",se.plot=TRUE)   
  # internal function for contouring 2-d smooths with 1 s.e. limits
  { gap<-median(zse)  
    zr<-max(z+zse)-min(z-zse) # plotting range  
    n<-10  
    while (n>1 && zr/n<2.5*gap) n<-n-1    
    zrange<-c(min(z-zse),max(z+zse))  
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
    contour(x,y,z+zse,levels=zlev,add=T,lty=2,col=2,labcex=cs*0.5)  
    contour(x,y,z-zse,levels=zlev,add=T,lty=3,col=3,labcex=cs*0.5)  
    
    xpos<-xrange[1]  
    xl<-c(xpos,xpos+xr/10);yl<-c(ypos,ypos)  
    lines(xl,yl,xpd=TRUE,lty=3,col=3)  
    text(xpos+xr/10,ypos,"-se",xpd=TRUE,pos=4,cex=cs,off=0.5*cs)  
  
    xpos<-xpos+3*xr/10  
    xl<-c(xpos,xpos+xr/10);  
    lines(xl,yl,xpd=TRUE,lwd=2)  
    text(xpos+xr/10,ypos,zlab,xpd=TRUE,pos=4,cex=cs,off=0.5*cs)  
    
    xpos<-xrange[2]-xr/5  
    xl<-c(xpos,xpos+xr/10);  
    lines(xl,yl,xpd=TRUE,lty=2,col=2)  
    text(xpos+xr/10,ypos,"+se",xpd=TRUE,pos=4,cex=cs,off=0.5*cs)  
  }   



  
  fix.form<-function(formula)
  # internal function to replace "$" on rhs of formula with ".", to allow construction of 
  # data frame for prediction.....  
  { ft<-as.character(formula)
    work.around<-as.character(deparse(formula[3]))
    ft[3]<-""
    for (i in 1:length(work.around)) ft[3]<-paste(ft[3],work.around[i],sep="")
    form<-paste(ft[2],ft[1])
    for (i in 1:(nchar(ft[3])-2))
    { c<-substring(ft[3],i,i)
      if (c=="$") c<-"."
      form<-paste(form,c,sep="")
    }
    as.formula(form)
  }
  
  rename<-function (x)
  # internal function for renaming rows of x$x
  { nm<-rownames(x)
    for (j in 1:length(nm))
    { tmp<-nm[j]
      nm[j]<-""
      for (i in 1:nchar(tmp))
      { c<-substring(tmp,i,i)
        if (c=="$") c<-"."
        nm[j]<-paste(nm[j],c,sep="")
      }
    }
    nm
  }

  # start of main function

  if (x$dim==0) stop("Model has no smooth terms - nothing for plot.gam() to do.")
  
  x$formula<-fix.form(x$formula)  # replace all "$" characters on rhs of formula
  rownames(x$x)<-rename(x$x)     # replace all "$" characters in row names of x$x
  
  m<-length(x$df) # number of smooth terms
  if (se && x$Vp[1,1]<=0) 
  { se<-FALSE
    warning("No variance estimates available")
  }
  # sort out number of pages and plots per page
  if (pages>m) pages<-m
  if (pages<0) pages<-0
  if (pages!=0)    # figure out how to display things
  { ppp<-m%/%pages
    if (m%%pages!=0) 
    { ppp<-ppp+1
      while (ppp*(pages-1)>=m) pages<-pages-1
      if (m%%pages) last.pages<-0 else last.ppp<-m-ppp*pages
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
	ylim<-c(r,c)
  } else
  { ppp<-1}
  
  # now create dataframe for 1-d terms to send to predict.gam() to return smooth terms
  xx1<-data.frame(1:n)
  if (x$nsdf>0) 
  for (i in 1:x$nsdf) xx1[i]<-seq(1:n)*0
  j<-1;md2<-FALSE
  for (i in 1:m)
  { #x0<-x$xp[i,1]       # min x
    #x1<-x$xp[i,x$df[i]] # max x
    if (x$dim[i]==1)
    { x0<-min(x$x[j+x$nsdf,])
      x1<-max(x$x[j+x$nsdf,])
      dx<-(x1-x0)/(n-1) 
      xx1[j+x$nsdf]<-seq(x0,x1,dx)+x$covariate.shift[j]
      j<-j+1
    } else
    { for (k in 1:x$dim[i])
      { xx1[j+x$nsdf]<-rep(0,n);j<-j+1 # multi-dim terms padded with zeroes
        if (x$dim[i]==2) md2<-TRUE
      }  
    }  
  }
  names(xx1)<-rownames(x$x)

  pl1<-predict.gam(x,xx1,type="terms",se,plot.call=TRUE) 

  if (md2) # then create data frames for 2-d plots
  { if (n2<10) n2<-10
    xx2<-data.frame(1:(n2*n2))
    xm<-data.frame(1:n2);ym<-data.frame(1:n2)
    if (x$nsdf>0) 
    for (i in 1:x$nsdf) xx2[i]<-rep(0,n2*n2)
    j<-1
    for (i in 1:m)
    { if (x$dim[i]==2)
      { x0<-min(x$x[j+x$nsdf,]);x1<-max(x$x[j+x$nsdf,])
        dx<-(x1-x0)/(n2-1);
        y0<-min(x$x[j+x$nsdf+1,]);y1<-max(x$x[j+x$nsdf+1,])
        dy<-(y1-y0)/(n2-1)
        
        xm[i]<-seq(x0,x1,dx)+x$covariate.shift[j]
        ym[i]<-seq(y0,y1,dy)+x$covariate.shift[j+1]
        xx2[j+x$nsdf]<-rep(xm[[i]],n2)
        xx2[j+x$nsdf+1]<-rep(ym[[i]],rep(n2,n2))
        j<-j+2
      } else
      {  for (k in 1:x$dim[i])
         { xx2[j+x$nsdf]<-rep(0,n2*n2);j<-j+1
           xm[i]<-rep(0,n2);ym[i]<-xm[i];
         }

      }
    } 
    names(xx2)<-rownames(x$x)
    pl2<-predict.gam(x,xx2,type="terms",se,plot.call=TRUE)
  }
  if (se)   # pl$fit and pl$se.fit
  { k<-0
    if (scale==-1) # getting common scale for 1-d terms
    for (i in 1:m)
    { if (x$dim[i]==1)
      { ul<-pl1$fit[x$nsdf+i,]+2*pl1$se.fit[x$nsdf+i,]
        ll<-pl1$fit[x$nsdf+i,]-2*pl1$se.fit[x$nsdf+i,]
        if (k==0) { ylim<-c(min(ll),max(ul));k<-1;}
        else
        { if (min(ll)<ylim[1]) ylim[1]<-min(ll)
	  if (max(ul)>ylim[2]) ylim[2]<-max(ul)
        }
      }
    }
    j<-1
    for (i in 1:m)
    { if (interactive() && x$dim[i]<3 && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
      if (x$dim[i]==1)
      { ul<-pl1$fit[x$nsdf+i,]+2*pl1$se.fit[x$nsdf+i,]
        ll<-pl1$fit[x$nsdf+i,]-2*pl1$se.fit[x$nsdf+i,]
        if (scale==0) { ylim<-c(min(ll),max(ul))}
        title<-paste("s(",rownames(x$x)[j+x$nsdf],",",as.character(round(x$edf[i],2)),")",sep="")
        plot(xx1[[x$nsdf+j]],pl1$fit[x$nsdf+i,],type="l",xlab=rownames(x$x)[j+x$nsdf],ylim=ylim,ylab=title)
	lines(xx1[[x$nsdf+j]],ul,lty=2)
        lines(xx1[[x$nsdf+j]],ll,lty=2)
	if (rug) rug(as.numeric(x$x[x$nsdf+j,]+x$covariate.shift[j]))
      } else if (x$dim[i]==2)
      { xla<-rownames(x$x)[j+x$nsdf];yla<-rownames(x$x)[j+x$nsdf+1]
        title<-paste("s(",xla,",",yla,",",as.character(round(x$edf[i],2)),")",sep="")
        if (pers) persp(xm[[i]],ym[[i]],matrix(pl2$fit[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,zlab=title,theta=theta,phi=phi)
        else
        { sp.contour(xm[[i]],ym[[i]],matrix(pl2$fit[x$nsdf+i,],n2,n2),matrix(pl2$se.fit[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,zlab=title)
          if (rug) points(x$x[j+x$nsdf,]+x$covariate.shift[j],x$x[j+1+x$nsdf,]+x$covariate.shift[j+1],pch=".")
        } 
      } else
      { warning("no automatic plotting for smooths of more than one variable")
      }  
      j<-j+x$dim[i]
    }
  } else # don't plot confidence limits
  { k<-0
    if (scale==-1)
    for (i in 1:m)
    { if (x$dim[i]==1)
      { if (k==0) { ylim<-range(pl1[x$nsdf+i,]);k<-1 }
	else
        { if (min(pl1[x$nsdf+i,])<ylim[1]) ylim[1]<-min(pl1[x$nsdf+i,])
	  if (max(pl1[x$nsdf+i,])>ylim[2]) ylim[2]<-max(pl1[x$nsdf+i,])
	}
      }
    }
    j<-1
    for (i in 1:m)
    { if (interactive() && x$dim[i]<3 && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
      if (x$dim[i]==1)
      { title<-paste("s(",rownames(x$x)[j+x$nsdf],",",as.character(round(x$edf[i],2)),")",sep="")
        if (scale==0) ylim<-range(pl1[x$nsdf+i,])
        plot(xx1[[x$nsdf+j]],pl1[x$nsdf+i,],type="l",,xlab=rownames(x$x)[j+x$nsdf],ylab=title,ylim=ylim)
        rug(as.numeric(x$x[x$nsdf+j,]+x$covariate.shift[j]))
      } else if (x$dim[i]==2)
      { xla<-rownames(x$x)[j+x$nsdf];yla<-rownames(x$x)[j+x$nsdf+1]
        title<-paste("s(",xla,",",yla,",",as.character(round(x$edf[i],2)),")",sep="")
       # sp.contour(xm[[i]],ym[[i]],matrix(pl2$fit[x$nsdf+i,],n2,n2),matrix(pl2$se.fit[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,zlab=title,se=FALSE)
        if (pers) persp(xm[[i]],ym[[i]],matrix(pl2[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,zlab=title,theta=theta,phi=phi)
        else
        { contour(xm[[i]],ym[[i]],matrix(pl2[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,main=title)
          if (rug) points(x$x[j+x$nsdf,]+x$covariate.shift[j],x$x[j+1+x$nsdf,]+x$covariate.shift[j+1],pch=".")
        }  

      } else
      { warning("no automatic plotting for smooths of more than one variable")}
      j<-j+x$dim[i]
    } 
  }
  if (pages>0) par(oldpar)
}


residuals.gam <-function(object, type = c("deviance", "pearson", "working", "response"),...)
# calculates residuals for gam object - defualt for glm (from which this is adapted) seems to be buggy
{ type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  family <- object$family
  wts <- object$prior.weights
  switch(type,working = object$residuals,
              pearson = (y-mu)/sqrt(object$sig2*object$family$variance(mu)),
              deviance = { d.res<-sqrt(pmax(object$family$dev.resids(y,mu,wts),0))
                           ifelse(y>mu , d.res, -d.res)             
                         },
              response = y - mu)
}

#####################################################################################################
# Code by Mike Lonergan, from here
#####################################################################################################

persp.gam<-function(x, view=NULL, slice=list(), sizes=c(20,20), mask=FALSE,
se=2,theta=0,phi=15, r = sqrt(3), d = 1, scale = TRUE, 
                    expand = 1, col = NULL,border = NULL,  ltheta = -135, lphi = 45, 
                    shade = 0.75, box = TRUE, axes = TRUE,  nticks = 5,
                    ticktype = "detailed",...)
# GAM visualization routine. Author: Mike Lonergan.  
{

   se<-as.numeric(se)
   variables<-rownames(x$x)[-1]  
   which.slice<-""

   if (length(variables)<2)  
       plot(x)
   else   {
      if (is.null(view) )
         view<-c(variables[1],variables[2])

      view.index<-pmatch(view,variables)

      for(v in 1:length(variables))   {
         if (eval(parse(text=paste("is.null(slice$",variables[v],")",sep=""))))
            eval(parse(text=paste("slice$",variables[v],
                    "<-mean(x$x[v+1,])+x$covariate.shift[v]",sep="")))
      }

      view.data1<-x$x[view.index[1]+1,]+x$covariate.shift[view.index[1]]
      view.data2<-x$x[view.index[2]+1,]+x$covariate.shift[view.index[2]]

      gridmaker<-paste("expand.grid(", view[1], "=seq(", min(view.data1), ",",
                           max(view.data1), ", length=", sizes[1], "),",
                       view[2], "=seq(", min(view.data2), ",",
                           max(view.data2), ", length=", sizes[2], "))", sep="")
      persp.dat<-eval(parse(text=gridmaker))

      attach(persp.dat)

      for (v in variables)   {
          if(eval(parse(text=paste("length(persp.dat$",v,")==0",sep=""))))   {
             this.slice<-eval(parse(text=slice[[v]]))
             if (length(this.slice)==1)
                this.slice<-rep(this.slice,dim(persp.dat)[1])
             eval(parse(text=paste("persp.dat$",v,"<-this.slice",sep="")))
             sl<-slice[[v]]
             if (is.numeric(sl)) sl<-signif(sl,3) 
             which.slice<-paste(which.slice, v, "=", sl , ", ",sep="") 
          }
      }

      if (se>0)
         which.slice<-paste(which.slice," red/green are +/-",se," se",sep="")
      else
         which.slice<-substr(which.slice,1,nchar(which.slice)-1)

      detach(persp.dat)

      if (se<=0)
        
persp.dat$response<-as.vector(predict(x,newdata=persp.dat,type="response"))   
   else   {          
response<-predict(x,newdata=persp.dat,type="response",se.fit=TRUE)          
persp.dat$response<-as.vector(response[[1]])          
persp.dat$se<-as.vector(response[[2]])       }

      if (mask!=FALSE)   {   # crudely removes parts of surface far from data points
                                            # in view dimensions
         if (is.numeric(mask))   {
            if (length(mask)<2) mask<-c(mask,mask)
         } else  {
            xd<-sort(view.data1)
            xd<-abs(xd[-1]-xd[-length(xd)])
            mask<-min(pmax(xd[-1]-xd[-length(xd)]))
            xd<-sort(view.data2)
            xd<-abs(xd[-1]-xd[-length(xd)])
            mask<-c(mask,min(pmax(xd[-1]-xd[-length(xd)])))
         }
         for (i in 1:length(persp.dat$response))   {
            if (min(abs(view.data1-persp.dat[i,1])/mask[1] 
                + abs(view.data2-persp.dat[i,2])/mask[2]) >1)
               persp.dat$response[i]<-NA
         }
      }
      fit.mat<-matrix(persp.dat$response,sizes[1],sizes[2])
      if (se>0)   {
         fit.mat.u<-matrix((persp.dat$response+se*persp.dat$se),sizes[1],sizes[2])
         fit.mat.l<-matrix((persp.dat$response-se*persp.dat$se),sizes[1],sizes[2])
         zlim<-c(min(persp.dat$response-se*persp.dat$se,na.rm=TRUE),
                           max(persp.dat$response+se*persp.dat$se,na.rm=TRUE))

        persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]),
               y = seq(min(view.data2), max(view.data2), len = sizes[2]), z=fit.mat,
               xlab = view[1], ylab = view[2], zlab = "", zlim=zlim,
               sub=which.slice, 
               theta=theta, phi=phi, r=r, d=d, scale=scale, expand=expand, col=col, 
               border=NA, ltheta=ltheta, lphi=lphi, shade=shade,
               box=box, axes=axes, nticks=nticks, ticktype=ticktype) 
 
        par(new=TRUE)
        persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]),
               y = seq(min(view.data2), max(view.data2), len = sizes[2]), z=fit.mat.l,
               xlab = view[1], ylab = view[2], zlab = "", zlim=zlim,
               sub=which.slice, 
               theta=theta, phi=phi, r=r, d=d, scale=scale, expand=expand, col=NA, 
               border="green", shade=NA,
               box=box, axes=axes, nticks=nticks, ticktype=ticktype) 

         par(new=TRUE)
         persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]),
               y = seq(min(view.data2), max(view.data2), len = sizes[2]), z=fit.mat,
               xlab = view[1], ylab = view[2], zlab = "", zlim=zlim,
               sub=which.slice, 
               theta=theta, phi=phi, r=r, d=d, scale=scale, expand=expand, col=NA, 
               border=border, shade=NA,
               box=box, axes=axes, nticks=nticks, ticktype=ticktype) 
         par(new=TRUE)
         persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]),
               y = seq(min(view.data2), max(view.data2), len = sizes[2]), z=fit.mat.u,
               xlab = view[1], ylab = view[2], zlab = "", zlim=zlim,
               sub=which.slice, 
               theta=theta, phi=phi, r=r, d=d, scale=scale, expand=expand, col=NA, 
               border="red", shade=NA,
               box=box, axes=axes, nticks=nticks, ticktype=ticktype) 

         invisible(list(fit.mat,fit.mat.l,fit.mat.u))
      }      
      else   {
         persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]),
               y = seq(min(view.data2), max(view.data2), len = sizes[2]), z=fit.mat,
               xlab = view[1], ylab = view[2], zlab = "",sub=which.slice, 
               theta=theta, phi=phi, r=r, d=d, scale=scale, expand=expand, col=col, 
               border=border, ltheta=ltheta, lphi=lphi, shade=shade,
               box=box, axes=axes, nticks=nticks, ticktype=ticktype) 

        invisible(fit.mat)
     }
   }
}

gam.nbut<-
function (G, link = log, control = gam.control(), scale = 0) 
{                                               
# Basically a modified version of glm.nb() from MASS (c) Venables and Ripley,
# suitable  for use with gam(). Modifications by Mike Lonergan. 
# This version does entire gam fitting each time between theta estimating

   loglik <- function(n, th, mu, y) {
           sum(lgamma(th + y) - lgamma(th) + th * log(th) + y * 
               log(mu + (y == 0)) - (th + y) * log(th + mu))
    }

   fam<- do.call("poisson", list(link = link))
   object <- gam.fit(G, family = fam, control = control)

   th <- as.vector(theta.maxl(object$y, object$fitted, limit = control$maxit, 
        trace = control$trace > 2))
   if (control$trace > 1) 
     cat("First value for theta:", signif(th), "\n")
   iter <- 0
   d1 <- sqrt(2 * max(1, object$df.residual))
   d2 <- del <- 1
   Lm <- loglik(length(object$y), th, object$fitted, object$y)
   Lm0 <- Lm + 2 * d1

   while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - 
         Lm)/d1 + abs(del)/d2) > control$epsilon) {
      
      family <- do.call("neg.bin", list(theta = th, link = link))
      if (iter==1) etastart<-family$linkfun(object$fitted)
      object <- gam.fit(G, etastart = etastart, 
                           family = family, control = control)
      t0 <- th
      th <- theta.maxl(object$y, object$fitted, limit = control$maxit,
            trace = control$trace > 2)
      del <- t0 - th
      df.resid<-object$df.null-object$nsdf-sum(object$edf)
      d1 <- sqrt(2 * max(1, df.resid))
      Lm0 <- Lm
      Lm <- loglik(length(object$y), th, object$fitted, object$y)
      if (iter==1) Lmin<-Lm
      if (Lm<Lmin)
        { Lmin<-Lm
          etastart<-family$linkfun(object$fitted)
        }  
      if (control$trace) {
          Ls <- loglik(length(object$y), th, object$y, object$y)
          Dev <- 2 * (Ls - Lm)
          cat("Theta(", iter, ") =", signif(th), ", 2(Ls - Lm) =", 
              signif(Dev), "\n")
      }
  }
  object$theta <- as.vector(th)
  object$SE.theta <- attr(th, "SE")   
  if (!is.null(attr(th, "warn"))) 
     object$th.warn <- attr(th, "warn")
  if (iter > control$maxit) {
     warning("n.b. not converged: it's likely that this is because your data are not distinguishable from Poisson, or that your model has too many spurious covariates. It's usually the case that the model fit is never-the-less quite adequate for practical purposes.")
     object$th.warn <- "n.b. not converged"
  }
  return(object) 
}


   

get.family <- function(family)
# routine bike Mike Lonergan, separating out deling with family objects from gam()
{  if (is.character(family)) 
     family<-eval(parse(text=family))
   if (is.function(family)) 
      family <- family()
   if (is.null(family$family)) {
      print(family)
      stop("`family' not recognized")
   }
   family
}


neg.bin <-function(theta = NA, link = "log")
# Provides a negative binomial family for use with gam()
# Routine is slight modification of negative.binomial family provided
# in MASS library (c) Venables and Ripley. The modification (M. Lonergan) 
# is to make it painless to use in a call to gam() even when theta=NA
{    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link")
            linktemp <- eval(link)
    }
    if (any(linktemp == c("log", "identity", "sqrt")))
        stats <- make.link(linktemp)
    else stop(paste(linktemp, "link not available for negative binomial",
                    "family; available links are", "\"identity\", \"log\" and \"sqrt\""))
    .Theta <- theta
    stats <- make.link("log")
    variance <- function(mu)
        mu + mu^2/.Theta
    validmu <- function(mu)
        all(mu > 0)
    dev.resids <- function(y, mu, wt)
        2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) *
                  log((y + .Theta)/ (mu + .Theta)))
    aic <- function(y, n, mu, wt, dev) {
        term <- (y + .Theta) * log((y + .Theta)/ (mu + .Theta)) - y * log(mu) +
            lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) - lgamma(.Theta+y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0)) stop(paste("Negative values not allowed for",
                                   "the Poisson family"))
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
    famname <- paste("Negative Binomial(", format(round(theta, 4)), ")",
                     sep = "")
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                   aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                   validmu = validmu, valideta = stats$valideta), class = "family")
}


theta.maxl<-function (y, mu, n = length(y), limit = 10, eps =
.Machine$double.eps^0.25,      trace = FALSE) 
# Modification of MASS theta.ml() [(c) Venables and Ripley] for use with gam()
# The modification stops \theta -> \infty when data close to Poisson  
{
    score <- function(n, th, mu, y) sum(digamma(th + y) - digamma(th) + 
        log(th) + 1 - log(th + mu) - (y + th)/(mu + th))
    info <- function(n, th, mu, y) sum(-trigamma(th + y) + trigamma(th) - 
        1/th + 2/(mu + th) - (y + th)/(mu + th)^2)
    if (inherits(y, "lm")) {
        mu <- y$fitted
        y <- if (is.null(y$y)) 
            mu + residuals(y)
        else y$y
    }
    t0 <- n/sum((y/mu - 1)^2)
    it <- 0
    del <- 1
    max.t0<-max(y+1)^2*100
    if (trace) 
        cat("theta.maxl: initial theta =", signif(t0), "\n")
    while ((it <- it + 1) < limit && abs(del) > eps) {
        t0 <- abs(t0)
        del <- score(n, t0, mu, y)/(i <- info(n, t0, mu, y))
        t1 <- t0 + del
        if (t1>max.t0) t1<-max.t0
        del<-t1-t0;t0<-t1
        if (trace) 
            cat("theta.maxl: iter", it, " theta =", signif(t0), 
                "\n")
    }
    if (t0 < 0) {
        t0 <- 0
        warning("estimator truncated at zero")
        attr(t0, "warn") <- "estimate truncated at zero"
    }
    if (it == limit) {
        warning("iteration limit reached")
        attr(t0, "warn") <- "iteration limit reached"
    }
    attr(t0, "SE") <- sqrt(1/i)
    t0
}




#####################################################################################################
# end of Mike Lonergan code
#####################################################################################################

.First.lib <- function(lib, pkg) {
    library.dynam("mgcv", pkg, lib)
    cat("This is mgcv 0.6.2\n")
}











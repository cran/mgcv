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

oo<-.C("RQT",as.double(A),as.integer(nrow(A)),as.integer(ncol(A)),PACKAGE="mgcv")
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
  oo<-.C("RMonoCon",as.double(A),as.double(b),as.double(x),as.integer(control),as.double(lower),
         as.double(upper),as.integer(n),PACKAGE="mgcv")
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
# vectorized function for calculating null space dimension for penalties of order m
# for dimension d data M=(m+d+1)!/(d!(m-d)!). Any m not satisfying 2m>d is reset so 
# that 2m>d+1 (assuring "visual" smoothness) 
{ ind<-2*m<d+1
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
{ M<-null.space.dimension(d=d,m=m)
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
      for (j in (i+1):n) if (sum(vn==c(v.names[j,],by.names[j]))==(m+1)) return(TRUE) # there is a repeat in this list of terms
    }
    return(FALSE) # no repeats
  }
  # end of local function declaration, start of main code.....
  if (G$m==0) return(G$C) # nothing to do
  n.vars<-length(G$vnames)
  if (length(unique(G$vnames[(G$nsdf+1):n.vars]))==(n.vars-G$nsdf))
  return(G$C) # nothing to do - there are no repeated variable names in the spline part
  # first assemble the array of term labels
  term.names<-array("",c(G$m,max(G$dim)))
  k<-G$nsdf   
  for (i in 1:G$m)   # term[i,j] holds jth variable name for ith term
  { term.names[i,1:G$dim[i]]<-G$vnames[(k+1):(k+G$dim[i])]
    k<-k+G$dim[i]
  }  
  # now work through the term dimensions
  all.ind<-1:G$m # index vector
  var.stack<-array("",0)  # stores names for null space basis vectors used so far - checked to avoid repeats
  stack.size<-0  # size of var.stack
  #by.index<-rep(0,G$m)
  by.names<-rep("NA",G$m)
  #by.index[G$by.exists]<-1:length(G$by.exists) # indexes existing rows of G$by
  if (sum(G$by.exists)) by.names[G$by.exists]<-rownames(G$by) # names of by variables
  for (d in 1:max(G$dim)) # note if d=1 were not first then would have to deal with "cr" terms specially 
  { # assemble all terms of dimension d, need variable names, p.order and parameter vector offsets
    ind<-all.ind[G$dim==d]
    n.d<-length(ind)
    if (n.d>0) # then there is something to do
    { d.vnames<-term.names[ind,]
      dim(d.vnames)<-c(n.d,dim(term.names)[2]) # otherwise vectors lose dimension attribute!
      d.by.names<-by.names[ind]
      d.off<-G$off[ind]
      d.dim<-G$dim[ind]
      d.p.order<-G$p.order[ind] # order of penalty 
      d.k<-G$df[ind] # basis dimension
      # check these terms for repeats within and between terms
      if(gam.repeat.check(v.names=d.vnames,by.names=d.by.names)) 
      stop("The smooth part of this model is not identifiable")
      # work through the terms imposing constraints on redundant null basis  
      # parameters and updating the null basis stack.....
      for (j in 1:n.d) # work through all n.d d-dimensional terms 
      { #cat(d.vnames[j,]," ",d.p.order[j]," ")
        by<-by.names[j]
        bs.label<-null.space.basis.labels(d.vnames[j,],d.p.order[j],by=by) # get unique null basis vector names
        bs.label<-bs.label[bs.label!="1"] # constants already constrained by sum to zero constraints
        #print(bs.label)
        for (i in 1:length(bs.label)) 
        if (length(var.stack[var.stack==bs.label[i]])>0) # already on the stack
        { # locate the parameter for this term 
          k<-1+d.off[j]+d.k[j]-null.space.dimension(d.dim[j],d.p.order[j])+i
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
# G$by.exists is an array of logicals indicating whether each term has a corresponding "by" 
#             variable
# G$by a compact array of by variables - only has rows for the smooths with by variables
# G$knots a compact array of knot locations for each smooth, in the order corresponding to the 
#         row order in G$x. There are G$dim[i] arrays of length G$n.knots[i] for the ith
#         smnooth - all these arrays are packed end to end in 1-d array G$knots - zero length 1 for no knots.
# G$n.knots - array giving number of knots of basis for each smooth term 0's for none
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
    if (is.null(G$knots)) { G$n.knots<-rep(0,G$m);G$knots<-0;}
    k<-1;for (i in 1:G$m)  # apply shift to basis, as well  
    { n<-G$n.knots[i]
      if (n>0) for (j in 1:G$dim[i]) 
      { G$knots[k:(k+n-1)]<-G$knots[k:(k+n-1)]-G$covariate.shift[i];k<-k+n;}
    }
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
        as.integer(G$nsdf),as.integer(G$dim),as.integer(G$s.type),as.integer(G$p.order),
        as.integer(G$by.exists),as.double(G$by),as.double(G$knots),as.integer(G$n.knots),
        PACKAGE="mgcv") # compiled code to set up matrices for model
  G$X<-matrix(o[[1]],G$n,q);
  G$C<-matrix(o[[2]],G$m,q);
  G$S<-array(o[[3]],sum(G$df^2));             #dim=c(G$m,mdf,mdf));
  G$UZ<-array(o[[4]],UZ.length)  #dim=c(m.type1,G$n+maxM,mdf.type1))
  G$xu.length<-array(o[[6]],m.type1)
  if (m.type1>0) Xu.length<-sum(G$xu.length*G$dim[G$s.type==1]) # Xu.length<-sum(G$xu.length[G$s.type==1]*G$dim[G$s.type==1])
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
        as.integer(length(M$off)),as.integer(nar),PACKAGE="mgcv")
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
# M$min.edf - minimum possible estimated degrees of freedom for model - useful for setting limits
#             on overall smoothing parameter. Set to zero or negative to ignore.
# M$target.edf - set to negative to ignore. This should only be used if cautious optimization
#                is to be used in mgcv searching. If this is non-negative then the local 
#                minimum closest to the target edf will be returned (which can be the global
#                optimum). Designed for use with non-convergent gams.
#
# The routine returns M with the following elements added (or reset):
#
# M$p    - the best fit parameter vector given the selected smoothing parameters
# M$sp   - the vector of smoothing parameters (\theta_i/\rho) estimated by the method
# M$sig2 - the estimate of the error variance (if GCV used)
# M$Vp   - the estimated covariance matrix of the parameters set Vp[1,1] <0 to not calculate
# M$edf  - the estimated degrees of freedom for the ith smooth term if Vp calculated
# M$conv - a list of convergence diagnostics
#          g - gradients of gcv/ubre score at termination, h - leading diagonal of Hessian
#          e - eigenvalues of Hessian, iter - iterations taken, init.ok - TRUE if second 
#          autonitialization guess ok (or intial values supplied), step.fail - TRUE
#          if algorithm terminated on step failure rather than convergence. 
#          edf - array of model edf's from final grid search for overall s.p.
#          score - array of gcv/ubre scores corresponding to edf.
# M$gcv.ubre - the gcv/ubre score at the minimum.
#  
  if (is.null(M$sig2)) M$sig2<-0
  C.r<-nrow(M$C)          # number of equality constraints
  if (is.null(C.r)) C.r<-0
  q<-ncol(M$X)            # number of parameters
  
  n<-nrow(M$X)            # number of data
  # need to strip out penalties for fixed df smooths..... 
  # k<-length(M$df)
  k<-M$m  # needed to allow models with no smooths!
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
  ddiag<-array(0,3*m)   # array for diagonostics
  idiag<-array(0,3)     # array for diagnostics
  Vp[1,1]<-1.0
  M$gcv.ubre<-1.0;
  direct.mesh<-100      # number of points for overall s.p. initial direct search
  sdiag<-array(0.0,2*direct.mesh) # array for gcv/ubre vs edf diagnostics
  if (is.null(M$target.edf)) M$target.edf<- -1 # set to signal no target edf

  oo<-.C("mgcv",as.double(M$y),as.double(M$X),as.double(M$C),as.double(M$w),as.double(S),
         as.double(p),as.double(M$sp),as.integer(off),as.integer(df),as.integer(m),
         as.integer(n),as.integer(q),as.integer(C.r),as.double(M$sig2),as.double(Vp),
		 as.double(edf),as.double(M$conv.tol),as.integer(M$max.half),as.double(ddiag),
                 as.integer(idiag),as.double(sdiag),as.integer(direct.mesh),as.double(M$min.edf),
                 as.double(M$gcv.ubre),as.double(M$target.edf),PACKAGE="mgcv")
   
  p<-matrix(oo[[6]],q,1);
  sig2<-oo[[14]]
  Vp<-matrix(oo[[15]],q,q)
  sp<-matrix(oo[[7]])
  edf<-oo[[16]]
  ddiag<-oo[[19]]
  idiag<-oo[[20]]
  sdiag<-oo[[21]]
  M$gcv.ubre<-oo[[24]]
  conv<-list(edf=sdiag[1:direct.mesh],score=sdiag[direct.mesh+1:direct.mesh],g=ddiag[1:m],h=ddiag[(m+1):(2*m)],
             e=ddiag[(2*m+1):(3*m)],iter=idiag[1],init.ok=as.logical(idiag[2]),step.fail=as.logical(idiag[3]))
  # unpack results back to correct place in output (which includes fixed d.f. and free d.f. terms)
  M$conv<-conv # the convergence diagnostics
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
  G<-list(m=4,n=n,nsdf=1,df=c(15,15,15,15),dim=c(1,1,1,1),s.type=c(0,0,0,0),
          by=0,by.exists=c(FALSE,FALSE,FALSE,FALSE),p.order=c(0,0,0,0),x=x,n.knots=rep(0,4)) # input list for GAMsetup
  H<-GAMsetup(G)  
  H$fix<-array(FALSE,H$m)
  H$y<-y;H$sig2<-sig2;H$w<-w # add data, variance and weights to mgcv input list
  H$sp<-array(-1,H$m)
  H$conv.tol<-1e-6;H$max.half<-15
  H$min.edf<-5
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
  term<-deparse(vars[[d]]) # last term in the ... arguments
  by.var<-deparse(substitute(by)) #getting the name of the by variable
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
    { if (i>ns) stop("Syntax error in s() term: missing f")
      c<-" ";while(c==" "&&i<=ns) { c<-substring(term,i,i);i<-i+1}
      if (c=="f"||c=="F") fx<-TRUE 
      else warning("Strange argument in s() call ignored")
    }
    d<-d-1 # reduce m since last ... argument specified basis dimension not variable
    if (d==1) bs<-"cr" # maintain strict backwards compatibility of code
  } else options(old.warn) # turn on regular warnings again
  term<-deparse(vars[[1]]) # first covariate
  if (d>1) # then deal with further covariates
  for (i in 2:d)
  { term[i]<-deparse(vars[[i]])
  }
  # term now contains the names of the covariates for this model term
  # now evaluate all the other 
  if (k==-1) k<-10*3^(d-1) # auto-initialize basis dimension
  if (bs=="cr") # set basis types
  { bs.type<-0
    if (d>1) { warning("cr basis only works with 1-d smooths!");bs.type<-1;}
  } else 
  { if (bs!="tp") warning("unrecognised basis specifier - set to t.p.r.s. basis.");
    bs.type<-1
  }
  if (m<0) m<-0
  # check for repeated variables in function argument list
  if (length(unique(term))!=d) stop("Repeated variables as arguments of a smooth are not permitted")
  # assemble version of call with all options expanded as text
  full.call<-paste("s(",term[1],sep="")
  if (d>1) for (i in 2:d) full.call<-paste(full.call,",",term[i],sep="")
  full.call<-paste(full.call,",k=",deparse(k),",fx=",deparse(fx),",bs=",deparse(bs),",m=",deparse(m),
                   ",by=",by.var,")",sep="")
  ret<-list(term=term,bs.dim=k,bs.type=bs.type,fixed=fx,dim=d,p.order=m,by=by.var,full.call=full.call)
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
# arguments of s() in GAM formula - changed 4/5/02 environment now taken from gf.
{ p.env<-environment(gf) # environment of formula
  tf<-terms.formula(gf,specials=c("s")) # specials attribute indicates which terms are smooth
  terms<-attr(tf,"term.labels") # labels of the model terms
  nt<-length(terms) # how many terms?
  pf<-"~";          # start of parametric model formula
  if (attr(tf,"response")>0)  # start the replacement formula
  rf<-paste(as.character(attr(tf,"variables")[2]),"~",sep="")
  else rf<-"~"
  sp<-attr(tf,"specials")$s     # array of indeces of smooth terms 
  off<-attr(tf,"offset") # location of offset in formula
  if (!is.null(off)) sp[sp>off]<-sp[sp>off]-1 # have to remove the offset from this index list 
  ns<-length(sp) # number of smooths
  ks<-1;kp<-1 # counters for terms in the 2 formulae
  df<- -1     # max d.f. array 
  fix<-FALSE  # fixed d.f. array 
  bs.type<-1  # basis type array
  dim<-0      # dimension of smoother array
  p.order<-0  # order of the penalties
  v.names<-as.character(attr(tf,"variables")[2])  # names of covariates for smooths starting with response
  n.cov<-1     # total number of covariates for smooths
  by.names<-"" # array of names of "by" variables for each smooth
  if (nt)
  for (i in 1:nt) # work through all terms
  { if (ks<=ns&&sp[ks]==i+1) # it's a smooth
    { #stxt<-paste(substring(terms[i],1,nchar(terms[i])-1),",parent.level=",deparse(parent.level+1),")",sep="")
      st<-eval(parse(text=terms[i]),envir=p.env)
      #sys.frame(sys.parent(n=parent.level))) # get smooth term information
      if (ks>1||kp>1) rf<-paste(rf,"+",st$full.call,sep="") # add to smooth formula
      else rf<-paste(rf,st$full.call,sep="")
      for (i in 1:st$dim) v.names[n.cov+i]<-st$term[i]
      by.names[ks]<-st$by
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
   
  if (!is.null(off)) # deal with offset
  { if (kp>1) pf<-paste(pf,"+",sep="")
    if (kp>1||ks>1) rf<-paste(rf,"+",sep="")
    pf<-paste(pf,as.character(attr(tf,"variables")[1+off]),sep="")
    rf<-paste(rf,as.character(attr(tf,"variables")[1+off]),sep="")
    kp<-kp+1          
  }
  if (attr(tf,"intercept")==0) {pf<-paste(pf,"-1",sep="");rf<-paste(rf,"-1",sep="");if (kp>1) pfok<-1 else pfok<-0}
  else { pfok<-1;if (kp==1) { pf<-"~1"; if (ks==1) rf<-paste(rf,"1",sep="");}}
  sfok<-0;if (ks>1) sfok<-1;
    ret<-list(pftext=pf,pf=as.formula(pf),pfok=pfok,v.names=v.names,by.names=by.names,fix=fix,df=df,
            bs.type=bs.type,s.dim=dim,p.order=p.order,full.formula=as.formula(rf))
  ret
}



gam.setup<-function(formula,data=list(),gam.call=NULL,predict=TRUE,parent.level=1,nsdf=-1,knots=NULL)

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
# 4/5/02 default environment now taken from environment of formula - no need for parent.level
#        this will match behaviour of glm()/lm(), but means that if you create a formula 
#        in a different environment from the one you call it in then things can go wrong
{ # now split the formula
  split<-gam.parser(formula,parent.level=parent.level+1) 
  dmiss<-missing(data)  
  p.env<-environment(formula)
  if (dmiss) data<-p.env
  #sys.frame(sys.parent(n=parent.level))
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
      { environment(split$pf)<-p.env
        #sys.frame(sys.parent(n = parent.level)) # associate appropriate environment with split formula 
        mf<-model.frame(split$pf,data,drop.unused.levels=TRUE) # evaluates in data picking up rest from environment associated with split$pf 
        G$offset <- model.offset(mf)   # get the model offset (if any)
        if (!is.null(G$offset) && length(attr(terms(split$pf),"variables"))<=2) # offest exists, but no more terms except "+1" or "-1"
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
  if (!predict) # obtain the "y variable"  
  { # evaluates rhs of model using data and failing that the calling environment...
    G$y<-eval(parse(text=split$v.names[1]),data,p.env)
    #sys.frame(sys.parent(n = parent.level)))
    # checks that this has worked.....
    if (is.null(G$y)) stop(paste("Failed to find variable for",split$v.names[1]))
  }
  if (m>0) 
  { # now find number of data without looking at y (needed for predicting)
    v.names<-split$v.names          #all.vars(split$sf)  # getting the list of variable names for the smooths
    # need to know number of rows in design matrix .....
    z<-eval(parse(text=v.names[2]),data,p.env)
    #sys.frame(sys.parent(n = parent.level))) 
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
  # count the "by" variables, and fill in logical structure indicating presence
  G$by.exists<-array(FALSE,m);n.by<-0
  if (m>0) for (jj in 1:m) if (split$by.names[jj]!="NA") 
  { G$by.exists[jj]<-TRUE;n.by<-n.by+1}
  G$by<-array(0,c(n.by,G$n)) # array for "by" variables
  n.by<-0 # use to count through by's
  kk<-k<-G$nsdf+1
  v.off<-2  # must be 2 to avoid response!!
  if (m>0)
  for (jj in 1:m)  # loop through smooths
  { if (G$by.exists[jj]) # first deal with by-variables 
    { z<-eval(parse(text=split$by.names[jj]),data,p.env)
      #sys.frame(sys.parent(n = parent.level)))
      if (is.null(z)) stop(paste("Failed to find by-variable",split$by.names[jj]))  
      if (length(z)!=G$n) stop("variable lengths don't all match!")
      n.by<-n.by+1;G$by[n.by,]<-z;
    } 
    # now deal with variables that smooth depends on
    vlist<-v.names[v.off:(v.off+G$dim[jj]-1)];v.off<-v.off+G$dim[jj];
    if (jj==1) 
    {  if (G$nsdf>0) G$vnames<-c(colnames(X),vlist) else
       G$vnames<-vlist 
    } else G$vnames<-c(G$vnames,vlist)
    G$names[kk]<-"s("
    for (i in 1:G$dim[jj]) # produces a column for each variable in this smooth
    { z<-eval(parse(text=vlist[i]),data,p.env)
      #sys.frame(sys.parent(n = parent.level)))
      if (is.null(z)) stop(paste("Failed to find variable",vlist[i]))  
      if (length(z)!=G$n) stop("variable lengths don't all match")
      G$x[k,]<-z
      if (i<length(vlist)) G$names[kk]<-paste(G$names[kk],vlist[i],",",sep="")
      else G$names[kk]<-paste(G$names[kk],vlist[i],sep="")
      k<-k+1
    }
    G$s.type[jj]<-split$bs.type[jj] 
    G$names[kk]<-paste(G$names[kk],")",sep="")
    kk<-kk+1
  }
  if (n.by>0) rownames(G$by)<-split$by.names[split$by.names!="NA"]
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
  # evaluate weights......
  if (is.null(gam.call$weights)) G$w<-rep(1,G$n)
  else  G$w<-eval(gam.call$weights,data)
  
  # now create information on knot locations, if it's been provided....
  kk<-G$nsdf;k<-1;G$knots<-0;G$n.knots<-array(0,G$m)
  if (!is.null(knots)) # then user has supplied a knot location list
  { for (i in 1:G$m)
    { for (j in 1:G$dim[i])
      { kk<-kk+1
        knot.seq<-eval(parse(text=paste("knots$",G$vnames[kk],sep="")))
        if (is.null(knot.seq)) stop(paste("knot sequence variable",G$vnames[kk],"not found"))
        if (G$s.type[i]==0)  knot.seq<-sort(unique(knot.seq))
        n<-length(knot.seq)
        if (j>1&&n!=G$n.knots[i]) stop("Knot location vectors of unequal length within a term")
        else G$n.knots[i]<-n
        G$knots[k:(k+n-1)]<-knot.seq
        k<-k+n
        if (G$s.type[i]==0&&n!=G$df[i]) stop("With a \"cr\" basis you must supply the same number of knots as the basis dimension!")
        if (n<G$df[i]) stop("You can't have fewer basis knots than the basis dimension!")
      }
    }
  }
  G
}


gam<-function(formula,family=gaussian(),data=list(),weights=NULL,control=gam.control(),scale=0,knots=NULL)

# Routine to fit a GAM to some data. The model is stated in the formula, which is then 
# parsed to figure out which bits relate to smooth terms and which to parametric terms.
# -ve binomial stuff by Mike Lonergan, rest by Simon Wood.

{ gam.call<-match.call()  # store the call to facilitate searching in gam.setup()

  family<-get.family(family) # deals with -ve binomial handling as well as any coercion

  if (missing(data)) G<-gam.setup(formula,gam.call=gam.call,parent.level=2,predict=FALSE,knots=knots)
  else G<-gam.setup(formula,data,gam.call=gam.call,parent.level=2,predict=FALSE,knots=knots)
  # ... note weights evaluated using gam.call in gam.setup
  G<-GAMsetup(G) 
  
  G$C<-gam.side.conditions(G) # check for identifiability of smooth part, constrain if poss.
 
  if (is.null(G$offset)) G$offset<-rep(0,G$n)
  
  if (scale==0) 
  { if (family$family=="binomial"||family$family=="poisson"
         || substr(family$family,1,17)=="Negative Binomial") scale<-1 #ubre
    else scale <- -1 #gcv
  }
  
  G$sig2<-scale
  G$sp<-array(-1,G$m) # set up smoothing parameters for autoinitialization at first run
  G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
  G$max.half<-control$mgcv.half # max step halving in Newton update mgcv
  G$min.edf<-G$nsdf+sum(null.space.dimension(G$dim,G$p.order))-dim(G$C)[1]

  if(family$family=="Negative Binomial(NA)") object<-gam.nbut(G, family$link, control, scale)
  else object<-gam.fit(G,family=family,control=control)

  object$covariate.shift<-G$covariate.shift # basis stabilizing linear translations  
  names(object$covariate.shift)<-G$vnames[(dim(G$x)[1]-length(G$covariate.shift)+1):dim(G$x)[1]]
  if (scale<0) object$gcv.used<-TRUE else object$gcv.used<-FALSE
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
  object$full.formula<-as.formula(G$full.formula)
  environment(object$full.formula)<-environment(formula) 
  object$formula<-formula
  object$x<-G$x
  object$by<-G$by;object$by.exists<-G$by.exists
  object$s.type<-G$s.type
  object$p.order<-G$p.order
  object$dim<-G$dim
  object$call<-gam.call
  object$min.edf<-G$min.edf
  rownames(object$x)<-G$vnames
  names(object$x)<-NULL
  class(object)<-"gam"
  object
}

gam.check<-function(b)
# takes a fitted gam object and produces some standard diagnostic plots
{ old.par<-par(mfrow=c(2,2))
  if (b$gcv.used) sc.name<-"GCV" else sc.name<-"UBRE"
  plot(b$mgcv.conv$edf,b$mgcv.conv$score,xlab="Estimated Degrees of Freedom",ylab=paste(sc.name,"Score"),main=paste(sc.name,"w.r.t. model EDF"),type="l")
  points(b$nsdf+sum(b$edf),b$gcv.ubre,col=2,pch=20)
  plot(fitted(b),residuals(b),main="Residuals vs. Fitted",xlab="Fitted Values",ylab="Residuals");
  hist(residuals(b),xlab="Residuals",main="Histogram of residuals");
  plot(fitted(b),b$y,xlab="Fitted Values",ylab="Response",main="Response vs. Fitted Values")
  cat("\nSmoothing parameter selection converged after",b$mgcv.conv$iter,"iteration")
  if (b$mgcv.conv$iter>1) cat("s")
  if (b$mgcv.conv$step.fail) cat(" by steepest\ndescent step failure.\n") else cat(".\n")
  if (length(b$df)>1)
  { cat("The mean absolute",sc.name,"score gradient at convergence was ",mean(abs(b$mgcv.conv$g)),".\n")
    if (sum(b$mgcv.conv$e<0)) cat("The Hessian of the",sc.name ,"score at convergence was not positive definite.\n")
    else cat("The Hessian of the",sc.name,"score at convergence was positive definite.\n")
  }
  if (!b$mgcv.conv$init.ok) cat("Note: the default second smoothing parameter guess failed.\n")
  cat("\n")
  par(old.par)
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
  if (is.null(x$gcv.used)) # then it's an old object
  { gcv<-x$df.null*x$sig2/(x$df.null-sum(x$edf)-x$nsdf)
    cat("\nGCV score: ",gcv," (old gam object: < 0.8)\n")
  } else
  { if (x$gcv.used)  
    cat("\nGCV score: ",x$gcv.ubre,"\n")
    else
    cat("\nUBRE score: ",x$gcv.ubre,"\n")
  }
}

gam.control<-function (epsilon = 1e-04, maxit = 20,globit = 20,mgcv.tol=1e-6,mgcv.half=15, trace = FALSE) 
# control structure for a gam
{   if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(globit) || globit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit,globit = globit, trace = trace, mgcv.tol=mgcv.tol,mgcv.half=mgcv.half)
}

gam.fit<-function (G, start = NULL, etastart = NULL, 
    mustart = NULL, family = gaussian(), 
    control = gam.control(),nb.iter=NULL) 
# fitting function for a gam, modified from glm.fit.
# note that smoothing parameter estimates from one irls iterate are carried over to the next irls iterate
# unless the range of s.p.s is large enough that numerical problems might be encountered (want to avoid 
# completely flat parts of gcv/ubre score). In the latter case autoinitialization is requested.

{
   
    conv <- FALSE
    nobs <- NROW(G$y)
    nvars <- NCOL(G$X) # check this needed
    y<-G$y # original data
    X<-G$X # original design matrix
    if (nvars == 0) stop("Model seems to contain no terms")
    if (family$family=="gaussian" && family$link=="identity") olm<-TRUE # i.e. only one iteration needed
    else olm<-FALSE 

    # obtain average element sizes for the penalties
    if (G$m>0)
    { k1<-0;
      S.size<-array(0,G$m)
      for (i in 1:G$m)
      { k0<-k1
        k1<-k0+G$df[i]^2
        k0<-k0+1
        S.size[i]<-mean(abs(G$S[k0:k1]))
      }
    }  

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
  #  if (is.null(mustart)) 
  #      eval(family$initialize, sys.frame(sys.nframe()))
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

    #eta <- if (!is.null(etastart) && valideta(etastart)) 
    #    etastart
    else if (!is.null(start)) 
    if (length(start) != nvars) 
    stop(paste("Length of start should equal", nvars,
        "and correspond to initial coefs.")) # 1.5.0
   # stop(paste("Length of start should equal", nvars ))
        else 
    #  as.vector(if (NCOL(G$X) == 1) 
    #        G$X * start
    #    else G$X %*% start)
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
   # coefold <- start
    boundary <- FALSE
    scale<-G$sig2
    
    for (iter in 1:(control$maxit+control$globit)) {
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
		if (dim(X)[2]==1) dim(G$X)<-c(length(X[good,]),1) # otherwise dim(G$X)==NULL !!
        ngoodobs <- as.integer(nobs - sum(!good))
        ncols <- as.integer(1)
        # must set G$sig2 to scale parameter or -1 here....
        G$sig2<-scale
        if (G$m>0&&sum(!G$fix>0)) # check that smoothing parameters haven't drifted too far apart
        { temp.sp<-G$sp[!G$fix];temp.S.size<-S.size[!G$fix]*temp.sp
          # check if there is a danger of getting stuck on a flat section of gcv/ubre score...
          if (min(temp.sp)>0 && min(temp.S.size)<Machine()$double.eps^0.5*max(temp.S.size)) 
          G$sp[!G$fix]<- -1.0 # .... if so use use auto-initialization in mgcv
          if (control$trace) cat("Re-initializing smoothing parameters\n")
        } 
        if (iter>control$globit) # solution could be cycling - use more cautious optimization approach
        { G$target.edf<-G$nsdf+sum(G$edf)
        } else
        G$target.edf<- -1 # want less cautious optimization - better at local minimum avoidance
        G<-mgcv(G) 
        if (control$trace)
        { cat("sp: ",G$sp,"\n")
          plot(G$conv$edf,G$conv$score,xlab="EDF",ylab="GCV/UBRE score",type="l");
          points(G$nsdf+sum(G$edf),G$gcv.ubre,pch=20,col=2)
        }
        if (any(!is.finite(G$p))) {
            conv <- FALSE   
            warning(paste("Non-finite coefficients at iteration",
                iter))
            break
        }

		
        start <- G$p
        eta <- drop(X %*% start) # 1.5.0
        #eta[good] <- drop(X[good, , drop = FALSE] %*% start)
        mu <- linkinv(eta <- eta + offset)
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
                #eta[good] <- drop(X[good, , drop = FALSE] %*%  start)
                eta<-drop(X %*% start)
                mu <- linkinv(eta <- eta + offset)
                dev <- sum(dev.resids(y, mu, weights))
            }
            boundary <- TRUE
          #  coef <- start
            if (control$trace) 
                cat("New Deviance =", dev, "\n")
        }
        if (!(valideta(eta) && validmu(mu))) {
            warning("Step size truncated: out of bounds.",call.=FALSE)
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit) 
                  stop("inner loop 2; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                #eta[good] <- drop(X[good, , drop = FALSE] %*% start)
                eta<-drop(X %*% start)
                mu <- linkinv(eta <- eta + offset)
            }
            boundary <- TRUE
            #coef <- start
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("New Deviance =", dev, "\n")
        }
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon || olm) {
            conv <- TRUE
            coef<-start #1.5.0
            break
        }
        else {
            devold <- dev
            coefold <- coef<-start
        }
    }
    if (!conv) 
    { if (is.null(nb.iter)) warning("Algorithm did not converge") 
      else warning("gam.fit didn't converge at nb iteration ",nb.iter)
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
        boundary = boundary,sp = G$sp,df=G$df,nsdf=G$nsdf,Vp=G$Vp,mgcv.conv=G$conv,gcv.ubre=G$gcv.ubre)
}


predict.gam<-function(object,newdata,type="link",se.fit=FALSE,...) {

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
    nc<-0;for (i in 1:m) nc<-nc+object$dim[i]
    for (i in 1:nc) x[i+object$nsdf,] <- x[i+object$nsdf,] + object$covariate.shift[i]
    G<-list(x=x,nsdf=object$nsdf,m=m,n=n,dim=object$dim,by=object$by)
    no.data<-TRUE # used to signal that no data provided, later
  }
  else 
  { G<-gam.setup(object$full.formula,newdata,gam.call=object$call,parent.level=2,predict=TRUE)    
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
            as.integer(control),as.integer(object$by.exists),as.double(G$by),PACKAGE="mgcv")
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
      if (no.data) # reconstrunct original offset
      { link<-object$family$linkfun
        fv<-fitted(object) # original fitted predictor
        offset<-link(fv)-eta # offset is original l.p. minus l.p. calculated without offset 
      } else
      offset<-G$offset
      if (!is.null(offset)) eta<-eta+offset # add offset to l.p.
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
   
  H # ... and return
}

plot.gam<-function(x,rug=TRUE,se=TRUE,pages=0,select=NULL,scale=-1,n=100,n2=40,pers=FALSE,theta=30,phi=30,jit=FALSE,...)

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
# Calls RGAMpredict() directly rather than go via predict.gam() - this is more
# efficient and avoids needing multiple prediction calls to deal with repeated 
# covariates. 
{ sp.contour<-function (x,y,z,zse,xlab="",ylab="",zlab="main",se.plot=TRUE,se.mult=1)   
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

  if (x$dim==0) stop("Model has no smooth terms - nothing for plot.gam() to do.")
  if (se)
  { if (is.numeric(se)) se2.mult<-se1.mult<-se else { se1.mult<-2;se2.mult<-1} 
    if (se1.mult<0) se1.mult<-0;if (se2.mult<0) se2.mult<-0
  } else se1.mult<-se2.mult<-1
  m<-length(x$df) # number of smooth terms
  if (se && x$Vp[1,1]<=0) 
  { se<-FALSE
    warning("No variance estimates available")
  }
  # plot should ignore all "by" variables
  by.exists<-rep(FALSE,m);by<-0
  
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
  
  # now array for 1-d terms to send to predict.gam() to return smooth terms
  
  xx1<-array(0,c(dim(x$x)[1],n))

  j<-1;md2<-FALSE
  for (i in 1:m)
  { if (x$dim[i]==1)
    { x0<-min(x$x[j+x$nsdf,])
      x1<-max(x$x[j+x$nsdf,])
      dx<-(x1-x0)/(n-1) 
      xx1[j+x$nsdf,]<-seq(x0,x1,dx)
      j<-j+1
    } else
    { for (k in 1:x$dim[i])
      { xx1[j+x$nsdf,]<-rep(0,n);j<-j+1 # multi-dim terms padded with zeroes
        if (x$dim[i]==2) md2<-TRUE
      }  
    }  
  }
  rownames(xx1)<-rownames(x$x)

  X<-0;
  if (se) control<-3 else control<-2

  eta<-array(0,dim=c(x$nsdf+m,n))
  se.fit<-array(0,dim=c(x$nsdf+m,n))

  o<-.C("RGAMpredict", as.integer(x$xu.length),as.double(x$Xu),as.double(x$UZ),
             as.double(x$xp),as.integer(x$nsdf),as.integer(x$dim),as.integer(x$s.type),
            as.integer(x$df),as.integer(x$p.order), as.integer(m),as.integer(length(x$y)),as.double(xx1),
            as.integer(n),as.double(x$coefficients),as.double(x$Vp),as.double(eta),as.double(se.fit),as.double(X),
            as.integer(control),as.integer(by.exists),as.double(by),PACKAGE="mgcv")
  # now eta contains linear predictor (terms) and se may contain corresponding standard errors
  eta<-array(o[[16]],c(x$nsdf+m,n));
  se.fit<-array(o[[17]],c(x$nsdf+m,n))*se1.mult;
  if (se) pl1<-list(fit=eta,se.fit=se.fit) else pl1<-eta
  rm(o)
  # shift covariates back
  j<-1;
  for (i in 1:m)
  { if (x$dim[i]==1)
    { xx1[j+x$nsdf,]<-xx1[j+x$nsdf,]+x$covariate.shift[j]
      j<-j+1
    } else j<-j+x$dim[i]
  }
 
  if (md2) # then create data frames for 2-d plots
  { if (n2<10) n2<-10
    xx2<-array(0,c(dim(x$x)[1],(n2*n2)))   # new x array for prediction
    xm<-data.frame(1:n2);ym<-data.frame(1:n2) # grid points for plotting
    if (x$nsdf>0) 
    for (i in 1:x$nsdf) xx2[i,]<-rep(0,n2*n2)
    j<-1
    for (i in 1:m)
    { if (x$dim[i]==2)
      { x0<-min(x$x[j+x$nsdf,]);x1<-max(x$x[j+x$nsdf,])
        dx<-(x1-x0)/(n2-1);
        y0<-min(x$x[j+x$nsdf+1,]);y1<-max(x$x[j+x$nsdf+1,])
        dy<-(y1-y0)/(n2-1)
        
        xm[i]<-seq(x0,x1,dx) 
        ym[i]<-seq(y0,y1,dy)
        xx2[j+x$nsdf,]<-rep(xm[[i]],n2)
        xx2[j+x$nsdf+1,]<-rep(ym[[i]],rep(n2,n2))
        xm[i]<-xm[i]+x$covariate.shift[j]
        ym[i]<-ym[i]+x$covariate.shift[j+1]
        j<-j+2
      } else
      {  for (k in 1:x$dim[i])
         { xx2[j+x$nsdf,]<-rep(0,n2*n2);j<-j+1
           xm[i]<-rep(0,n2);ym[i]<-xm[i];
         }

      }
    } 
    rownames(xx2)<-rownames(x$x)
    X<-0;
    if (se) control<-3 else control<-2

    eta<-array(0,dim=c(x$nsdf+m,n2*n2))
    se.fit<-array(0,dim=c(x$nsdf+m,n2*n2))

    o<-.C("RGAMpredict", as.integer(x$xu.length),as.double(x$Xu),as.double(x$UZ),
             as.double(x$xp),as.integer(x$nsdf),as.integer(x$dim),as.integer(x$s.type),
            as.integer(x$df),as.integer(x$p.order), as.integer(m),as.integer(length(x$y)),as.double(xx2),
            as.integer(n2*n2),as.double(x$coefficients),as.double(x$Vp),as.double(eta),as.double(se.fit),as.double(X),
            as.integer(control),as.integer(by.exists),as.double(by),PACKAGE="mgcv")
    # now eta contains linear predictor (terms) and se may contain corresponding standard errors
    eta<-array(o[[16]],c(x$nsdf+m,n2*n2));
    se.fit<-array(o[[17]],c(x$nsdf+m,n2*n2))*se2.mult;
    if (se) pl2<-list(fit=eta,se.fit=se.fit) else pl2<-eta
    rm(o)
    # shift covariates back
    j<-1;
    for (i in 1:m)
    { if (x$dim[i]==2)
      { xx1[j+x$nsdf,]<-xx1[j+x$nsdf,]+x$covariate.shift[j]
        j<-j+1
        xx1[j+x$nsdf,]<-xx1[j+x$nsdf,]+x$covariate.shift[j]
        j<-j+1
      } else j<-j+x$dim[i]
    }

  }
  if (se)   # pl$fit and pl$se.fit
  { k<-0
    if (scale==-1) # getting common scale for 1-d terms
    for (i in 1:m)
    { if (x$dim[i]==1)
      { ul<-pl1$fit[x$nsdf+i,]+pl1$se.fit[x$nsdf+i,]
        ll<-pl1$fit[x$nsdf+i,]-pl1$se.fit[x$nsdf+i,]
        if (k==0) { ylim<-c(min(ll),max(ul));k<-1;}
        else
        { if (min(ll)<ylim[1]) ylim[1]<-min(ll)
	  if (max(ul)>ylim[2]) ylim[2]<-max(ul)
        }
      }
    }
    j<-1
    for (i in 1:m)
    if (is.null(select)||i==select)
	{ if (interactive() && x$dim[i]<3 && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
      if (x$dim[i]==1)
      { ul<-pl1$fit[x$nsdf+i,]+pl1$se.fit[x$nsdf+i,]
        ll<-pl1$fit[x$nsdf+i,]-pl1$se.fit[x$nsdf+i,]
        if (scale==0) { ylim<-c(min(ll),max(ul))}
        title<-paste("s(",rownames(x$x)[j+x$nsdf],",",as.character(round(x$edf[i],2)),")",sep="")
        plot(xx1[x$nsdf+j,],pl1$fit[x$nsdf+i,],type="l",xlab=rownames(x$x)[j+x$nsdf],ylim=ylim,ylab=title)
	lines(xx1[x$nsdf+j,],ul,lty=2)
        lines(xx1[x$nsdf+j,],ll,lty=2)
	if (rug) 
        { if (jit) rug(jitter(as.numeric(x$x[x$nsdf+j,]+x$covariate.shift[j])))
           else rug(as.numeric(x$x[x$nsdf+j,]+x$covariate.shift[j]))
	}
      } else if (x$dim[i]==2)
      { xla<-rownames(x$x)[j+x$nsdf];yla<-rownames(x$x)[j+x$nsdf+1]
        title<-paste("s(",xla,",",yla,",",as.character(round(x$edf[i],2)),")",sep="")
        if (pers) 
        { persp(xm[[i]],ym[[i]],matrix(pl2$fit[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,zlab=title,theta=theta,phi=phi)
        } else
        { sp.contour(xm[[i]],ym[[i]],matrix(pl2$fit[x$nsdf+i,],n2,n2),matrix(pl2$se.fit[x$nsdf+i,],n2,n2),
                     xlab=xla,ylab=yla,zlab=title,se.mult=se2.mult)
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
	if (is.null(select)||i==select)
    { if (interactive() && x$dim[i]<3 && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
      if (x$dim[i]==1)
      { title<-paste("s(",rownames(x$x)[j+x$nsdf],",",as.character(round(x$edf[i],2)),")",sep="")
        if (scale==0) ylim<-range(pl1[x$nsdf+i,])
        plot(xx1[x$nsdf+j,],pl1[x$nsdf+i,],type="l",,xlab=rownames(x$x)[j+x$nsdf],ylab=title,ylim=ylim)
        if (rug) 
		{ if (jit) rug(jitter(as.numeric(x$x[x$nsdf+j,]+x$covariate.shift[j])))
          else rug(as.numeric(x$x[x$nsdf+j,]+x$covariate.shift[j]))
		}
      } else if (x$dim[i]==2)
      { xla<-rownames(x$x)[j+x$nsdf];yla<-rownames(x$x)[j+x$nsdf+1]
        title<-paste("s(",xla,",",yla,",",as.character(round(x$edf[i],2)),")",sep="")
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

summary.gam<-function (object,...) 
# summary method for gam object - provides approximate p values for terms + other diagnostics
{ pinv<-function(V,M)
  { D<-svd(V)
    M1<-length(D$d[D$d>1e-12*D$d[1]]);if (M>M1) M<-M1 # avoid problems with zero eigen-values
    D$d[(M+1):length(D$d)]<-1
    D$d<- 1/D$d
    D$d[(M+1):length(D$d)]<-0
    D$u%*%diag(D$d)%*%t(D$v)
  } 
  se<-0;for (i in 1:length(object$coefficients)) se[i]<-object$Vp[i,i]^0.5
  residual.df<-length(object$y)-object$nsdf-sum(object$edf)
  if (object$nsdf>0)
  { p.coeff<-object$coefficients[1:object$nsdf]
    p.t<-p.coeff/se[1:object$nsdf]
    p.pv<-2*pt(abs(p.t),df=round(residual.df),lower.tail=FALSE)
  } 
  else {p.coeff<-p.t<-p.pv<-array(0,0)}
  m<-length(object$edf)
  s.pv<-chi.sq<-array(0,0)
  if (m>0) # form test statistics for each smooth
  { stop<-object$nsdf
    for (i in 1:m)
    { start<-stop+1;stop<-stop+object$df[i]
      V<-object$Vp[start:stop,start:stop] # cov matrix for smooth
      p<-object$coefficients[start:stop]  # params for smooth
      # now get null space dimension for this term
      M1<-object$df[i]-1
      M<-round(object$edf[i])
      V<-pinv(V,M1) # get rank M pseudoinverse of V
      chi.sq[i]<-t(p)%*%V%*%p
      er<-names(object$coefficients)[start]
      er<-substring(er,1,nchar(er)-2)
      names(chi.sq)[i]<-er
      s.pv[i]<-pchisq(chi.sq[i],df=max(1,object$edf[i]),lower.tail=FALSE) 
    }
  }
  r.sq<- 1 - var(object$y-object$fitted.values)*(object$df.null-1)/(var(object$y)*residual.df) 
  ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,edf=object$edf,
       s.pv=s.pv,scale=object$sig2,r.sq=r.sq,family=object$family,formula=object$formula,n=object$df.null)
  if (is.null(object$gcv.used))  # then it's an old object
  {  ret$gcv<-object$df.null*object$sig2/residual.df
  } else # new object
  if (object$gcv.used) ret$gcv<-object$gcv.ubre else ret$ubre<-object$gcv.ubre
  class(ret)<-"summary.gam"
  ret
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
  cat("\nAdjusted r-sq. = ",formatC(x$r.sq,digits=3,width=5))
  if (is.null(x$ubre)) cat("    GCV score = ",formatC(x$gcv,digits=5))
  else cat("   UBRE score = ",formatC(x$ubre,digits=5))
   cat("\nScale estimate = ",formatC(x$scale,digits=5,width=8,flag="-"),"         n = ",x$n,"\n",sep="")
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
    se <- as.numeric(se)
    variables <- rownames(x$x)
    variables <-variables[variables!="(Intercept)"]
    variables <-variables[variables!="constant"]
    which.slice <- ""
    if (length(variables) < 2) 
        plot(x)
    else {
        if (is.null(view)) 
            view <- variables[1:2]
        for (v in variables) {
            if (eval(parse(text = paste("is.null(slice$",v,")", sep = "")))) {
                if (v %in% dimnames(x$covariate.shift)[[1]])
                   eval(parse(text = paste("slice$", v, "<-x$covariate.shift[v]",
                                                    sep = "")))
                else
                   eval(parse(text = paste("slice$", v, "<-mean(x$x[v,])", sep = "")))
            }
        }
        view.data1 <- x$x[view[1], ] 
        if (view[1] %in% dimnames(x$covariate.shift)[[1]])
           view.data1 <- view.data1 + x$covariate.shift[view[1]]
        view.data2 <- x$x[view[2], ]
        if (view[2]%in% dimnames(x$covariate.shift)[[1]])
           view.data2 <- view.data2 + x$covariate.shift[view[2]]
        gridmaker <- paste("expand.grid(", view[1], "=seq(", 
            min(view.data1), ",", max(view.data1), ", length=", 
            sizes[1], "),", view[2], "=seq(", min(view.data2), 
            ",", max(view.data2), ", length=", sizes[2], "))", 
            sep = "")
        persp.dat <- eval(parse(text = gridmaker))
        for (v in variables) {
            if (eval(parse(text = paste("length(persp.dat$", 
                v, ")==0", sep = "")))) {
                this.slice <- eval(parse(text = slice[[v]]), persp.dat)
                if (length(this.slice) == 1) 
                  this.slice <- rep(this.slice, dim(persp.dat)[1])
                eval(parse(text = paste("persp.dat$", v, "<-this.slice", 
                                                 sep = "")))
                sl <- slice[[v]]
                if (is.numeric(sl)) 
                  sl <- signif(sl, 3)
                which.slice <- paste(which.slice, v, "=", sl, 
                  ", ", sep = "")
            }
        }
        if (se > 0) 
            which.slice <- paste(which.slice, " red/green are +/-", 
                se, " se", sep = "")
        else which.slice <- substr(which.slice, 1, nchar(which.slice) - 
            1)
        if (se <= 0) 
            persp.dat$response <- as.vector(predict(x, newdata = persp.dat, 
                type = "response"))
        else {
            response <- predict(x, newdata = persp.dat, type = "response", 
                se.fit = TRUE)
            persp.dat$response <- as.vector(response[[1]])
            persp.dat$se <- as.vector(response[[2]])
        }
        if (mask != FALSE) {
            if (is.numeric(mask)) {
                if (length(mask) < 2) 
                  mask <- c(mask, mask)
            }
            else {
                xd <- sort(view.data1)
                xd <- abs(xd[-1] - xd[-length(xd)])
                mask <- min(pmax(xd[-1] - xd[-length(xd)]))
                xd <- sort(view.data2)
                xd <- abs(xd[-1] - xd[-length(xd)])
                mask <- c(mask, min(pmax(xd[-1] - xd[-length(xd)])))
            }
            for (i in 1:length(persp.dat$response)) {
                if (min(abs(view.data1 - persp.dat[i, 1])/mask[1] + 
                  abs(view.data2 - persp.dat[i, 2])/mask[2]) > 
                  1) 
                  persp.dat$response[i] <- NA
            }
        }
        fit.mat <- matrix(persp.dat$response, sizes[1], sizes[2])
        if (se > 0) {
            fit.mat.u <- matrix((persp.dat$response + se * persp.dat$se), 
                sizes[1], sizes[2])
            fit.mat.l <- matrix((persp.dat$response - se * persp.dat$se), 
                sizes[1], sizes[2])
            zlim <- c(min(persp.dat$response - se * persp.dat$se, 
                na.rm = TRUE), max(persp.dat$response + se * 
                persp.dat$se, na.rm = TRUE))
            persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]), 
                y = seq(min(view.data2), max(view.data2), len = sizes[2]), 
                z = fit.mat, xlab = view[1], ylab = view[2], 
                zlab = "", zlim = zlim, sub = which.slice, theta = theta, 
                phi = phi, r = r, d = d, scale = scale, expand = expand, 
                col = col, border = NA, ltheta = ltheta, lphi = lphi, 
                shade = shade, box = box, axes = axes, nticks = nticks, 
                ticktype = ticktype)
            par(new = TRUE)
            persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]), 
                y = seq(min(view.data2), max(view.data2), len = sizes[2]), 
                z = fit.mat.l, xlab = view[1], ylab = view[2], 
                zlab = "", zlim = zlim, sub = which.slice, theta = theta, 
                phi = phi, r = r, d = d, scale = scale, expand = expand, 
                col = NA, border = "green", shade = NA, box = box, 
                axes = axes, nticks = nticks, ticktype = ticktype)
            par(new = TRUE)
            persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]), 
                y = seq(min(view.data2), max(view.data2), len = sizes[2]), 
                z = fit.mat, xlab = view[1], ylab = view[2], 
                zlab = "", zlim = zlim, sub = which.slice, theta = theta, 
                phi = phi, r = r, d = d, scale = scale, expand = expand, 
                col = NA, border = border, shade = NA, box = box, 
                axes = axes, nticks = nticks, ticktype = ticktype)
            par(new = TRUE)
            persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]), 
                y = seq(min(view.data2), max(view.data2), len = sizes[2]), 
                z = fit.mat.u, xlab = view[1], ylab = view[2], 
                zlab = "", zlim = zlim, sub = which.slice, theta = theta, 
                phi = phi, r = r, d = d, scale = scale, expand = expand, 
                col = NA, border = "red", shade = NA, box = box, 
                axes = axes, nticks = nticks, ticktype = ticktype)
            invisible(list(fit.mat, fit.mat.l, fit.mat.u))
        }
        else {
            persp(x = seq(min(view.data1), max(view.data1), len = sizes[1]), 
                y = seq(min(view.data2), max(view.data2), len = sizes[2]), 
                z = fit.mat, xlab = view[1], ylab = view[2], 
                zlab = "", sub = which.slice, theta = theta, 
                phi = phi, r = r, d = d, scale = scale, expand = expand, 
                col = col, border = border, ltheta = ltheta, 
                lphi = lphi, shade = shade, box = box, axes = axes, 
                nticks = nticks, ticktype = ticktype)
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
      
      family <- do.call("neg.binom", list(theta = th, link = link))
      if (iter==1) etastart<-family$linkfun(object$fitted)
      object <- gam.fit(G, etastart = etastart, 
                           family = family, control = control,nb.iter=iter)
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


neg.binom <-function(theta = NA, link = "log")
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
    cat("This is mgcv 0.8.1\n")
}











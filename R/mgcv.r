# These are the R routines for the package mgcv (c) Simon Wood 2000-2003

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
if (!is.matrix(A)) stop("Argument to QT must be a matrix.")
if (nrow(A)>ncol(A)) stop("Matrix argument to QT has more rows than columns.")
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
        by<-d.by.names[j]
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
# G$sp an array of supplied smoothing parameters
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
# G$fit.method - one of "mgcv" or "magic", which determines the exact form of H$S and H$off
# G$min.sp - minimum values for the smoothing parameters, only used with fit.method=="magic"
# G$H - offset penalty matrix, only used with fit.method=="magic"
# G$fix array of logicals indicating whether or not smooth should have no penalty
# The function returns a list H containing all the above named elements plus the following:
#
# H$X the full design matrix.
# H$S contains the elements of the G$m penalty matrices. let start.k=sum(G$df[1:(k-1)]^2) and start_1=0
#     Element i,j of the kth penalty matrix is S[start.k+i+G$df[k]*(j-1)]. Note however that this
#     matrix is infact the smallest matrix containing all the non-zero elements of the full
#     penalty. Paste it into location starting at M$off[k], M$off[k] of an appropriate
#     matrix of zeroes to get the full penalty. If fit.method=="magic" this is a list of matrices! 
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
# H$rank - the ranks of the penalty matrices
# H$m.free - number of free penalties (magic only)
# H$m.off - offeste for free penalties (magic only)
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
    G$rank<- G$df-null.space.dimension(G$dim,G$p.order) # ranks of penalty matrices
  } else mdf<-G$covariate.shift<-xp<-S<-C<-m.type1<-M<-UZ<-Xu<-xu.length<-off<-UZ.length<-Xu.length<-0;
  o<-.C("RGAMsetup",as.double(X),as.double(C),as.double(S),as.double(UZ),as.double(Xu),as.integer(xu.length),as.double(xp),
        as.integer(off),as.double(G$x),as.integer(G$m),as.integer(G$n),as.integer(G$df),
        as.integer(G$nsdf),as.integer(G$dim),as.integer(G$s.type),as.integer(G$p.order),
        as.integer(G$by.exists),as.double(G$by),as.double(G$knots),as.integer(G$n.knots),
        PACKAGE="mgcv") # compiled code to set up matrices for model
  G$X<-matrix(o[[1]],G$n,q);
  G$C<-matrix(o[[2]],G$m,q);
  G$off<-array(o[[8]],c(G$m));
  if (G$fit.method=="magic") # magic requires a different penalty matrix format to mgcv
  { G$S<-list();
    if (!is.null(G$H)&&(dim(G$H)[1]!=q||dim(G$H)[2]!=q)) stop("dimensions of H incorrect.")
    k<-1;G$m.free<-0;G$m.off<-0
    if (G$m>0)
    for (i in 1:G$m)
    { j<-G$df[i];Si<-matrix(o[[3]][k:(k+j*j-1)],j,j);k<-k+j*j;
      if (G$sp[i]<0&&!G$fix[i]) # then s.p. to be estimated
      { G$m.free<-G$m.free+1;G$S[[G$m.free]]<-Si;
        G$m.off[G$m.free]<-G$off[i]+1
        G$rank[G$m.free]<-G$rank[i]
      } else
      if (G$sp[i]>0&&!G$fix[i]) # then s.p. supplied by user
      { if (is.null(G$H)) G$H<-matrix(0,q,q)
        off1<-G$off[i]+1;off2<-off1+j-1
        G$H[off1:off2,off1:off2]<-G$H[off1:off2,off1:off2]+G$sp[i]*Si
      } else G$sp[i]=0; # s.p. fixed at zero 
    }
    # if minimum smoothing parameters supplied then penalties times these must be added to H
    if (!is.null(G$min.sp)&&G$m>0)
    { if (is.null(G$H)) G$H<-matrix(0,q,q)
      for (i in 1:G$m)
      { off1<-G$off[i]+1;off2<-off1+G$df[i]-1
        G$H[off1:off2,off1:off2]<-G$H[off1:off2,off1:off2]+G$min.sp[i]*G$S[[i]]
      }
    }
  } else # method is "mgcv"
  { G$S<-array(o[[3]],sum(G$df^2));             #dim=c(G$m,mdf,mdf));
    if (!is.null(G$H)) warning("H ignored for fit method mgcv.")
    if (!is.null(G$min.sp)) warning("min.sp ignored for fit method mgcv.")
  }
  G$UZ<-array(o[[4]],UZ.length)  #dim=c(m.type1,G$n+maxM,mdf.type1))
  G$xu.length<-array(o[[6]],m.type1)
  if (m.type1>0) Xu.length<-sum(G$xu.length*G$dim[G$s.type==1]) # Xu.length<-sum(G$xu.length[G$s.type==1]*G$dim[G$s.type==1])
  else Xu.length<-0
  G$Xu<-array(o[[5]],Xu.length)    #dim=c(m.type1,G$n,mdim))
  G$xp<-matrix(o[[7]],G$m,mdf);
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
# M$m  - number of penalties
# M$X  - the design matrix (declared as a matrix, with correct numbers of rows and columns,
#        as number of data and number of parameters are read from this)
# M$C  - matrix defining linear equality constraints on parameters (Cp=0). 
#        Number of rows is number of constraints.
# M$w  - weight vector (often proportional to inverse of standard deviation)
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
# M$fixed.sp - set to 1 to treat smoothing parameters as fixed.
#
# The routine returns M with the following elements added (or reset):
#
# M$p    - the best fit parameter vector given the selected smoothing parameters
# M$sp   - the vector of smoothing parameters (\theta_i/\rho) estimated by the method
# M$sig2 - the estimate of the error variance (if GCV used)
# M$Vp   - the estimated covariance matrix of the parameters set Vp[1,1] <0 to not calculate
# M$edf  - the estimated degrees of freedom for the ith smooth term if Vp calculated
# M$hat  - array of same length as M$y with elements from leading diagonal of hat/influence matrix
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
      M$sp[j]<-M$sp[i] 
      j<-j+1
    }  
  }
  # check smoothing parameters if fixed....
  if (M$fixed.sp)
  { if (sum(M$sp<0)) stop("mgcv called with fixed negative smoothing parameter(s).")
    if (sum(is.na(M$sp))) stop("mgcv called with fixed NA smoothing parameter(s).")
  }

  # deal with quantities that will be estimated
  p<-matrix(0,q,1)      # set up parameter vector
  Vp<-matrix(0.0,q,q)   # estimated covariance matrix
  edf<-array(0,m)       # estimated degrees of freedom
  hat<-array(0,n)       # elements on leading diagonal of hat matrix
  ddiag<-array(0,3*m)   # array for diagonostics
  idiag<-array(0,3)     # array for diagnostics
  Vp[1,1]<-1.0
  M$gcv.ubre<-1.0;
  direct.mesh<-100      # number of points for overall s.p. initial direct search
  sdiag<-array(0.0,2*direct.mesh) # array for gcv/ubre vs edf diagnostics
  if (is.null(M$target.edf)) M$target.edf<- -1 # set to signal no target edf

  oo<-.C("mgcv",as.double(M$y),as.double(M$X),as.double(M$C),as.double(M$w^2),as.double(S),
         as.double(p),as.double(M$sp),as.integer(off),as.integer(df),as.integer(m),
         as.integer(n),as.integer(q),as.integer(C.r),as.double(M$sig2),as.double(Vp),
		 as.double(edf),as.double(M$conv.tol),as.integer(M$max.half),as.double(ddiag),
                 as.integer(idiag),as.double(sdiag),as.integer(direct.mesh),as.double(M$min.edf),
                 as.double(M$gcv.ubre),as.double(M$target.edf),as.integer(M$fixed.sp),
                 as.double(hat),PACKAGE="mgcv")
   
  p<-matrix(oo[[6]],q,1);
  sig2<-oo[[14]]
  Vp<-matrix(oo[[15]],q,q)
  sp<-matrix(oo[[7]])
  edf<-oo[[16]]
  ddiag<-oo[[19]]
  idiag<-oo[[20]]
  sdiag<-oo[[21]]
  M$gcv.ubre<-oo[[24]]
  M$hat<-oo[[27]]
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
    if (c=="|")  # check whether term specified as fixed d.f.
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
  k.new<-round(k) # in case user has supplied non-integer basis dimension
  if (!all.equal(k.new,k)) {warning("argument k of s() should be integer and has been rounded")}
  k<-k.new
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
  by.names<-array("",0) # array of names of "by" variables for each smooth
  if (nt)
  for (i in 1:nt) # work through all terms
  { if (ks<=ns&&sp[ks]==i+1) # it's a smooth
    { st<-eval(parse(text=terms[i]),envir=p.env)
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
  class(ret)<-"split.gam.formula"
  ret
}



gam.setup<-function(formula,data=stop("No data supplied to gam.setup"),predict=TRUE,nsdf=-1,knots=NULL,sp=NULL,
                    min.sp=NULL,H=NULL,fit.method="magic")
# This gets the data referred to in the model formula and sets up 
# G$names[i] contains the names associated with each column of the design matrix,
# followed by the name of each smooth term - the constant is un-named! 
# if nsdf==-1 then it assumed that the design matrix for the non-spline part of the 
# model is to be found here. Otherwise it is assumed that nsdf is the known 
# number of degrees of freedom for the non-spline part (including intercept), but
# the columns corresponding to the non-spline part of the model and the offset 
#  are set to zero.

{ # split the formula if the object being passed is a formula, otherwise it's already split
  if (class(formula)=="formula") split<-gam.parser(formula) 
  else if (class(formula)=="split.gam.formula") split<-formula
  else stop("First argument is no sort of formula!") 
  if (split$df[1]==-1)
  { if (split$pfok==0) stop("You've got no model....")
    m<-0
  }  
  else  m<-length(split$df) # number of smooth terms
  G<-list(m=m,df=split$df,full.formula=split$full.formula,min.sp=min.sp,H=H)
  G$fix<-split$fix
  if (fit.method=="fastest") 
  { if (G$m==1) G$fit.method<-"mgcv" else G$fit.method<-"magic"
  } else G$fit.method<-fit.method
  if (!is.null(sp)) # then user has supplied fixed smoothing parameters
  { ok<-TRUE
    if (length(sp)!=m) { ok<-FALSE;warning("Fixed smoothing parameter vector is too short - ignored.")}
    if (sum(is.na(sp))) { ok<-FALSE;warning("NA's in fixed smoothing parameter vector - ignoring.")}
    if (ok) { G$sp<-sp; G$fixed.sp<-1} else { G$fixed.sp<-0;G$sp<-rep(-1,m)}
  } else # set up for auto-initialization
  { G$fixed.sp<-0;G$sp<-rep(-1,m)}
  if (!is.null(min.sp)) # then minimum s.p.'s supplied
  { if (length(min.sp)!=m) stop("length of min.sp is wrong.")
    if (sum(is.na(min.sp))) stop("NA's in min.sp.")
    if (sum(min.sp<0)) stop("elements of min.sp must be non negative.")
  }
  add.constant<-FALSE
  if (split$pfok)   # deal with the strictly parametric part of the model 
  { if (length(attr(terms(split$pf),"variables"))==1) # then the formula is ~ 1 and model.frame(), etc can't cope 
    { G$nsdf<-1
      add.constant<-TRUE  # must construct model matrix later when length of response is known
    } else
    { if (nsdf>=0)  # set G$nsdf to supplied value, but deal with X later
      { G$nsdf<-nsdf
      } else 
      { if (is.null(attr(data,"terms"))) # then data is not a model frame
        mf<-model.frame(split$pf,data,drop.unused.levels=FALSE) # mudt be false or can end up with wrong prediction matrix!
        else mf<-data # data is already a model frame
        # ... evaluates in data  
        G$offset <- model.offset(mf)   # get the model offset (if any)
        if (!is.null(G$offset) && length(attr(terms(split$pf),"variables"))<=2) # offset exists, but no more terms except "+1" or "-1"
        { if (length(grep("-",as.character(split$pf[2])))>=1) # then there is a "-1" term in formula
          X<-matrix(0,length(G$offset),0)
          else X <- model.matrix(split$pf,mf)  # then there is a constant in the model and model.matrix() can cope 
        }
        else X <- model.matrix(split$pf,mf)       # getting the model matrix
        G$nsdf <- dim(X)[2]                  # extract the total non- smooth df
      }        
    }	
  } else
  { if (m==0) stop("Model appears to have no terms") 
    G$nsdf<-0
  }
  if (!predict) # obtain the "y variable"  
  { # evaluates rhs of model using data and failing that the calling environment...
    G$y<-eval(parse(text=paste("data$\"",split$v.names[1],"\"",sep=""))) # works if data is model frame
    if (is.null(G$y)) G$y<-eval(parse(text=split$v.names[1]),data)       # otherwise evaluate in data
    # checks that this has worked.....
    if (is.null(G$y)) stop(paste("Failed to find variable for",split$v.names[1]))
  }
  if (m>0) v.names<-split$v.names    # getting the list of variable names for the smooths

  G$n<-dim(data)[1]
    
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
  if (G$nsdf>0) G$vnames<-colnames(X) 
  if (m>0)
  for (jj in 1:m)  # loop through smooths
  { if (G$by.exists[jj]) # first deal with by-variables 
    { z<-eval(parse(text=split$by.names[jj]),data)#,p.env)
      if (is.null(z)) stop(paste("Failed to find by-variable",split$by.names[jj]))  
      if (length(z)!=G$n) stop("variable lengths don't all match!")
      n.by<-n.by+1;G$by[n.by,]<-z;
    } 
    # now deal with variables that smooth depends on
    vlist<-v.names[v.off:(v.off+G$dim[jj]-1)];v.off<-v.off+G$dim[jj];
    if (jj==1) 
    {  if (G$nsdf>0) G$vnames<-c(G$vnames,vlist) else
       G$vnames<-vlist 
    } else G$vnames<-c(G$vnames,vlist)
    G$names[kk]<-"s("
    for (i in 1:G$dim[jj]) # produces a column for each variable in this smooth
    { z<-eval(parse(text=paste("data$\"",vlist[i],"\"",sep=""))) # is term in model frame ?
      if (is.null(z)) z<-eval(parse(text=vlist[i]),data)#,p.env) # if not, evaluate in data frame
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
      { warning("Max d.f. for term must be > than penalty null space dimension - max. d.f. has been raised for term");
        G$df[i]<-M+1
      }
      k<-k+G$df[i]
      off<-off+G$dim[i] 
    }
    xx<-G$x[(G$nsdf+1):off,]
    if (off-G$nsdf==1) n<-length(unique(xx))
    else n<-dim(uniquecombs(xx))[2]
    if (k>n) 
    stop("Total max model d.f. must be <= unique covariate combinations - use uniquecombs() to investigate further") 
  }
  # evaluate weights......
 
  if (is.null(data$"(weights)")) G$w<-rep(1,G$n)
  else G$w<-data$"(weights)"  

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
        if (G$s.type[i]==0&&n!=G$df[i]) 
        stop("With a \"cr\" basis you must supply the same number of knots as the basis dimension!")
        if (n<G$df[i]) stop("You can't have fewer basis knots than the basis dimension!")
      }
    }
  }
  rownames(G$x)<-G$vnames
  G
}




gam<-function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,control=gam.control(),
              scale=0,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,...)

# Routine to fit a GAM to some data. The model is stated in the formula, which is then 
# parsed and interpreted to figure out which bits relate to smooth terms and which to parametric terms.

{  # create model frame.....
  gp<-gam.parser(formula) # interpret the formula 
  cl<-match.call() # call needed in gam object for update to work
  mf<-match.call(expand.dots=FALSE)
  ff<-paste(gp$v.names[1],gp$pftext) # fake formula to collect necessary data
  n<-length(gp$v.names) # pick up arguments of smooths
  if (n>1) { ff1<-paste(gp$v.names[2:n],collapse="+");ff<-paste(ff,"+",ff1)}
  if (sum(gp$by.names!="NA")) # ... and "by" variables
  { ff1<-paste(gp$by.names[gp$by.names!="NA"],collapse="+")
    ff<-paste(ff,"+",ff1)
  }
  mf$formula<-as.formula(ff,env=environment(formula))
  mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-mf$gamma<-mf$...<-NULL
  mf$drop.unused.levels<-TRUE
  mf[[1]]<-as.name("model.frame")
  mf <- eval(mf, parent.frame()) # the model frame now contains all the data 

  if (is.character(family)) family<-eval(parse(text=family))
  if (is.function(family)) family <- family()
  if (is.null(family$family)) stop("family not recognized")
  
  if (sum(control$fit.method==c("mgcv","magic","fastest"))==0) stop("Unknown fit method.") 
  G<-gam.setup(gp,data=mf,predict=FALSE,knots=knots,sp=sp,min.sp=min.sp,H=H,fit.method=control$fit.method)
  
  G<-GAMsetup(G) 
  
  G$C<-gam.side.conditions(G) # check for identifiability of smooth part, constrain if poss.
 
  if (is.null(G$offset)) G$offset<-rep(0,G$n)
  
  if (scale==0) 
  { if (family$family=="binomial"||family$family=="poisson") scale<-1 #ubre
    else scale <- -1 #gcv
  }
  
  G$sig2<-scale

  G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
  G$max.half<-control$mgcv.half # max step halving in Newton update mgcv
  G$min.edf<-G$nsdf+sum(null.space.dimension(G$dim,G$p.order))-dim(G$C)[1]

  object<-gam.fit(G,family=family,control=control,gamma=gamma,...)

 
  mgcv.conv<-object$mgcv.conv
  # need to check that i) following is wanted; (ii) there are free s.p.s to estimate...
  if (!control$perf.iter&&((G$fit.method=="magic"&&G$m.free>0)||(G$fit.method=="mgcv"&&G$fixed.sp<0.5&&sum(!G$fix)!=0)))
  { lsp<-log(object$sp[G$sp<0&!G$fix]) # make sure only free s.p.s are optimized!
    um<-nlm(full.score,lsp,typsize=lsp,fscale=abs(object$gcv.ubre),stepmax=1,
            ndigit=12,gradtol=1e-4,steptol=0.01,G=G,family=family,control=control,gamma=gamma)
    lsp<-um$estimate
    
    object<-attr(full.score(lsp,G,family,control,gamma=gamma),"full.gam.object")
    object$mgcv.conv<-mgcv.conv # want info on power iteration, not single evaluation calls used by nlm
    object$mgcv.conv$nlm.iterations<-um$iterations
  }
  object$fit.method<-G$fit.method
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
  object$model<-mf # store the model frame
  object$by<-G$by;object$by.exists<-G$by.exists
  object$s.type<-G$s.type
  object$p.order<-G$p.order
  object$dim<-G$dim
  object$min.edf<-G$min.edf
  object$call<-cl # needed for update() to work
  class(object)<-"gam"
  object
}

gam.check<-function(b)
# takes a fitted gam object and produces some standard diagnostic plots
{ { old.par<-par(mfrow=c(2,2))
    if (b$gcv.used) sc.name<-"GCV" else sc.name<-"UBRE"
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
    { if (length(b$df)>1&&b$mgcv.conv$iter>0)
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
}

print.gam<-function (x,...) 
# default print function for gam objects
{ print(x$family)
  cat("Formula:\n")
  print(x$formula)
  if (x$dim[1]==0)
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

gam.control<-function (irls.reg=0.0,epsilon = 1e-04, maxit = 20,globit = 20,mgcv.tol=1e-6,mgcv.half=15, 
                       nb.theta.mult=10000,trace = FALSE,fit.method="magic",perf.iter=TRUE,rank.tol=.Machine$double.eps^0.5) 
# Control structure for a gam. 
# irls.reg is the regularization parameter to use in the GAM fitting IRLS loop.
# epsilon is the tolerance to use in the IRLS MLE loop. maxit is the number 
# of IRLS iterations to use with local search for optimal s.p. after globit iterations have used global 
# searches. mgcv.tol is the tolerance to use in the mgcv call within each IRLS. mgcv.half is the 
# number of step halvings to employ in the mgcv search for the optimal GCV score, before giving up 
# on a search direction. trace turns on or off some de-bugging information.
# nb.theta.mult controls the upper and lower limits on theta estimates - for use with negative binomial  
# fit.method can be "magic" for QR/SVD method or "mgcv" for Wood (2000) method, or "fastest" to use "mgcv" for single s.p. 
# case and "magic" otherwise. 
# perf.iter TRUE to use Gu's performance iteration, FALSE to follow it up with O'Sullivan's slower approach.
# rank.tol is the tolerance to use for rank determination
{   if (!is.numeric(irls.reg) || irls.reg <0.0) stop("IRLS regularizing parameter must be a non-negative number.")
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(globit) || globit <= 0) 
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(nb.theta.mult)||nb.theta.mult<2) 
        stop("nb.theta.mult must be >= 2")
    if (!is.logical(perf.iter)) stop("power.iter must be one of TRUE or FALSE.")
    if (rank.tol<0||rank.tol>1) 
    { rank.tol=.Machine$double.eps^0.5
      warning("silly value supplied for rank.tol: reset to square root of machine precision.")
    }
    list(irls.reg=irls.reg,epsilon = epsilon, maxit = maxit,globit = globit, trace = trace, mgcv.tol=mgcv.tol,
         mgcv.half=mgcv.half,nb.theta.mult=nb.theta.mult,fit.method=fit.method,perf.iter=perf.iter,rank.tol=rank.tol)
    
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
{ free.sp<-G$sp<0&!G$fix
  G$sp[free.sp]<-exp(sp);
  G$sp[G$fix]<-0
  G$fixed.sp<-TRUE
  if (G$fit.method=="magic") # magic requires a different penalty matrix format to mgcv
  { k<-1;G$m.free<-0;G$m.off<-0;q<-NCOL(G$X)
    if (is.null(G$H)) G$H<-matrix(0,q,q)
    for (i in 1:G$m)
    { j<-G$df[i]
      off1<-G$off[i]+1;off2<-off1+j-1
      if (free.sp[i]) 
      { G$H[off1:off2,off1:off2]<-G$H[off1:off2,off1:off2]+G$sp[i]*G$S[[k]];k<-k+1}
    }
    G$S<-list() # have to reset since magic uses length of this as number of penalties
  }
  xx<-gam.fit(G,family=family,control=control,gamma=gamma)
  res<-xx$gcv.ubre
  attr(res,"full.gam.object")<-xx
  res
}

gam.fit<-function (G, start = NULL, etastart = NULL, 
    mustart = NULL, family = gaussian(), 
    control = gam.control(),gamma=1) 
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
    if (G$m>0)
    { if (G$fit.method=="mgcv")
      { k1<-0;
        S.size<-array(0,G$m)
        for (i in 1:G$m)
        { k0<-k1
          k1<-k0+G$df[i]^2
          k0<-k0+1
          S.size[i]<-mean(abs(G$S[k0:k1]))
        }
      } else # method is "magic"
      if (G$m.free>0)
      { S.size<-0
        for (i in 1:G$m.free) S.size[i]<-mean(abs(G$S[[i]])) 
      }
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
    { msp<-rep(-1,G$m.free) # free smoothing parameter vector for magic
      magic.control<-list(tol=G$conv.tol,step.half=G$max.half,maxit=control$maxit+control$globit,rank.tol=control$rank.tol)
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
        #w<-(w+control$irls.reg*mean(w))/(1+control$irls.reg) # optional weight regularization         
        
        G$w<-w
        G$X<-X[good,]  # truncated design matrix       
		if (dim(X)[2]==1) dim(G$X)<-c(length(X[good,]),1) # otherwise dim(G$X)==NULL !!
        ngoodobs <- as.integer(nobs - sum(!good))
        ncols <- as.integer(1)
        # must set G$sig2 to scale parameter or -1 here....
        G$sig2<-scale

        if (G$fit.method=="mgcv"&&G$m>0&&sum(!G$fix)>0&&!G$fixed.sp) # check that s.p.'s haven't drifted too far apart
        { temp.sp<-G$sp[!G$fix];temp.S.size<-S.size[!G$fix]*temp.sp
          # check if there is a danger of getting stuck on a flat section of gcv/ubre score...
          if (min(temp.sp)>0 && min(temp.S.size)<.Machine$double.eps^0.5*max(temp.S.size)) 
          G$sp[!G$fix]<- -1.0 # .... if so use use auto-initialization in mgcv
          if (control$trace) cat("Re-initializing smoothing parameters\n")
        } 
        if (iter>control$globit) # solution could be cycling - use more cautious optimization approach
        { G$target.edf<-G$nsdf+sum(G$edf)
        } else
        G$target.edf<- -1 # want less cautious optimization - better at local minimum avoidance
        
        if (sum(!is.finite(G$y))+sum(!is.finite(G$w))>0) 
        stop("iterative weights or data non-finite in gam.fit - regularization may help. See ?gam.control.")

        if (G$fit.method=="mgcv") G<-mgcv(G) 
        else
        { mr<-magic(G$y,G$X,msp,G$S,G$m.off,G$rank,G$H,G$C,G$w,gamma=gamma,G$sig2,G$sig2<0,
                    ridge.parameter=control$irls.reg,control=magic.control)
          G$p<-mr$b;msp<-mr$sp;G$sig2<-mr$scale;G$gcv.ubre<-mr$score
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
   
    wtdmu <- sum(weights * y)/sum(weights)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok 
    if (G$fit.method=="magic") # then some post processing is needed to extract covariance matrix etc...
    { mv<-magic.post.proc(G$X,mr,w=G$w)
      G$Vp<-mv$Vb;G$hat<-mv$hat;
      G$edf<-array(0,0) # Now some edf's for each term....
      if (G$m) for (i in 1:G$m) G$edf[i]<-sum(mv$edf[G$off[i]:(G$off[i]+G$df[i]-1)])
      G$conv<-mr$gcv.info
      G$sp[G$sp<0]<-msp
    }
	list(coefficients = as.vector(coef), residuals = residuals, fitted.values = mu, 
        family = family,linear.predictor = eta, deviance = dev,
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,  
        df.null = nulldf, y = y, converged = conv,sig2=G$sig2,edf=G$edf,hat=G$hat,
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
  { if (object$dim[1]==0) m<-0
    else m<-length(object$sp)
    n<-length(object$y)
    x<-gam.setup(object$full.formula,object$model)$x
    nc<-0;for (i in 1:m) nc<-nc+object$dim[i]
    G<-list(x=x,nsdf=object$nsdf,m=m,n=n,dim=object$dim,by=object$by)
     no.data<-TRUE # used to signal that no data provided, later
  }
  else 
  { # check that factor levels match for prediction and original fit 
    names(newdata)->nn # new data names
    for (i in 1:dim(newdata)[2]) 
    if (is.factor(object$model[,nn[i]])) # then so should newdata[[i]] be 
    { newdata[[i]]<-factor(newdata[[i]],levels=levels(object$model[,nn[i]])) # set prediction levels to fit levels
    }
    G<-gam.setup(object$full.formula,newdata,predict=TRUE)    
    if (G$nsdf!=object$nsdf) stop("Problem in predict.gam: number of model parameters don't match between fit and prediction??")
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
    { s.name<-array("",0)
      if (G$m>0) # get names for each smooth
      { stop<-object$nsdf
        by.count<-1
        for (i in 1:G$m)
        { start<-stop+1;stop<-stop+object$df[i]
          er<-names(object$coefficients)[start]
          er<-substring(er,1,nchar(er)-2)
          if (object$by.exists[i]) 
          { er<-paste(er,":",rownames(object$by)[by.count],sep="");by.count<-by.count+1} 
          s.name[i]<-er 
        }
      }
      all.names<-c(names(object$coefficients)[1:object$nsdf],s.name)
      eta<-array(o[[16]],c(G$nsdf+G$m,np));
      rownames(eta)<-all.names
      se<-array(o[[17]],c(G$nsdf+G$m,np));
      rownames(se)<-all.names 
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

  if (x$dim[1]==0) stop("Model has no smooth terms - nothing for plot.gam() to do.")
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
    oldpar<-par(mfrow=c(r,c),...)
	ylim<-c(r,c)
  } else
  { ppp<-1;oldpar<-par(...)}
  
  # now array for 1-d terms to send to predict.gam() to return smooth terms
  
  x.x<-gam.setup(x$full.formula,x$model)$x  
  xx1<-array(0,c(dim(x.x)[1],n))
  if (m>0) # should be able to do following by recycling, but it fails (non conformable array error?)
  for (i in (x$nsdf+1):(x$nsdf+length(x$covariate.shift) ) ) x.x[i,]<-x.x[i,] -x$covariate.shift[i-x$nsdf]
  j<-1;md2<-FALSE
  for (i in 1:m)
  { if (x$dim[i]==1)
    { x0<-min(x.x[j+x$nsdf,])
      x1<-max(x.x[j+x$nsdf,])
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
  rownames(xx1)<-rownames(x.x)

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
    xx2<-array(0,c(dim(x.x)[1],(n2*n2)))   # new x array for prediction
    xm<-data.frame(1:n2);ym<-data.frame(1:n2) # grid points for plotting
    if (x$nsdf>0) 
    for (i in 1:x$nsdf) xx2[i,]<-rep(0,n2*n2)
    j<-1
    for (i in 1:m)
    { if (x$dim[i]==2)
      { x0<-min(x.x[j+x$nsdf,]);x1<-max(x.x[j+x$nsdf,])
        dx<-(x1-x0)/(n2-1);
        y0<-min(x.x[j+x$nsdf+1,]);y1<-max(x.x[j+x$nsdf+1,])
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
    rownames(xx2)<-rownames(x.x)
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
    { if (is.null(select)||i==select)
      { if (interactive() && x$dim[i]<3 && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
        if (x$dim[i]==1)
        { ul<-pl1$fit[x$nsdf+i,]+pl1$se.fit[x$nsdf+i,]
          ll<-pl1$fit[x$nsdf+i,]-pl1$se.fit[x$nsdf+i,]
          if (scale==0) { ylim<-c(min(ll),max(ul))}
          title<-paste("s(",rownames(x.x)[j+x$nsdf],",",as.character(round(x$edf[i],2)),")",sep="")
          plot(xx1[x$nsdf+j,],pl1$fit[x$nsdf+i,],type="l",xlab=rownames(x.x)[j+x$nsdf],ylim=ylim,ylab=title)
	  lines(xx1[x$nsdf+j,],ul,lty=2)
          lines(xx1[x$nsdf+j,],ll,lty=2)
	  if (rug) 
          { if (jit) rug(jitter(as.numeric(x.x[x$nsdf+j,]+x$covariate.shift[j])))
             else rug(as.numeric(x.x[x$nsdf+j,]+x$covariate.shift[j]))
	  }
        } else if (x$dim[i]==2)
        { xla<-rownames(x.x)[j+x$nsdf];yla<-rownames(x.x)[j+x$nsdf+1]
          title<-paste("s(",xla,",",yla,",",as.character(round(x$edf[i],2)),")",sep="")
          if (pers) 
          { persp(xm[[i]],ym[[i]],matrix(pl2$fit[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,zlab=title,theta=theta,phi=phi)
          } else
          { sp.contour(xm[[i]],ym[[i]],matrix(pl2$fit[x$nsdf+i,],n2,n2),matrix(pl2$se.fit[x$nsdf+i,],n2,n2),
                     xlab=xla,ylab=yla,zlab=title,se.mult=se2.mult)
            if (rug) points(x.x[j+x$nsdf,]+x$covariate.shift[j],x.x[j+1+x$nsdf,]+x$covariate.shift[j+1],pch=".")
          } 
        } else
        { warning("no automatic plotting for smooths of more than one variable")
        }
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
    { if (is.null(select)||i==select)
      { if (interactive() && x$dim[i]<3 && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
        if (x$dim[i]==1)
        { title<-paste("s(",rownames(x.x)[j+x$nsdf],",",as.character(round(x$edf[i],2)),")",sep="")
          if (scale==0) ylim<-range(pl1[x$nsdf+i,])
          plot(xx1[x$nsdf+j,],pl1[x$nsdf+i,],type="l",,xlab=rownames(x.x)[j+x$nsdf],ylab=title,ylim=ylim)
          if (rug) 
	  	{ if (jit) rug(jitter(as.numeric(x.x[x$nsdf+j,]+x$covariate.shift[j])))
          else rug(as.numeric(x.x[x$nsdf+j,]+x$covariate.shift[j]))
		}
        } else if (x$dim[i]==2)
        { xla<-rownames(x.x)[j+x$nsdf];yla<-rownames(x.x)[j+x$nsdf+1]
          title<-paste("s(",xla,",",yla,",",as.character(round(x$edf[i],2)),")",sep="")
          if (pers) persp(xm[[i]],ym[[i]],matrix(pl2[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,zlab=title,theta=theta,phi=phi)
          else
          { contour(xm[[i]],ym[[i]],matrix(pl2[x$nsdf+i,],n2,n2),xlab=xla,ylab=yla,main=title)
            if (rug) points(x.x[j+x$nsdf,]+x$covariate.shift[j],x.x[j+1+x$nsdf,]+x$covariate.shift[j+1],pch=".")
          }  

        } else
        { warning("no automatic plotting for smooths of more than one variable")}
      }
      j<-j+x$dim[i]
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
  residual.df<-length(object$y)-object$nsdf-sum(object$edf)
  if (object$nsdf>0)
  { p.coeff<-object$coefficients[1:object$nsdf]
    p.t<-p.coeff/se[1:object$nsdf]
    if (object$gcv.used)
    p.pv<-2*pt(abs(p.t),df=residual.df,lower.tail=FALSE)
    else
    p.pv<-2*pnorm(abs(p.t),lower.tail=FALSE)
  } 
  else {p.coeff<-p.t<-p.pv<-array(0,0)}
  m<-length(object$edf)
  s.pv<-chi.sq<-array(0,0)
  if (m>0) # form test statistics for each smooth
  { stop<-object$nsdf
    by.count<-1
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
      if (object$by.exists[i]) 
      { er<-paste(er,":",rownames(object$by)[by.count],sep="");by.count<-by.count+1} 
      names(chi.sq)[i]<-er
      if (object$gcv.used)
      s.pv[i]<-pf(chi.sq[i]/object$edf[i],df1=max(1,object$edf[i]),df2=residual.df,lower.tail=FALSE) 
      else
      s.pv[i]<-pchisq(chi.sq[i],df=max(1,object$edf[i]),lower.tail=FALSE)
    }
  }
  r.sq<- 1 - var(object$y-object$fitted.values)*(object$df.null-1)/(var(object$y)*residual.df) 
  dev.expl<-(object$null.deviance-object$deviance)/object$null.deviance
  ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,edf=object$edf,
       s.pv=s.pv,scale=object$sig2,r.sq=r.sq,family=object$family,formula=object$formula,n=object$df.null,
       dev.expl=dev.expl)
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
  cat("\nR-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5),
      "   Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%",sep="")
  if (is.null(x$ubre)) cat("\nGCV score = ",formatC(x$gcv,digits=5)," ",sep="")
  else cat("\nUBRE score = ",formatC(x$ubre,digits=5),sep="")
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
  if (length(g2)!=n) stop("grid vectors are different lengths")
  if (length(d1)!=length(d2)) stop("data vectors are of different lengths")
  if (dist<0) stop("supplied dist negative")
  res<-array(FALSE,n)
  for (i in 1:n)
  { md<-min(((d1-g1[i])^2+(d2-g2[i])^2)^0.5)
    if (md>dist) res[i]<-TRUE # grid is too far from data
  }
  res
}

vis.gam<-function(x,view=NULL,cond=list(),n.grid=30,too.far=0,col=NA,color="topo",se=-1,type="link",zlim=NULL,...)
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
  
  if (is.null(view)) # get default view if none supplied
  { gp<-gam.parser(x$full.formula)
    view<-c(all.vars(gp$pf),gp$v.names[-1],gp$by.names[gp$by.names!="NA"])
    if (length(view)>1) view<-view[1:2]
    else stop("Model doesn't seem to have enough terms to do anything useful")
  }
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
  names(newd)<-m.name
  newd[[view[1]]]<-v1
  newd[[view[2]]]<-v2
  # call predict.gam to get predictions.....
  if (type=="link") zlab<-paste("linear predictor")
  else if (type=="response") zlab<-type
  else stop("type must be \"link\" or \"response\"")
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
    surf.col<-round(surf.col*50)
    if (color=="heat") pal<-heat.colors(50)
    else if (color=="topo") pal<-topo.colors(50)
    else if (color=="cm") pal<-cm.colors(50)
    else if (color=="terrain") pal<-terrain.colors(50)
    else stop("color scheme not recognised")
    surf.col[surf.col<1]<-1;surf.col[surf.col>50]<-50 # otherwise NA tiles can get e.g. -ve index
    if (is.na(col)) col<-pal[as.array(surf.col)]
    z<-matrix(fv$fit,n.grid,n.grid)
    persp(m1,m2,z,xlab=view[1],ylab=view[2],zlab=zlab,col=col,zlim=c(min.z,max.z),...)
  } else # add standard error surfaces
  { subs<-paste("red/green are +/-",se,"s.e.")
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
    persp(m1,m2,z,xlab=view[1],ylab=view[2],col=col,zlab=zlab,zlim=zlim,border="green",sub=subs,...)
    par(new=TRUE) # don't clean device
    z<-fv$fit;z<-matrix(z,n.grid,n.grid)
    persp(m1,m2,z,xlab=view[1],ylab=view[2],col=col,zlim=zlim,zlab=zlab,border="black",sub=subs,...)
    par(new=TRUE) # don't clean device
    z<-fv$fit+se*fv$se.fit;z<-matrix(z,n.grid,n.grid)
    persp(m1,m2,z,xlab=view[1],ylab=view[2],col=col,zlim=zlim,zlab=zlab,border="red",sub=subs,...)
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
  if (!is.null(w)) X<-w*X # use recycling rule to form diag(w)%*%X cheaply 
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






.First.lib <- function(lib, pkg) {
    library.dynam("mgcv", pkg, lib)
    cat("This is mgcv 0.9-5 \n")
}


###############################################################################
### ISSUES.....






/* library of routines designed for unstructured GCV problems:
   */
#define ANSI   
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"
#include "gcv.h"

/* routines based on singular value decomposition.... They are not as well
   structured as they could be and have been largely superceeded by the
   more efficient tridiagonalisation based routines that follow */

double svdopt(B,S,X,W,y,f,r,method,outcv) matrix B,S,X,W,y,f;double r,*outcv;int method;

/* Let A be an influence matrix:

   A = X inv(B+rS)X'W

   and assume f contains the fitted values corresponding to the B and r
   supplied:
   The routine finds r to minimise either:
   1: GCV score
   2: OCV score
   using a singular value decomposition algorithm.
   B must be positive definite.
   There is O(n^3) start up cost for this routine, but thereafter the
   evaluations are (O(n^2)) (about 4n^2, actually ). The routine returns the
   optimal r value.
   X in this formulation is a design matrix. So that usually
   B = X'WX

*/

{ matrix e,XU,U,E;
  int gridpoints=25,i,ok;
  double range=100000.0,rmin,cv,cvmin,dr,c1,c2,c3,cnew,r1,r2,r3,rnew,tau;
  tau=2.0/(1.0+sqrt(5.0));

  /* Now find the decomposition of (B+rS) that allows quick inversion */
  U=initmat(B.r,B.r);E=initmat(B.c,1L);
  suminvert(S,B,U,E);
  /* The influence matrix is now: A=X U Diag(r/(r E.V[i]+1)) U'X'W/r. Accordingly
     obtain XU..... */
  XU=initmat(X.r,U.c);multi(2,XU,X,U,0,0,0); /* added 9/11/96 - since previous version lacked design matrix. */
  /* form f=Ay */
  e=initmat(y.r,1L); /* required as offset by svdoptcv */


  r=m1norm(B)/m1norm(S);
  freemat(U);
  rmin=r;cvmin=svdoptcv(r,XU,E,W,e,y,method);
  dr=exp(2.0/(gridpoints-1)*log(range));r/=range;
  /* gridsearch for minimiser */
  for (i=0;i<gridpoints;i++)
  { cv=svdoptcv(r,XU,E,W,e,y,method);
    if (cv<cvmin) { cvmin=cv;rmin=r; }
    r*=dr;
  }
  /* Now polish the minimum, by golden section */
  r3=rmin*dr;r1=rmin/dr;c2=cvmin;
  c3=svdoptcv(r3,XU,E,W,e,y,method);
  c1=svdoptcv(r1,XU,E,W,e,y,method);
  if ((c2<c3)&&(c2<c1))
  { r2=tau*(r3-r1);
    c2=svdoptcv(r2,XU,E,W,e,y,method);
    /* perform golden section search for smoothing parameter */
    ok=1;
    if (ok)
    while ((r3-r1)>0.00001*r2)
    { if ((r2-r1)<(r3-r2))
      { rnew=r1+(r3-r1)*tau;
        cnew=svdoptcv(rnew,XU,E,W,e,y,method);
	     if (cnew<c2)
	     { r1=r2;r2=rnew;c2=cnew;}
	     else { r3=rnew;}
      } else
      { rnew=r1+(r3-r2)*(1.0-tau);
	     cnew=svdoptcv(rnew,XU,E,W,e,y,method);
	     if (cnew<c2)
	     { r3=r2;r2=rnew;c2=cnew;}
	     else { r1=rnew;}
      }
    }
    rmin=rnew;cvmin=cnew;
  }
  r=rmin;
  for (i=0;i<f.r;i++) e.V[i]=y.V[i]*W.V[i];
  e=vecmult(XU,e,1);
  for (i=0;i<f.r;i++) e.V[i]/=r*E.V[i]+1.0;
  e=vecmult(XU,e,0);
  for (i=0;i<f.r;i++) f.V[i]=e.V[i];
  freemat(e);freemat(XU);freemat(E);
  *outcv=cvmin;
  return(rmin);
}



double svdoptcv(r,ZU,E,W,e,y,method) double r;matrix ZU,E,W,e,y;int method;

/* service routine for svdopt() - calculates CV scores */

{ matrix EUZ,TrA,z;
  int i,j;
  double tr,cv;
  EUZ=initmat(ZU.c,ZU.r);TrA=initmat(ZU.r,1L);
  /* get elements of trace of A given r */
  for (i=0;i<EUZ.r;i++) for (j=0;j<EUZ.c;j++)
  EUZ.M[i][j]=ZU.M[j][i]/(r*E.V[i]+1.0);
  for (i=0;i<ZU.r;i++)
  { TrA.V[i]=0.0;
    for (j=0;j<ZU.c;j++) TrA.V[i]+=ZU.M[i][j]*EUZ.M[j][i];TrA.V[i]*=W.V[i];
  }
  /* Now form the residual vector */
  z=initmat(y.r,1L);
  for (i=0;i<y.r;i++) z.V[i]=W.V[i]*y.V[i];
  z=vecmult(EUZ,z,0);
  z=vecmult(ZU,z,0);
  for (i=0;i<y.r;i++) z.V[i]=y.V[i]-e.V[i]-z.V[i];
  cv=0.0;tr=0.0;
  if (method==1) /* then its GCV */
  { for (i=0;i<z.r;i++) { cv+=z.V[i]*W.V[i]*z.V[i];tr+=1.0-TrA.V[i];}
    cv/=tr*tr;
  } else
  if (method==2) /* OCV */
  for (i=0;i<z.r;i++)
  { cv+=z.V[i]*W.V[i]*z.V[i]/((1.0-TrA.V[i])*(1.0-TrA.V[i]));}
  freemat(z);freemat(EUZ);freemat(TrA);return(cv);
}


double zsvdopt(Z,B,S,X,W,f,y,r,method,outcv) matrix Z,B,S,X,W,f,y;double r,*outcv;int method;

/* Let A be an influence matrix:

   A = XZ inv(Z'(B+rS)Z) Z'X'W

   in general f must contain the fitted values corresponding to the B and r
   supplied, unless it is true that Ay gives the expected values according
   to the model in which case supply f as a vector with f.r=0;
   (this condition is not in general true for inequlity constrained problems -
    but is true for thin plate splines subjected only to the usual equality
    conditions).
   The routine finds r to minimise either:
   1: GCV score
   2: OCV score
   using a singular value decomposition algorithm.
   Z'BZ must be positive definite.
   There is O(n^3) start up cost for this routine, but thereafter the
   evaluations are (O(n^2)) (about 4n^2, actually ). The routine returns the
   optimal r value.
   X in this formulation is a design matrix.

*/

{ matrix e,ZU,BZ,ZSZ,ZBZ,U,E;
  int gridpoints=25,i,j,k,ok;
  double range=100000.0,rmin,cv,cvmin,dr,c1,c2,c3,cnew,r1,r2,r3,rnew,tau;
  tau=2.0/(1.0+sqrt(5.0));
  /* form Z'BZ and Z'SZ these are symmetric matrices, so the op count isn't too bad */
  BZ=initmat(B.r,Z.c);
  for (i=0;i<BZ.r;i++) for (j=0;j<BZ.c;j++) for (k=0;k<B.c;k++)
  BZ.M[i][j]+=B.M[i][k]*Z.M[k][j];ZBZ=initmat(Z.c,Z.c);
  for (i=0;i<ZBZ.r;i++) for (j=0;j<=i;j++)
  { for (k=0;k<Z.r;k++) ZBZ.M[i][j]+=Z.M[k][i]*BZ.M[k][j];
    ZBZ.M[j][i]=ZBZ.M[j][i];
  }
  freemat(BZ);BZ=initmat(S.r,Z.c);
  for (i=0;i<BZ.r;i++) for (j=0;j<BZ.c;j++) for (k=0;k<S.c;k++)
  BZ.M[i][j]+=S.M[i][k]*Z.M[k][j];ZSZ=initmat(Z.c,Z.c);
  for (i=0;i<ZSZ.r;i++) for (j=0;j<=i;j++)
  { for (k=0;k<Z.r;k++) ZSZ.M[i][j]+=Z.M[k][i]*BZ.M[k][j];
    ZSZ.M[j][i]=ZSZ.M[j][i];
  }
  /* Now find the decomposition of Z'(B+rS)Z that allows quick inversion */
  U=initmat(Z.c,Z.c);E=initmat(Z.c,1L);
  suminvert(ZSZ,ZBZ,U,E);
  /* The influence matrix is now: A=ZU Diag(r/(r E.V[i]+1)) Z'U'W/r. Accordingly
     obtain ZU..... */
  ZU=initmat(X.r,U.c);multi(3,ZU,X,Z,U,0,0,0); /* added 9/11/96 - since previous version lacked design matrix. */
  /* form e=f-Ay */
  e=initmat(y.r,1L);
  if (f.r) /* calculate offset for fitting problem */
  { for (i=0;i<e.r;i++) e.V[i]=y.V[i]*W.V[i];
    e=vecmult(ZU,e,1);
    for (i=0;i<e.r;i++) e.V[i]/=r*E.V[i]+1.0;
    e=vecmult(ZU,e,0);
    for (i=0;i<e.r;i++) e.V[i]=f.V[i]-e.V[i]; /* as A changes f=e + Ay */
  }
  r=m1norm(ZBZ)/m1norm(ZSZ);
  freemat(BZ);freemat(ZSZ);freemat(ZBZ);
  freemat(U);
  rmin=r;cvmin=svdoptcv(r,ZU,E,W,e,y,method);
  dr=exp(2.0/(gridpoints-1)*log(range));r/=range;
  /* gridsearch for minimiser */
  for (i=0;i<gridpoints;i++)
  { cv=svdoptcv(r,ZU,E,W,e,y,method);
    if (cv<cvmin) { cvmin=cv;rmin=r; }
    r*=dr;
  }
  /* Now polish the minimum, by golden section */
  r3=rmin*dr;r1=rmin/dr;c2=cvmin;
  c3=svdoptcv(r3,ZU,E,W,e,y,method);
  c1=svdoptcv(r1,ZU,E,W,e,y,method);
  if ((c2<c3)&&(c2<c1))
  { r2=tau*(r3-r1);
    c2=svdoptcv(r2,ZU,E,W,e,y,method);
    /* perform golden section search for smoothing parameter */
    ok=1;
    if (ok)
    while ((r3-r1)>0.00001*r2)
    { if ((r2-r1)<(r3-r2))
      { rnew=r1+(r3-r1)*tau;
	     cnew=svdoptcv(rnew,ZU,E,W,e,y,method);
	     if (cnew<c2)
	     { r1=r2;r2=rnew;c2=cnew;}
	     else { r3=rnew;}
      } else
      { rnew=r1+(r3-r2)*(1.0-tau);
	     cnew=svdoptcv(rnew,ZU,E,W,e,y,method);
	     if (cnew<c2)
	     { r3=r2;r2=rnew;c2=cnew;}
	     else { r1=rnew;}
      }
    }
    rmin=rnew;cvmin=cnew;
  }
  freemat(e);freemat(ZU);freemat(E);  *outcv=cvmin;
  return(rmin);
}

/* Routines based on tridiagonalisation of a transformed smoothness matrix
   See Gu, as reported in Wahba 1990. These routines are better than the svd
   routines, furthermore multiple smoothing parameters are dealt with. */

double TrInf(matrix *X,matrix *Z,matrix *w,matrix *S,double rho)

/* Finds the trace of the influence matrix:

   rho*XZ[Z'X'WXZ*rho + Z'SZ]^{-1}Z'X'W

   Z is the column basis for a null space of constraints.
   It is assumed that these are stored efficiently as householder
   transformations in the manner of routine QT. For example to impose
   Cp=0: QT(Z,C,0);Z.r=C.r;
   written and tested 26/11/97

*/
{ double *rw,zz;
  matrix R,Q,U,T,l0,l1;
  long n,i,j,k;
  n=X->r;
  /* Get Q'R=W^{0.5} X Z */
  rw=(double *)calloc((size_t)n,sizeof(double));
  for (i=0;i<n;i++) rw[i]=sqrt(w->V[i]);   /* rw contains l.d. of W^{0.5} */
  if (Z->r)
  { R=initmat(n,Z->c);
    mcopy (X,&R);HQmult(R,*Z,0,0);R.c-=Z->r;  /* R=XZ */
  } else
  { R=initmat(n,X->c);mcopy(X,&R);}          /* case when Z=I */
  for (i=0;i<n;i++) for (j=0;j<R.c;j++) R.M[i][j]*=rw[i];  /* R=W^{0.5} XZ */
  Q=initmat(n,n);
  QR(&Q,&R);  /* done in up to O(n^3) */
  /* invert R - it is this step that requires n>=nz*/
  freemat(Q);
  R.r=R.c; /* getting rid of zero part */
  InvertTriangular(&R);   /* R now contains L */
  T=initmat(S->r,S->c);
  mcopy(S,&T);
  if (Z->r)
  { HQmult(T,*Z,1,1);HQmult(T,*Z,0,0);
    T.r=T.c=Z->c-Z->r;
  }
  U=initmat(T.r,T.c);
  multi(3,U,R,T,R,1,0,0);
  for (j=T.c-1;j>=0;j--) for (i=0;i<T.r;i++)
  { zz=0.0; for (k=0;k<=j;k++) zz+=T.M[i][k]*R.M[k][j];T.M[i][j]=zz;}
  for (i=T.r-1;i>=0;i--) for (j=0;j<=i;j++)
  { zz=0.0; for (k=0;k<=i;k++) zz+=R.M[k][i]*T.M[k][j];T.M[i][j]=zz;}
  for (i=T.r-1;i>=0;i--) for (j=0;j<=T.c;j++) T.M[j][i]=T.M[i][j];
  zz=0.0;
  for (i=0;i<T.r;i++) for (j=0;j<T.c;j++) zz+=fabs(T.M[i][j]-U.M[i][j]);
  freemat(U);
  freemat(R);U=initmat(R.c,R.c);
  UTU(&T,&U);
  l0=initmat(T.r,1L);l1=initmat(T.r-1,1L);
  for (i=0;i<T.r;i++) T.M[i][i] += rho;  /* forming I*r + T */
  tricholeski(&T,&l0,&l1);                  /* get LL'=I*r + T */
  zz=triTrInvLL(&l0,&l1)*rho;
  freemat(l0);freemat(l1);freemat(U);freemat(T);
  free(rw);
  return(zz);
}

double EScv(matrix *T,matrix *l0,matrix *l1,matrix *x,double nx, matrix *z,double r,long n,
            double *trace,double *ress,double *sig2)

/* service routine for EasySmooth - obtains GCV score (*sig<=0.0) or UBRE score */

{ long i;
  int ubre=0;
  double rss=0.0,el,tr;
  if (*sig2>0.0) ubre=1;
  for (i=0;i<T->r;i++) T->M[i][i] += r;  /* forming I*r + T */
  tricholeski(T,l0,l1);                  /* get LL'=I*r + T */
  tr=1.0-triTrInvLL(l0,l1)*r/n;

  z->r=x->r;
  bicholeskisolve(x,z,l0,l1);
  for (i=0;i<x->r;i++)
  { el=z->V[i]- r * x->V[i];rss+=el*el;T->M[i][i] -= r;
  }
  rss+=nx;
  if (!ubre) *sig2=rss/(n*tr); /* variance estimate - GCV only */
  z->r=n;rss/=n;
  if (ubre) /* then use UBRE */
  { el = rss - 2*(*sig2)*tr+ *sig2;tr=tr*tr;}
  else /* use GCV */
  { tr=tr*tr;el=rss/tr;}
  *ress=rss;*trace=tr;
  /* debugging */
  return(el)/*(n-tr*n-5)*(n-tr*n-5))*/;
}

double EasySmooth(matrix *T,matrix *z,double *v,double *df,long n,double *sig2)

/* routine that minimises
   (||(I-r*(I*r + T)^{-1})z ||^2+||x||^2) / (n - r*Tr(I*r+T)^{-1})^2
   This is a gcv score for a problem with an influence matrix:
   r JZ[Z'J'WJZ r + Z'SZ]^{-1}Z'J'W
   This will have been re-written using the decomposition:
   W^{0.5}JZ=QR as....
   r W^{-0.5}QR(R'R r + Z'SZ)^{-1} R'Q'W^{0.5}
   which becomes....
   r W^{-0.5}(I*r+L'Z'SZL)^{-1}W^{0.5} where L=R^{-1}
   The decomposition UTU' = L'Z'SZL is then formed. As is z=UQ'W^{0.5} y.

   Alternatively the routine minimises an equivalent UBRE score if *sig2>0.0
   on entry.

   returns smoothing parameter and minimum gcv/ubre score.
*/

{ matrix l0,l1,x;
  double r,V,rm,tr,maxr,minr,minV,rss,r0,r1,rt,r1t,ft,f1t,tau,nx;
  long k,mesh=50L,i;
  int gcv=1;
  minV=0.0;
  if (*sig2>0.0) gcv=0;
  /* Initialise work space matrices */
  l0=initmat(T->r,1L);l1=initmat(T->r-1,1L);x=initmat(l0.r,1L);
  /* get initial smooth estimates */
  tr=0.0;for (i=0;i<T->r;i++) tr+=T->M[i][i];tr/=T->r;
  minr=tr*0.0000001;maxr=tr*100000.0;rm=exp(log(maxr/minr)/mesh);
  r=minr/rm;
  nx=0.0;for (i=x.r;i<n;i++) nx+=z->V[i]*z->V[i]; /* ||x||^2 */
  for (k=0;k<mesh;k++)
  { r*=rm;
    if (gcv) *sig2=-1.0;
    V=EScv(T,&l0,&l1,&x,nx,z,r,n,&tr,&rss,sig2);
    if (V<minV||k==0L) { minV=V;minr=r;}
    if (n*tr<1.0)
    break;
  }
  /* golden section search to polish minimisation */
  r0=minr/rm;r1=rm*minr;
  tau=2.0/(1.0+sqrt(5.0));
  rt=r0+(r1-r0)*tau;
  if (gcv) *sig2=-1.0;
  ft=EScv(T,&l0,&l1,&x,nx,z,rt,n,&tr,&rss,sig2);
  r1t=r0+(r1-r0)*(1.0-tau);
  if (gcv) *sig2=-1.0;
  f1t=EScv(T,&l0,&l1,&x,nx,z,r1t,n,&tr,&rss,sig2);
  while ((rt-r1t)>1e-5*fabs(rt+r1t))
  { if (ft<f1t)
    { r0=r1t;r1t=rt;f1t=ft;rt=r0+(r1-r0)*tau;
      if (gcv) *sig2=-1.0;
      ft=EScv(T,&l0,&l1,&x,nx,z,rt,n,&tr,&rss,sig2);
    } else
    { r1=rt;rt=r1t;ft=f1t;r1t=r0+(r1-r0)*(1.0-tau);
      if (gcv) *sig2=-1.0;
      f1t=EScv(T,&l0,&l1,&x,nx,z,r1t,n,&tr,&rss,sig2);
    }
  }
  minr=rt;
  minV=ft;

  *v = minV;
  *df=(1.0-sqrt(tr))*n;
  if (gcv) *sig2=-1.0;
  EScv(T,&l0,&l1,&x,nx,z,minr,n,&tr,&rss,sig2); /* here for debugging purposes */
  freemat(l0);freemat(l1);freemat(x);
  return(minr);
}

double SingleSmooth(matrix *y,matrix *X,matrix *Z,matrix *w,matrix *S,matrix *p,
                    double *df,double *sig2)

/* NOTE: Does not work well with rank deficient penalties.
   solves smoothing problems with influence matrix
   rho*XZ[Z'X'WXZ*rho + Z'SZ]^{-1}Z'X'W
   smoothed values are returned in y. The algorithm is a generalisation of Gu's
   (1989), the generalisation widens the class of design matrices usable.
   returns rho, and parameter values in p.
   Z is the column basis for a null space of constraints.
   It is assumed that these are stored efficiently as householder
   transformations in the manner of routine QT. For example to impose
   Cp=0: QT(Z,C,0);Z.r=C.r;

   if *sig2>0.0 on entry then UBRE is used (which assumes that the variance
   sig2 is known), otherwise GCV is used.
*/
{ double *rw,v,rho,zz;
  matrix Wy,R,Q,U,T,z,l0,l1;
  long n,i,j,k;
  n=y->r;
  /* Get Q'R=W^{0.5} X Z */
  rw=(double *)calloc((size_t)n,sizeof(double));
  for (i=0;i<n;i++) rw[i]=sqrt(w->V[i]);   /* rw contains l.d. of W^{0.5} */
  Wy=initmat(n,1L);
  for (i=0;i<n;i++) Wy.V[i]=rw[i]*y->V[i];
  if (Z->r)
  { R=initmat(n,Z->c);
    mcopy (X,&R);HQmult(R,*Z,0,0);R.c-=Z->r;  /* R=XZ */
    /*matmult(R,*X,*Z,0,0); FZ *//* R=XZ - up to O(n^3) */
  } else
  { R=initmat(n,X->c);mcopy(X,&R);}          /* case when Z=I */
  for (i=0;i<n;i++) for (j=0;j<R.c;j++) R.M[i][j]*=rw[i];  /* R=W^{0.5} XZ */
  Q=initmat(n,n);
  QR(&Q,&R);  /* done in up to O(n^3) */
  /* invert R - it is this step that requires n>=nz*/
  R.r=R.c; /* getting rid of zero part */
  InvertTriangular(&R);   /* R now contains L */
  /* following code is inefficient and should be optimized - DONE 26/11/97 and
     tested (quickly with matest)*/
  T=initmat(R.c,R.c);
  U=initmat(S->r,S->r);mcopy(S,&U);
  if (Z->r)
  { HQmult(U,*Z,1,1);HQmult(U,*Z,0,0);
    U.r=U.c=Z->c-Z->r;
  }
  //multi(3,T,R,U,R,1,0,0);
  for (j=U.c-1;j>=0;j--) for (i=0;i<U.r;i++)
  { zz=0.0; for (k=0;k<=j;k++) zz+=U.M[i][k]*R.M[k][j];U.M[i][j]=zz;}
  for (i=U.r-1;i>=0;i--) for (j=0;j<=i;j++)
  { zz=0.0; for (k=0;k<=i;k++) zz+=R.M[k][i]*U.M[k][j];U.M[i][j]=zz;}
  for (i=U.r-1;i>=0;i--) for (j=0;j<=U.c;j++) U.M[j][i]=U.M[i][j];
  mcopy(&U,&T);
  freemat(U);
  U=initmat(R.c,R.c);
  UTU(&T,&U);
  z=initmat(n,1L);
  for (i=0;i<n;i++) z.V[i]=rw[i]*y->V[i]; /* z=W^{1/2}y */
  OrthoMult(&Q,&z,0,Q.r,0,1,1);           /* z=QW^{1/2}y */
  z.r=T.r;                                /* z=[I,0]QW^{1/2}y */
  OrthoMult(&U,&z,1,T.r-2,1,1,0);         /* z=U'[I,0]QW^{1/2}y */
  z.r=n;
  rho=EasySmooth(&T,&z,&v,df,n,sig2);
  l0=initmat(T.r,1L);l1=initmat(T.r-1,1L);
  for (i=0;i<T.r;i++) T.M[i][i] += rho;
  tricholeski(&T,&l0,&l1);
  for (i=0;i<T.r;i++) T.M[i][i] -= rho;
  z.r=T.r;
  y->r=z.r;
  bicholeskisolve(y,&z,&l0,&l1);        /* y=(I*rho+T)^{-1} z */
  OrthoMult(&U,y,1,T.r-2,0,1,0);        /* y=U(I*rho+T)^{-1} z */
  p->r=R.r;
  matmult(*p,R,*y,0,0);
  if (Z->r)
  { p->r=Z->c;
    for (i=R.r;i<Z->c;i++) p->V[i]=0.0;
    HQmult(*p,*Z,1,0);
    /**p=vecmult(*Z,*p,0); FZ */
  }
  for (i=0;i<p->r;i++) p->V[i]*=rho;    /* p=rho*ZLU(I+rho+T)^{-1} z */
  y->r=n;for (i=T.r;i<y->r;i++) y->V[i]=0.0;
  OrthoMult(&Q,y,0,Q.r,1,1,1);   /* smoothed values now in y */
  for (i=0;i<n;i++) y->V[i]*=rho/rw[i];
  freemat(Q);freemat(U);freemat(R);freemat(T);freemat(z);freemat(Wy);
  freemat(l0);freemat(l1);free(rw);
  return(rho);
}

/* tediouscv and boringHg are debugging routines for MultiSmooth() */

double tediouscv(matrix R,matrix Q,matrix *LZSZL,matrix *y,double *rw,
                 double *trial,double rho,int m,double *tr,double *rss,double sig2)


{ long i,j,l,n;
  matrix T,U,z,l0,l1,x;
  double v,nx;
  n=y->r;
  T=initmat(LZSZL[0].r,LZSZL[0].r);
  U=initmat(T.r,T.r);
  z=initmat(n,1L);
  for (i=0;i<T.r;i++) for (j=0;j<T.c;j++)
  T.M[i][j]=exp(trial[0])*LZSZL[0].M[i][j];
  for (l=1;l<m;l++) for (i=0;i<T.r;i++) for (j=0;j<T.c;j++)
  T.M[i][j] += exp(trial[l])*LZSZL[l].M[i][j];
  UTU(&T,&U);    /* Form S=UTU' */
  z.r=n;
  for (i=0;i<n;i++) z.V[i]=rw[i]*y->V[i]; /* z=W^{1/2}y */
  OrthoMult(&Q,&z,0,Q.r,0,1,1);           /* z=QW^{1/2}y */
  nx=0.0;for (i=R.r;i<n;i++) nx+=z.V[i]*z.V[i];
  z.r=R.r;                                /* z=[I,0]QW^{1/2}y */
  OrthoMult(&U,&z,1,T.r-2,1,1,0);         /* z=U'[I,0]QW^{1/2}y */
  z.r=n;
  l0=initmat(T.r,1L);l1=initmat(T.r-1,1L);x=initmat(T.r,1L);
  v=EScv(&T,&l0,&l1,&x,nx,&z,rho,n,tr,rss,&sig2);
  freemat(l0);freemat(l1);freemat(x);freemat(T);freemat(U);freemat(z);
  return(v);
}

void boringHg(matrix R,matrix Q,matrix *LZSZL,matrix *y,double *rw,
                 double *trial,double rho,int m,double sig2,double dt1)
/* Does f.d. estimation of gradient and Hessian dt is interval to use for
   differencing  */

{ double f,v,v1,v2,v3,v4,tr,rss,tr1,rss1,
         r1,r2,r3,r4,t1,t2,t3,t4,t,r;
  int i,j,k;
  matrix a,M,p;
  printf("\nHit Return ... ");getc(stdin);
  v=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&tr,&rss,sig2);t=tr;r=rss;
  printf("\ntedious cv = %g\n",v);
  for (i=0;i<m;i++)
  { trial[i]+=dt1;
    v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&tr1,&rss1,sig2);
    trial[i] -= dt1;
    tr1=(tr1-tr)/(dt1);
    rss1=(rss1-rss)/(dt1);
    printf("\ng%d = %g drss=%g  dtr=%g",i,(v1-v)/dt1,rss1,tr1);
  }
  printf("\n");
  for (i=0;i<m;i++) for (j=0;j<=i;j++)
  { if (i!=j)
    { /*trial[i] += dt1/2;trial[j]+=dt2/2;
      v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t1,&r1,sig2);
      trial[i]-=dt1;
      v2=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t2,&r2,sig2);
      trial[j] -= dt1;
      v3=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t3,&r3,sig2);
      trial[i] +=dt1;
      v4=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t4,&r4,sig2);
      trial[i] -=dt1/2; trial[j]+=dt1/2;
      f=(v4+v2)-(v3+v1);f/= -dt1*dt1;
      f=(r4+r2)-(r3+r1);f/= -dt1*dt1;
      f=(t4+t2)-(t3+t1);f/= -dt1*dt1;
      printf("%8.4g  ",f);

*/
      M=initmat(6L,6L);a=initmat(6L,1L);p=initmat(6L,1L);
      trial[i]+=dt1/2;
      v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t1,&r1,sig2);
      k=0;
      M.M[k][0]=1.0;M.M[k][1]=trial[i];M.M[k][2]=trial[j];M.M[k][3]=trial[i]*trial[j];
      M.M[k][4]=trial[i]*trial[i];M.M[k][5]=trial[j]*trial[j];a.V[k]=v1;
      trial[i]-=dt1;
      v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t1,&r1,sig2);
      k=1;
      M.M[k][0]=1.0;M.M[k][1]=trial[i];M.M[k][2]=trial[j];M.M[k][3]=trial[i]*trial[j];
      M.M[k][4]=trial[i]*trial[i];M.M[k][5]=trial[j]*trial[j];a.V[k]=v1;
      trial[i]-=dt1/2;trial[j]-=dt1;
      v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t1,&r1,sig2);
      k=2;
      M.M[k][0]=1.0;M.M[k][1]=trial[i];M.M[k][2]=trial[j];M.M[k][3]=trial[i]*trial[j];
      M.M[k][4]=trial[i]*trial[i];M.M[k][5]=trial[j]*trial[j];a.V[k]=v1;
      trial[j]+=2*dt1;
      v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t1,&r1,sig2);
      k=3;
      M.M[k][0]=1.0;M.M[k][1]=trial[i];M.M[k][2]=trial[j];M.M[k][3]=trial[i]*trial[j];
      M.M[k][4]=trial[i]*trial[i];M.M[k][5]=trial[j]*trial[j];a.V[k]=v1;
      trial[i]+=2*dt1;
      v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t1,&r1,sig2);
      k=4;
      M.M[k][0]=1.0;M.M[k][1]=trial[i];M.M[k][2]=trial[j];M.M[k][3]=trial[i]*trial[j];
      M.M[k][4]=trial[i]*trial[i];M.M[k][5]=trial[j]*trial[j];a.V[k]=v1;
      trial[j]-=2*dt1;
      v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t1,&r1,sig2);
      k=5;
      M.M[k][0]=1.0;M.M[k][1]=trial[i];M.M[k][2]=trial[j];M.M[k][3]=trial[i]*trial[j];
      M.M[k][4]=trial[i]*trial[i];M.M[k][5]=trial[j]*trial[j];a.V[k]=v1;
      trial[i]-=dt1;trial[j]+=dt1;
      svdLS(M,p,a,1e-10);
      printf("%8.4g  ",p.V[3]);
      freemat(p);freemat(M);freemat(a);
    } else
    { trial[i] += dt1;
      v1=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t1,&r1,sig2);
      trial[i] -= 2*dt1;
      v2=tediouscv(R,Q,LZSZL,y,rw,trial,rho,m,&t2,&r2,sig2);
      trial[i] +=dt1;
      f=(v1-2*v+v2);f/=dt1*dt1;
    //  f=(r1-2*r+r2);f/=dt1*dt1;
     // f=(t1-2*t+t2);f/=dt1*dt1;

      printf("%8.4g\n",f);
    }
  }
}

void MultiSmooth(matrix *y,matrix *J,matrix *Z,matrix *w,matrix *S,matrix *p,
                 double *theta,long *off,int m,double *sig2)

/* routine for multiple smoothing parameter problems, with influence matrix:
   A=r*JZ(Z'J'WJZ*r + \sum_{i=1}^m \theta_i Z'S_iZ)^{-1} Z'J'W
   The routine uses Gu and Wahba's approach of alternating exact searches
   for overall smoothing parameter r, with Newton updates of the theta_i's
   The actual calculations are rather different because the method doesn't
   assume the special structure needed by Gu and Wahba's method.
   The S_i matrices should be supplied in their smallest possible form:
   if p is the total parameter vector then in principle a smoothness measure
   will be of the form p'S_ip, but much of S_i will be zeroes. Suppose that S_i
   contains an s by s non-zero sub matrix S[i], and that q is the s-vector
   consisting of the s elements of p starting from element off[i]. Then the
   measure is given by q'S[i]q. It is the S[i] matrices and array off[] that
   should be supplied to the routine. This way of proceeding means that the
   smoothing parameter vector updates require only around n^3 operations not
   m n^3.
   Problem must have datapoints >= to null space dimension (parameter vector
   length, or columns of Z) - I'm not sure if this restriction is fundamental.

   If any of the m elements of theta[] are -ve then automatic initialisation
   of smoothing parameters takes place. Otherwise the supplied theta[i]s are
   used as initial estimates. On exit theta[i] contains best estimates of
   theta[i]/rho.

   NOTE that for some problems it may make sense to initialise by 
   running first with very low s.p.s and then using the Wahba suggested 2nd
   guesses - not implemented yet.

   Z is actually supplied as a matrix containing a series of householder
   transformations, as produced by QT when fullQ==0... this improves efficiency
   considerably.... old, full Z code in brackets labelled FZ.
   To obtain Z for Cp=0: call QT(Q,C,0). Then set Q.r=C.r; Q is now used as Z.

   Supplying *sig2 as a positive number signals the routine to use UBRE, rather
   than GCV, since *sig2 is known. Supplying *sig2 as zero or negative causes
   the routine to use GCV. (added 15/1/99)
*/

{ double *rw,tr,*eta,*del,*trial,trA,a,b,*trdA,**trd2A,*db,**d2b,*da,**d2a,
          rho,v,vmin=0.0,x,**RM,*pp,**MM,**M1M,**M2M,*pp1,tdf,xx1,xx2,tol,
          *pinf,*ninf;
  long i,j,k,l,n,np,nz;
  int iter,reject,ok=1,autoinit=0,op=0,trials,ubre=0;
  matrix z,l0,l1,Q,R,C,ZC,*LZrS,*LZSZL,T,U,A,Wy,Hess,g,c,*H,*ULZSZLU,
         *ULZrS,*dAy,*dpAy,d,Ay;
  if (*sig2>0.0) ubre=1; /* use UBRE rather than GCV */
  n=y->r;  /* number of datapoints */
  np=J->c; /* number of parameters */
  for (i=0;i<m;i++) if (theta[i]<=0.0) autoinit=1;
  if (Z->r) /*nz=Z->c FZ */ nz=np-Z->r;else nz=np; /* dimension of null space */
  A=initmat(np,np); /* workspace matrix */
  c=initmat(n,1L); /*     "     vector  */
  d=initmat(n,1L);
  Ay=initmat(n,1L);
  Hess=initmat((long)m,(long)m);g=initmat((long)m,1L);
  /*for (l=0;l<m;l++) ... appears to be totally redundant matrix allocation */
  { trdA=(double *)calloc((size_t)m,sizeof(double));
    da=(double *)calloc((size_t)m,sizeof(double));
    db=(double *)calloc((size_t)m,sizeof(double));
    trd2A=(double **)calloc((size_t)m,sizeof(double));
    d2a=(double **)calloc((size_t)m,sizeof(double));
    d2b=(double **)calloc((size_t)m,sizeof(double));
    ninf=(double *)calloc((size_t)m,sizeof(double));
    pinf=(double *)calloc((size_t)m,sizeof(double));
    for (k=0;k<m;k++)
    { trd2A[k]=(double *)calloc((size_t)m,sizeof(double));
      d2a[k]=(double *)calloc((size_t)m,sizeof(double));
      d2b[k]=(double *)calloc((size_t)m,sizeof(double));
    }
  }
  /* Get Q'R=W^{0.5} J Z */
  rw=(double *)calloc((size_t)n,sizeof(double));
  for (i=0;i<n;i++) rw[i]=sqrt(w->V[i]);   /* rw contains l.d. of W^{0.5} */
  Wy=initmat(n,1L);
  for (i=0;i<n;i++) Wy.V[i]=rw[i]*y->V[i];
  if (Z->r)
  { R=initmat(n,np);mcopy(J,&R);
    HQmult(R,*Z,0,0);R.c=nz;
    /* matmult(R,*J,*Z,0,0); FZ */ /* R=JZ - up to O(n^3) */
  } else
  { R=initmat(n,np);mcopy(J,&R);}          /* case when Z=I */
  RM=R.M;
  for (i=0;i<n;i++)
  { x=rw[i];for (pp=RM[i];pp<RM[i]+R.c;pp++) *pp *= x;}  /* R=W^{0.5} JZ */
  Q=initmat(np,n);/* altered for efficient storage */
  QR(&Q,&R);  /* done in up to O(n^3) */
  /* invert R - it is this step that requires n>=nz*/
  R.r=R.c; /* getting rid of zero part */
  InvertTriangular(&R);   /* R now contains L - explicit R formation isn't very efficient */
  /* Form the matrices L'Z'S^{0.5} and L'Z'S_iZL */
  ZC=initmat(np,np);
  LZrS=(matrix *)calloc((size_t)m,sizeof(matrix));
  LZSZL=(matrix *)calloc((size_t)m,sizeof(matrix));
  ULZSZLU=(matrix *)calloc((size_t)m,sizeof(matrix));
  H= (matrix *)calloc((size_t)m,sizeof(matrix));
  ULZrS=(matrix *)calloc((size_t)m,sizeof(matrix));
  for (l=0;l<m;l++)
  { root(S+l,&C,1e-14);    /* S[l]=CC' */
    if (Z->r)
    { ZC.c=C.c;ZC.r=np;
      /* set ZC = [0,C',0]' */
      for (i=0;i<off[l];i++) for (j=0;j<C.c;j++) ZC.M[i][j]=0.0;
      for (i=off[l];i<off[l]+C.r;i++)
      for (j=0;j<C.c;j++) ZC.M[i][j]=C.M[i-off[l]][j];
      for (i=off[l]+C.r;i<np;i++)
      for (j=0;j<C.c;j++) ZC.M[i][j]=0.0;
      /* ...and apply Z efficiently, to get Z'C... */
      HQmult(ZC,*Z,1,1);ZC.r=nz;  /* bug fixed here 17/8/99 - Z was untransposed! */

      /* MM=ZC.M;M1M=Z->M;M2M=C.M;
      ZC.c=C.c;ZC.r=nz;
      for (i=0;i<ZC.r;i++) for (j=0;j<ZC.c;j++)
      { MM[i][j]=0.0;
        for (k=0;k<C.r;k++)
        MM[i][j]+=M1M[k+off[l]][i]*M2M[k][j];
      }  FZ*/
    } else
    { ZC.c=C.c;ZC.r=np;
      for (i=0;i<ZC.r;i++) for (j=0;j<ZC.c;j++) ZC.M[i][j]=0.0;
      for (i=0;i<C.r;i++) for (j=0;j<C.c;j++) ZC.M[i+off[l]][j]=C.M[i][j];
    }
    freemat(C);
    LZrS[l]=initmat(R.c,ZC.c);
    MM=LZrS[l].M;M1M=R.M;M2M=ZC.M;
    for (i=0;i<LZrS[l].r;i++) for (j=0;j<LZrS[l].c;j++)
    for (k=0;k<=i;k++) MM[i][j]+=M1M[k][i]*M2M[k][j];
    LZSZL[l]=initmat(R.c,R.c); // this memory requirement could be halved
    matmult(LZSZL[l],LZrS[l],LZrS[l],0,1);
    H[l]=initmat(R.c,R.c);
    ULZSZLU[l]=initmat(R.c,R.c); // memory requirement could be halved, but would need own choleskisolve
    //ULZrS[l]=initmat(R.c,ZC.c);
    ULZrS[l].M=H[l].M;ULZrS[l].r=R.c;ULZrS[l].c=ZC.c; // sharing memory with H[l]
  }
  /* Start the main loop */
  freemat(ZC);
  eta=(double *)calloc((size_t)m,sizeof(double));
  del=(double *)calloc((size_t)m,sizeof(double)); /* change in s.p.s */
  trial=(double *)calloc((size_t)m,sizeof(double));
  /* get initial estimates for theta_i and eta_i=log(theta_i) */
  if (autoinit)
  for (l=0;l<m;l++)
  { tr=0.0;for (i=0;i<LZSZL[l].r;i++) tr+=LZSZL[l].M[i][i];
    eta[l]=log(1.0/(n*tr));//eta[l]= -40.0;
  } else
  { x=0.0;for (i=0;i<m;i++) { eta[i]=log(theta[i]);x+=eta[i];}
    x/=m;
    for (i=0;i<m;i++) { ninf[i]=eta[i]-x-300.0;pinf[i]=ninf[i]+600.0;}
  }
  T=initmat(R.c,R.c);
  U=initmat(T.r,T.r);
  z=initmat(n,1L);
  l0=initmat(T.r,1L);l1=initmat(T.r-1,1L);
  dAy=(matrix *)calloc((size_t)m,sizeof(matrix));
  dpAy=(matrix *)calloc((size_t)m,sizeof(matrix));
  for (l=0;l<m;l++) { dAy[l]=initmat(n,1L);dpAy[l]=initmat(T.r,1L);}
  iter=0;  /* iteration counter */
  while (ok)
  { /* form combined smoothness measure */
    reject=1;
    while (reject) /* try current search direction (unless first step) */
    { x=0.0;for (i=0;i<m;i++) { trial[i]=eta[i]+del[i];x+=trial[i];}
      x/=m;
      for (i=0;i<m;i++) trial[i] -= x; /* normalising smooths */
      if (iter>1||!autoinit)   /* check smoothing parameters won't lead to overflow */
      for (i=0;i<m;i++)
      if (trial[i]<ninf[i])
      trial[i]=ninf[i];
      else if (trial[i]>pinf[i])
      trial[i]=pinf[i];
      /* form S the combined smooth measure */
      for (i=0;i<T.r;i++) for (j=0;j<T.c;j++)
      T.M[i][j]=exp(trial[0])*LZSZL[0].M[i][j];
      for (l=1;l<m;l++) for (i=0;i<T.r;i++) for (j=0;j<T.c;j++)
      T.M[i][j] += exp(trial[l])*LZSZL[l].M[i][j];
      UTU(&T,&U);    /* Form S=UTU' */
      z.r=n;
      for (i=0;i<n;i++) z.V[i]=rw[i]*y->V[i]; /* z=W^{1/2}y */
     // matrixintegritycheck();
      OrthoMult(&Q,&z,0,Q.r,0,1,1);           /* z=QW^{1/2}y */
      z.r=R.r;                                /* z=[I,0]QW^{1/2}y */
      OrthoMult(&U,&z,1,T.r-2,1,1,0);         /* z=U'[I,0]QW^{1/2}y */
      z.r=n;                                  /* z->[z,x_1]' */
      if (!ubre) *sig2=-1.0; /* setting to signal GCV rather than ubre */
      rho=EasySmooth(&T,&z,&v,&tdf,n,sig2);    /* do a cheap minimisation in rho */
      z.r=R.r;
      if (!iter||v<vmin) /* accept current step */
      { reject=0;
        /* test for convergence */
        tol=1e-4;ok=0;
        if (vmin-v>tol*(1+v)) ok=1;
        xx1=0.0;for (i=0;i<m;i++) { xx2=eta[i]-trial[i];xx1+=xx2*xx2;}
        xx1=sqrt(xx1);
        xx2=0.0;for (i=0;i<m;i++) xx2+=trial[i]*trial[i];
        xx2=sqrt(xx2);xx2=(1+xx2)*sqrt(tol);
        if (xx1>xx2) ok=1;
        xx1=0.0;for (i=0;i<m;i++) xx1+=g.V[i]*g.V[i];xx1=sqrt(xx1);
        if (xx1>pow(tol,1.0/3)*(1+v)) ok=1;
        for (i=0;i<m;i++) eta[i]=trial[i];
        vmin=v;
      } else   /* contract step */
      { reject++;
        for (i=0;i<m;i++) del[i]*=0.5;
        if (reject==8) for (i=0;i<m;i++) del[i]=0.0;
        if (reject==9)
        reject=0;
        if (!reject&&iter>3) ok=0;
      }
      if (op) printf("\n%12.6g  %12.6g",v,vmin);
    }
    /* get choleski decomposition of (I*rho+T) */
    for (i=0;i<T.r;i++) T.M[i][i] += rho;
    tricholeski(&T,&l0,&l1);
    for (i=0;i<T.r;i++) T.M[i][i] -= rho;
    if ((!iter&&autoinit)||!ok) /* do second update guess at start parameters */
    { /* get the current best parameter vector ZLU(I*r+T)^{-1}z */
      c.r=z.r;
      bicholeskisolve(&c,&z,&l0,&l1);
      OrthoMult(&U,&c,1,T.r-2,0,1,0);
      for (i=0;i<c.r;i++)
      { x=0.0;for (j=i;j<c.r;j++) x+=c.V[j]*R.M[i][j];c.V[i]=x*rho;}
      if (ok) /* then it's the second parameter guess */
      { for (i=0;i<c.r;i++) p->V[i]=c.V[i]; /* use p_z not Zp_z */
        c.r=np;
        for (l=0;l<m;l++)
        { for (i=0;i<LZrS[l].c;i++)
          { c.V[i]=0.0;for (j=0;j<LZrS[l].r;j++) c.V[i]+=p->V[j]*LZrS[l].M[j][i];
          }
          tr=0.0;for (i=0;i<LZrS[l].c;i++) tr+=c.V[i]*c.V[i];
          trial[l]=log(exp(eta[l])*exp(eta[l])*n*tr);
          del[l]=trial[l]-eta[l];
        }
        /* now estimate effective -ve and +ve infinities */
        x=0.0;for (l=0;l<m;l++) x+=trial[l];x/=m;
        for (l=0;l<m;l++) ninf[l]=trial[l]-x-300.0;
        for (l=0;l<m;l++) pinf[l]=trial[l]-x+300.0;

      }  else /* it's the final calculation */
      { if (Z->r)
        { p->r=np;for (i=0;i<nz;i++) p->V[i]=c.V[i];
          for (i=nz;i<np;i++) p->V[i]=0.0;
          HQmult(*p,*Z,1,0);  /* Q [c,0]'= [Z,Y][c,0]' = Zc */
         /* p->r=Z->r;matmult(*p,*Z,c,0,0); FZ */
        } else
        { p->r=np;for (i=0;i<np;i++) p->V[i]=c.V[i];}
        for (l=0;l<m;l++) theta[l]=exp(eta[l])/rho; /* return smoothing parameters */
      }
    } else /* full  Newton update */
    { /* form U'L'Z'S_i^{0.5} and U'L'Z'S_iZLU - This is the expensive part -
         <= O(n^3) worth further optimization of address calculation later.... */
      for (l=0;l<m;l++)
      { mcopy(LZrS+l,ULZrS+l);
        OrthoMult(&U,&ULZrS[l],1,T.r-2,1,1,0);
        MM=ULZrS[l].M;M1M=ULZSZLU[l].M;
        for (i=0;i<ULZrS[l].r;i++) for (j=0;j<=i;j++)
        { x=0.0;pp=MM[i];pp1=MM[j];
          for (k=0;k<ULZrS[l].c;k++) {x+= *pp * *pp1;pp++;pp1++;}
          M1M[i][j]=M1M[j][i]=x;
        }
      }
      /* now calculate the various traces needed in the calculation.
         These will be stored in trA, trdA[], trd2A[][]. These calculations
         are up to O(n^2)*/
      trA=triTrInvLL(&l0,&l1)*rho;
      for (l=0;l<m;l++)
      { A.r=ULZrS[l].r;A.c=ULZrS[l].c;
        bicholeskisolve(&A,ULZrS+l,&l0,&l1);
        trdA[l]=0.0;
        for (i=0;i<A.r;i++) for (j=0;j<A.c;j++)
        trdA[l] -= A.M[i][j]*A.M[i][j];
        trdA[l]*=rho;
      }
      /* first form (Ir+T)^{-1}U'L'Z'S_iZLU(Ir+T)^{-0.5} for all i */
      for (l=0;l<m;l++)
      { A.r=ULZSZLU[l].r;A.c=ULZSZLU[l].c;
        bicholeskisolve(&A,ULZSZLU+l,&l0,&l1);
        for (i=0;i<A.r;i++) H[l].M[i][0]=A.M[i][0]/l0.V[0];
        for (j=1;j<A.c;j++) for (i=0;i<A.r;i++)
        H[l].M[i][j]=(A.M[i][j]-H[l].M[i][j-1]*l1.V[j-1])/l0.V[j];
        /* H algorithm checked directly - ok */
      }
      for (l=0;l<m;l++) for (k=0;k<=l;k++)
      { x=0.0;
        for (i=0;i<H[l].r;i++) for (j=0;j<H[k].c;j++)
        x+=H[l].M[i][j]*H[k].M[i][j];
        trd2A[l][k]=trd2A[k][l]=2*x*rho;
        if (l==k) trd2A[l][k]+=trdA[l]/exp(eta[l]);
      }
      /* Now form b, db[], d2b[][] ...... */
      b=1-trA/n;b=b*b;
      for (l=0;l<m;l++) db[l]= exp(eta[l])*2.0/n*(trA/n-1)*trdA[l];
      for (l=0;l<m;l++) for (k=0;k<=l;k++)
      d2b[k][l]=d2b[l][k]=
      exp(eta[l]+eta[k])*2.0/n*(trdA[l]*trdA[k]/n-(1-trA/n)*trd2A[l][k]);
      /* Need the derivatives of the rss term next. Use W^{0.5}y in place of y.
         Form and store vector Ay = Q'[I,0]'U(I*r+T)^{-1}U'[I,0]Qy */
      d.r=n;for (i=0;i<n;i++) d.V[i]=Wy.V[i];  /* W^{0.5}y */
      OrthoMult(&Q,&d,0,Q.r,0,1,1);
      Ay.r=d.r=T.r;
      OrthoMult(&U,&d,1,T.r-2,1,1,0);
      bicholeskisolve(&Ay,&d,&l0,&l1); /* Ay = (I*r+T)^{-1}U'[I,0]Qy */
      /* now for the dAy[] and dpAy[] terms, before finishing off Ay :
         dpAy[i] = (I*r+T)^{-1}U'L'Z'S_iZLU(I*r+T)^{-1}U'[I,0]Qy
          dAy[i] = -Q'[I,0]'U dpAy[i]
         the dpAy[] terms save time in calculating d^2A/dt_idt_j
      */
      for (l=0;l<m;l++)
      { A.r=T.r;A.c=1L;
        matmult(A,ULZSZLU[l],Ay,0,0);
        bicholeskisolve(dpAy+l,&A,&l0,&l1);
        for (i=0;i<dpAy[l].r;i++) dAy[l].V[i]= -dpAy[l].V[i];
        for (i=dpAy[l].r;i<dAy[l].r;i++) dAy[l].V[i]=0.0;
        dAy[l].r= dpAy[l].r;
        OrthoMult(&U,dAy+l,1,T.r-2,0,1,0);
        dAy[l].r=y->r;
        OrthoMult(&Q,dAy+l,0,Q.r,1,1,1);
        for (i=0;i<dAy[l].r;i++) dAy[l].V[i]*=rho;
      }
      /* now finish Ay.... */
      OrthoMult(&U,&Ay,1,T.r-2,0,1,0);
      Ay.r=y->r;for (i=T.r;i<y->r;i++) Ay.V[i]=0.0;
      OrthoMult(&Q,&Ay,0,Q.r,1,1,1);
      for (i=0;i<Ay.r;i++) Ay.V[i]*=rho;
      /* form a, da[] & d2a[][]..... */
      for (l=0;l<m;l++) for (k=0;k<=l;k++)  // starting with d2a[][]...
      { matmult(A,ULZSZLU[k],dpAy[l],0,0); // forming (A=) U'L'Z'S_kZLU(I*r+T)^{-1}U'L'Z'S_lZLU(I*r+T)^{-1}U'[I,0]Qy
        c.r=T.r;
        bicholeskisolve(&c,&A,&l0,&l1);  // forming (c=) (I*r+T)^{-1}U'L'Z'S_kZLU(I*r+T)^{-1}U'L'Z'S_lZLU(I*r+T)^{-1}U'[I,0]Qy
        matmult(A,ULZSZLU[l],dpAy[k],0,0);  // This line and next 3 are a bug
        d.r=T.r;                            // fix for incorrect original
        bicholeskisolve(&d,&A,&l0,&l1);     // derivation - 18/8/99
        for (i=0;i<c.r;i++) c.V[i]+=d.V[i];
        OrthoMult(&U,&c,1,T.r-2,0,1,0);  // forming (c=) U(I*r+T)^{-1}U'L'Z'S_kZLU(I*r+T)^{-1}U'L'Z'S_lZLU(I*r+T)^{-1}U'[I,0]Qy
        c.r=y->r;
        for (i=T.r;i<y->r;i++) c.V[i]=0.0; // and premutiplying result by (0',I')'
        OrthoMult(&Q,&c,0,Q.r,1,1,1);      // premultiplying by Q'
        for (i=0;i<c.r;i++) c.V[i]*=rho; /* c=d^2A/dt_idt_j y */
        if (l==k) /* then operator needs additional term dA/dt_i y / e^eta_i*/
        { for (i=0;i<c.r;i++)
          c.V[i]+=dAy[l].V[i]/exp(eta[l]);
        }

        x=0.0;
        for (i=0;i<n;i++)
        x+= dAy[l].V[i]*dAy[k].V[i]+(Ay.V[i]-Wy.V[i])*c.V[i];
        x*=2.0/n;d2a[l][k]=d2a[k][l]=exp(eta[l]+eta[k])*x;
      }  // form da[]
      for (l=0;l<m;l++)
      { x=0.0;
        for (i=0;i<n;i++) x+=(Ay.V[i]-Wy.V[i])*dAy[l].V[i];
        x*=2.0/n;
        da[l]=x*exp(eta[l]);
      }
      a=0.0;
      for (i=0;i<n;i++) a+=Wy.V[i]*Wy.V[i]+(Ay.V[i]-2*Wy.V[i])*Ay.V[i];
      a/=n;
      /* with luck and a fair wind we now have the ingredients for the gradient
         and Hessian */
      for (i=0;i<m;i++)
      { if (ubre)
        { g.V[i]=da[i]- *sig2/sqrt(b)*db[i];
          for (j=0;j<=i;j++)
          Hess.M[j][i]=Hess.M[i][j]=d2a[i][j]- *sig2/sqrt(b)*d2b[i][j] +
                                    *sig2/(2*sqrt(b)*b)*db[i]*db[j];
        } else /* it's GCV */
        { g.V[i]=da[i]/b - a*db[i]/(b*b);
          for (j=0;j<=i;j++)
          Hess.M[j][i]=Hess.M[i][j]=d2a[i][j]/b-(da[i]*db[j]+da[j]*db[i])/(b*b)
                              +2*a*db[i]*db[j]/(b*b*b)-a*d2b[i][j]/(b*b);
        }
      }
      // DEBUGGING CODE Checking Hessian and other 2nd derivative results
      //boringHg(R,Q,LZSZL,y,rw,trial,rho,m,ubre*(*sig2),1e-3);
      if (op)
      { printf("\n");
        for (i=0;i<m;i++)
        { for (j=0;j<=i;j++)
          printf("%8.4g  ",Hess.M[i][j]);
      //    printf("%8.4g  ",d2a[i][j]);
          printf("\n");
        }
        for (i=0;i<m;i++)
        printf("\n%g",g.V[i]);
        //printf("\n%g",da[i]);
      }
      /* and finally the update ........ */
      A.c=A.r=Hess.r;
      x=0.0;for (i=0;i<m;i++) x+=Hess.M[i][i];x/=m;x*=0.0001;
      x=fabs(x);trials=0;
      while(!chol(Hess,A,0,0)&&(trials<500))
      { for (i=0;i<m;i++) Hess.M[i][i]+=x;
        x*=2.0;trials++;
      }
      if (trials==500)
      ok=0;
      else
      { c.r=g.r;
        choleskisolve(A,c,g);
        for (i=0;i<m;i++) del[i]= -c.V[i];
      }
    }
    iter++;
  }
  freemat(A);freemat(c);freemat(Ay);freemat(d);
  freemat(Wy);freemat(Q);freemat(R);
  freemat(T);freemat(U);freemat(z);freemat(l0);
  freemat(l1);freemat(g);freemat(Hess);
  for (i=0;i<m;i++)
  { freemat(LZrS[i]);freemat(LZSZL[i]);
    freemat(H[i]);
    //freemat(ULZrS[i]); not freed - memory shared with H[i]
    freemat(ULZSZLU[i]);
    freemat(dAy[i]);freemat(dpAy[i]);
    free(trd2A[i]);free(d2b[i]);free(d2a[i]);
  }
  free(LZrS);free(LZSZL);free(H);free(ULZrS);free(ULZSZLU);
  free(dAy);free(dpAy);free(ninf);free(pinf);
  free(trd2A);free(d2b);free(d2a);free(trdA);free(db);free(da);
  free(rw);free(eta);free(del);free(trial);
}


void MSmooth(double ft(int,int,int,double*,double*,int,int,int),
             matrix *y,matrix *J,matrix *Z,matrix *w,matrix *S,matrix *p,
             double *theta,long *off,int m,int mp,double *sig2,int transform)

/* This routine is a generalization of multismooth(). The generalization is that
   the smoothing parameters applying to each penalty (eta) can be functions of a
   smaller set of smoothing parameters (lam).
   m is the number of penalties and mp is the number of smoothing parameters.
   ft(int e,int m,int mp, double *eta, double *lam,int i,int j,int k)
        is the function defining the transformation from log true smoothing
        parameters to log parameters multiplying each penalty. Note the log
        transformation. e is a control code: -1 to transform initial eta
        guesses to lam guesses; 0 to transform lam to eta; 1 first partial
        derivative of eta_i w.r.t. lam_j; 2 for 2nd partial of eta_i wrt
        lam_j and lam_k; -2 is used to transform lam so as to achieve a
        uniform reduction in the etas - i.e. the transform of lam required to
        reduce all the eta values by a fixed amount. Works by passing in the
        transformed eta values and the original lam values, passing out
        transformed lam.  

   To reproduce Multismooth ft() would apply identity transforms, return 1
   for 1st partials and zero for second.

   routine for multiple smoothing parameter problems, with influence matrix:
   A=r*JZ(Z'J'WJZ*r + \sum_{i=1}^m \theta_i Z'S_iZ)^{-1} Z'J'W
   The routine is more or less Gu and Wahba's method (although it should
   work with alot less structure).
   The S_i matrices should be supplied in their smallest possible form:
   if p is the total parameter vector then in principle a smoothness measure
   will be of the form p'S_ip, but much of S_i will be zeroes. Suppose that S_i
   contains an s by s non-zero sub matrix S[i], and that q is the s-vector
   consisting of the s elements of p starting from element off[i]. Then the
   measure is given by q'S[i]q. It is the S[i] matrices and array off[] that
   should be supplied to the routine. This way of proceeding means that the
   smoothing parameter vector updates require only around n^3 operations not
   m n^3.
   Problem must have datapoints >= to null space dimension (parameter vector
   length, or columns of Z) - I'm not sure if this restriction is fundamental.

   If any of the m elements of theta[] are -ve then automatic initialisation
   of smoothing parameters takes place. Otherwise the supplied theta[i]s are
   used as initial estimates. On exit theta[i] contains best estimates of
   theta[i]/rho.

   NOTE that for some problems it may make sense to initialise by setting
   running first with very low s.p.s and then using the Wahba suggested 2nd
   guesses - not implemented yet.

   Z is actually supplied as a matrix containing a series of householder
   transformations, as produced by QT when fullQ==0... this improves efficiency
   considerably.... old, full Z code in brackets labelled FZ.
   To obtain Z for Cp=0: call QT(Q,C,0). Then set Q.r=C.r; Q is now used as Z.

   If *sig2>0.0 on entry then UBRE replaces GCV as the s.p. selection criterion.
*/

{ double *rw,tr,*eta,*del,*trial,trA,a,b,*trdA,**trd2A,*db,**d2b,*da,**d2a,
          rho,v,vmin=0.0,x,**RM,*pp,**MM,**M1M,**M2M,*pp1,tdf,xx1,xx2,tol,
          *pinf,*ninf,*lam;
  long i,j,k,l,n,np,nz;
  int iter,reject,ok=1,autoinit=0,op=0,ubre=0;
  matrix z,l0,l1,Q,R,C,ZC,*LZrS,*LZSZL,T,U,A,Wy,Hess,g,c,*H,*ULZSZLU,
         *ULZrS,*dAy,*dpAy,d,Ay,Ht;
  if (*sig2>0.0) ubre=1; /* signals UBRE rather than GCV */
  n=y->r;  /* number of datapoints */
  np=J->c; /* number of parameters */
  for (i=0;i<mp;i++) if (theta[i]<=0.0) autoinit=1; //** m ->mp
  if (Z->r) /*nz=Z->c FZ */ nz=np-Z->r;else nz=np; /* dimension of null space */
  A=initmat(np,np); /* workspace matrix */
  c=initmat(n,1L); /*     "     vector  */
  d=initmat(n,1L);
  Ay=initmat(n,1L);
  Hess=initmat((long)m,(long)m);g=initmat((long)m,1L);
  /*for (l=0;l<m;l++) ... appears to be totally redundant matrix allocation */
  { trdA=(double *)calloc((size_t)m,sizeof(double));
    da=(double *)calloc((size_t)m,sizeof(double));
    db=(double *)calloc((size_t)m,sizeof(double));
    trd2A=(double **)calloc((size_t)m,sizeof(double));
    d2a=(double **)calloc((size_t)m,sizeof(double));
    d2b=(double **)calloc((size_t)m,sizeof(double));
    ninf=(double *)calloc((size_t)m,sizeof(double));
    pinf=(double *)calloc((size_t)m,sizeof(double));
    for (k=0;k<m;k++)
    { trd2A[k]=(double *)calloc((size_t)m,sizeof(double));
      d2a[k]=(double *)calloc((size_t)m,sizeof(double));
      d2b[k]=(double *)calloc((size_t)m,sizeof(double));
    }
  }
  /* Get Q'R=W^{0.5} J Z */
  rw=(double *)calloc((size_t)n,sizeof(double));
  for (i=0;i<n;i++) rw[i]=sqrt(w->V[i]);   /* rw contains l.d. of W^{0.5} */
  Wy=initmat(n,1L);
  for (i=0;i<n;i++) Wy.V[i]=rw[i]*y->V[i];
  if (Z->r)
  { R=initmat(n,np);mcopy(J,&R);
    HQmult(R,*Z,0,0);R.c=nz;
    /* matmult(R,*J,*Z,0,0); FZ */ /* R=JZ - up to O(n^3) */
  } else
  { R=initmat(n,np);mcopy(J,&R);}          /* case when Z=I */
  RM=R.M;
  for (i=0;i<n;i++)
  { x=rw[i];for (pp=RM[i];pp<RM[i]+R.c;pp++) *pp *= x;}  /* R=W^{0.5} JZ */
  Q=initmat(np,n);/* altered for efficient storage */
  QR(&Q,&R);  /* done in up to O(n^3) */
  /* invert R - it is this step that requires n>=nz*/
  R.r=R.c; /* getting rid of zero part */
  InvertTriangular(&R);   /* R now contains L */
  /* Form the matrices L'Z'S^{0.5} and L'Z'S_iZL */
  ZC=initmat(np,np);
  LZrS=(matrix *)calloc((size_t)m,sizeof(matrix));
  LZSZL=(matrix *)calloc((size_t)m,sizeof(matrix));
  ULZSZLU=(matrix *)calloc((size_t)m,sizeof(matrix));
  H= (matrix *)calloc((size_t)m,sizeof(matrix));
  ULZrS=(matrix *)calloc((size_t)m,sizeof(matrix));
  for (l=0;l<m;l++)
  { root(S+l,&C,1e-14);    /* S[l]=CC' */
    if (Z->r)
    { ZC.c=C.c;ZC.r=np;
      /* set ZC = [0,C',0]' */
      for (i=0;i<off[l];i++) for (j=0;j<C.c;j++) ZC.M[i][j]=0.0;
      for (i=off[l];i<off[l]+C.r;i++)
      for (j=0;j<C.c;j++) ZC.M[i][j]=C.M[i-off[l]][j];
      for (i=off[l]+C.r;i<np;i++)
      for (j=0;j<C.c;j++) ZC.M[i][j]=0.0;
      /* ...and apply Z efficiently... */
      HQmult(ZC,*Z,1,1);ZC.r=nz; /* bug here fixed 17/8/99 - should be Z'C but was attempting ZC*/ 
    } else
    { ZC.c=C.c;ZC.r=np;
      for (i=0;i<ZC.r;i++) for (j=0;j<ZC.c;j++) ZC.M[i][j]=0.0;
      for (i=0;i<C.r;i++) for (j=0;j<C.c;j++) ZC.M[i+off[l]][j]=C.M[i][j];
    }
    freemat(C);
    LZrS[l]=initmat(R.c,ZC.c);
    MM=LZrS[l].M;M1M=R.M;M2M=ZC.M;
    for (i=0;i<LZrS[l].r;i++) for (j=0;j<LZrS[l].c;j++)
    for (k=0;k<=i;k++) MM[i][j]+=M1M[k][i]*M2M[k][j];
    LZSZL[l]=initmat(R.c,R.c); // this memory requirement could be halved
    matmult(LZSZL[l],LZrS[l],LZrS[l],0,1);
    H[l]=initmat(R.c,R.c);
    ULZSZLU[l]=initmat(R.c,R.c); // memory requirement could be halved, but would need own choleskisolve
    //ULZrS[l]=initmat(R.c,ZC.c);
    ULZrS[l].M=H[l].M;ULZrS[l].r=R.c;ULZrS[l].c=ZC.c; // sharing memory with H[l]
  }
  /* Start the main loop */
  freemat(ZC);
  eta=(double *)calloc((size_t)m,sizeof(double));
  lam=(double *)calloc((size_t)mp,sizeof(double)); //** added
  del=(double *)calloc((size_t)mp,sizeof(double)); /* change in s.p.s */
  trial=(double *)calloc((size_t)mp,sizeof(double));
  /* get initial estimates for theta_i and eta_i=log(theta_i) */
  if (autoinit)
  for (l=0;l<m;l++)
  { tr=0.0;for (i=0;i<LZSZL[l].r;i++) tr+=LZSZL[l].M[i][i];
    eta[l]=log(1.0/(n*tr));//eta[l]= -40.0;
  } else
  { for (i=0;i<mp;i++) theta[i]=log(theta[i]);
    if (transform) ft(0,m,mp,eta,theta,0,0,0);else
    for (i=0;i<m;i++) eta[i]=theta[i];
    x=0.0;for (i=0;i<m;i++) { x+=eta[i];}x/=m;
    for (i=0;i<m;i++) { ninf[i]=eta[i]-x-300.0;pinf[i]=ninf[i]+600.0;}
  }
  if (transform)               //**
  { ft(-1,m,mp,eta,lam,0,0,0); // get initial lam estimates
    ft(0,m,mp,eta,lam,0,0,0);  // transform these back to initial eta estimates
  } else
  for (i=0;i<m;i++) lam[i]=eta[i];
  T=initmat(R.c,R.c);
  U=initmat(T.r,T.r);
  z=initmat(n,1L);
  l0=initmat(T.r,1L);l1=initmat(T.r-1,1L);
  dAy=(matrix *)calloc((size_t)m,sizeof(matrix));
  dpAy=(matrix *)calloc((size_t)m,sizeof(matrix));
  for (l=0;l<m;l++) { dAy[l]=initmat(n,1L);dpAy[l]=initmat(T.r,1L);}
  iter=0;  /* iteration counter */
  while (ok)
  { /* form combined smoothness measure */
    reject=1;
    while (reject) /* try current search direction (unless first step) */
    { for (i=0;i<mp;i++) { trial[i]=lam[i]+del[i];x+=trial[i];}
      if (transform) ft(0,m,mp,eta,trial,0,0,0); else
      for (i=0;i<m;i++) eta[i]=trial[i];
      x=0.0;for (i=0;i<m;i++) x+=eta[i];x/=m;
      for (i=0;i<m;i++) eta[i] -= x; /* normalising smooths */
      if (transform) ft(-2,m,mp,eta,trial,0,0,0); // making trial consistent with eta
      else for (i=0;i<m;i++) trial[i] -= x;
      if (iter>1||!autoinit)   /* check smoothing parameters won't lead to overflow */
      for (i=0;i<m;i++)
      if (eta[i]<ninf[i])
      eta[i]=ninf[i];
      else if (eta[i]>pinf[i])
      eta[i]=pinf[i];
      /* form S the combined smooth measure */
      for (i=0;i<T.r;i++) for (j=0;j<T.c;j++)
      T.M[i][j]=exp(eta[0])*LZSZL[0].M[i][j];
      for (l=1;l<m;l++) for (i=0;i<T.r;i++) for (j=0;j<T.c;j++)
      T.M[i][j] += exp(eta[l])*LZSZL[l].M[i][j];
      UTU(&T,&U);    /* Form S=UTU' */
      z.r=n;
      for (i=0;i<n;i++) z.V[i]=rw[i]*y->V[i]; /* z=W^{1/2}y */
     // matrixintegritycheck();
      OrthoMult(&Q,&z,0,Q.r,0,1,1);           /* z=QW^{1/2}y */
      z.r=R.r;                                /* z=[I,0]QW^{1/2}y */
      OrthoMult(&U,&z,1,T.r-2,1,1,0);         /* z=U'[I,0]QW^{1/2}y */
      z.r=n;                                  /* z->[z,x_1]' */
      if (!ubre) *sig2=-1.0; /* signalling use of GCV */
      rho=EasySmooth(&T,&z,&v,&tdf,n,sig2);    /* do a cheap minimisation in rho */
      z.r=R.r;
      if (!iter||v<vmin) /* accept current step */
      { reject=0;
        /* test for convergence */
        tol=1e-4;ok=0;
        if (vmin-v>tol*(1+v)) ok=1;
        xx1=0.0;for (i=0;i<mp;i++) { xx2=lam[i]-trial[i];xx1+=xx2*xx2;}
        xx1=sqrt(xx1);
        xx2=0.0;for (i=0;i<mp;i++) xx2+=trial[i]*trial[i];
        xx2=sqrt(xx2);xx2=(1+xx2)*sqrt(tol);
        if (xx1>xx2) ok=1;
        xx1=0.0;for (i=0;i<mp;i++) xx1+=g.V[i]*g.V[i];xx1=sqrt(xx1);
        if (xx1>pow(tol,1.0/3)*(1+v)) ok=1;
        for (i=0;i<mp;i++) lam[i]=trial[i];
        vmin=v;
      } else   /* contract step */
      { reject++;
        for (i=0;i<mp;i++) del[i]*=0.5;
        if (reject==8) for (i=0;i<mp;i++) del[i]=0.0;
        if (reject==9) reject=0;
        if (!reject&&iter>3) ok=0;
      }
      if (op) printf("\n%12.6g  %12.6g",v,vmin);
    }
    /* get choleski decomposition of (I*rho+T) */
    for (i=0;i<T.r;i++) T.M[i][i] += rho;
    tricholeski(&T,&l0,&l1);
    for (i=0;i<T.r;i++) T.M[i][i] -= rho;
    if ((!iter&&autoinit)||!ok) /* do second update guess at start parameters */
    { /* get the current best parameter vector ZLU(I*r+T)^{-1}z */
      c.r=z.r;
      bicholeskisolve(&c,&z,&l0,&l1);
      OrthoMult(&U,&c,1,T.r-2,0,1,0);
      for (i=0;i<c.r;i++)
      { x=0.0;for (j=i;j<c.r;j++) x+=c.V[j]*R.M[i][j];c.V[i]=x*rho;}
      if (ok) /* then it's the second parameter guess */
      { for (i=0;i<c.r;i++) p->V[i]=c.V[i]; /* use p_z not Zp_z */
        c.r=np;
        for (l=0;l<m;l++)
        { for (i=0;i<LZrS[l].c;i++)
          { c.V[i]=0.0;for (j=0;j<LZrS[l].r;j++) c.V[i]+=p->V[j]*LZrS[l].M[j][i];
          }
          tr=0.0;for (i=0;i<LZrS[l].c;i++) tr+=c.V[i]*c.V[i];
          eta[l]=log(exp(eta[l])*exp(eta[l])*n*tr);
        }
        if (transform)
        { ft(-1,m,mp,eta,trial,0,0,0);
          ft(0,m,mp,eta,trial,0,0,0);
        } else
        for (i=0;i<mp;i++) trial[i]=eta[i];
        for (l=0;l<mp;l++) del[l]=trial[l]-lam[l];
        /* now estimate effective -ve and +ve infinities */
        x=0.0;for (l=0;l<m;l++) x+=eta[l];x/=m;
        for (l=0;l<m;l++) ninf[l]=eta[l]-x-300.0;
        for (l=0;l<m;l++) pinf[l]=eta[l]-x+300.0;
      }  else /* it's the final calculation */
      { if (Z->r)
        { p->r=np;for (i=0;i<nz;i++) p->V[i]=c.V[i];
          for (i=nz;i<np;i++) p->V[i]=0.0;
          HQmult(*p,*Z,1,0);  /* Q [c,0]'= [Z,Y][c,0]' = Zc */
         /* p->r=Z->r;matmult(*p,*Z,c,0,0); FZ */
        } else
        { p->r=np;for (i=0;i<np;i++) p->V[i]=c.V[i];}
        for (l=0;l<m;l++) eta[l]-=log(rho);
        if (transform) ft(-2,m,mp,eta,lam,0,0,0);
        else for (i=0;i<m;i++) lam[i]=eta[i];
        for (l=0;l<mp;l++) theta[l]=exp(lam[l]); /* return smoothing parameters */
      }
    } else /* full  Newton update */
    { /* form U'L'Z'S_i^{0.5} and U'L'Z'S_iZLU - This is the expensive part -
         <= O(n^3) worth further optimization of address calculation later.... */
      for (l=0;l<m;l++)
      { mcopy(LZrS+l,ULZrS+l);
        OrthoMult(&U,&ULZrS[l],1,T.r-2,1,1,0);
        MM=ULZrS[l].M;M1M=ULZSZLU[l].M;
        for (i=0;i<ULZrS[l].r;i++) for (j=0;j<=i;j++)
        { x=0.0;pp=MM[i];pp1=MM[j];
          for (k=0;k<ULZrS[l].c;k++) {x+= *pp * *pp1;pp++;pp1++;}
          M1M[i][j]=M1M[j][i]=x;
        }
      }
      /* now calculate the various traces needed in the calculation.
         These will be stored in trA, trdA[], trd2A[][]. These calculations
         are up to O(n^2)*/
      trA=triTrInvLL(&l0,&l1)*rho;
      for (l=0;l<m;l++)
      { A.r=ULZrS[l].r;A.c=ULZrS[l].c;
        bicholeskisolve(&A,ULZrS+l,&l0,&l1);
        trdA[l]=0.0;
        for (i=0;i<A.r;i++) for (j=0;j<A.c;j++)
        trdA[l] -= A.M[i][j]*A.M[i][j];
        trdA[l]*=rho;
      }
      /* first form (Ir+T)^{-1}U'L'Z'S_iZLU(Ir+T)^{-0.5} for all i */
      for (l=0;l<m;l++)
      { A.r=ULZSZLU[l].r;A.c=ULZSZLU[l].c;
        bicholeskisolve(&A,ULZSZLU+l,&l0,&l1);
        for (i=0;i<A.r;i++) H[l].M[i][0]=A.M[i][0]/l0.V[0];
        for (j=1;j<A.c;j++) for (i=0;i<A.r;i++)
        H[l].M[i][j]=(A.M[i][j]-H[l].M[i][j-1]*l1.V[j-1])/l0.V[j];
        /* H algorithm checked directly - ok */
      }
      for (l=0;l<m;l++) for (k=0;k<=l;k++)
      { x=0.0;
        for (i=0;i<H[l].r;i++) for (j=0;j<H[k].c;j++)
        x+=H[l].M[i][j]*H[k].M[i][j];
        trd2A[l][k]=trd2A[k][l]=2*x*rho;
        if (l==k) trd2A[l][k]+=trdA[l]/exp(eta[l]);
      }
      /* Now form b, db[], d2b[][] ...... */
      b=1-trA/n;b=b*b;
      for (l=0;l<m;l++) db[l]= exp(eta[l])*2.0/n*(trA/n-1)*trdA[l];
      for (l=0;l<m;l++) for (k=0;k<=l;k++)
      d2b[k][l]=d2b[l][k]=
      exp(eta[l]+eta[k])*2.0/n*(trdA[l]*trdA[k]/n-(1-trA/n)*trd2A[l][k]);
      /* Need the derivatives of the rss term next. Use W^{0.5}y in place of y.
         Form and store vector Ay = Q'[I,0]'U(I*r+T)^{-1}U'[I,0]Qy */
      d.r=n;for (i=0;i<n;i++) d.V[i]=Wy.V[i];  /* W^{0.5}y */
      OrthoMult(&Q,&d,0,Q.r,0,1,1);
      Ay.r=d.r=T.r;
      OrthoMult(&U,&d,1,T.r-2,1,1,0);
      bicholeskisolve(&Ay,&d,&l0,&l1); /* Ay = (I*r+T)^{-1}U'[I,0]Qy */
      /* now for the dAy[] and dpAy[] terms, before finishing off Ay :
         dpAy[i] = (I*r+T)^{-1}U'L'Z'S_iZLU(I*r+T)^{-1}U'[I,0]Qy
          dAy[i] = -Q'[I,0]'U dpAy[i]
         the dpAy[] terms save time in calculating d^2A/dt_idt_j
      */
      for (l=0;l<m;l++)
      { A.r=T.r;A.c=1L;
        matmult(A,ULZSZLU[l],Ay,0,0);
        bicholeskisolve(dpAy+l,&A,&l0,&l1);
        for (i=0;i<dpAy[l].r;i++) dAy[l].V[i]= -dpAy[l].V[i];
        for (i=dpAy[l].r;i<dAy[l].r;i++) dAy[l].V[i]=0.0;
        dAy[l].r= dpAy[l].r;
        OrthoMult(&U,dAy+l,1,T.r-2,0,1,0);
        dAy[l].r=y->r;
        OrthoMult(&Q,dAy+l,0,Q.r,1,1,1);
        for (i=0;i<dAy[l].r;i++) dAy[l].V[i]*=rho;
      }
      /* now finish Ay.... */
      OrthoMult(&U,&Ay,1,T.r-2,0,1,0);
      Ay.r=y->r;for (i=T.r;i<y->r;i++) Ay.V[i]=0.0;
      OrthoMult(&Q,&Ay,0,Q.r,1,1,1);
      for (i=0;i<Ay.r;i++) Ay.V[i]*=rho;
      /* form a, da[] & d2a[][]..... */
      for (l=0;l<m;l++) for (k=0;k<=l;k++)
      { matmult(A,ULZSZLU[k],dpAy[l],0,0);
        c.r=T.r;
        bicholeskisolve(&c,&A,&l0,&l1);
        matmult(A,ULZSZLU[l],dpAy[k],0,0);
        d.r=T.r;
        bicholeskisolve(&d,&A,&l0,&l1);
        for (i=0;i<c.r;i++) c.V[i]+=d.V[i];
        OrthoMult(&U,&c,1,T.r-2,0,1,0);
        c.r=y->r;
        for (i=T.r;i<y->r;i++) c.V[i]=0.0;
        OrthoMult(&Q,&c,0,Q.r,1,1,1);
        for (i=0;i<c.r;i++) c.V[i]*=rho; /* b=d^2A/dt_idt_j y */
        if (l==k) /* then operator needs additional term dA/dt_i y / e^eta_i*/
        { for (i=0;i<c.r;i++)
          c.V[i]+=dAy[l].V[i]/exp(eta[l]);
        }

        x=0.0;
        for (i=0;i<n;i++)
        x+= dAy[l].V[i]*dAy[k].V[i]+(Ay.V[i]-Wy.V[i])*c.V[i];
        x*=2.0/n;d2a[l][k]=d2a[k][l]=exp(eta[l]+eta[k])*x;
      }
      for (l=0;l<m;l++)
      { x=0.0;
        for (i=0;i<n;i++) x+=(Ay.V[i]-Wy.V[i])*dAy[l].V[i];
        x*=2.0/n;
        da[l]=x*exp(eta[l]);
      }
      a=0.0;
      for (i=0;i<n;i++) a+=Wy.V[i]*Wy.V[i]+(Ay.V[i]-2*Wy.V[i])*Ay.V[i];
      a/=n;
      /* with luck and a fair wind we now have the ingredients for the gradient
         and Hessian */
      Hess.r=Hess.c=g.r=(long)m;
      for (i=0;i<m;i++)
      { if (ubre)
        { g.V[i]=da[i]- *sig2/sqrt(b)*db[i];
          for (j=0;j<=i;j++)
          Hess.M[j][i]=Hess.M[i][j]=d2a[i][j]- *sig2/sqrt(b)*d2b[i][j] +
                                    *sig2/(2*sqrt(b)*b)*db[i]*db[j];
        } else /* it's GCV */
        { g.V[i]=da[i]/b - a*db[i]/(b*b);
          for (j=0;j<=i;j++)
          Hess.M[j][i]=Hess.M[i][j]=d2a[i][j]/b-(da[i]*db[j]+da[j]*db[i])/(b*b)
                              +2*a*db[i]*db[j]/(b*b*b)-a*d2b[i][j]/(b*b);
        }
      }
    /*  boringHg(R,Q,LZSZL,y,rw,trial,rho,m);*/
      if (op)
      { printf("\n");
        for (i=0;i<m;i++)
        { for (j=0;j<=i;j++) printf("%8.4g  ",Hess.M[i][j]);printf("\n");}
        for (i=0;i<m;i++) printf("\n%g",g.V[i]);
      }
      /* Now transform the Hessian and gradient if necessary */
      if (transform)
      { Ht=initmat((long)mp,(long)mp);
        for (i=0;i<mp;i++) for (j=0;j<mp;j++)
        for (k=0;k<m;k++)
        { Ht.M[i][j]+=g.V[k]*ft(2,m,mp,eta,lam,k,i,j);
          for (l=0;l<m;l++)
          Ht.M[i][j]+=Hess.M[k][l]*ft(1,m,mp,eta,lam,k,i,0)*ft(1,m,mp,eta,lam,l,j,0);
        }
        for (i=0;i<mp;i++) for (j=0;j<mp;j++) Hess.M[i][j]=Ht.M[i][j];
        freemat(Ht);
        Ht=initmat((long)mp,1L);
        for (i=0;i<mp;i++) for (k=0;k<m;k++)
        Ht.V[i]+=g.V[k]*ft(1,m,mp,eta,lam,k,i,0);
        for (i=0;i<mp;i++) g.V[i]=Ht.V[i];
        freemat(Ht);
        Hess.r=Hess.c=g.r=(long)mp;
      }

      /* and finally the update ........ */

      A.c=A.r=Hess.r;
      x=0.0;for (i=0;i<mp;i++) x+=Hess.M[i][i];x/=m;x*=0.0001;
      x=fabs(x);
      while(!chol(Hess,A,0,0))
      { for (i=0;i<mp;i++) Hess.M[i][i]+=x;
        x*=2.0;
      }
      c.r=g.r;
      choleskisolve(A,c,g);
      for (i=0;i<mp;i++) del[i]= -c.V[i];
    }
    iter++;
  }
  freemat(A);freemat(c);freemat(Ay);freemat(d);
  freemat(Wy);freemat(Q);freemat(R);
  freemat(T);freemat(U);freemat(z);freemat(l0);
  freemat(l1);freemat(g);freemat(Hess);
  for (i=0;i<m;i++)
  { freemat(LZrS[i]);freemat(LZSZL[i]);
    freemat(H[i]);
    //freemat(ULZrS[i]); not freed - memory shared with H[i]
    freemat(ULZSZLU[i]);
    freemat(dAy[i]);freemat(dpAy[i]);
    free(trd2A[i]);free(d2b[i]);free(d2a[i]);
  }
  free(LZrS);free(LZSZL);free(H);free(ULZrS);free(ULZSZLU);
  free(dAy);free(dpAy);free(ninf);free(pinf);
  free(trd2A);free(d2b);free(d2a);free(trdA);free(db);free(da);
  free(rw);
  free(eta);
  free(del);
  free(lam);
  free(trial);
}

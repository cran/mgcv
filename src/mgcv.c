/* Source code for mgcv.dll multiple smoothing parameter estimation code,
suitable for interfacing to R or S-PLUS

Copyright (C) 2000 Simon N. Wood  snw@st-and.ac.uk

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */


//#include <windows.h>   // useful to turn this on and alter infobox for Windows debugging
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gcv.h"
#include "mgcv.h"
#include "matrix.h"
#include <R.h>
#define round(a) ((a)-floor(a) <0.5 ? (int)floor(a):(int) floor(a)+1)


void ErrorMessage(char *msg,int fatal)

{ //MessageBox(HWND_DESKTOP,msg,"Info!",MB_ICONEXCLAMATION|MB_OK);
  if (fatal) error("%s",msg);
  else warning("%s",msg);
}

void infobox(char *msg)

{ //MessageBox(HWND_DESKTOP,msg,"Info!",MB_ICONEXCLAMATION|MB_OK);
  warning("%s",msg);
}

/* The following are some rather ancient routines used to set up an example
   additive model using regression (cubic) splines, via RGAMsetup(). */


/********** The following are from spline.c (rather plodding coding) ***********/

/* The next 4 functions are basis functions for the 1st derivative
   representation of a cubic spline */


double b0(x0,x1,x) double x0,x1,x;

/* multiplies function value at x0 */

{ double res,h,xx1;
  h=x1-x0;xx1=x-x1;
  res=2.0*(x-x0+0.5*h)*xx1*xx1/(h*h*h);
  return(res);
}

double b1(x0,x1,x) double x0,x1,x;

/* multiplies function value at x1 */

{ double res,h,xx0;
  h=x1-x0;xx0=x-x0;
  res= -2.0*(x-x1-0.5*h)*xx0*xx0/(h*h*h);
  return(res);
}

double d0(x0,x1,x) double x0,x1,x;

/* multiplies gradient at x0 */

{ double res,h,xx1;
  h=x1-x0;xx1=x-x1;
  res=(x-x0)*xx1*xx1/(h*h);
  return(res);
}

double d1(x0,x1,x) double x0,x1,x;

/* multiplies gradient at x1 */

{ double res,h,xx0;
  h=x1-x0;xx0=x-x0;
  res=xx0*xx0*(x-x1)/(h*h);
  return(res);
}


/* The next 4 functions are derivatives of the spline basis functions used
   above. */

double db0(x0,x1,x) double x0,x1,x;

{ double res,h,xx1;
  h=x1-x0;xx1=x-x1;
  res=2.0*(2.0*(x-x0+0.5*h)*xx1+xx1*xx1)/(h*h*h);
  return(res);
}
double db1(x0,x1,x) double x0,x1,x;

{ double res,h,xx0;
  h=x1-x0;xx0=x-x0;
  res= -2.0*(2.0*(x-x1-0.5*h)*xx0+xx0*xx0)/(h*h*h);
  return(res);
}

double dd0(x0,x1,x) double x0,x1,x;

{ double res,h,xx1;
  h=x1-x0;xx1=x-x1;
  res=(xx1*xx1+2.0*(x-x0)*xx1)/(h*h);
  return(res);
}

double dd1(x0,x1,x) double x0,x1,x;

{ double res,h,xx0;
  h=x1-x0;xx0=x-x0;
  res=(xx0*xx0+2.0*xx0*(x-x1))/(h*h);
  return(res);
}

matrix getD(h,nak) matrix h;int nak;

/* the matrix mapping the value of the spline to the gradients at the knots.
   nak is true for 'not-a-knot' end conditions at the early end, otherwise
   'natural' end conditions are used. If there are only 2 knots then the spline
   is taken as a straight line if only 1 a constant. */

{ long i,j,n;
  matrix T,D,Res;
  n=h.r+1;
  T=initmat(n,n);D=initmat(n,n);Res=initmat(n,n);
  for (i=0;i<n;i++) for (j=0;j<n;j++)
  { T.M[i][j]=0.0;D.M[i][j]=0.0;}
  if (n==1L)
  { Res.M[0][0]=0.0;
  } else
  if (n==2L)
  { Res.M[0][0]=Res.M[1][0]=-1.0/h.V[0];
    Res.M[0][1]=Res.M[1][1]=1.0/h.V[0];
  } else
  { for (i=0;i<n;i++) T.M[i][i]=2.0;
    for (i=1;i<n-1;i++)
    { T.M[i][i-1]=h.V[i]/(h.V[i]+h.V[i-1]);
      T.M[i][i+1]=1.0-T.M[i][i-1];
      D.M[i][i-1]= -3.0*T.M[i][i-1]/h.V[i-1];
      D.M[i][i+1]=3.0*T.M[i][i+1]/h.V[i];
      D.M[i][i]= -(D.M[i][i+1]+D.M[i][i-1]);
    }
    if (!nak)
    { T.M[0][1]=1.0;D.M[0][0]= -3.0/h.V[0];D.M[0][1]= -D.M[0][0];}
    else
    { T.M[0][1]=2.0*(h.V[0]+h.V[1])/h.V[1];
      D.M[0][0]= -2.0*(3.0*h.V[0]+2.0*h.V[1])/
		(h.V[0]*(h.V[0]+h.V[1]));
      D.M[0][2]=2.0*h.V[0]*h.V[0]/
      (h.V[1]*h.V[1]*(h.V[0]+h.V[1]));
      D.M[0][1]= -D.M[0][0]-D.M[0][2];
    }
    T.M[n-1][n-2]=1.0;D.M[n-1][n-2]= -3.0/h.V[n-2];
    D.M[n-1][n-1]= -D.M[n-1][n-2];
    invert(&T);
    matmult(Res,T,D,0,0);
  }
  freemat(T);freemat(D);
  return(Res);
}

void tmap(tm,tgm,t,time,kill)
matrix tm,tgm,t;
double time;
int kill; /* to release static matrix allocation set to 1 otherwise 0 and
	     prepare for a new sequence of knot positions in t*/

/* tm maps values of a function at the t values contained in vector t to
   the value of a spline through those points at 'time' ;tgm does the same
   for the gradient of the spline */

{ static matrix D;static char first=1;
  matrix h;
  long i,k;
  if (first)
  { first=0;h=initmat(t.r-1,1L);
    for (i=0L;i<t.r-1;i++) h.V[i]=t.V[i+1]-t.V[i];
    D=getD(h,0); /* time trajectories always have natural end conditions */
    freemat(h);
  }
  if (t.r==1L)
  { tm.V[0]=1.0;tgm.V[0]=0.0;}
  else
  { i=0L;while((time>t.V[i+1])&&(i<t.r-2)) i++;
    for (k=0;k<t.r;k++)
    tm.V[k]=D.M[i][k]*d0(t.V[i],t.V[i+1],time)+
	    D.M[i+1][k]*d1(t.V[i],t.V[i+1],time);
    tm.V[i]+=b0(t.V[i],t.V[i+1],time);
    tm.V[i+1]+=b1(t.V[i],t.V[i+1],time);
    for (k=0;k<t.r;k++)
    tgm.V[k]=D.M[i][k]*dd0(t.V[i],t.V[i+1],time)+
	     D.M[i+1][k]*dd1(t.V[i],t.V[i+1],time);
    tgm.V[i]+=db0(t.V[i],t.V[i+1],time);
    tgm.V[i+1]+=db1(t.V[i],t.V[i+1],time);
  }
  if (kill)
  { first=1;
    freemat(D);
  }
}

void getHBH(HBH,h,nak,rescale) matrix *HBH,h;int nak,rescale;

/* Generates the wiggliness measure matrix for vector h; nak=0 for natural
   end conditions or nak=1 to use the not a knot condition at the lower end;
   set rescale=1 to produce a measure rescaled for the unit interval, set to
   zero otherwise */

{ long n,i,j;
  matrix C,B,BI,H,hn;
  double interval=0.0;
  n=h.r;
  if (rescale)
  { for (i=0;i<h.r;i++) interval+=h.V[i];
    hn=initmat(h.r,1L);
    for (i=0;i<h.r;i++) hn.V[i]=h.V[i]/interval;
  } else hn=h;
  (*HBH)=initmat(n+1,n+1);
  if (!nak)
  { C=initmat(n-1,n+1);
    B=initmat(n-1,n-1);
    H=initmat(n-1,n+1);
    for (i=0;i<n-1;i++)
    { for (j=0;j<n-1;j++)
      { B.M[i][j]=0.0;
	     H.M[i][j]=0.0;
      }
      H.M[i][n-1]=0.0;
      H.M[i][n]=0.0;
    }
    for (i=0;i<n-1;i++)
    { B.M[i][i]=(hn.V[i]+hn.V[i+1])/3.0;
      H.M[i][i]=1.0/hn.V[i];
      H.M[i][i+1]= -1.0/hn.V[i]-1.0/hn.V[i+1];
      H.M[i][i+2]=1.0/hn.V[i+1];
    }
    for (i=0;i<n-2;i++)
    { B.M[i][i+1]=hn.V[i+1]/6.0;
      B.M[i+1][i]=hn.V[i+1]/6.0;
    }
    invert(&B);
    matmult(C,B,H,0,0);
    matmult((*HBH),H,C,1,0);
    freemat(C);freemat(B);freemat(H);
  } else
  { H=initmat(n,n+1);
    BI=initmat(n,n);B=initmat(n,n);
    for (i=0;i<H.r;i++) for (j=0;j<H.c;j++) H.M[i][j]=0.0;
    for (i=1;i<n;i++)
    { H.M[i][i-1]=1.0/hn.V[i-1];H.M[i][i]= -1.0/hn.V[i-1]-1.0/hn.V[i];
      H.M[i][i+1]=1.0/hn.V[i];
    }
    for (i=0;i<n;i++) for (j=0;j<n;j++)
    { BI.M[i][j]=0.0;B.M[i][j]=0.0;}
    for (i=1;i<n;i++)
    { B.M[i][i-1]=hn.V[i-1]/6.0;B.M[i][i]=(hn.V[i-1]+hn.V[i])/3.0;
      if (i<(n-1))
      { B.M[i][i+1]=hn.V[i]/6.0;
	BI.M[i][i+1]=B.M[i][i+1];
      }
      for (j=0;j<2;j++) BI.M[i][j+i-1]=B.M[i][j+i-1];
    }
    B.M[0][0]= -hn.V[1];B.M[0][1]=hn.V[0]+hn.V[1];
    B.M[0][2]= -hn.V[0];
    BI.M[0][0]=hn.V[0]/3.0;BI.M[0][1]=hn.V[0]/6.0;
    C=initmat(n,n);
    invert(&B);
    matmult(C,BI,B,0,0);
    matmult(BI,B,C,1,0);
    freemat(B);freemat(C);
    C=initmat(n,n+1);
    matmult(C,BI,H,0,0);
    matmult((*HBH),H,C,1,0);
    freemat(C);freemat(BI);freemat(H);
  }
  if (rescale) freemat(hn);
}


void getSmooth(S,x,rescale) matrix *S,x;int rescale;

/* gets a natural wiggliness measure for a spline with knots at the elements
   of vector x. Set rescale not zero to pretend that the domain is the unit
   interval. */

{ matrix h;
  long i;
  h=initmat(x.r-1L,1L);
  for (i=0;i<x.r-1;i++) h.V[i]=x.V[i+1]-x.V[i];
  getHBH(S,h,0,rescale);
  freemat(h);
}






/***************** The following is from glamp.c **************************************/

void GAMsetup(matrix *X,matrix *Z,matrix *S,matrix *xp,long *off,double **x,
              int m,int n,long *df,int nsdf,int getZ)

/* Sets up the Design matrix and Smoothness constraints for a regression spline
   based GAM. There are m smoothers, n datapoints. The ith smoother has
   df[i] parameters starting at off[i], with smoothing constraint matrix
   S[i].
   x[i][j], contains the jth value of the ith covariate.
   The first nsdf x[i]s are columns of the design matrix for the non-spline
   part of the model (including the intercept term) .
   X & Z are initialised in the routine.
   Each S[i] is also initialised in the routine.
   xp[i] is returned containing the vector of knot positions for the ith spline.
   Each xp[i] is initialised within the routine.

   Z contains the null space for the GAM side constraints if getZ is not set to
     0 - if getZ==0 then Z is returned containing the constraint matrix itself.

   All smooths are constrained so that their parameters sum to zero - the first
   parameter in the model will then be a mean. Hence for a pure gam() x[1] should
   be a vector of 1's and nsdf should be equal to one, to allow for a constant.


*/

{ long np,anp,l,i,j,lp,k;
  double xx,dx;
  matrix my,mg,T,y;
//  char *msg[200];
  np=nsdf+df[0];for (i=1;i<m;i++) np+=df[i];
  T=initmat((long)m,np);
  *X=initmat((long)n,(long)np);
  off[0]=nsdf;
  y=initmat((long)n,1L);
  for (j=0;j<nsdf;j++) for (i=0;i<n;i++) X->M[i][j]=x[j][i];
  for (l=0;l<m;l++) /* work through the smooths */
  { if (l) off[l]=off[l-1]+df[l-1];else off[l]=nsdf;
    /* set up part of design matrix for this smooth */
    /* get knot sequence */
    lp=l+nsdf;
   
    // arrange knots by spreading them roughly every dx x values
    for (j=0;j<n;j++) y.V[j]=x[lp][j];
    y.r=(long)n;
    sort(y); // next reduce to list of unique values.....
    k=0;for (i=0;i<n;i++) if (y.V[k]!=y.V[i]) { k++;y.V[k]=y.V[i];} y.r=(long)k+1;
    dx=(y.r-1)/(df[l]-1.0);
    xp[l]=initmat((long)df[l],1L);   /* knot vector */
   
    xp[l].V[0]=y.V[0];
    for (i=1;i<df[l]-1;i++)  // place knots
    { xx=dx*i;
      k=(int)floor(xx);
      //if (xx-k>0.5) k++;
      xx -= k;
      xp[l].V[i]=(1-xx)*y.V[k]+xx*y.V[k+1];
    } 
    xp[l].V[xp[l].r-1]=y.V[y.r-1];

    getSmooth(S+l,xp[l],0);
    my=initmat(xp[l].r,1L);mg=initmat(xp[l].r,1L);
   /* tmap(my,mg,xp[l],x[lp][0],0);*/
    anp=my.r;
/*    for (i=0;i<anp;i++) X->M[0][off[l]+i]=my.V[i];*/
    for (j=0;j<n;j++)
    { tmap(my,mg,xp[l],x[lp][j],0);
      for (i=0;i<anp;i++) X->M[j][off[l]+i]=my.V[i];
    }
    /* contribute to Constraint matrix */
    for (i=0;i<df[l];i++)
    T.M[l][off[l]+i]=1.0;
    tmap(my,mg,xp[l],x[lp][0],1); /* kill matrix allocation in tmap */
    freemat(my);freemat(mg);
  }
 
  freemat(y);
  if (getZ)  // obtain null space of T
  { *Z=initmat(np,np);
    QT(*Z,T,0);Z->r = T.r;
    freemat(T);
  } else     // return T itself
  { *Z=T;}
}

/*************************** End of cut and paste code! ***********************************/


void RArrayFromMatrix(double *a,long r,matrix *M)

/* copies matrix *M into R array a where r is the number of rows of a treated as
  a matrix by R */

{ int i,j;
  for (i=0;i<M->r;i++) for (j=0;j<M->c;j++) a[i+r*j]=M->M[i][j];
}


matrix Rmatrix(double *A,long r,long c)

/* produces a matrix from the array containing a (default) R matrix stored:
   A[0,0], A[1,0], A[2,0] .... etc */

{ int i,j;
  matrix M;
  M=initmat(r,c);
  for (i=0;i<r;i++) for (j=0;j<c;j++) M.M[i][j]=A[i+j*r];
  return(M);
}


void mgcv(double *yd,double *Xd,double *Cd,double *wd,double *Sd,
          double *pd, double *sp,int *offd,int *dimd,int *md,
          int *nd,int *qd,int *rd,double *sig2d,double *Vpd,double *edf)


/* Solves :

   minimise || W^{0.5} (Xp-y) ||^2 + \sum_i sp[i] p'S_i p subject to Cp=b

   selecting sp[] by GCV (sig2<=0.0) or UBRE (sig2>0.0) - W=diag(w).

   y is n by 1 data vector
   w is n by 1 weight vector (W above)
   X is n by q design matrix
   p is q by 1 parameter vector
   C is r by q constraint matrix
   S is m dimensional array of square matrices S[i] is dim[i] by dim[i]: the
     physical dimensions of each S[i] is given by the largest dim[i]- any
     padding may be used, as it will be ignored.
   sp is m dimensional array of smoothing parameters
   off[i] is m dimensional array of offsets specifiying row and column of
          overall penalty matrix which will contain S[i][0][0]
   Vp the q by q cov matrix for p set Vpd[0]<= 0.0 for now cov matrix calc
      and to > 0.0 for cov matrix to be returned
   edf is an array of estiamted degrees of freedom for each smooth
   The routine is an interface routine suitable for calling from the
   R and Splus.
  
*/

{ long n,q,r,*dim,*off,dimmax;
  int m,i,j,k;
  matrix *S,y,X,p,C,Z,w,Vp,L;
  double sig2,xx,gcvubre;
  //char msg[100];
  sig2= *sig2d;
  m=(int)*md;
  n=(long)*nd;
  q=(long)*qd;
  r=(long)*rd;
  if (m)
  { dim=(long *)calloc((size_t)m,sizeof(long));
    off=(long *)calloc((size_t)m,sizeof(long));
    S=(matrix *)calloc((size_t)m,sizeof(matrix));
  } 
  for (i=0;i<m;i++) off[i]=(long)offd[i];
  for (i=0;i<m;i++) dim[i]=(long)dimd[i];
  // set up matrices for MultiSmooth

  X=Rmatrix(Xd,n,q);p=Rmatrix(pd,q,1L);y=Rmatrix(yd,n,1L);w=Rmatrix(wd,n,1L);
  C=Rmatrix(Cd,r,q);
  dimmax=0;for (i=0;i<m;i++) if (dim[i]>dimmax) dimmax=dim[i];
  
  for (k=0;k<m;k++)
  { S[k]=initmat(dim[k],dim[k]);
    for (i=0;i<dim[k];i++) for (j=0;j<dim[k];j++) S[k].M[i][j]=Sd[k+i*m+j*m*dimmax];
  }

  Z=initmat(q,q);QT(Z,C,0);Z.r=C.r; // finding null space of constraints

  if (m)
  { 
    gcvubre=MultiSmooth(&y,&X,&Z,&w,S,&p,sp,off,m,&sig2);
   
    *sig2d=sig2;
  } else // no penalties, just solve the least squares problem
  { L=initmat(X.r,X.c); 
    mcopy(&X,&L);
    HQmult(L,Z,0,0);  // L=XZ
    L.c -= Z.r;
    for (i=0;i<w.r;i++) w.V[i]=sqrt(w.V[i]);
    p.r -= Z.r;        // solving in null space
    leastsq(L,p,y,w);  // Solve min ||w(Lp-y)||^2
    Vp=initmat(y.r,1L);// temp for fitted values
    matmult(Vp,L,p,0,0);
    sig2=0.0;for (i=0;i<y.r;i++) { xx=(Vp.V[i]-y.V[i])*w.V[i];sig2+=xx*xx;}
    
    if (n==p.r) sig2=0.0; else sig2/=(n-p.r);*sig2d=sig2;
    gcvubre=sig2/(n-p.r);
    for (i=0;i<w.r;i++) w.V[i] *=w.V[i];
    p.r+=Z.r;
    HQmult(p,Z,1,0); // back out of null space
    freemat(L);freemat(Vp);
  }
  for (i=0;i<q;i++) pd[i]=p.V[i];
  if (Vpd[0]>0.0)
  // calculate an estimate of the cov. matrix for the parameters:
  // Z[Z'(X'WX + \sum_i sp_i S[i])Z]^{-1}Z'
  { Vp=initmat(p.r,p.r);
    for (i=0;i<p.r;i++) for (j=0;j<=i;j++) // form X'WX
    { xx=0.0;
      for (k=0;k<X.r;k++) xx+=X.M[k][i]*w.V[k]*X.M[k][j];
      Vp.M[i][j]=Vp.M[j][i]=xx;
    }
    for (k=0;k<m;k++) // add on penalty terms
    { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++)
      Vp.M[i+off[k]][j+off[k]]+=sp[k]*S[k].M[i][j];
    }
    // now project into null space
    HQmult(Vp,Z,1,1);
    HQmult(Vp,Z,0,0);  
    for (j=0;j<Z.r;j++)
    for (i=0;i<p.r;i++) Vp.M[p.r-j-1][i]=Vp.M[i][p.r-j-1]=0.0;
    Vp.r -=Z.r;Vp.c -= Z.r;
    L=initmat(Vp.r,Vp.c);
    if (chol(Vp,L,1,1)) // invert
    { Vp.r+=Z.r;Vp.c+=Z.r;    // image in full space
      HQmult(Vp,Z,1,0);
      HQmult(Vp,Z,0,1);  
    } else
    { Vp.M[0][0]=-1.0; // signals failure to find cov - matrix - co-linearity problem
    }
    if (Vp.M[0][0]>-0.5) // work out edf per term
    { freemat(L);
      L=initmat(Vp.r,X.r);
      matmult(L,Vp,X,0,1);
      for (i=0;i<L.r;i++) for (j=0;j<L.c;j++) L.M[i][j]*=w.V[j];
      for (i=0;i<m;i++)
      { edf[i]=0.0;
        for (j=0;j<X.r;j++) for (k=off[i];k<off[i]+S[i].r;k++) 
        edf[i]+=X.M[j][k]*L.M[k][j];
      }
    
      for (i=0;i<Vp.r;i++) for (j=0;j<Vp.c;j++) Vp.M[i][j]*=sig2;
      freemat(L);
    }
    RArrayFromMatrix(Vpd,Vp.r,&Vp); // convert to R format
    freemat(Vp);
  } else Vpd[0]=0.0;

  // tidy up

  freemat(y);freemat(X);freemat(p);freemat(w);freemat(C);

  for (k=0;k<m;k++) freemat(S[k]);

  freemat(Z);
  if (m) {free(dim);free(off);free(S);}
}

void RQT(double *A,int *r,int*c)

/* Obtains the QT decomposition of matrix A (stored according to R conventions)
   AQ=[0,T] where T is reverse lower triangular (upper left is zero). r<c and 
   first c-r columns of Q are basis vectors for the null space of A 
   (Q orthogonal). Let this null space basis be Z. It is actually stored as 
   a series of r Householder rotations over the rows of A. Let u_i be the ith
   row of A (A[i,], i>=1) then the last i-1 elements of u_i are zero, while if 
   H_i=(I-u_i u_i') then Q=H_1 H_2 H_3 ...H_r.
   
   The main purpose of this routine *was* to provide a suitable representation 
   of the null space of any equality constraints on the problem addressed by 
   mgcv(). So if the constraints are Cp=0, RQT() was called to get an 
   appropriate null space basis in A. The non-obvious representation usually 
   saves much computing, since there are usually few constraints, resulting in 
   a high dimensional null space - in this case the Householder representation
   is very efficient. 

   However, the current version of mgcv() expects to get the constraint matrix    itself and not the null space. 
*/

{ matrix Q,B;
  B=Rmatrix(A,(long)(*r),(long)(*c));
  Q=initmat(B.r,B.c);
  QT(Q,B,0);
  RArrayFromMatrix(A,(long)(*r),&Q);
  freemat(Q);freemat(B);
}


void RGAMsetup(double *Xd,double *Cd,double *Sd,double *xpd,
             int *offd,double *xd,int  *md,int  *nd,int *dfd,int *nsdfd)

/* Interface routine to GAMsetup from R (and hopefully Splus)
   The arrays pointed to by all arguments must be initialised to the
   correct size prior to calling this routine.
   Inputs are:
   1. x - in R terms x[i][j] is the jth obs of the ith covariate (corresponding
          to the jth datapoint being modelled i.e. y[j]). Note however that the
          first nsdf x[i]'s are the first nsdf columns of the design matrix
          corresponding to the parametric part of the model (including the intercept); 
          the remaining x[i]'s are covariates that will be treated using smooths.
          Dimension (m+nsdf) by n.
   2. nsdf - the number of parametric terms (including the mean/intercept term)
   3. df - m terms; df[i] is the max. d.f. for the ith smooth
   4. m - number of smooths
   5. n - number of data to be modelled.
   Ouputs are:
   1. X - the n by q design matrix where q = \sum_i df[i]+nsdf
   2. xp - xp[i][j] is the jth knot position for the ith smooth. Dimension must
           be set to m by (max_i df[i]).
   3. S - array containing smoothness matrices. S[k][i][j] is ith row jth col of
          the kth smoothness penalty matrix. Dimension m by (max_i df[i]) by
          (max_i df[i]).
   4. C - Constraint matrix for problem. Input dimension m by q (see 1).
   5. off - array of offsets locating S[k]'s within overall constraint matrix
            S[k][0][0] goes in row off[k] col off[k]

*/

{ //char msg[200];
  matrix X,C,*S,*xp;
  long mdf,*off,*df;
  int m,n,nsdf,i,j,k;
  double **x;
  /* setup x[][], df[], off[], m, n, nsdf, S[], xp[], for calling GAMsetup.
     X, C, S[i]'s and xp[i]'s are initialised in the function */
  m= *md;
  n= *nd;
  nsdf= *nsdfd;
  S=(matrix *)calloc((size_t)m,sizeof(matrix));
  xp=(matrix *)calloc((size_t)m,sizeof(matrix));
  x=(double **)calloc((size_t)m+nsdf,sizeof(double *));
  for (i=0;i<m+nsdf;i++) 
  { x[i]=(double *)calloc((size_t)n,sizeof(double));
    for (j=0;j<n;j++) x[i][j]=xd[i+(m+nsdf)*j]; // loading data passed by R into x[][]
  }
  df=(long *)calloc((size_t)m,sizeof(long));
  for (i=0;i<m;i++) df[i]=(long)dfd[i];
  off=(long *)calloc((size_t)m,sizeof(long));
  /* now run GAMsetup to get X, C, S[], xp, off */
  GAMsetup(&X,&C,S,xp,off,x,m,n,df,nsdf,0);
  /* unload returned matrices into R arrays: X, C, S[k], xp */
  RArrayFromMatrix(Xd,(long)n,&X);
  RArrayFromMatrix(Cd,(long)m,&C);
  mdf=0;for (i=0;i<m;i++) if (mdf<df[i]) mdf=df[i];
  for (k=0;k<m;k++)
  { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++)
    Sd[k+i*m+j*m*mdf]=S[k].M[i][j];
    for (i=0;i<df[k];i++) xpd[k+m*i]=xp[k].V[i];
  }
  for (i=0;i<m;i++) offd[i]=(int)off[i];
  /* tidy up */

  freemat(X);freemat(C);
  for (i=0;i<m;i++) {freemat(S[i]);freemat(xp[i]);}
  free(S);free(xp);free(off);free(df);
}

void gam_map(matrix tm, matrix *t, double *x, int m,int nsdf, int kill)

/* Consider a gam with a set of parametric terms given by the first nsdf terms of 
   vector x and m smooth terms    
   represented by splines, each with df[i] parameters. It is assumed the the splines
   are parameterized by the parameter values at the knots, and that the model parameters
   are arranged: constant, parametric terms, spline terms_1, splineterms_2 etc...
   If p is the parameter vector then this routine returns tm, such that tm'p = f(x).
   
   To save flops this routine saves some matrices statically - to free them, and set 
   things up for a new knot sequence call with kill==1.
*/


{ static matrix *D;static char first=1;
  static int terms=0;
  matrix h;
  double xx;
  long i,j,k,l;
  if (first)
  { first=0;
    D=(matrix *)calloc((size_t)m,sizeof(matrix));
    for (i=0;i<m;i++)
    { h=initmat(t[i].r-1,1L);
      for (j=0L;j<t[i].r-1;j++) h.V[j]=t[i].V[j+1]-t[i].V[j];
      D[i]=getD(h,0); 
      freemat(h);
    } 
    terms=m;
  }
  // deal with the  parametric terms first....
  for (k=0;k<nsdf;k++) tm.V[k]=x[k];
  // now the splines.....
  for (j=0;j<m;j++) 
  { xx=x[nsdf+j];
    i=0L;while((xx>t[j].V[i+1])&&(i<t[j].r-2)) i++; 
    for (l=0;l<t[j].r;l++)
    { tm.V[k]=D[j].M[i][l]*d0(t[j].V[i],t[j].V[i+1],xx)+
	          D[j].M[i+1][l]*d1(t[j].V[i],t[j].V[i+1],xx);
      if (l==i) tm.V[k]+=b0(t[j].V[i],t[j].V[i+1],xx);
      if (l==(i+1)) tm.V[k]+=b1(t[j].V[i],t[j].V[i+1],xx);
      k++;
    }
  }
  if (kill)
  { first=1;
    for (i=0;i<terms;i++) freemat(D[i]);free(D);
    terms=0;
  }
}

void RGAMpredict(double *xpd,int *nsdf,int *df,int *m,double *xd,int *np,double *p,
                 double *Vpd,double *etad,double *sed,int *control)

/* Routine for prediction from GAMs made up of penalized cubic regression splines.
   
   xp[i][j] is the jth knot position for the ith smooth (xp[i] dimension must be set to 
            max_i(df[i])...)
   nsdf     is the number of non spline terms including the constant.
   df[i]    is the number of parameters for the ith spline
   m        is the number of splines
   x[i][j]  is the ith observation of jth covariate, for which predictions are required - 
            first nsdf are none spline terms
   np       is the number of predictions required (number of obs. per cov.)
   p[i]     is the ith parameter: 0 to nsdf-1 are the none spline terms
                                  nsdf to nsdf+df[0]-1 is the first spline
                                  nsdf+df[i-1] to nsdf+df[i]-1 is the ith spline 
   Vp[i][j] is covariance of p[i] and p[j]
   
   mu is the o/p linear predictor
   se is the o/p standard error on l.p.
   
   control sets the type of output: 
           0: l.p. no s.e.
           1: l.p. with s.e.
           2: predictor for each term, no s.e.
           3: predictor for each term with s.e.

   note  that the dimensions of mu and se must be appropriate for the control option selected.
*/
           
{ matrix tm,*xp,Vp,eta,se;
  int i,j,k,nb,kk,l; 
  double **x,z,Vt;
  char info[200];
  // perform bits of unpacking ......
  x=(double **)calloc((size_t) *np,sizeof(double *));
  for (i=0;i<*np;i++) 
  { x[i]=(double *)calloc((size_t)*m+*nsdf,sizeof(double));
    for (j=0;j< *m + *nsdf;j++) 
    { x[i][j]=xd[j+(*m + *nsdf)*i]; // loading data passed by R into x[][]
    } 
  }
  nb =  *nsdf;  // number of parameters
  for (i=0;i < *m;i++) nb += df[i];
  Vp=Rmatrix(Vpd,(long)nb,(long)nb); // param cov matrix
  // initialise a couple of storage matrices
  if (*control<2) k=1;else k=  *nsdf + *m;  // don't allocate more than is needed
  eta=initmat((long)k,(long)*np); 
  se=initmat((long)k,(long)*np);
  // need to unpack xpd into xp here .....
  xp=(matrix *)calloc((size_t) *m,sizeof(matrix));
  for (k=0;k< *m;k++)
  { xp[k]=initmat((long)df[k],1L);
    for (i=0;i<df[k];i++) xp[k].V[i]=xpd[k+ *m *i];
  }  

  tm=initmat((long)nb,1L);

  // the constant...... 
  for (k=0;k< *np;k++) // loop through the predictions
  { gam_map(tm,xp,x[k],*m,*nsdf,0);
    if (*control<2) // then linear predictor required
    { z=0.0;for (i=0;i<tm.r;i++) z+=tm.V[i]*p[i];eta.M[0][k]=z;
      if (*control==1) // get s.e.
      { z=0.0;
        for (i=0;i<tm.r;i++) 
        { Vt=0.0; for (j=0;j<tm.r;j++) Vt+=Vp.M[i][j]*tm.V[j];
          z+=tm.V[i]*Vt; 
        }
        se.M[0][k]=sqrt(z);
      }
    } else // individual predictors wanted
    { // the constant .....
      //eta.M[0][k]=p[0];
      //if (*control==3) se.M[0][k]=sqrt(Vp.M[0][0]); 
      // other parametric terms .....
      for (i=0;i< *nsdf;i++) 
      { eta.M[i][k]=p[i]*x[k][i];
      if (*control==3) se.M[i][k]=sqrt(x[k][i]*x[k][i]*Vp.M[i][i]); 
      }
      // and now the splines .....
      kk= *nsdf;
      for (l=0;l< *m;l++)
      { z=0.0;for (i=0;i<df[l];i++) { z+=tm.V[kk+i]*p[kk+i];}
        eta.M[ *nsdf + l][k] = z;
        if (*control==3)
        { z=0.0;
          for (i=0;i<df[l];i++) 
          { Vt=0.0; 
            for (j=0;j<df[l];j++) Vt+=Vp.M[kk+i][kk+j]*tm.V[kk+j];
            z+=Vt*tm.V[kk+i];              
          } 
          if (z<0.0) {sprintf(info,"z = %g in RGAMpredict",z);
          infobox(info);}     
          se.M[ *nsdf + l][k]=sqrt(z);
        }
        kk+=df[l];
      }
    } 
  } 
  // convert results for o/p
  gam_map(tm,xp,x[0],*m,*nsdf,1);  // free memory in gam_map()
  freemat(tm);
  RArrayFromMatrix(etad,eta.r,&eta);
  if ((*control)%2) RArrayFromMatrix(sed,se.r,&se);
  // tidy up....
  for (i=0;i< *np;i++) free(x[i]);free(x);
  freemat(Vp);freemat(eta);freemat(se);
  for (k=0;k< *m;k++) freemat(xp[k]);free(xp);
}    

/*********************************************************************************************/
/* Bug fix record:

1. 20/10/00: Knot placement method in GAMsetup() modified. Previous method had an error, so 
   that when df for a term was close to the number of data, a non-existent covariate value
   (i.e. out of array bound). New code also yields more regular placement, and now deals with 
   repeat values of covariates.
2. 20/10/00: Modified mgcv() to cope with problems with no penalties, by call to leastsq() -
   this is needed to allow gam() to fit models with fixed degrees of freedom.
3. 5/1/01: Modified RGAMsetup(), GAMsetup(), gam_map() and RGAMpredict() so that nsdf is now
   total number of non-spline parameters including any constant. Hence R code must now provide 
   column for constant explicitly.
4. 5/1/01: fixed bug in RGAMpredict - standard errors of parametric components of linear predictor
   were wrongly calculated.
*/


/* Source code for mgcv.dll/.so multiple smoothing parameter estimation code,
suitable for interfacing to R 

Copyright (C) 2000-2002 Simon N. Wood  snw@st-and.ac.uk

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

    /*#define WINDOWS*/ 

#ifdef WINDOWS
#include <windows.h>   /* For easy self window porting, without neading everything R.h needs */
#include <stdarg.h>
#else
#include <R.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gcv.h"
#include "tprs.h"
#include "mgcv.h"
#include "matrix.h"
#include "qp.h"
#include "general.h"

#define round(a) ((a)-floor(a) <0.5 ? (int)floor(a):(int) floor(a)+1)


#ifdef WINDOWS
/* following routine provides a version of Rprintf that outputs to a file instead 
   of console, to facilitate debugging under windows.....
*/
typedef struct{
double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20;
} _dodgy_print_storage;

void Rprintf(char *str,...)
{ _dodgy_print_storage ug; 
  FILE *f;
  va_list argptr;
  f=fopen("d:/simon/Rdump.txt","at");
  va_start(argptr,str);
  ug=va_arg(argptr,_dodgy_print_storage);
  fprintf(f,str,ug);
  va_end(argptr);
  fclose(f);
}
#endif

void ErrorMessage(char *msg,int fatal)

{ 
#ifdef WINDOWS
  MessageBox(HWND_DESKTOP,msg,"Info!",MB_ICONEXCLAMATION|MB_OK); 
#else  
if (fatal) error("%s",msg);else warning("%s",msg);
#endif
}

void infobox(char *msg)

{ 
#ifdef WINDOWS
MessageBox(HWND_DESKTOP,msg,"Info!",MB_ICONEXCLAMATION|MB_OK); 
#else
warning("%s",msg);
#endif
}

/* The following are some rather ancient routines used to set up an example
   additive model using regression (cubic) splines, via RGAMsetup(). */


/********** The following are from spline.c (rather plodding coding) ***********/

/* The next 4 functions are basis functions for the 1st derivative
   representation of a cubic spline */


double b0(x0,x1,x) double x0,x1,x;

/* multiplies function value at x0 */

{ double res,h,xx1;
  if (x<x0) return(1.0);
  if (x>x1) return(0.0);
  h=x1-x0;xx1=x-x1;
  res=2.0*(x-x0+0.5*h)*xx1*xx1/(h*h*h);
  return(res);
}

double b1(x0,x1,x) double x0,x1,x;

/* multiplies function value at x1 */

{ double res,h,xx0;
  if (x<x0) return(0.0);
  if (x>x1) return(1.0);
  h=x1-x0;xx0=x-x0;
  res= -2.0*(x-x1-0.5*h)*xx0*xx0/(h*h*h);
  return(res);
}

double d0(x0,x1,x) double x0,x1,x;

/* multiplies gradient at x0 */

{ double res,h,xx1;
  if (x<x0) res=x-x0;  /* before start of interval - use linear extrapolation */
  else 
  if (x>x1) res=0.0;  /* after end of interval - d0 plays no part */
  else                /* within interval */
  { h=x1-x0;xx1=x-x1;
    res=(x-x0)*xx1*xx1/(h*h);
  }
  return(res);
}

double d1(x0,x1,x) double x0,x1,x;

/* multiplies gradient at x1 */

{ double res,h,xx0;
  if (x<x0) res=0.0; /* d1 plays no part */
  else
  if (x>x1) res=x-x1; /* linear extrapolation after end of interval */
  else                /* in interval */
  { h=x1-x0;xx0=x-x0;
    res=xx0*xx0*(x-x1)/(h*h);
  }
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

void MonoCon(matrix *A,matrix *b,matrix *x,int control,double lower,double upper ) 

/* gets matrices A and b for constraints of the form Ay>=b ensuring monotonic
   change  of the cubic spline interpolating (x_i,y_i) where h_i=x_{i+1}-x_i 
   
   control indicates type of constraints:
   up=control/4 - 0 for decrease, 1 for increase
   lo=(control-up*4)/2 - 1 for lower bound, 0 no lower bound
   hi=(control-up*4-lo*2) - 1 for upper bound, 0 no upper bound
   control = 4*up+2*lo+hi
*/

{ long i,j,n;
  int up,lo,hi;
  double m;
  matrix h,D;
  h=initmat(x->r-1,1L);
  n=h.r;
  for (i=0;i<n;i++) h.V[i]=x->V[i+1]-x->V[i];
  D=getD(h,0);
  up=control/4;control=control%4;
  lo=control/2;control=control%2;
  hi=control;
  if (up) m= -1.0; else m=1.0; 
  (*A)=initmat(4*n+hi+lo,n+1);
  for (i=0;i<n;i++)
  { for (j=0;j<n+1;j++)
    { if (j==i)
      { A->M[i][j]=(D.M[i][j]+3.0/h.V[i])*m;   /**not certain of d.M update**/
	    A->M[i+n][j]=(D.M[i+1][j]+3.0/h.V[i])*m;
	    A->M[i+2*n][j]=m;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      } else
      if (j==(i+1))
      { A->M[i][j]=(D.M[i][j]-3.0/h.V[i])*m;
	    A->M[i+n][j]=(D.M[i+1][j]-3.0/h.V[i])*m;
	    A->M[i+2*n][j]= -m;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      } else
      { A->M[i][j]=D.M[i][j]*m;
	    A->M[i+n][j]=D.M[i+1][j]*m;
	    A->M[i+2*n][j]=0.0;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      }
    }
  }
  *b = initmat(A->r,1L);
  if (lo)
  { for (j=0;j<n+1;j++) A->M[4*n][j]=0.0;
    if (up) A->M[4*n][0]=1.0; else A->M[4*n][n]=1.0;
    b->V[4*n]=lower;
  }
  if (hi)
  { for (j=0;j<n+1;j++) A->M[4*n][j]=0.0;
    if (up) A->M[4*n+lo][n]=-1.0; else A->M[4*n+lo][0]=-1.0;
    b->V[4*n+lo]=upper;
  }
  freemat(D);
  freemat(h);
}


void tmap(matrix tm,matrix t,double time,int kill)
 /* to release static matrix allocation set kill to 1 otherwise 0 and
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
  { tm.V[0]=1.0;}
  else
  { i=0L;while((time>t.V[i+1])&&(i<t.r-2)) i++;
    for (k=0;k<t.r;k++)
    tm.V[k]=D.M[i][k]*d0(t.V[i],t.V[i+1],time)+
	    D.M[i+1][k]*d1(t.V[i],t.V[i+1],time);
    tm.V[i]+=b0(t.V[i],t.V[i+1],time);
    tm.V[i+1]+=b1(t.V[i],t.V[i+1],time);
   
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


void crspline(double *x,int n,int knots,matrix *X,matrix *S, matrix *C, matrix *xp)

/* Sets up a cubic regression spline, by specifying a set of knot locations 
   spread evenly throughout the covariate values, and using cubic hermite
   polynomials to represent the spline.

   This is a very fast way of setting up 1-d regression splines, but the
   basis provided by tprs() will be better. This code gives the basis used
   in mgcv prior to version 0.6.

   The inputs are:
   x - the array of covariate values
   n - the length of x
   knots - the number of knots

   The outputs are:

   X - design matrix 
   S - wiggliness penalty matrix
   C - constraint matrix 
   xp - knot position vector 

   so that a penalized regression spline could be fitted by minimising:

       ||Xb-y||^2 + \lambda b'Sb

   subject to Cb=0 if the sum of the values of the spline at the xp values 
   is to be zero - a constraint that is useful in GAM modelling, to ensure
   model identifiability.
*/

{ int j,i,k;
  matrix y,my;
  double dx,xx;
  /* sort x values into order and get list of unique x values .... */
  y=initmat((long)n,1L);
  for (j=0;j<n;j++) y.V[j]=x[j];
  y.r=(long)n;
  sort(y); /* next reduce to list of unique values..... */
  k=0;for (i=0;i<n;i++) if (y.V[k]!=y.V[i]) { k++;y.V[k]=y.V[i];} y.r=(long)k+1;
  dx=(y.r-1)/(knots-1.0);
  /* now place the knots..... */
  if (xp->V[0]==xp->V[1]) /* then knot positions have not been supplied */
  { xp->V[0]=y.V[0];
    for (i=1;i<knots-1;i++)  /* place knots */
    { xx=dx*i;
      k=(int)floor(xx);
      xx -= k;
      xp->V[i]=(1-xx)*y.V[k]+xx*y.V[k+1];
    } 
    xp->V[knots-1]=y.V[y.r-1];
  }
  freemat(y);
 
  /* create the wiggliness measure matrix...... */
  getSmooth(S,*xp,0);
 
  /* create the constraint matrix ...... */
  *C=initmat(1L,(long)knots);
  for (i=0;i<knots;i++) C->M[0][i]=1.0;
  /* and finally, create the design matrix ...... */
  *X=initmat((long)n,xp->r);
  my=initmat(xp->r,1L);
  for (j=0;j<n;j++)
  { tmap(my,*xp,x[j],0);
    for (i=0;i<my.r;i++) X->M[j][i]=my.V[i];
  }   
  tmap(my,*xp,x[0],1); /* kill matrix allocation in tmap */
 
  freemat(my);
}

void GAMsetup(matrix *X,matrix *Z,matrix *S,matrix *UZ,matrix *Xu,matrix *xp,
              long *off,double **x,double **by,double *knots,int m,int n,int *dim,int *s_type, long *df,
              int *p_order,int *by_exists,int *n_knots,int nsdf,int getZ)

/* Sets up the Design matrix and Smoothness constraints for a regression spline
   based GAM. There are m smoothers, n datapoints. The ith smoother has
   df[i] parameters starting at off[i], with smoothing constraint matrix
   S[i], and is dim[i] dimensional, of type s_type[i]. s_type 0 is crs.
   x[i][j], contains the jth value of the ith covariate.
   by[i][j] is the jth value for the ith "by" variable.
   by_exists[] is an array of 0/1's indicating existence of otherwise of by
   variables for each smooth
   knots[] may contain information specifying knot locations for use when 
           setting up a basis. Arrays of covariate values are stacked end 
           on end in this array in an order corresponding to the ordering 
           of the rows of x. 
   n_knots[] - number of knots for each smooth term.
   The first nsdf x[i]s are columns of the design matrix for the non-spline
   part of the model (including the intercept term) .
   X & Z are initialised in the routine.
   Each S[i], UZ[i] and Xu[i] is also initialised in the routine, but note that 
   no memory is allocated for UZ and Xu for cubic regression splines: their rows
   and columns are just set to 0.
   xp[i] is returned containing the vector of knot positions for the ith spline.
   Each xp[i] is initialised within the routine.

   Z contains the null space for the GAM side constraints if getZ is not set to
     0 - if getZ==0 then Z is returned containing the constraint matrix itself.

   All smooths are constrained so that their parameters sum to zero - the first
   parameter in the model will then be a mean. Hence for a pure gam() x[1] should
   be a vector of 1's and nsdf should be equal to one, to allow for a constant.
*/

{ long np,l,i,j,k,lp,by_counter=0,nsm;
  matrix T,Xl,Cl;
  double *XMj,*XlMj,bylj,**knt;
  np=nsdf;nsm=0;for (i=0;i<m;i++) {np+=df[i];nsm+=dim[i];}
  if (m) T=initmat((long)m,np);T.r=0L;
  *X=initmat((long)n,(long)np);
  knt=(double **)calloc((size_t)nsm+1,sizeof(double *)); /* knt[i][j] is jth knot value for ith variable */
  k=1;knt[0]=knots;for (i=0;i<m;i++) for (j=0;j<dim[i];j++) { knt[k]=knt[k-1]+n_knots[i];k++;}
  off[0]=nsdf;
  for (j=0;j<nsdf;j++) for (i=0;i<n;i++) X->M[i][j]=x[j][i];
  lp=nsdf;  /* where we are in the x[i]'s */
  for (l=0;l<m;l++) /* work through the smooths */
  { if (l) off[l]=off[l-1]+df[l-1];else off[l]=nsdf;
    
    /* get design matrix etc for this particular term...... */
    if (s_type[l]==0)
    { xp[l]=initmat((long)df[l],1L);
      if (n_knots[l]==df[l]) /* knots provided */
      { for (i=0;i<df[l];i++) xp[l].V[i]=knt[lp-nsdf][i];
      }
      crspline(x[lp],n,df[l],&Xl,S+l, &Cl, xp+l);
      UZ[l].r=UZ[l].c=Xu[l].r=Xu[l].c=0L;
    }
    else
    { tprs_setup(x+lp,knt+lp-nsdf,p_order[l],dim[l],n,df[l],1,&Xl,S+l,UZ+l,Xu+l,n_knots[l]);
    }
    /* copy Cl into overall constraint matrix T..... */
    if (s_type[l]==0)
    { for (i=0;i<df[l];i++)
      T.M[T.r][off[l]+i]=Cl.M[0][i];
      freemat(Cl);T.r++;
    } else /* impose a sum to zero side condition on the spline term */
    { for (i=0;i<df[l];i++)
      for (j=0;j<Xl.r;j++)
      T.M[T.r][off[l]+i]+=Xl.M[j][i];
      T.r++;
    }
    /* multiply design matrix component by "by" variable if there is one */
    if (by_exists[l])
    { for (j=0;j<n;j++)
      { bylj=by[by_counter][j];
        for (XlMj=Xl.M[j];XlMj<Xl.M[j]+Xl.c;XlMj++) *XlMj *= bylj;
      } 
      by_counter++;  
    }
    /* copy Xl into overall design matrix X  */
    for (j=0;j<n;j++)
    { XMj=X->M[j]+off[l];
      XlMj=Xl.M[j];
      for (i=0;i<Xl.c;i++) XMj[i]=XlMj[i];
    }
    lp+=dim[l];
    freemat(Xl);
   
  }
  if (getZ)  /* obtain null space of T */
  { if (T.r)
    { *Z=initmat(np,np);
      QT(*Z,T,0);
    } 
    Z->r = T.r;
    if (T.r) freemat(T);
  } else     /* return T itself */
  { *Z=T;}
  free(knt);
  
}

/*************************** End of cut and paste code! ***********************************/


void RUnpackSarray(int m,matrix *S,double *RS)
/* unpacks the R array RS into an array of matrices initialized to the correct dimensions 
   let kk = sum_{i=0}^k S[i].r*S[i].c
   Then the kth matrix starts at element kk of RS and stops at element k(k+1)
   ... let this extracted array be M. S[k].M[i][j]=M[i+S[k].r*j] - in this way we ensure that 
   M can be extracted straight to a matrix in R with 
   A<-matrix(M,S[k].r,S[k].c) 
*/ 
{ int start,i,j,k;
  start=0;
  for (k=0;k<m;k++)
  { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++) S[k].M[i][j]=RS[start+i+S[k].r*j];
    start += S[k].r*S[k].c;
  }
}

void RPackSarray(int m,matrix *S,double *RS)
/* Packs an array of matrices S[] into an R array RS in the manner described in RUnpackSarray
*/
{ int start,i,j,k;
  start=0;
  for (k=0;k<m;k++)
  { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++) RS[start+i+S[k].r*j]=S[k].M[i][j];
    start += S[k].r*S[k].c;
  }

}



void RArrayFromMatrix(double *a,long r,matrix *M)

/* copies matrix *M into R array a where r is the number of rows of A treated as
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
          int *nd,int *qd,int *rd,double *sig2d,double *Vpd,double *edf,
          double *conv_tol,int *ms_max_half,double *ddiag, int *idiag,double *sdiag,
          int *direct_mesh,double *min_edf,double *gcvubre,double *target_edf,
          int *fixed_sp,double *hat)


/* Solves :

   minimise || W^{0.5} (Xp-y) ||^2 + \sum_i sp[i] p'S_i p subject to Cp=b

   selecting sp[] by GCV (sig2<=0.0) or UBRE (sig2>0.0) - W=diag(w).

   y is n by 1 data vector
   w is n by 1 weight vector (W above)
   X is n by q design matrix
   p is q by 1 parameter vector
   C is r by q constraint matrix
   S is contains m square matrices S[i] each of dimension dim[i] by dim[i]: 
     S is actually a 1-d array, let start_k=sum_{i=0}^{k-1} S[i].r*S[i].c
     S[k].M[i][j]= S[start_k+i+S[k].r*j] - use RUnpackSarray() above
   sp is m dimensional array of smoothing parameters
   off[i] is m dimensional array of offsets specifiying row and column of
          overall penalty matrix which will contain S[i][0][0]
   Vp the q by q cov matrix for p set Vpd[0]<= 0.0 for now cov matrix calc
      and to > 0.0 for cov matrix to be returned
   edf is an array of estimated degrees of freedom for each smooth - need 
       Vp to be calculated to work these out.

   The routine is an interface routine suitable for calling from the
   R and Splus.

  
   direct_mesh controls the number of steps in the direct search for the 
   overall smoothing parameter: 100 is reasonable.

   min_edf is the minimum possible estimated degrees of freedom - useful
   for setting limits on overall smoothing parameter: <=0 => ignore.

   fixed_sp should be set to 1 if the supplied smoothing parameters are to be 
            treated as fixed, in which case model is fitted as an augmented 
            least squares problem. 
   
   hat is an array of elements from the leading diagonal of the hat/influence matrix

   gcvubre is the minimum GCV or UBRE score achieved.

   There are 3 arrays of convergence diagnostics returned:

   ddiag[] is a 3*m array of doubles - first m are gradients of gcv/ubre
           score at convergence, next m is leading diagonal of Hessian, last
		   m are eigenvalues of Hessian.
   idiag[] is a length 3 array of integers. idiag[0] is iterations taken to 
           converg, idiag[1] is 1/0 if second guess was/wasn't good, idiag[2]
		   is 1/0 if terminated/not on step failure.
   sdiag[] is a 2*direct_mesh array of doubles. The first half are model edf's, 
           the second half corresponding GCV/UBRE scores. Useful for checking 
           for multiple minima (but not conclusive, of course). 
  
*/

{ long n,q,r,*dim,*off,dimmax;
  int m,i,j,k,gcv;
  msctrl_type msctrl;
  msrep_type msrep;
  matrix *S,y,X,p,C,rS,Z,w,Vp,L,ya,wa;
  double sig2,xx,trA;
  /*char msg[100]; */
  sig2= *sig2d;
  if (sig2<=0) gcv=1; else gcv=0;
  m=(int)*md;
  n=(long)*nd;
  q=(long)*qd;
  r=(long)*rd;
  if (m)
  { dim=(long *)calloc((size_t)m,sizeof(long));
    off=(long *)calloc((size_t)m,sizeof(long));
    S=(matrix *)calloc((size_t)m,sizeof(matrix));
  } else { dim=off=&n;S=&y;} /* purely to avoid compiler warning */
  for (i=0;i<m;i++) off[i]=(long)offd[i];
  for (i=0;i<m;i++) dim[i]=(long)dimd[i];
  /* set up matrices for MultiSmooth */

  X=Rmatrix(Xd,n,q);p=Rmatrix(pd,q,1L);y=Rmatrix(yd,n,1L);w=Rmatrix(wd,n,1L);
  if (r) C=Rmatrix(Cd,r,q); else C.r=C.c=0L;
  dimmax=0;for (i=0;i<m;i++) if (dim[i]>dimmax) dimmax=dim[i];
  
  for (k=0;k<m;k++) S[k]=initmat(dim[k],dim[k]);
  RUnpackSarray(m,S,Sd);  

  if (C.r) {Z=initmat(q,q);QT(Z,C,0);}
  Z.r=C.r; /* finding null space of constraints */
  if (m&&!*fixed_sp) /* then there are some penalties with unknown smoothing parameters */
  { init_msrep(&msrep,m,*direct_mesh);
    msctrl.conv_tol= *conv_tol;
    msctrl.max_step_half= *ms_max_half;
    msctrl.min_edf= *min_edf;
    msctrl.target_edf= *target_edf;
    *gcvubre=MultiSmooth(&y,&X,&Z,&w,S,&p,sp,off,m,&sig2,&msctrl,&msrep,*direct_mesh);
    for (i=0;i<m;i++) /* copy out diagnostics for return */
    { ddiag[i]=msrep.g[i];
      ddiag[m+i]=msrep.h[i];
      ddiag[2*m+i]=msrep.e[i];
    }
    for (i=0;i< *direct_mesh;i++)
    { sdiag[i]=msrep.edf[i];
      sdiag[i+ *direct_mesh]=msrep.score[i];
    }
    idiag[0]=msrep.iter;idiag[1]=msrep.inok;idiag[2]=msrep.step_fail;
    free_msrep(&msrep);
    if (gcv) *sig2d=sig2;
  } else /* no penalties, just solve the least squares problem */
  { if (m&& *fixed_sp) /* then there are penalties with fixed smoothing parameters to contend with */
    { Vp=initmat(p.r,p.r); /* temporary storage for sum of smooths */     
      for (k=0;k<m;k++) /* add up penalty terms */
      { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++)
        Vp.M[i+off[k]][j+off[k]]+=sp[k]*S[k].M[i][j];
      }
      /* now project into null space */
      HQmult(Vp,Z,1,1);
      HQmult(Vp,Z,0,0);  
      Vp.r -=Z.r;Vp.c -= Z.r;
      /* Now find square root of total penalty: rS */
      root(&Vp,&rS,8*DOUBLE_EPS);
      freemat(Vp);
      L=initmat(X.r+rS.c,X.c); /* augmented design matrix */
      L.r=X.r;
      mcopy(&X,&L);     /* design matrix part */ 
      HQmult(L,Z,0,0);  /* L=XZ */
      L.c -= Z.r;
      /* now add on penalty */
      for (i=X.r;i<X.r+rS.c;i++) for (j=0;j<L.c;j++) L.M[i][j]=rS.M[j][i-X.r];
      L.r=X.r+rS.c;freemat(rS);
      /* create augmented y and w vectors */
      ya=initmat(L.r,1L);wa=initmat(ya.r,1L);
      for (i=0;i<y.r;i++) { ya.V[i]=y.V[i];wa.V[i]=sqrt(w.V[i]);}
      for (i=y.r;i<L.r;i++) { wa.V[i]=1.0;ya.V[i]=0.0;}
      /* get tr(A) */
      Vp=initmat(X.c,X.c);Vp.c=Vp.r=L.c;
      matmult(Vp,L,L,1,0); /* X'X+S, basically */
      rS=initmat(Vp.r,Vp.c); /* dummy */
      if (!chol(Vp,rS,1,1)) /* try cheap inversion by choleski */  
      { pinv(&Vp,0.0);     /* but use svd if this fails */
      }
      Vp.r+=Z.r;Vp.c+=Z.r;    /* get image in full space */
      HQmult(Vp,Z,1,0);HQmult(Vp,Z,0,1);  
      freemat(rS);
      rS=initmat(Vp.r,X.r);
      matmult(rS,Vp,X,0,1); /* basically (X'X+S)^{-1}X' */
      freemat(Vp); 
      for (i=0;i<rS.r;i++) for (j=0;j<rS.c;j++) rS.M[i][j]*=w.V[j]; /* (X'X+S)^{-1}X'W */
      trA=0.0;for (j=0;j<X.r;j++) for (k=0;k<X.c;k++) trA+=X.M[j][k]*rS.M[k][j];
      freemat(rS);
    } else /* no penalties to consider */
    { L=initmat(X.r,X.c); 
      mcopy(&X,&L);
      HQmult(L,Z,0,0);  /* L=XZ */
      L.c -= Z.r;
      wa=initmat(w.r,1L);
      for (i=0;i<w.r;i++) wa.V[i]=sqrt(w.V[i]);
      ya=y;
      trA=p.r-Z.r;
    }
    p.r -= Z.r;        /* solving in null space */
    leastsq(L,p,ya,wa);  /* Solve min ||wa(Lp-ya)||^2 */
    Vp=initmat(ya.r,1L); /* temp for fitted values */
    matmult(Vp,L,p,0,0); /* oversized if penalties used, but latter elements ignored below */
    sig2=0.0;for (i=0;i<y.r;i++) { xx=(Vp.V[i]-y.V[i])*wa.V[i];sig2+=xx*xx;}
    if (n==p.r) sig2=0.0; else sig2/=(n-trA);
    if (gcv) 
    { *sig2d=sig2;
      *gcvubre=n*sig2/(n-trA);
    } else /* UBRE */
    { *gcvubre = sig2 / n * (n-trA)-2.0 / n * *sig2d  * ( n - trA) + *sig2d;  
    } 
    p.r+=Z.r;
    HQmult(p,Z,1,0); /* back out of null space */
    freemat(L);freemat(Vp);freemat(wa);
    if (m&& *fixed_sp) freemat(ya);
  }
  for (i=0;i<q;i++) pd[i]=p.V[i];
  if (Vpd[0]>0.0)
  /* calculate an estimate of the cov. matrix for the parameters: */
  /* Z[Z'(X'WX + \sum_i sp_i S[i])Z]^{-1}Z' */
  { Vp=initmat(p.r,p.r);
    for (i=0;i<p.r;i++) for (j=0;j<=i;j++) /* form X'WX */
    { xx=0.0;
      for (k=0;k<X.r;k++) xx+=X.M[k][i]*w.V[k]*X.M[k][j];
      Vp.M[i][j]=Vp.M[j][i]=xx;
    }
    for (k=0;k<m;k++) /* add on penalty terms */
    { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++)
      Vp.M[i+off[k]][j+off[k]]+=sp[k]*S[k].M[i][j];
    }
    /* now project into null space */
    HQmult(Vp,Z,1,1);
    HQmult(Vp,Z,0,0);  
    for (j=0;j<Z.r;j++)
    for (i=0;i<p.r;i++) Vp.M[p.r-j-1][i]=Vp.M[i][p.r-j-1]=0.0;
    Vp.r -=Z.r;Vp.c -= Z.r;
    L=initmat(Vp.r,Vp.c);
    if (!chol(Vp,L,1,1)) /* try cheap inversion by choleski */  
    { Rprintf("rank of Vp = %ld\n",pinv(&Vp,DOUBLE_EPS));     /* but use svd if this fails */
      Rprintf("using pinv\n");
    }
    Vp.r+=Z.r;Vp.c+=Z.r;    /* get image in full space */
    HQmult(Vp,Z,1,0);
    HQmult(Vp,Z,0,1);  
    freemat(L);
    /* now work out edf per term - meaningful only if penalties are non-overlapping */
    L=initmat(Vp.r,X.r);
    matmult(L,Vp,X,0,1);
    for (i=0;i<L.r;i++) for (j=0;j<L.c;j++) L.M[i][j]*=w.V[j];
    for (i=0;i<m;i++)
    { edf[i]=0.0;
      for (j=0;j<X.r;j++) for (k=off[i];k<off[i]+S[i].r;k++) 
      edf[i]+=X.M[j][k]*L.M[k][j];
    } 
    /* work out elements on leading diagonal of hat matrix */
    for (i=0;i<X.r;i++)
    { hat[i]=0.0;
      for (j=0;j<X.c;j++) hat[i]+=X.M[i][j]*L.M[j][i]; 
    }
    /* multiply Vp by estimated scale parameter so that it is proper covariance matrix estimate */
    for (i=0;i<Vp.r;i++) for (j=0;j<Vp.c;j++) Vp.M[i][j] *= *sig2d;
    freemat(L);
    RArrayFromMatrix(Vpd,Vp.r,&Vp); /* convert to R format */
    freemat(Vp);
  } else { Vpd[0]=0.0;hat[0]=0.0;}

  /* tidy up */

  freemat(y);freemat(X);freemat(p);freemat(w);
  if (C.r) freemat(C);

  for (k=0;k<m;k++) freemat(S[k]);

  if (Z.r) freemat(Z);
  if (m) {free(dim);free(off);free(S);}
#ifdef MEM_CHECK
  dmalloc_log_unfreed(); dmalloc_verify(NULL);
#endif
}


void RuniqueCombs(double *X,int *r, int *c)

/* X is a matrix. This routine finds its unique rows and strips out the 
   duplicates. This is useful for finding out the number of unique covariate
   combinations present in a set of data. */

{ matrix B,Xd;
  int i,*ind;
  B=Rmatrix(X,(long)(*r),(long)(*c));
  Xd=initmat(B.r,B.c+1);
  Xd.c--;mcopy(&B,&Xd);freemat(B);Xd.c++;
  for (i=0;i<Xd.r;i++) Xd.M[i][Xd.c-1]=(double)i;
  ind=Xd_strip(&Xd);
  Xd.c--; /* hide index array  */
  RArrayFromMatrix(X,Xd.r,&Xd);  /* NOTE: not sure about rows here!!!! */
  *r = (int)Xd.r; 
  freemat(Xd);free(ind);
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif 
}

void RMonoCon(double *Ad,double *bd,double *xd,int *control,double *lower,double *upper,int *n)
/* obtains coefficient matrices for imposing monotoicity (and optionally bounds) on a 
   cubic regression spline with n knots located as specified in xd. 
   
   control indicates type of constraints:
   up=control/4 - 0 for decrease, 1 for increase
   lo=(control-up*4)/2 - 1 for lower bound, 0 no lower bound
   hi=(control-up*4-lo*2) - 1 for upper bound, 0 no upper bound
   control = 4*up+2*lo+hi
   
   lower and upper are the bounds to impose (ignored if control doesn't
   indicate that they should be used).

   Ad will have 4(n-1)+lo+hi rows and n columns
   bd will have 4(n-1)+lo+hi rows

*/ 

{ int i;
  matrix x,A,b;
  x=initmat((long)*n,1L);
  for (i=0;i<x.r;i++) x.V[i]=xd[i];
  MonoCon(&A,&b,&x,*control,*lower,*upper); 
  RArrayFromMatrix(Ad,A.r,&A);
  RArrayFromMatrix(bd,b.r,&b);
 
  freemat(x);freemat(A);freemat(b);  
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif 
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
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif
}


void RGAMsetup(double *Xd,double *Cd,double *Sd,double *UZd,double *Xud,int *xu,double *xpd,
             int *offd,double *xd,int  *md,int  *nd,int *dfd,int *nsdfd, int *dim, int *s_type, 
             int *p_order,int *by_exists,double *byd,double *knots,int *n_knots)

/* Interface routine to GAMsetup from R (and hopefully Splus)
   The arrays pointed to by all arguments must be initialised to the
   correct size prior to calling this routine.
   Inputs are:
   1. x - in R terms x[i,j] is the jth obs of the ith covariate (corresponding
          to the jth datapoint being modelled i.e. y[j]). Note however that the
          first nsdf x[i]'s are the first nsdf columns of the design matrix
          corresponding to the parametric part of the model (including the intercept); 
          the remaining x[i]'s are covariates that will be treated using smooths.
          Dimension (m+nsdf) by n.
   2. nsdf - the number of parametric terms (including the mean/intercept term)
   3. df - m terms; df[i] is the max. d.f. for the ith smooth
   4. m - number of smooths
   5. n - number of data to be modelled.
   6. dim - dim[i] is the dimension of smooth term i
   7. s_type - s_type[i] is the smooth type for smooth i - currently only 0 
               for crs allowed
   8. p_order - p_order[i] is the order of the penalty for the ith term
   9. by_exists - 0/1 m-array indicating whether or not there is a "by" variable
      for each smooth.
   10. by - in R terms by[i,j] is the jth observation for the ith existing by-var
   11. knots - a densely packed array of covariate values specifying knot locations
               sub-arrays for each covariate are packed within this 1-d array
               in order corresponding to the ordering for the smooth related parts of x
   12. n_knots - array of number of knots per smooth.     
   Ouputs are:
   1. X - the n by q design matrix where q = \sum_i df[i]+nsdf
   2. xp - xp[i][j] is the jth knot position for the ith smooth. Dimension must
           be set to m by (max_i df[i]).
   3. S - array containing smoothness matrices. let start_k=sum_{i=0}^{k-1} S[i].r*S[i].c
          S[k].M[i][j]= S[start_k+i+S[k].r*j] - see RPackSarray() above
   4. UZd - Array Containing matrices UZ[k] such that \delta = UZ[k] \delta_k for 
           smooth k.  i.e. UZ[k] is basis of truncated parameter space.
           UZ[k].M[i][j] is row i, col j of kth basis. These matrices are actually
           stored in UZd in a packed form as follows:
           let start[k]=\sum_{k=0}^{k-1} UZ[k].c*UZ[k].r if k>0 and start[0]=0 
           UZ[k].M[i][j]=UZd[start[k]+i+UZ[k].r*j]
   5. Xud - array containing matrices the rows of which give all the unique covariate 
            combinations for each term. row i col j of the kth matrix is given by:
            Xud[start[k]+i+j*xu[k]] 
            where start[k+1]=start[k]+ xu[k]*dim[k], but start[0]=0
            xu[k] is number of unique covariate combinations, dim[k] is number of 
            covariates.      
   6. xu - array of numbers of unique cov.values for each s_type=1 term 
   4. C - Constraint matrix for problem. Input dimension h by q (see 1), where
          h = sum_i (1-s_type[i]).
   5. off - array of offsets locating S[k]'s within overall constraint matrix
            S[k][0][0] goes in row off[k] col off[k]

*/

{ matrix X,C,*S,*xp,*UZ,*Xu;
  long *off,*df,dumb;
  int m,n,nsdf,i,j,k,nx,k1,*M,nby;
  double **x,**by=(double **)NULL; /* avoid compiler warning */
  /* setup x[][], df[], off[], m, n, nsdf, S[], xp[], for calling GAMsetup.
     X, C, S[i]'s and xp[i]'s are initialised in the function */
  m= *md;
  n= *nd;
  nsdf= *nsdfd;
  /* deal with any "by" variables */
  nby=0;for (i=0;i<m;i++) nby+=by_exists[i];
  if (nby) by=(double **)calloc((size_t)nby,sizeof(double *));
  for (i=0;i<nby;i++) 
  { by[i]=(double *)calloc((size_t)n,sizeof(double));
    for (j=0;j<n;j++) by[i][j]=byd[i+nby*j];
  }
  /* and now the covariates + parametric design matrix columns */ 
  nx=nsdf;for (i=0;i<m;i++) nx+=dim[i];
  x=(double **)calloc((size_t)nx,sizeof(double *));
  for (i=0;i<nx;i++) 
  { x[i]=(double *)calloc((size_t)n,sizeof(double));
    for (j=0;j<n;j++) x[i][j]=xd[i+nx*j]; /* loading data passed by R into x[][] */
  }
  if (m>0)
  { S=(matrix *)calloc((size_t)m,sizeof(matrix));
    xp=(matrix *)calloc((size_t)m,sizeof(matrix));
    UZ=(matrix *)calloc((size_t)m,sizeof(matrix));
    Xu=(matrix *)calloc((size_t)m,sizeof(matrix));
    M=(int *)calloc((size_t) m,sizeof(int));
    for (i=0;i<m;i++) M[i]=null_space_dimension(dim[i],p_order[i]);  
    df=(long *)calloc((size_t)m,sizeof(long));
    for (i=0;i<m;i++) df[i]=(long)dfd[i];
    off=(long *)calloc((size_t)m,sizeof(long));
  } else {S=xp=UZ=Xu=&X;off=df=&dumb;M=&m;} /* purely to avoid compiler warnings */
  /* now run GAMsetup to get X, C, S[],UZ[],Xu[], xp, off */
  GAMsetup(&X,&C,S,UZ,Xu,xp,off,x,by,knots,m,n,dim,s_type,df,p_order,by_exists,n_knots,nsdf,0); /* Initializes X and C */
  /* unload returned matrices into R arrays: X, C, S[k], xp */
  RArrayFromMatrix(Xd,(long)n,&X);
  RArrayFromMatrix(Cd,(long)m,&C);
  /*mdf=0;for (i=0;i<m;i++) if (mdf<df[i]) mdf=df[i];*/
  /*m1=0; for (i=0;i<m;i++) if (s_type[i]) m1++;*/   /* number of tprs terms */
  /*mdf1=0;for (i=0;i<m;i++) if (s_type[i]) if (mdf1<df[i]) mdf1=df[i];*/ /* max df over tprs terms */
  /*mdim=0;for (i=0;i<m;i++) if (dim[i]>mdim) mdim=dim[i];*/   /* max dimension of terms   */
 
  RPackSarray(m,S,Sd);
  RPackSarray(m,UZ,UZd);
  RPackSarray(m,Xu,Xud);
  k1=0;  /*counter for tprs terms */  
  for (k=0;k<m;k++)
  { /*for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++)
      Sd[k+i*m+j*m*mdf]=S[k].M[i][j]; OLD CODE INEFFICIENT - AS IS CODE BELOW!!*/
    if (s_type[k]==1) /* then save the basis UZ and reduced covariate set Xu */
    { /*for (i=0;i<UZ[k].r;i++) for (j=0;j<UZ[k].c;j++)
	UZd[k1+i*m1+j*(n+maxM)*m1]=UZ[k].M[i][j];*/
    
      /*for (i=0;i<Xu[k].r;i++) for (j=0;j<Xu[k].c;j++)
	Xud[k1+i*m1+j*n*m1]=Xu[k].M[i][j];*/
     
      xu[k1]=(int)Xu[k].r;k1++;
    } else
    {  for (i=0;i<df[k];i++) xpd[k+m*i]=xp[k].V[i];
    } 
   
  } 
  for (i=0;i<m;i++) offd[i]=(int)off[i];
  /* tidy up */
  freemat(X);if (C.r) freemat(C);
  if (m>0)
  { for (i=0;i<m;i++) if (s_type[i]) { freemat(UZ[i]); freemat(Xu[i]);}
    free(UZ);free(Xu);
    for (i=0;i<m;i++) {freemat(S[i]); if (s_type[i]==0) freemat(xp[i]);}
    free(S);free(xp);free(off);free(df);free(M);
  }
  for (i=0;i<nx;i++) free(x[i]);free(x);
  for (i=0;i<nby;i++) free(by[i]); if (nby) free(by);
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif 
}

void gam_map(matrix tm, matrix *t, double *x,double *by,matrix *UZ,matrix *Xu,int *s_type, int *dim, int *p_order,
             int *by_exists,int m,int nsdf, int kill)

/* Consider a gam with a set of parametric terms given by the first nsdf terms of 
   vector x and m smooth terms each of dimension dim[i]. 

   The splines are one of two types:
   s_type[i]==0 => a cubic regression spline  parameterized by the function values at the knots, 
   s_type[i]==1 => thin plate regression spline                   

   the model parameters are assumed arranged: constant, parametric terms, spline terms_1, spline terms_2 etc...

   each smooth may also be multiplied by a "by"-variable... see by_exists.

   If p is the parameter vector then this routine returns tm, such that tm'p = f(x).
   
   To save time this routine saves some matrices statically - to free them, and set 
   things up for a new knot sequence call with kill==1.

*/


{ static matrix *D,b;static char first=1;
  static int terms=0,bl=0;
  matrix h,p;
  double xx,by_mult,**UZjM,*bV,*tmVk;
  int offset,by_counter=0;
  long i,j,k,l;
  
  if (kill)
  { first=1;
    for (i=0;i<terms;i++) 
    if (s_type[i]==0) 
    { freemat(D[i]);
    }
    free(D);
    terms=0;
    tps_g(UZ,&p,x,0,0,&h,1); /* clear up here */ 
    if (bl) freemat(b);
    return;
  }
  if (first)
  { first=0;
    D=(matrix *)calloc((size_t)m,sizeof(matrix));
    bl=0;
    for (i=0;i<m;i++)
    if (s_type[i]==0)
    { h=initmat(t[i].r-1,1L);
      for (j=0L;j<t[i].r-1;j++) h.V[j]=t[i].V[j+1]-t[i].V[j];
      D[i]=getD(h,0); 
      freemat(h);
    } else
    { if (UZ[i].r>bl) bl=UZ[i].r;}
    if (bl) b=initmat((long)bl,1L);
    terms=m;
  }
  /* deal with the  parametric terms first.... */
  for (k=0;k<nsdf;k++) tm.V[k]=x[k];
  /* now the splines..... */
  p.r=0;
  offset=nsdf; /* how far through x's we are */
  for (j=0;j<m;j++) 
  { if (by_exists[j]) { by_mult=by[by_counter];by_counter++;} else by_mult=1.0;
    if (s_type[j]==0) /* it's a cubic regression spline  */
    { xx=x[offset];
      i=0L;while((xx>t[j].V[i+1])&&(i<t[j].r-2)) i++; 
      for (l=0;l<t[j].r;l++)
      { tm.V[k]=D[j].M[i][l]*d0(t[j].V[i],t[j].V[i+1],xx)+
	          D[j].M[i+1][l]*d1(t[j].V[i],t[j].V[i+1],xx);
        if (l==i) tm.V[k]+=b0(t[j].V[i],t[j].V[i+1],xx);
        if (l==(i+1)) tm.V[k]+=b1(t[j].V[i],t[j].V[i+1],xx);
        tm.V[k]*=by_mult;
        k++;
      }
      offset++;
    }
    else  /* it's a thin plate regression spline  */
    { if (by_mult==0.0)         /* then don't waste flops on calculating stuff that will only be zeroed */
      { for (i=0;i<UZ[j].c;i++) {tm.V[k]=0.0;k++;}
      } else                    /* proceed as normal */
      { b.r=UZ[j].r; /*  - don't forget null space dimension */
        tps_g(Xu+j,&p,x+offset,dim[j],p_order[j],&b,1);
        UZjM=UZ[j].M;bV=b.V;     
        for (i=0;i<UZ[j].c;i++) 
        { tmVk=tm.V+k;*tmVk=0.0; 
          for (l=0;l<b.r;l++) *tmVk += bV[l]*UZjM[l][i]; /* forming b'UZ */
          *tmVk *= by_mult;
          k++;
        }
      }
      offset+=dim[j];
    }
  } 
 
}

void RGAMpredict(int *xu,double *Xud,double *UZd,double *xpd,int *nsdf,int *dim,int *s_type,int *df,int *p_order,int *m,
                int *n,double *xd,int *np,double *p, double *Vpd,double *etad,double *sed,double *X,int *control,
                int *by_exists, double *by)

/* Routine for prediction from GAMs made up of penalized cubic regression splines.
   
   xu[i]    is the number of unique combinations of covariates for the ith tprs term
   Xud -    array containing matrices the rows of which give all the unique covariate 
            combinations for each term. row i col j of the kth matrix is given by:
            Xud[start[k]+i+j*xu[k]] 
            where start[k+1]=start[k]+ xu[k]*dim[k], but start[0]=0
            xu[k] is number of unique covariate combinations, dim[k] is number of 
            covariates. 
   UZd      Array Containing matrices UZ[k] such that \delta = UZ[k] \delta_k for 
            smooth k.  i.e. UZ[k] is basis of truncated parameter space.
            UZ[k].M[i][j] is row i, col j of kth basis. These matrices are actually
            stored in UZd in a packed form as follows:
            let start[k]=\sum_{k=0}^{k-1} UZ[k].c*UZ[k].r if k>0 and start[0]=0 
            UZ[k].M[i][j]=UZd[start[k]+i+UZ[k].r*j]   
    
   xp[i][j] is the jth knot position for the ith smooth (xp[i] dimension must be set to 
            max_i(df[i])... only applicable for s_type=0 terms)
   nsdf     is the number of non spline terms including the constant.
   dim[i]   is the dimension of the ith term
   s_type[i] is the type of the ith term - cubic regression spline or tprs
   p_order[i] is the order of the ith penalty
   df[i]    is the number of parameters for the ith spline
   m        is the number of splines
   n        is the number of data that the model were originally fitted to - needed for unpacking 
            UZ[i]'s
   x[i][j]  is the ith observation of jth covariate, for which predictions are required - 
            first nsdf are none spline terms
   np       is the number of predictions required (number of obs. per cov.)
   p[i]     is the ith parameter: 0 to nsdf-1 are the none spline terms
                                  nsdf to nsdf+df[0]-1 is the first spline
                                  nsdf+df[i-1] to nsdf+df[i]-1 is the ith spline 
   Vp[i][j] is covariance of p[i] and p[j]
   by_exists[i] is 1 if there is a by-variable to multiply the whole of smooth term i - 0 if not
   by       in R terms by[i,j] is jth value for ith by-variable  
 
   eta is the o/p linear predictor
   se is the o/p standard error on l.p.
   X is the o/p matrix mapping the parameters to the l.p. vector
   

   control sets the type of output: 
           0: l.p. no s.e.
           1: l.p. with s.e.
           2: predictor for each term, no s.e.
           3: predictor for each term with s.e.
           4: matrix mapping parameters to l.p.

   note  that the dimensions of X, eta and se must be appropriate for the control option selected.
*/
           
{ matrix tm,*xp,Vp,eta,se,*UZ,*Xu;
  int i,j,k,nb,kk,l,nx,*M,nby; 
  double **x,z,Vt;
  char info[200];
  /* perform bits of unpacking ...... */
  nby=0;for (i=0;i< *m;i++) nby+=by_exists[i];
  /* by[nby*j+i] is jth entry for ith by variable */
  x=(double **)calloc((size_t) *np,sizeof(double *));
  nx= *nsdf;for (i=0;i<*m;i++) nx+=dim[i]; /* the number of covariates (inc. none spline terms) */
  for (i=0;i<*np;i++) 
  { x[i]=(double *)calloc((size_t)nx,sizeof(double));
    for (j=0;j<nx;j++) 
    { x[i][j]=xd[j+nx*i]; /* loading data passed by R into x[][] */
    } 
  }
  /* now unpack any UZ[i], Xu[i], terms */
  if (*m>0)
  { M=(int *)calloc((size_t) *m,sizeof(int));
    for (i=0;i< *m;i++) M[i]=null_space_dimension(dim[i],p_order[i]);  
    UZ=(matrix *)calloc((size_t) *m,sizeof(matrix));
    Xu=(matrix *)calloc((size_t) *m,sizeof(matrix)); 
    /* m1=0;for (i=0;i< *m;i++) m1+=s_type[i];*/ /* number of tprs terms */
  } else { Xu=UZ=&Vp;M=&i;} /* just avoids spurious fussy compiler warnings */
  kk=0; /* counter for tprs terms */
 
  for (k=0;k<*m;k++)
  if (s_type[k]==1)
  { UZ[k]=initmat((long)xu[kk]+M[k],(long)df[k]);
    /*for (i=0;i<UZ[k].r;i++) for (j=0;j<UZ[k].c;j++)
      UZ[k].M[i][j]=UZd[kk+i * m1+j * (*n+maxM) * m1];*/
    Xu[k]=initmat((long)xu[kk],(long)dim[k]);
    /*for (i=0;i<Xu[k].r;i++) for (j=0;j<Xu[k].c;j++)
      Xu[k].M[i][j]= Xud[kk+i * m1+j * *n * m1];*/  
    kk++; 
  } else 
  { Xu[k].r=Xu[k].c=UZ[k].r=UZ[k].c=0; /* needs to be set to allow correct unpacking */ 
  }
  RUnpackSarray(*m,Xu,Xud);
  RUnpackSarray(*m,UZ,UZd);
  nb =  *nsdf;  /* number of parameters */
  for (i=0;i < *m;i++) nb += df[i];
  Vp=Rmatrix(Vpd,(long)nb,(long)nb); /* param cov matrix */
 
  /* initialise a couple of storage matrices */
  if (*control<2) k=1;else k=  *nsdf + *m;  /* don't allocate more than is needed */
  eta=initmat((long)k,(long)*np); 
  se=initmat((long)k,(long)*np);
  /* need to unpack xpd into xp here ..... */
  if (*m>0) xp=(matrix *)calloc((size_t) *m,sizeof(matrix));
  else xp=&Vp; /* avoids spurious compiler warning */
  for (k=0;k< *m;k++)
  if (s_type[k]==0)
  { xp[k]=initmat((long)df[k],1L);
    for (i=0;i<df[k];i++) xp[k].V[i]=xpd[k+ *m *i];
  }  
  tm=initmat((long)nb,1L);
 
  for (k=0;k< *np;k++) /* loop through the predictions */
  { gam_map(tm,xp,x[k],by+nby*k,UZ,Xu,s_type,dim,p_order,by_exists,*m,*nsdf,0);
    if (*control==4) /* matrix X such that eta=Xp, is required */
    { for (j=0;j<tm.r;j++) X[k + *np * j]=tm.V[j];
    }
    else if (*control<2) /* then linear predictor required */
    { z=0.0;
     for (i=0;i<tm.r;i++) 
     { z+=tm.V[i]*p[i];eta.M[0][k]=z;
       
     }
      if (*control==1) /* get s.e. */
      { z=0.0;
        for (i=0;i<tm.r;i++) 
        { Vt=0.0; for (j=0;j<tm.r;j++) Vt+=Vp.M[i][j]*tm.V[j];
          z+=tm.V[i]*Vt; 
        }
        se.M[0][k]=sqrt(z);
      }
    } else /* individual predictors wanted */
      { /* parametric terms..... */
      for (i=0;i< *nsdf;i++) 
      { eta.M[i][k]=p[i]*x[k][i];
      if (*control==3) se.M[i][k]=sqrt(x[k][i]*x[k][i]*Vp.M[i][i]); 
      }
      /* and now the splines ..... */
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
  /* convert results for o/p */
 
  gam_map(tm,xp,x[0],by,UZ,Xu,s_type,dim,p_order,by_exists,*m,*nsdf,1);  /* free memory in gam_map() */
 
  freemat(tm);
  if (*control!=4)
  { RArrayFromMatrix(etad,eta.r,&eta);
    if ((*control)%2) RArrayFromMatrix(sed,se.r,&se);
  }
  /* tidy up.... */
  for (i=0;i< *np;i++) free(x[i]);free(x);  
  if (*m>0)
  { for (i=0;i< *m;i++) if (s_type[i]==1) { freemat(UZ[i]);freemat(Xu[i]);} 
    free(UZ);free(Xu);free(M);
  }
  freemat(Vp);freemat(eta);freemat(se);

  for (k=0;k< *m;k++) if (s_type[k]==0) freemat(xp[k]); if (*m>0) free(xp);
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif 
}    


void RprintM(matrix *A)

{ 
#ifdef WINDOWS
#else
int i,j;
  if (A->c==1L) 
  { for (i=0;i<A->r;i++) 
   Rprintf("%8.3g ",A->V[i]);Rprintf("\n");
  } 
else
 for (i=0;i<A->r;i++)
     { for (j=0;j<A->c;j++) Rprintf("%8.3g ",A->M[i][j]);Rprintf("\n");}
#endif
}

void  RPCLS(double *Xd,double *pd,double *yd, double *wd,double *Aind,double *bd,double *Afd,double *Hd,double *Sd,
      int *off,int *dim,double *theta, int *m,int *nar)

/* Interface routine for PCLS the constrained penalized weighted least squares solver.
   nar is an array of dimensions. Let:
   n=nar[0] - number of data
   np=nar[1] - number of parameters
   nai=nar[2] - number of inequality constraints
   naf=nar[3] - number of fixed constraints
   getH=nar[4] - 0 for no hat matrix, 1 to produce one. 
   
   Problem to be solved is:

   minimise      ||W^0.5 (y - Xp)||^2 + p'Bp
   subject to    Ain p >= b  & Af p = "constant"

   where B = \sum_{i=1}^m \theta_i S_i and W=diag(w)

   - in fact S_i are not stored whole - rather the smallest non-zero sub-matrix of each S_i is 
   stored in a densely packed form in S[]: see routines RpackSarray() and RUnpackSarray() for 
   details of the sub-matrix packing. off[i],off[i] is the location within the full S_i to
   insert the sub-matrix actually stored which is of dimension dim[i] by dim[i].

   W = diag(w) 

   on exit p contains the best fit parameter vector. 

*/
{ matrix y,X,p,w,Ain,Af,b,H,*S;
  int n,np,i,*active;
 
  np=nar[1];n=nar[0];
  /* unpack from R into matrices */
  X=Rmatrix(Xd,(long)n,(long)np);
  p=Rmatrix(pd,(long)np,1L);
  y=Rmatrix(yd,(long)n,1L);
  w=Rmatrix(wd,(long)n,1L);
  if (nar[2]>0) Ain=Rmatrix(Aind,(long)nar[2],(long)np); else Ain.r=0L;
  if (nar[3]>0) Af=Rmatrix(Afd,(long)nar[3],(long)np); else Af.r=0L;
  if (nar[2]>0) b=Rmatrix(bd,(long)nar[2],1L);else b.r=0L;
 
  if (*m) S=(matrix *)calloc((size_t) *m,sizeof(matrix));
  else S=&H; /* avoid spurious compiler warning */
  for (i=0;i< *m;i++) S[i]=initmat((long)dim[i],(long)dim[i]);
  RUnpackSarray(*m,S,Sd);
  
  if (nar[4]) H=initmat(y.r,y.r); else H.r=H.c=0L;
  active=(int *)calloc((size_t)p.r+1,sizeof(int)); /* array for active constraints at best fit active[0] will be  number of them */
  /* call routine that actually does the work */
 
  PCLS(&X,&p,&y,&w,&Ain,&b,&Af,&H,S,off,theta,*m,active);

  /* copy results back into R arrays */ 
  for (i=0;i<p.r;i++) pd[i]=p.V[i];
 
  if (H.r) RArrayFromMatrix(Hd,H.r,&H);
  /* clear up .... */
  free(active);
 
  for (i=0;i< *m;i++) freemat(S[i]);
  if (*m) free(S);
 
  freemat(X);freemat(p);freemat(y);freemat(w);
  if (H.r) freemat(H);
  if (Ain.r) freemat(Ain);
  if (Af.r) freemat(Af);
  if (b.r) freemat(b);
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif
}

/*********************************************************************************************/
/* Bug fix and revision record:

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
5. 30/5/01: GAMsetup re-organised to ease introduction of new tprs basis and multi-dimensional 
   smooths
6. 10/2001: b0,b1,d0,d1 modified so that extrapolation is linear beyond ends of spline as
   it should be for a natural spline.
7. 31/10/01: mgcv.c covariance and edf calculations made more robust - for poorly conditioned 
   cases choleski can fail in calculation of covariance matrix: in these cases use svd instead.
8. 2/11/01: RGAMpredict, RGAMsetup, GAMsetup, and gam_map modified to take array of penalty
   orders, p_order. This allows user explicit control of the order of penalty in spline terms, 
   while still supporting autoselection when p_order[i]==0. 
9. 9/11/01: RGAMpredict modified to allow 5th control option, which returns a matrix mapping 
            params to l.p. vector
10. 9/11/01: New routine RPackSArray and RUnpackSarray so that storage of arrays of penalty 
             matrices is not so wasteful.
11. 12/11/01: UZ and Xu now packed efficiently using above 2 routines.
12. 12/11/01: Routine RPCLS added for solving linearly constrained penalized least squares problems 
              by quadratic programming
13. 13/11/01: Routine RMonoCon added for finding monotonic constraint matrices for cubic regression 
              splines.
14. 5/2/02: GAMsetup and RGAMsetup modified to deal with "by" variables - covariates that multiply a 
            whole smooth term. The centering conditions have not been changed.
15. 6/9/02: Slight modification to gam_map() - terms are not calculated if corresponding by variable 
            is zero. This can save flops in fairly advanced use (e.g. posum package)

*/


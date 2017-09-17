/* Copyright (C) 1991-2005 Simon N. Wood  simon.wood@r-project.org

14/9/17 --- cleaned out routines no longer needed

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
USA.*/

/* Routines for basic matrix manipulation creation, destruction and file i/o
   for matrices. See end of file for update log */

#include <R.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "mgcv.h"
#include "matrix.h"
#include "general.h"
#define RANGECHECK
#define PAD 1

#define ROUND(a) ((a)-(int)floor(a)>0.5) ? ((int)floor(a)+1):((int)floor(a))


matrix null_mat; /* matrix for passing when you don't actually need to */
#define PADCON (-1.234565433647588392902028934e270)

/* counter for memory used */


long memused=0L,matrallocd=0L;

/* the routines */

struct mrec
{ matrix mat;
  struct mrec *fp,*bp;
};
typedef struct mrec MREC;

matrix null_mat;
MREC *top,*bottom;

matrix initmat(int rows,int cols)
/* Don't alter this without altering freemat() as well !! */
{ matrix A;int i,j,pad;
#ifdef RANGECHECK
  pad=PAD;
#else
  pad=0;
#endif
  A.vec=0;
  A.M=(double **)CALLOC((size_t)(rows+2*pad),sizeof(double *));
  if ((cols==1)||(rows==1))
  { if (A.M)
    A.M[0]=(double *)CALLOC((size_t)(cols*rows+2*pad),sizeof(double));
    for (i=1;i<rows+2*pad;i++) A.M[i]=A.M[0]+i*cols;
    A.vec=1;
  } else
  { if (A.M)
    for (i=0;i<rows+2*pad;i++)
    A.M[i]=(double *)CALLOC((size_t)(cols+2*pad),sizeof(double));
  }
  A.mem=(long)(rows*cols*sizeof(double));
  memused+=A.mem;matrallocd++;
  A.original_r=A.r=rows;A.original_c=A.c=cols;
  if (((!A.M)||(!A.M[rows-1+2*pad]))&&(rows*cols>0))
  { error(_("Failed to initialize memory for matrix."));}
  if (pad)  /* This lot is debugging code that checks out matrix errors
		       on allocation and release */
  { if (A.vec)
    { A.V=A.M[0];for (i=0;i<pad;i++) { A.V[i]=PADCON;A.V[i+pad+A.r*A.c]=PADCON;}
    } else
    { for (i=0;i<A.r+2*pad;i++)
      { for (j=0;j<pad;j++) A.M[i][j]=PADCON;
	     for (j=A.c+pad;j<A.c+2*pad;j++) A.M[i][j]=PADCON;
      }
      for (i=0;i<A.c+2*pad;i++)
      { for (j=0;j<pad;j++) A.M[j][i]=PADCON;
	     for (j=A.r+pad;j<A.r+2*pad;j++) A.M[j][i]=PADCON;
      }
    }
    for (i=0;i<A.r+2*pad;i++)
    for (j=0;j<pad;j++) A.M[i]++;  /* shifting pointers forward past padding */
    if (!A.vec) for (j=0;j<pad;j++) A.M++;
    A.V=A.M[0];
  /* putting a record of the matrix on the linked list of all extant matrices */
    if (matrallocd==1) /*new list*/
    { top=bottom=(MREC *)CALLOC(1,sizeof(MREC));
      bottom->mat=top->mat=A;top->bp=bottom;bottom->fp=top;
    } else  /* expanding the linked list by one */
    { top->fp=(MREC *)CALLOC(1,sizeof(MREC));
      top->fp->mat=A;top->fp->bp=top;top=top->fp; /* crystal clear, no? */
    }
  }
  A.V=A.M[0];/* This allows vectors to be accessed using A.V[i] */
  return(A);
} /* initmat */

matrix initvec(int rows)

{ return(initmat(1,rows));}

void freemat(matrix A)

{ int i,j,pad;int ok=1;
  MREC *delet;
#ifdef RANGECHECK
  pad=PAD;
#else
  pad=0;
#endif
/*  if (A.original_r*A.original_c!=0L) */
  { if (pad)
    { if (A.vec)
      { for (i=-pad;i<0;i++)
	     if ((A.V[i]!=PADCON)||(A.V[i+A.original_r*A.original_c+pad]!=PADCON))
	     ok=0;
      } else
      { for (i=-pad;i<A.original_r+pad;i++)
	     { for (j=A.original_c;j<A.original_c+pad;j++) if (A.M[i][j]!=PADCON) ok=0;
	       for (j=-pad;j<0;j++) if (A.M[i][j]!=PADCON) ok=0;
	     }
	     for (i=-pad;i<A.original_c+pad;i++)
	     { for (j=A.original_r;j<A.original_r+pad;j++) if (A.M[j][i]!=PADCON) ok=0;
	       for (j=-pad;j<0;j++) if (A.M[j][i]!=PADCON) ok=0;
	     }
      }
      if (!ok)
      { error(_("An out of bound write to matrix has occurred!"),1);
      }
      /* find the matrix being deleted in the linked list of extant matrices */
      i=0;delet=bottom;
      while ((i<matrallocd)&&(delet->mat.M!=A.M)) { i++;delet=delet->fp;}
      if (i==matrallocd)
      { error(_("INTEGRITY PROBLEM in the extant matrix list."));
      } else
      { if (i)
	     delet->bp->fp=delet->fp;
	     else bottom=delet->fp;
	     if (i!=matrallocd-1)
	     delet->fp->bp=delet->bp;
	     else top=delet->bp;
	     FREE(delet);
      }
      /* repositioning pointers so that what was allocated gets freed */
      if (!A.vec) for (i=0;i<pad;i++) A.M--;
      for (i=0;i<A.original_r+2*pad;i++)
      for (j=0;j<pad;j++) A.M[i]--;
    }
    if (A.vec) FREE(A.M[0]); else
    for (i=0;i<A.original_r+2*pad;i++) if (A.M[i]) FREE(A.M[i]);
    if (A.M) FREE(A.M);
    memused -= A.mem;matrallocd--;
  }
} /* freemat */

void matrixintegritycheck()

/* iff RANGECHECK is defined above then you can call this routine to check
   on the integrity of the matrix system. The routine looks for writing out
   of bounds from the matrix */

{ MREC *B;
  int ok=1,pad=PAD,i,j,k=0;
  matrix A;
#ifndef RANGECHECK
  error(_("You are trying to check matrix integrity without defining RANGECHECK."));
#endif
  B=bottom;
  while (k<matrallocd) {
    A=B->mat;
    if (A.vec) {
      for (i=-pad;i<0;i++)
      if ((A.V[i]!=PADCON)||(A.V[i+A.original_r*A.original_c+pad]!=PADCON))
      ok=0;
    } else {
      for (i=-pad;i<A.original_r+pad;i++) {
        for (j=A.original_c;j<A.original_c+pad;j++) if (A.M[i][j]!=PADCON) ok=0;
	for (j=-pad;j<0;j++) if (A.M[i][j]!=PADCON) ok=0;
      }
      for (i=-pad;i<A.original_c+pad;i++) {
	for (j=A.original_r;j<A.original_r+pad;j++) if (A.M[j][i]!=PADCON) ok=0;
	for (j=-pad;j<0;j++) if (A.M[j][i]!=PADCON) ok=0;
      }
    }
    if (!ok) {
      error(_("An out of bound write to matrix has occurred!"));
    }
    k++;B=B->fp;
  }
} /* matrixintegritycheck */




void vmult(matrix *A,matrix *b,matrix *c,int t)

/* fast multiplication of vector by matrix c=Ab if t==0 c=A'b otherwise*/

{ double **AM,*bV,*cV,*p;
  int i,j,cr,br;
  cr=c->r;br=b->r;
  AM=A->M;bV=b->V;cV=c->V;
  if (t) /* then A transposed */
  for (i=0;i<cr;i++)
  { *cV=0.0;
    for (j=0;j<br;j++) *cV += AM[j][i]*bV[j];
    cV++;
  } else
  for (i=0;i<cr;i++)
  { *cV=0.0;
    p=AM[i];
    for (j=0;j<br;j++)
    *cV += p[j]*bV[j];
    cV++;
  }
} /* vmult */

void mcopy(matrix *A,matrix *B)

/* copies A into B */

{ int Ac;
  double *pA,*pB,**AM,**BM;
  if (A->r>B->r||A->c>B->c) error(_("Target matrix too small in mcopy"));
  BM=B->M;Ac=A->c;
  for (AM=A->M;AM<A->M+A->r;AM++)
  { pB= *BM;
    for (pA= *AM;pA< *AM+Ac; pA++) *(pB++) = *pA;
    BM++;
  }
} /* mcopy */

void matmult(C,A,B,tA,tB) matrix C,A,B;int tA,tB;

/* Puts A*B in C. A will be transposed in this calculation if tA is not zero.
   B will be transposed if tB is not zero */

{ int i,j,k;
  double temp,*p,*p1,*p2,**CM,**AM,**BM;
  AM=A.M;BM=B.M;CM=C.M; /* Saves address calculation involved in C.M */
  if (tA)
  { if (tB)
    { if ((A.r!=B.c)||(A.c!=C.r)||(B.r!=C.c))
      { error(_("Incompatible matrices in matmult."));}
      for (i=0;i<A.c;i++) for (j=0;j<B.r;j++)
      { p2=CM[i]+j;(*p2)=0.0;p=BM[j];
	for (k=0;k<A.r;k++)
	(*p2)+=AM[k][i]*(*p++);
      }
    } else
    { if ((A.r!=B.r)||(A.c!=C.r)||(B.c!=C.c))
      { error(_("Incompatible matrices in matmult."));}
      for (i=0;i<A.c;i++)
      for (p=CM[i];p<(CM[i]+C.c);p++)
      (*p)=0.0;
      for (k=0;k<A.r;k++) for (i=0;i<A.c;i++)
      { temp=AM[k][i];p1=BM[k];
	for (p=CM[i];p<(CM[i]+B.c);p++)
	(*p)+=temp*(*p1++);
      }
    }
  } else
  { if (tB)
    { if ((A.c!=B.c)||(A.r!=C.r)||(B.r!=C.c))
      { error(_("Incompatible matrices in matmult."));}
      for (i=0;i<A.r;i++) for (j=0;j<B.r;j++)
      { p2=CM[i]+j;*p2=0.0;p1=BM[j];
	for (p=AM[i];p<(AM[i]+A.c);p++)
	(*p2)+=(*p)*(*p1++);
      }
    } else
    { if ((A.c!=B.r)||(C.r!=A.r)||(C.c!=B.c))
      { error(_("Incompatible matrices in matmult."));}
      for (i=0;i<A.r;i++) for (p=CM[i];p<(CM[i]+B.c);p++) *p=0.0;
      for (k=0;k<A.c;k++) for (i=0;i<A.r;i++)
      { p1=BM[k];temp=AM[i][k];
	     for (p=CM[i];p<(CM[i]+B.c);p++)
	     (*p)+=temp*(*p1++);
      }
    }
  }
} /* matmult */



void invert(matrix *A)

/* Matrix inversion by Guass-Jordan Elimination with full pivoting.
   See "Numerical Recipes", and Burden and Faires "Numerical Analysis", for basis
   of method (but not actual code). This version written as part of elimination of 
   "Numerical Recipes" routines from my code. Tested against Numerical Recipes code 
   with a variety of random matrices - fine on accuracy and speed. 13/1/2000
 */

{ double **AM,*p,*p1,max,x;
  int *c,*rp,*cp,i,j,k,pr=0,pc=0,*d,cj,ck;
  if (A->r!=A->c) error(_("Attempt to invert() non-square matrix"));
  c=(int *)CALLOC((size_t)A->c,sizeof(int)); /* index of columns, used for column pivoting */
  d=(int *)CALLOC((size_t)A->c,sizeof(int));
  rp=(int *)CALLOC((size_t)A->c,sizeof(int)); /* row changes */
  cp=(int *)CALLOC((size_t)A->c,sizeof(int)); /* row changes */
  for (i=0;i<A->c;i++) { c[i]=i;d[i]=i;}
  AM=A->M;            /* saving adress calculations*/
  for (j=0;j<A->c;j++) /* loop through columns to be reduced */
  { max=0.0; 
    for (i=j;i<A->r;i++) /* loop through rows to search for pivot */
    { p=AM[i];
      for (k=j;k<A->c;k++) /* loop through cols to search for pivot */
      { x=p[c[k]];if (fabs(x)>max) { max=fabs(x);pr=i;pc=k;}}
    }
    /* now move pivot to element j,j */
    p=AM[j];AM[j]=AM[pr];AM[pr]=p; /* rows exchanged */
    k=c[j];c[j]=c[pc];c[pc]=k;   /* columns exchanged */
    rp[j]=pr;  /* stores row pivoted with */
    cp[j]=pc;  /* stores column pivoted with */
    cj=c[j]; /* save time */
    /* Now reduce the column */
    x=AM[j][cj];
    if (x==0.0) error(_("Singular Matrix passed to invert()"));
    for (p=AM[j];p<AM[j]+A->c;p++) *p/=x; /* divide row j by pivot element */
    AM[j][cj]=1.0/x;
    for (i=0;i<A->r;i++) /* work down rows eliminating column j */
    { p=AM[i];p1=AM[j];
      if (i!=j)
      { x = -p[cj]; /* multiplier for this row */
        for (k=0;k<j;k++) /* work through columns of the inverse matrix */
        { ck=c[k];p[ck]+=x*p1[ck];}
        p[cj]=x*p1[cj]; /* new column for inverse (entries in A implicitly zeroed except jj) */
        for (k=j+1;k<A->c;k++)   /* cols of A */
        { ck=c[k];p[ck]+=x*p1[ck];}
      }
    }
  } 
 
  for (i=A->r-1;i>=0;i--) /*work down through column re-ordering  */
  { if (cp[i]!=i)
    { p=AM[i];AM[i]=AM[cp[i]];AM[cp[i]]=p; /* row exchange */
    }
  }
  
  for (j=0;j<A->c-1;j++) /* implement column exchange */
  if (c[j]!=j)
  { if (c[j]<j) k=c[c[j]]; else k=c[j]; 
    for (i=0;i<A->r;i++)
    { p=AM[i];x=p[j];p[j]=p[k];p[k]=x;}  
    d[k]=d[j];d[j]=c[j];
    c[d[k]]=k;
  } 
  
  for (i=A->r-1;i>=0;i--) /* column exchange implied by row re-ordering */
  if (rp[i]!=i)
  { for (k=0;k<A->r;k++) 
    { p=AM[k];x=p[i];p[i]=p[rp[i]];p[rp[i]]=x;} /* column exchange  */
  }
   
  FREE(c);FREE(rp);FREE(cp);FREE(d);
} /* invert */



double dot(a,b) matrix a,b;

{ int i,k=0;double c=0.0,*p,*p1;
  if (a.vec) { p1=b.V;for (p=a.V;p<a.V+a.c*a.r;p++) c+=(*p)*(*p1++);}
  else
  for (i=0;i<a.r;i++) for (p=a.M[i];p<(a.M[i]+a.c);p++)
  { c+=(*p)*b.M[k/b.c][k%b.c];k++;}
  return(c);
} /* dot */




double enorm(d) matrix d;

/* Euclidian norm of vector d, or rms for matrix */

{ double e=0.0,m=0.0,y,*p;
  int i;
  if (d.vec) for (p=d.V;p<d.V+d.r*d.c;p++) { y=fabs(*p); if (y>m) m=y; }
  else for (i=0;i<d.r;i++) for (p=d.M[i];p<d.M[i]+d.c;p++) 
  { y=fabs(*p);if (y>m) m=y;}/* m=max(m,fabs(*p)); */
  if (!m) return(0.0);
  if (d.vec) for (p=d.V;p<d.V+d.r*d.c;p++)
  { y= *p / m; e+=y*y;} else
  for (i=0;i<d.r;i++) for (p=d.M[i];p<(d.M[i]+d.c);p++)
  { y= *p / m;e+=y*y;}
  e=sqrt(e)*m;
  return(e);
} /* enorm */





void householder(u,a,b,t1) matrix *u,a,b;int t1;

/* transforms a to b, iff they are of equal Euclidian length. u is the
   (t1+1) vector such that the full post multiplying householder matrix is
   H' = [ I - vv' ]   where v'  = [u',0] where 0 is a vector of zeroes.   */

{ int i;double v,*aV,*bV,*uV;
  aV=a.V;bV=b.V;uV=u->V;
  u->r=t1+1;
  for (i=0;i<u->r;i++) uV[i]=aV[i]-bV[i];
  v=enorm((*u))/sqrt(2.0);
  for (i=0;i<u->r;i++) uV[i]/=v;
} /* householder */



void Hmult(C,u) matrix C,u;

/* This routine is for post multiplication by Housholder matrices only */

{ double temp,*p,*p1,*uV,**CuM,**CM;
  int i,j;
  matrix Cu;
  Cu=initmat(C.r,u.c);
  uV=u.V;CuM=Cu.M;CM=C.M;
  for (i=0;i<Cu.r;i++)
  { p1=CuM[i];(*p1)=0.0;p=CM[i];
    for (j=0;j<u.r;j++) (*p1)+=(*p++)*uV[j];
  }
  for (i=0;i<Cu.r;i++)
  { temp=Cu.V[i];p=CM[i];
    for (j=0;j<u.r;j++)
    (*p++) -= temp*uV[j];
  }
  freemat(Cu);
} /* Hmult */

void HQmult(C,U,p,t) matrix C,U;int p,t;

/* This routine is for multiplying by householder matrices, stored as a
   series of vectors of the type u defined in routine householder above.
   If ui is ith row of U, then Hi=(I-ui ui') = Hi'. Let Q = H1 H2 H3 ....
   then Q' = .... H3 H2 H1.
   Returned in C: p==0,t==0 => CQ; 
                  p==0,t==1 => CQ'; 
                  p==1,t==0 => QC;
		          p==1,t==1 => Q'C
   
   NOTE that the routine expects C to be compatible with the Hi's - if
   this routine is being used for projection in and out of Null spaces, then
   make sure that C is appropriately packed with zeroes.

   If appropriate zero packing conventions have been used then OrthMult() is
   more efficient....
   */

{ double *u,*CuV,**CM;
  matrix Cu;
  int i,j,k;
  if (p) Cu=initmat(C.c,1);else Cu=initmat(C.r,1);
  CuV=Cu.V;CM=C.M;
  if (p)
  { if (t)
    { for (k=0;k<U.r;k++) /* loop through the householder matrices */
      { u=U.M[k];
        for (i=0;i<C.c;i++)
	     { CuV[i]=0.0;
          for (j=0;j<C.r;j++) CuV[i]+=CM[j][i]*u[j];
	     }
	     for (i=0;i<C.r;i++) for (j=0;j<C.c;j++) CM[i][j] -= CuV[j]*u[i];
      }
    }  else
    { for (k=U.r-1;k>=0;k--) /* loop through the householder matrices */
      { u=U.M[k];
	     for (i=0;i<C.c;i++)
	     { CuV[i]=0.0;
	       for (j=0;j<C.r;j++) CuV[i]+=CM[j][i]*u[j];
	     }
	     for (i=0;i<C.r;i++) for (j=0;j<C.c;j++) CM[i][j] -= CuV[j]*u[i];
      }
    }
  } else  /* postmultiplication */
  { if (t)
    { for (k=U.r-1;k>=0;k--) /* loop through the householder matrices */
      { u=U.M[k];
	     for (i=0;i<C.r;i++)
	     { CuV[i]=0.0;
	       for (j=0;j<C.c;j++) CuV[i]+=CM[i][j]*u[j];
	     }
	     for (i=0;i<C.r;i++) for (j=0;j<C.c;j++) CM[i][j] -= CuV[i]*u[j];
      }
      } else
      { for (k=0;k<U.r;k++) /* loop through the householder matrices */
      { u=U.M[k];
	     for (i=0;i<C.r;i++)
	     { CuV[i]=0.0;
	       for (j=0;j<C.c;j++) CuV[i]+=CM[i][j]*u[j];
	     }
	     for (i=0;i<C.r;i++) for (j=0;j<C.c;j++) CM[i][j] -= CuV[i]*u[j];
      }
    }
  }
  freemat(Cu);
} /* HQmult */


void QT(Q,A,fullQ) matrix Q,A;int fullQ;

/* Uses householder matrices to perform the factorization of A (nxm),n<=m: */
/*                  AQ=[0,T]  where Tij=0 if i+j<n for T an (nxn) matrix   */
/*                    i.e. T reverse lower triangular - top left is zero   */
/* if fullQ!=0 then Q is formed explicitly (what a waste!), otherwise the  */
/* first A.r rows of Q will be used to store the vectors u used to define  */
/* the householder matrix. These can be used to multiply some other matrix */
/* using the routine HQmult. Revised 13/1/2000 - more efficient, over/under*/
/* flow protected, cancellation error free, but still behaves as before.   */
/* Tested using variety of random matrices.                                */

{ int i,j,Ar,Ac,k;
  double lsq,*p,*p1,**QM,**AM,g,x,m;
  QM=Q.M;AM=A.M;Ar=A.r;Ac=A.c;
  if (fullQ) for (i=0;i<Ac;i++) 
  { p=QM[i];
    for (j=0;j<Ac;j++) if (i==j)
    p[j]=1.0; else p[j]=0.0;
  }
  if (Ar>0)
  { for (i=0;i<Ar;i++)
    { /* rotate elements 0 to A.c-i-1 row i of A into element A.c-i-1 of that row  */
      p=AM[i];
      m=0.0;for (j=0;j<Ac-i;j++) { x=p[j];x=fabs(x);if (x>m) m=x;} /* scale factor */
      if (m) for (j=0;j<Ac-i;j++) p[j]/=m; /* avoid over/underflow */
      lsq=0.0;for (j=0;j<Ac-i;j++) lsq+=p[j]*p[j];
      lsq=sqrt(lsq);
      if (p[Ac-i-1]<0.0) lsq= -lsq;
      p[Ac-i-1]+=lsq;
      if (lsq)
      g=1/(lsq*p[Ac-i-1]); /* multiplier for HH rotation (I-g*uu') */
      else g=0.0;
      lsq*=m; /* Element to end up on A.M[i][A.c-i-1] */
      for (j=i+1;j<Ar;j++) /* Apply rotation down through the rows of A */
      { x=0.0;p1=AM[j];
        for (k=0;k<Ac-i;k++) x+=p[k]*p1[k];
        x*=g;
        for (k=0;k<Ac-i;k++) p1[k] += -x*p[k];
      } 
      if (fullQ)
      for (j=0;j<Q.r;j++) /* down through rows of Q */
      { x=0.0;p=AM[i];p1=QM[j];
        for (k=0;k<Ac-i;k++) x+=p[k]*p1[k];
        x*=g;
        for (k=0;k<Ac-i;k++) p1[k] += -x*p[k];
      } else
      { g=sqrt(g);
        p=QM[i];p1=AM[i]; /* address saving */
        for (j=0;j<Ac-i;j++) p[j]=p1[j]*g;
        for (j=Ac-i;j<Ac;j++) p[j]=0.0;  
      }
      AM[i][Ac-i-1]=-lsq;
      for (j=0;j<Ac-i-1;j++) AM[i][j]=0.0;
    }
  }
} /* QT */


void OrthoMult(matrix *Q,matrix *A,int off,int rows,int t,int pre,int o_pre)

/* The first `rows' of Q are vectors for a householder transformation. The ith
   vector starts with i+off zero elements (i starts at 0). The vectors are
   stored in the order that they should be applied. If o_pre==1 then they
   originally pre-multiplied, otherwise they originally post multiplied.
   if t==1 then the transform is transposed. If pre==1 then it is applied from
   the left.
   Each householder transform has the same form: (I-uu') (treating u as a column
   vector).
   The transformation is applied to A.
*/

{ double au,*u,*a,**AtM=NULL,**AM,**QM;
  int i,j,k,Ar,Qc,kk;
  matrix At;
  if (o_pre) t=1-t; /* default assumption is that creation was for post mult. */
  if (pre) /* use fact that QA=(A'Q')' and Q'A=(A'Q)' */
  { At=initmat(A->c,A->r);
    AM=A->M;AtM=At.M;
    for (i=0;i<A->r;i++) for (j=0;j<A->c;j++) AtM[j][i]=AM[i][j];
    t=1-t;
  } else At=*A;
  AM=At.M;QM=Q->M;Ar=At.r;Qc=Q->c;
  for (kk=0;kk<rows;kk++)
  { if (t) k=rows-1-kk; else k=kk;
    u=QM[k];
    for (i=0;i<Ar;i++)
    { a=AM[i];au=0.0;
      for (j=off+k;j<Qc;j++) au+=a[j]*u[j];
      for (j=off+k;j<Qc;j++) a[j] -= au*u[j];
    }
  }
  if (pre)
  { AM=A->M;
    for (i=0;i<At.r;i++) for (j=0;j<At.c;j++) AM[j][i]=AtM[i][j];
    freemat(At);
  }
} /* OrthoMult */


void Rsolv(matrix *R,matrix *p,matrix *y, int transpose)

/* Solves Rp=y for p when R is upper triangular - i.e. lower left 0.
   dimensions of p and y are not checked
   if transpose!=0 then solves: R'p=y
   
*/

{ int i,j,k;
  double x,*pV,*yV,*RMi,**RM,*dum,**pM,**yM;
  pV=p->V;yV=y->V;
  if (y->r==1) /* then p and y are vectors */
  { if (transpose) /* solve R'p=y for p */
    { RM=R->M;
      for (i=0;i<R->r;i++)
      { x=0.0;dum=pV;for (j=0;j<i;j++) { x+=RM[j][i] * *dum;dum++;}
        *dum=(yV[i]-x)/RM[i][i];
      } 
    } else /* solve Rp=y for p */
    for (i=R->r-1;i>=0;i--) 
    { RMi=R->M[i];
      x=0.0;for (j=i+1;j<R->r;j++) x+=RMi[j]*pV[j];
      pV[i]=(yV[i]-x)/RMi[i];
    }
  } else /* p and y are matrices */
  { pM=p->M;yM=y->M;
    if (transpose) /* solve R'p=y for p */
    { RM=R->M;
      for (k=0;k<p->c;k++)
      for (i=0;i<R->r;i++)
      { x=0.0;for (j=0;j<i;j++) x+=RM[j][i] * pM[j][k];
        pM[i][k]=(yM[i][k]-x)/RM[i][i];
      } 
    } else /* solve Rp=y for p */
    for (k=0;k<p->c;k++)
    for (i=R->r-1;i>=0;i--) 
    { RMi=R->M[i];
      x=0.0;for (j=i+1;j<R->r;j++) x+=RMi[j]*pM[j][k];
      pM[i][k]=(yM[i][k]-x)/RMi[i];
    }
  }
} /* Rsolv */


int QR(matrix *Q,matrix *R)

/* Does a QR factorisation of the matrix supplied in R. In Q the householder
   vectors are supplied to perform the transformation
   QR(in) -> R(out)
   R(out) is upper triangular (elements are 0 below leading diagonal).
   If Q->r is none zero then the vectors u are stored in successive rows of
   Q. The u vectors make up Q as a series of (stable) householder transformations.
   (I-uu'). The transformations are to be applied from the left in row order.
   The first i elements of the ith u are zero (i starting at zero).
   If A is the matrix input in R then QA=R, so that A=Q'R.
   Q can be used with OrthoMult(). 
   Under/overflow avoidance added 13/1/2000 along with more efficient calculation
   of length of u (modifications tested).
   
*/

{ int i,j,k,n,Rr;
  double *u,t,z,**RM,*p,m;
  RM=R->M;Rr=R->r;
  if (Rr<R->c) n=Rr; else n=R->c;
  u=(double *)CALLOC((size_t)Rr,sizeof(double));
  for (k=0;k<n;k++)
  { m=0.0;for (i=k;i<Rr;i++) { z=RM[i][k];z=fabs(z);if (z>m) m=z;}
    if (m) for (i=k;i<Rr;i++) RM[i][k]/=m; /* avoid over/underflow problems */
    t=0.0;for (i=k;i<Rr;i++) { z=RM[i][k];t+=z*z;} /* get euclidean length of column */
    if (RM[k][k]>0.0) t = -sqrt(t);else t= sqrt(t);  /* value of new RM[k][k] (stable) */
    for (i=k+1;i<Rr;i++) { u[i]=RM[i][k];RM[i][k]=0.0;}
    z=RM[k][k]; 
    u[k]=RM[k][k]-t;RM[k][k]=t*m;
    t=t*t;t+=u[k]*u[k]-z*z; /* efficient t calculation */
    /* t=0.0;for (p=u+k;p<u+Rr;p++) t+= *p * *p; - old inefficient calculation */
    t=sqrt(t/2);
    if (t==0.0) {FREE(u);return(0);} /* singular matrix */
    for (p=u+k;p<u+Rr;p++) *p /= t;
    for (j=k+1;j<R->c;j++)
    { t=0.0;for (i=k;i<Rr;i++) t+=u[i]*RM[i][j];
      for (i=k;i<Rr;i++) RM[i][j]-=u[i]*t;
    }
    if (Q->r) /* store vectors u for making Q */
    { p=Q->M[k];
      for (i=k;i<Rr;i++) p[i]=u[i];
    }
  }
  FREE(u);
  return(1);
} /* QR */




int real_elemcmp(const void *a,const void *b,int el)
/* declaring this inline static slows it down!! */
{ static int k=0;
  double *na,*nb,*nak;
  if (el>0) { k=el;return(0);}
  na=(*(double **)a);nb=(*(double **)b);
  nak = na + k;
  for (;na<nak;na++,nb++) { 
    if (*na < *nb) return(-1);
    if (*na > *nb) return(1);
  }
  return(0);
}

int melemcmp(const void *a,const void *b)

{ return(real_elemcmp(a,b,-1));
}


void msort(matrix a)

/* sorts a matrix, in situ, using standard routine qsort so 
   that its first col is in ascending order, its second col
   is in ascending order for any ties in the first col, and 
   so on.....
*/

{ double z=0.0;
  real_elemcmp(&z,&z,a.c); 
  qsort(a.M,(size_t)a.r,sizeof(a.M[0]),melemcmp);
}

void RArrayFromMatrix(double *a,int r,matrix *M)

/* copies matrix *M into R array a where r is the number of rows of A treated as
  a matrix by R */

{ int i,j;
  for (i=0;i<M->r;i++) for (j=0;j<M->c;j++) a[i+r*j]=M->M[i][j];
}


matrix Rmatrix(double *A,int r,int c)

/* produces a matrix from the array containing a (default) R matrix stored:
   A[0,0], A[1,0], A[2,0] .... etc */

{ int i,j;
  matrix M;
  M=initmat(r,c);
  for (i=0;i<r;i++) for (j=0;j<c;j++) M.M[i][j]=A[i+j*r];
  return(M);
}




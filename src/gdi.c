/* Copyright (C) 2007-2015 Simon N. Wood  simon.wood@r-project.org

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include "mgcv.h"

#ifdef OPENMP_ON
#include <omp.h>
#endif
#define ANSI
/*#define DEBUG*/



double trBtAB(double *A,double *B,int *n,int*m) 
/* form tr(B'AB) where A is n by n and B is n by m, m < n,
   basic point is that this is sum_ijk A_ik B_ij B_kj
 */
{ double tr=0.0,x,*p,*p1,*p2;
  int j,k;
  for (j=0;j<*m;j++)
    for (k=0;k<*n;k++) {
      p = A + *n * k;p2 = p + *n;
      p1 = B + *n * j;
      x = B[k + j * *n];
      for (;p<p2;p++,p1++) tr += *p * *p1 * x;   
    }
  return(tr);
}


void multSk(double *y,double *x,int *xcol,int k,double *rS,int *rSncol,int *q,double *work)
/* function to form y = Sk x, where a square root of Sk 
   is packed somewhere inside rS. x must be q by xcol. The 
   kth square root is q by rSncol[k]. The square roots are packed 
   one after another columnwise (R default).
   
   work and y must be the same dimension as x.
*/
{ int i,off,nc,bt,ct;
  double *rSk;
  off=0; /* the start of the kth square root */
  for (i=0;i<k;i++) off += *q * rSncol[i];
  rSk = rS + off; /* pointer to the kth square root */
  nc = rSncol[k];
  bt=1;ct=0;
  mgcv_mmult(work,rSk,x,&bt,&ct,&nc,xcol,q);
  bt=0;
  mgcv_mmult(y,rSk,work,&bt,&ct,q,xcol,&nc);
}

double diagABt(double *d,double *A,double *B,int *r,int *c)
/* obtain diag(AB') as efficiently as possible, and return tr(AB') A and B are
   r by c stored column-wise. 
   Thread safety: Modifies d. rest unmodified
*/
{ int j;
  double tr=0.0,*pa,*pb,*p1,*pd;
  if (*c>0) {
    for (pa=A,pb=B,p1=pa + *r,pd=d;pa<p1;pa++,pb++,pd++) *pd = *pa * *pb;
    for (j=1;j < *c;j++)
    for (p1=d + *r,pd=d;pd<p1;pa++,pb++,pd++) *pd += *pa * *pb;
    /* d now contains diag(AB') */
    for (pd=d,p1=d + *r;pd < p1;pd++) tr += *pd;
  }
  return(tr);
}


double trAB(double *A,double *B,int *n, int *m)
/* Get tr(AB) where A is n by m and B is m by n 
*/ 
{ double *p,*pa,*pb,tr=0.0;
  int i,j;
  for (pa=A,pb=B,j=0;j<*m;j++,pb++)
    for (p=pb,i=0;i<*n;i++,p += *m,pa++) tr+= *p * *pa;
  return(tr);
}


void get_bSb(double *bSb,double *bSb1, double *bSb2,double *sp,double *E,
             double *rS,int *rSncol,int *Enrow, int *q,int *M,int *M0,
             double *beta, double *b1, double *b2,int *deriv)
/*
  Routine to obtain beta'Sbeta and its derivatives w.r.t. the log smoothing 
   parameters, this is part of REML calculation... 
   
   sp is of length M, but b1 and b2 contain derivative w.r.t. to M0 + M 
   parameter [theta,log sp] where theta are parameters of the likelihood,
   but not the penalties.  The returned bSb1 and bSb2 should contain 
   derivatives w.r.t. theta as well as log sp (although these derivatives 
   will be zero). 

   b1 and b2 contain first and second derivatives of q-vector beta w.r.t. 
   \pho_k. They are packed as follows....

   * b1 will contain dbeta/d\rho_0, dbeta/d\rho_1 etc. So, for example, dbeta_i/d\rho_j
     (indices starting at zero) is located in b1[q*j+i].
   
   * b2 will contain d^2beta/d\rho_0d\rho_0, d^2beta/d\rho_1d\rho_0,... but rows will not be
     stored if they duplicate an existing row (e.g. d^2beta/d\rho_0d\rho_1 would not be 
     stored as it already exists and can be accessed by interchanging the sp indices).
     So to get d^2beta_k/d\rho_id\rho_j: 
     i)   if i<j interchange the indices
     ii)  off = (j*m-(j+1)*j/2+i)*q (m is number of sp's) 
     iii) v2[off+k] is the required derivative.       

*/
{ double *Sb,*Skb,*work,*work1,*p1,*p0,*p2,xx;
  int i,j,bt,ct,one=1,m,k,rSoff,mk,km,Mtot; 
  j = *q; for (i=0;i<*M;i++) if (rSncol[i]>j) j = rSncol[i];
  j += *M0; /* work space needed */
  
  work = (double *)CALLOC((size_t)j,sizeof(double)); 
  Sb = (double *)CALLOC((size_t)*q,sizeof(double));
  bt=0;ct=0;mgcv_mmult(work,E,beta,&bt,&ct,Enrow,&one,q);
  bt=1;ct=0;mgcv_mmult(Sb,E,work,&bt,&ct,q,&one,Enrow); /* S \hat \beta */

  for (*bSb=0.0,i=0;i<*q;i++) *bSb += beta[i] * Sb[i]; /* \hat \beta' S \hat \beta */

  if (*deriv <=0) {FREE(work);FREE(Sb);return;}

  work1 = (double *)CALLOC((size_t)j,sizeof(double));
  Skb = (double *)CALLOC((size_t)*M * *q,sizeof(double));
  for (p1=Skb,rSoff=0,i=0;i<*M;i++) { /* first part of first derivatives */
     /* form S_k \beta * sp[k]... */
     bt=1;ct=0;mgcv_mmult(work,rS + rSoff ,beta,&bt,&ct,rSncol+i,&one,q);
     for (j=0;j<rSncol[i];j++) work[j] *= sp[i]; 
     bt=0;ct=0;mgcv_mmult(p1,rS + rSoff ,work,&bt,&ct,q,&one,rSncol+i);
     rSoff += *q * rSncol[i];

     /* now the first part of the first derivative */
     for (xx=0.0,j=0;j<*q;j++,p1++) xx += beta[j] * *p1;
     bSb1[i + *M0] = xx; 
  }
  for (i=0;i<*M0;i++) bSb1[i] = 0.0;
  Mtot = *M + *M0;  

  if (*deriv>1)  for (m=0;m < Mtot;m++) { /* Hessian */
     bt=0;ct=0;mgcv_mmult(work1,E,b1 + m * *q,&bt,&ct,Enrow,&one,q);
     bt=1;ct=0;mgcv_mmult(work,E,work1,&bt,&ct,q,&one,Enrow);  /* S dbeta/drho_m */

    for (k=m;k < Mtot;k++) {
      km= k * Mtot + m ;
      mk= m * Mtot + k ;  /* second derivatives needed */
      /* d2beta'/drho_k drho_m S beta */
      for (xx=0.0,p0=Sb,p1=Sb + *q;p0<p1;p0++,b2++) xx += *b2 * *p0;
      bSb2[km] = 2*xx; 
       
      /* dbeta'/drho_k S dbeta/drho_m */
      for (xx=0.0,p0=b1+ k  * *q,p1=p0 + *q,p2=work;p0<p1;p0++,p2++) xx += *p2 * *p0;
      bSb2[km] += 2*xx;

      /* dbeta'/drho_m S_k beta sp[k] */
      if (k >= *M0) {
        for (xx=0.0,p0=Skb + (k- *M0) * *q,p1=p0 + *q,p2= b1+ m  * *q;p0<p1;p0++,p2++) 
             xx += *p2 * *p0;
        bSb2[km] += 2*xx;
      }
      /* dbeta'/drho_k S_m beta sp[m] */
      if (m >= *M0) {
        for (xx=0.0,p0=Skb + (m - *M0) * *q,p1=p0 + *q,p2= b1 + k  * *q;p0<p1;p0++,p2++) 
             xx += *p2 * *p0;
        bSb2[km] += 2*xx;
      }

      if (k==m) bSb2[km] += bSb1[k]; else bSb2[mk] = bSb2[km];
    }
  } /* done Hessian */

    
  /* Now finish off the first derivatives */
  bt=1;ct=0;mgcv_mmult(work,b1,Sb,&bt,&ct,&Mtot,&one,q);
  for (i=0;i<Mtot;i++) bSb1[i] += 2*work[i];
  
  FREE(Sb);FREE(work);FREE(Skb);FREE(work1);

} /* end get_bSb */

double frobenius_norm(double *X,int *r, int *c)
/* The Frobenius norm of r by c matrix X. Interestingly, this gives an 
   upper bound on the two norm (largest singular value). 
*/
{ double fnorm=0.0,*p1;
  int n;
  n = *r * *c;
  for (p1=X+n;X<p1;X++) fnorm += *X * *X;
  return(sqrt(fnorm));
}

void pivoter(double *x,int *r,int *c,int *pivot, int *col, int *reverse)
/* Routine for pivoting or unpivoting r by c matrix x 
   according to what's in `pivot'.

   The ith pivoted element comes from the original element pivot[i] 
   i.e. pivot[i] is the unpivoted element that pivoted element i 
   should end up in.

   If `reverse' is non-zero then x is unpivoted. Otherwise pivoting is 
   applied. 

   If `col' is non zero then columns are un/pivoted, otherwise rows.

   Typical applications are to pivot matrices in the same way that a 
   qr decomposition has been pivoted, or to reverse such a pivoting.
*/
{ double *dum,*px,*pd,*pd1,*p,*p1;
  int *pi,*pi1,i,j;
  if (*col) { /* pivot columns */ 
     dum = (double *) CALLOC((size_t)*c,sizeof(double));
     if (*reverse) /* unpivot x */
       for (i=0;i< *r;i++) {
	 for (px=x+i,pi=pivot,pi1=pi+*c;pi<pi1;pi++,px+=*r) dum[*pi]= *px; /*dum[pivot[j]] = x[j* *r + i] */ 
         for (px=x+i,pd=dum,pd1=dum+*c;pd<pd1;pd++,px += *r) *px = *pd;    /* x[j * *r + i] = dum[j]; */
       } else /* pivot x */
       for (i=0;i< *r;i++) {
	 for (px = x+i,pd=dum,pd1=dum + *c,j=0;pd < pd1;pd++,j++)  *pd = px[pivot[j] * *r];
	 for (px=x+i,pd=dum,pd1=dum+*c;pd<pd1;pd++,px += *r) *px = *pd;  /* x[j * *r + i] = dum[j]; */ 
     }
  } else { /* pivot rows */
    dum = (double *) CALLOC((size_t)*r,sizeof(double));
    if (*reverse) /* unpivot x */
    for (p=x,j=0;j<*c;j++,p += *r) { /* work column by column using dum as working storage */
      for (pi=pivot,pi1=pi+*r,p1=p;pi<pi1;pi++,p1++) dum[*pi] = *p1; /*dum[pivot[i]] = p[i]; ith row of pivoted -> pivot[i] row of unpivoted */
      for (pd=dum,pd1=dum+*r,p1=p;pd<pd1;pd++,p1++) *p1 = *pd;       /* store unpivoted column in x */

    } else  /* pivot x */
      for (p=x,j=0;j<*c;j++,p += *r) { /* work column by column using dum as working storage */
	for (pi=pivot,pi1=pi+*r,pd=dum;pi<pi1;pd++,pi++) *pd = p[*pi];  /* dum[i] = p[pivot[i]]; pivot[i] row of unpivoted -> ith row of pivoted */
        for (pd=dum,pd1=dum+*r,p1=p;pd<pd1;pd++,p1++) *p1 = *pd;        /* store pivoted column in x */
    }
  } 
  FREE(dum);
} /* end pivoter */


double qr_ldet_inv(double *X,int *r,double *Xi,int *get_inv) 
/* Obtains the log|X| and the inverse of X (r by r), by pivoted QR decomposition. 
   The inverse is returned (unpivoted) in Xi. 
   The function returns log|X| as its value.
   X is overwritten in the process
*/
{ double *tau,ldet,*p,*Qt;
  int *pivot,i,TRUE=1,j,FALSE=0;
  /* Allocated working storage ...*/
  pivot = (int *)CALLOC((size_t)*r,sizeof(int));
  tau = (double *)CALLOC((size_t)*r,sizeof(double));
  
  mgcv_qr(X,r,r,pivot,tau); /* get QR=X itself */

  /* evaluate log|X| = sum_i log(|R_ii|) ...*/
  for (ldet=0.0,p=X,i=0;i<*r;i++,p += *r+1) ldet += log(fabs(*p));
  
  if (*get_inv) {
  /* Now get the inverse of X. X^{-1} = R^{-1}Q' */
    Qt = (double *)CALLOC((size_t)*r * *r,sizeof(double));
    for (p=Qt,i=0;i<*r;i++,p += *r+1) *p = 1.0;
    mgcv_qrqy(Qt,X,tau,r,r,r,&TRUE,&TRUE); /* Extracting the orthogonal factor Q' */

    mgcv_backsolve(X,r,r,Qt,Xi,r,&FALSE); /* Now Xi contains the row pivoted inverse of X */

    /* Finally unpivot Xi. 
       pivot[i] is the unpivoted row that pivoted row i should end up in  
    */

    for (p=Xi,j=0;j<*r;j++,p += *r) { /* work column by column using tau as working storage */

      for (i=0;i<*r;i++) tau[pivot[i]] = p[i]; /* ith row of pivoted -> pivot[i] row of unpivoted */
      for (i=0;i<*r;i++) p[i] = tau[i];        /* store unpivoted column in Xi */

    }
    FREE(Qt);
  } /* end if (*get_inv) */
  FREE(pivot);FREE(tau);
  return(ldet);
} /* end qr_ldet_inv */



void get_detS2(double *sp,double *sqrtS, int *rSncol, int *q,int *M, int * deriv, 
               double *det, double *det1, double *det2, double *d_tol,
               double *r_tol,int *fixed_penalty)
/* Routine to evaluate log|S| and its derivatives wrt log(sp), in a stable manner, using 
   an orthogonal transformation strategy based on QR decomposition.

   Inputs are:
   `sp' the array of smoothing parameters.
   `sqrtS' the `M' square root penalty matrices. The ith is `q' by `rSncol[i]'. They are 
        packed one after the other. 
   `deriv' is the order of derivatives required. 0,1 or 2.
   `d_tol' is the tolerance to use for grouping dominant terms. 
   `r_tol' (<< d_tol) is the tolerance used for rank determination.
   `fixed_penalty' non-zero indicates that there is a fixed component of 
          total penalty matrix S, the square root of which is in the final 
          q * rSncol[M+1] elements of sqrtS.                 

   Outputs are:
   `det' the log determinant.
   `det1' M-array of derivatives of log det wrt log sp. 
   `det2' M by M Hessian of log det wrt log sp.   
   
*/
{ double *R,*work,*tau,*rS1,*rS2, *S,*Si,*Sb,*B,*Sg,*p,*p1,*p2,*p3,*p4,*frob,max_frob,x,*spf,Rcond;
  int *pivot,iter,i,j,k,bt,ct,rSoff,K,Q,Qr,*gamma,*gamma1,*alpha,r,max_col,Mf,tot_col=0,left,tp;

  if (*fixed_penalty) { 
    Mf = *M + 1;  /* total number of components, including fixed one */
    spf = (double *)CALLOC((size_t)Mf,sizeof(double));
    for (i=0;i<*M;i++) spf[i]=sp[i];
    spf[*M]=1.0; /* includes sp for fixed term */
  } 
  else {spf=sp;Mf = *M;} /* total number of components, including fixed one */

  /* Create working copies of sqrtS, which can be modified:
     rS1 is repeatedly orthogonally transformed, while rS2 is row pivoted. 
  */
  if (*deriv) { /* only need to modify if derivatives needed */
    for (j=i=0;i<Mf;i++) j += rSncol[i];
    tot_col=j;
    j *= *q;
    rS1 = (double *)CALLOC((size_t) j,sizeof(double));
    rS2 = (double *)CALLOC((size_t) j,sizeof(double));
    for (p=rS1,p3=rS2,p1=rS1+j,p2=sqrtS;p<p1;p++,p2++,p3++) *p3 = *p = *p2;
  } else {rS1=rS2=NULL;}
  /* Explicitly form the Si (stored in a single block), so S_i is stored
     in Si + i * q * q (starting i from 0). As iteration progresses,
     blocks are shrunk -- always q by Q */

  max_col = *q; /* need enough storage just in case square roots are over-sized */ 
  for (i=0;i<Mf;i++) if (rSncol[i]>max_col) max_col=rSncol[i];

  p = Si = (double *)CALLOC((size_t)*q * max_col * Mf,sizeof(double));
  
  for (rSoff=i=0;i<Mf;p+= *q * *q,rSoff+=rSncol[i],i++) {
    bt=0;ct=1;mgcv_mmult(p,sqrtS+rSoff * *q,sqrtS+rSoff * *q,&bt,&ct,q,q,rSncol+i);   
  }

 
  /* Initialize the sub-dominant set gamma and the counters */
  K = 0;Q = *q;
  frob =  (double *)CALLOC((size_t)Mf,sizeof(double)); 
  gamma = (int *)CALLOC((size_t)Mf,sizeof(int));  /* terms remaining to deal with */
  gamma1 = (int *)CALLOC((size_t)Mf,sizeof(int)); /* new gamma */
  alpha = (int *)CALLOC((size_t)Mf,sizeof(int));  /* dominant terms */
  for (i=0;i<Mf;i++) gamma[i] = 1; /* no terms dealt with initially */
  
  /* Other storage... */

  S = (double *) CALLOC((size_t) Q * Q,sizeof(double)); /* Transformed S (total) */
  Sb = (double *) CALLOC((size_t) Q * Q,sizeof(double)); /* summation storage */
  pivot = (int *)CALLOC((size_t) Q,sizeof(int)); /* pivot storage */
  tau = (double *) CALLOC((size_t) Q,sizeof(double)); /* working storage */  
  work = (double *)CALLOC((size_t)(4 * Q),sizeof(double));

  Sg = (double *) CALLOC((size_t) Q * Q,sizeof(double)); /* summation storage */
  B = (double *) CALLOC((size_t) Q * max_col,sizeof(double)); /* Intermediate storage */
  R = (double *) CALLOC((size_t) Q * Q,sizeof(double)); /* storage for unpivoted QR factor */
  /* Start the main orthogonal transform loop */
  iter =0;
  while(1) {
    iter ++;

  /* Find the Frobenius norms of the Si in set gamma */
    max_frob=0.0;
    for (p=Si,i=0;i<Mf;i++,p += *q * Q) 
      if (gamma[i]) { /* don't bother if already dealt with */ 
        frob[i] = frobenius_norm(p,q,&Q);
        if (frob[i] *spf[i] >max_frob) max_frob=frob[i]  * spf[i];
    }
  /* Find sets alpha and gamma' */
    for (i=0;i<Mf;i++) {
      if (gamma[i]) { /* term is still to be dealt with */
        if (frob[i]  * spf[i] > max_frob * *d_tol) { 
          alpha[i] = 1;gamma1[i] = 0; /* deal with it now */
        } else {
          alpha[i] = 0;gamma1[i] = 1; /* put it off */ 
        }
      } else { /* wasn't in gamma, so not in alpha or gamma1 */
        alpha[i] = gamma1[i] = 0;
      }
    }

  /* Form the scaled sum of the Si in alpha and get its rank by pivoted QR
     and condition estimation...
  */
    for (p=Sb,p1=p + *q * Q;p<p1;p++) *p=0.0; /* clear Sb */
    for (p=Si,i=0;i<Mf;i++,p += *q * Q) if (alpha[i]) { 
      x = frob[i];
      for (p1=p,p2=Sb,p3=p + *q * Q;p1<p3;p1++,p2++) *p2 += *p1 / x;
    } 
    for (i=0;i<*q;i++) pivot[i]=0; 
    mgcv_qr(Sb, &Q, q ,pivot,tau); /* obtain pivoted QR decomposition of Sb */
    /* Now obtain the rank, r, of Sb (see Golub and van Loan, 1996, p.129 & p.260)... */ 
    r = Q;
    R_cond(Sb,&Q,&r,work,&Rcond);
    while (*r_tol * Rcond > 1) { r--;R_cond(Sb,&Q,&r,work,&Rcond);}
    Qr = Q-r;  

    /* ...  r is the rank of Sb, or any other positively weighted sum over alpha */

    /*  printf("\n iter = %d,  rank = %d,   Q = %d",iter,r,Q);
    printf("\n gamma = ");for (i=0;i<Mf;i++) printf(" %d",gamma[i]);
    printf("\n alpha = ");for (i=0;i<Mf;i++) printf(" %d",alpha[i]);
    printf("\n gamma1 = ");for (i=0;i<Mf;i++) printf(" %d",gamma1[i]);*/
   

  /* If Q==r then terminate (form S first if it's the first iteration) */
    
    if (Q==r) { 
      if (iter==1 ) { /* form S */
        for (p=Si,i=0;i<Mf;i++,p += Q*Q) { 
          x = spf[i];
          for (p1=p,p2=S,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
        }
        break; 
      } else break; /* just use current S */ 
    } /* end if (Q==r) */

  /* Form the dominant term and QR-decompose it */
    for (p=Sb,p1=p + *q * Q;p<p1;p++) *p = 0.0; /* clear Sb */
    for (p=Si,i=0;i<Mf;i++,p += *q * Q) if (alpha[i]) { /* summing S[[i]]*sp[i] over i in alpha */
      x = spf[i];
      for (p1=p,p2=Sb,p3=p+ *q * Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
    }
    for (i=0;i<*q;i++) pivot[i]=0; 
    mgcv_qr(Sb, &Q, q ,pivot,tau); /* obtain pivoted QR decomposition of Sb */
    
  /* unpivot R, which means that no further pivoting is needed */
    for (p=R,p1=R + *q * r;p<p1;p++) *p=0.0; /* clear R */
    for (i=0;i<r;i++) for (j=i;j<*q;j++) R[i + pivot[j] * r] = Sb[i + j * Q]; 

   /* DEBUG ONLY... */
    /*  printf("\npivot = ");for (j=0;j<*q;j++) printf("%d ",pivot[j]);
    printf("Current R...\n");
    for (i=0;i<r;i++) { for (j=0;j<*q;j++) printf("%7.2g  ",Sb[i + Q *j]); printf("\n");} */

  /* Form the sum over the elements in gamma1, Sg */

    for (p=Sg,p1=p + *q * Q;p<p1;p++) *p=0.0; /* clear Sg */
    for (p=Si,i=0;i<Mf;i++,p += *q * Q) if (gamma1[i]) { /* summing S[[i]]*sp[i] over i in gamma1 */
      x = spf[i];
      for (p1=p,p2=Sg,p3=p+ *q * Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
    } 

  /* Form S' the orthogonal transform of S */
  
  /* Form Q'Sg... */ 
    left=1;tp=1; 
    mgcv_qrqy(Sg,Sb,tau,&Q,q,&Q,&left,&tp);

  /* copy transformed Sg into remainder of transformed S */
    for (i=0;i<Q;i++) for (j=0;j<*q;j++) S[i+K + j * *q] = Sg[i + j * Q];
  /* and add R in the appropriate place ... */ 
    for (i=0;i<r;i++) for (j=0;j<*q;j++) S[i+K + j * *q] += R[i + j * r];       


  /* transform remaining S_i in gamma1 */ 
    for (p1=p=Si,i=0;i<Mf;i++,p += *q * Q,p1 += *q *Qr) if (gamma1[i]) {
      left=1;tp=1;
      mgcv_qrqy(p,Sb,tau,&Q,q,&Q,&left,&tp); /* Q'Si */
      p2=p+r;p3=p1;
      for (j=0;j<*q;j++,p2+=r) for (k=0;k<Qr;k++,p3++,p2++) *p3 = *p2; /* copy to correct place */ 
    } 
 
  /* Transform the square roots of Si */
    if (*deriv) { /* transformed rS1 only needed for derivatives */
      /* copy last Q rows of rS1 into rS2 */
      for (i=0;i<Q;i++) for (j=0;j<tot_col;j++) rS2[i+Q*j] = rS1[K + i + *q * j];
      /* pre-multiply rS2 by Q */ 
      left=1;tp=1;
      mgcv_qrqy(rS2,Sb,tau,&Q,&tot_col,&Q,&left,&tp); /* Q'rS2 */
      
      /* copy rS2 into last Q rows of rS1 */
      for (i=0;i<Q;i++) for (j=0;j<tot_col;j++) rS1[K + i + *q * j] = rS2[i+Q*j];
    
      /* zero the last Qr rows of the rS1 in alpha */
      for (p=rS1,k=0;k<Mf;p +=rSncol[k] * *q, k++) if (alpha[k]) {
        for (i=K+r;i<*q;i++) for (j=0;j<rSncol[k];j++) p[i + j * *q] = 0.0;
      }
    }


    /* DEBUG ONLY... */
    /* printf("Current S...\n");
       for (i=0;i<*q;i++) { for (j=0;j<*q;j++) printf("%7.2g  ",S[i + *q *j]); printf("\n");}*/
 
  /* Update K, Q and gamma */   
    K = K + r; Q = Qr;
    for (i=0;i<Mf;i++) gamma[i] = gamma1[i];
  } /* end of Orthogonal Transform Loop */

  /* transpose S */
  for (i=0;i<*q;i++) for (j=0;j<*q;j++) R[i + *q * j] = S[j + *q * i]; 

  /* DEBUG ONLY... */
  /* printf("Final S...\n");
     for (i=0;i<*q;i++) { for (j=0;j<*q;j++) printf("%7.2g  ",S[i + *q *j]); printf("\n");}*/

  /* Now get the determinant and inverse of the transformed S (stored in B) */
  *det = qr_ldet_inv(R,q,B,deriv); /* R=S' here */
  /* finally, the derivatives, based on transformed S inverse and transformed square roots */  
  
  if (*deriv) { /* get the first derivatives */
    /* first accumulate S^{-T} sqrtS into Si */
    bt=0;ct=0;mgcv_mmult(Si,B,sqrtS,&bt,&ct,q,&tot_col,q);
    /* Now get the required derivatives */
    for (p=Si,p1=rS1,i=0;i<*M;p += *q *rSncol[i],p1+= *q *rSncol[i],i++) {
      for (x=0.0,p2=p1,p3=p,p4=p1+ *q*rSncol[i];p2<p4;p2++,p3++) x += *p2 * *p3; 
      det1[i] = x*sp[i]; /* tr(S^{-1}S_i) */
    }
  }
  
  if (*deriv==2) { /* get second derivatives, as well */
    for (p=Si,p1=rS2,p2 = p1 + *q * tot_col;p1<p2;p1++,p++) *p1 = *p; /* copy S^{-1} sqrtS into rS2 */   

    /* loop through creating S^{-1} S_i and storing in Si...*/
    for (p1=Si,p=rS2,p2=rS1,i=0;i<*M;p2+= *q * rSncol[i], p += *q *rSncol[i],i++,p1 += *q * *q) {
      bt=0;ct=1;mgcv_mmult(p1,p,p2,&bt,&ct,q,q,rSncol+i);
    }

    /* DEBUG ONLY...
       for (i=0;i<*M;i++) { for (x=0.0,j=0;j<*q;j++) x += Si[i* *q * *q + j + j* *q];det1[i]=x*sp[i];}*/

    for (i=0;i<*M;i++) for (j=i;j<*M;j++) 
      det2[i + *M * j] = det2[j + *M * i] = -sp[i]*sp[j]*trAB(Si + *q * *q *i,Si + *q * *q *j,q,q);
    for (i=0;i<*M;i++) det2[i + *M * i] += det1[i];
  }
  FREE(R);
  FREE(work);
  FREE(frob);
  FREE(gamma);
  FREE(gamma1);
  FREE(alpha);
  FREE(S);
  FREE(Sb);
  FREE(Sg);
  if (*deriv) { FREE(rS1);FREE(rS2);}
  if (*fixed_penalty) {FREE(spf);}
  FREE(Si);
  FREE(B);
  FREE(pivot);FREE(tau);
} /* end of get_detS2 */




void get_stableS(double *S,double *Qf,double *sp,double *sqrtS, int *rSncol, int *q,int *M, int * deriv, 
               double *det, double *det1, double *det2, double *d_tol,
               double *r_tol,int *fixed_penalty)
/* Routine to similarity transform S = \sum_i \lambda_i S_i, to produce an S which facilitates
   stable computation. 

   THEORETICAL NOTE:  If the square-root of S is found by choleski on the diagonally pre-conditioned 
                     S, then a well behaved root is obtained, with no `large-zero' leakage beyond 
                     the range space of each component of S. This should be compared with `mroot' 
                     (just pivoted Choleski, where some leakage does occur -- into the penalty null 
                      space, but within the non-zero block of the penalty). 

   Also evaluates log|S| and its derivatives wrt log(sp), in a stable manner, using 
   a similarity transform strategy.

   Inputs are:
   `sp' the array of smoothing parameters.
   `sqrtS' the `M' square root penalty matrices. The ith is `q' by `rSncol[i]'. They are 
        packed one after the other. 
   `deriv' is the order of derivatives required. 0,1 or 2.
   `d_tol' is the tolerance to use for grouping dominant terms. 
   `r_tol' (<< d_tol) is the tolerance used for rank determination.
   `fixed_penalty' non-zero indicates that there is a fixed component of 
          total penalty matrix S, the square root of which is in the final 
          q * rSncol[M+1] elements of sqrtS.                 

   Outputs are:
   `det' the log determinant.
   `det1' M-array of derivatives of log det wrt log sp. 
   `det2' M by M Hessian of log det wrt log sp.   
   `S' - the similarity transformed total penalty matrix
   `Qf' - the orthogonal factor of the similarity transform. If S0 is the 
          original total penalty then S = Qf' S0 Qf 
    sqrtS - the square roots of the components of S, transformed as S itself.        
   
*/
{ double *rS, *Un, *U, *Si,*Sb,*B,*C,*Sg,*p,*p1,*p2,*p3,*frob,*ev,max_frob,x,*spf;
  int iter,i,j,k,bt,ct,rSoff,K,Q,Qr,*gamma,*gamma1,*alpha,TRUE=1,FALSE=0,r,max_col,Mf,n_gamma1;

  if (*fixed_penalty) { 
    Mf = *M + 1;  /* total number of components, including fixed one */
    spf = (double *)CALLOC((size_t)Mf,sizeof(double));
    for (i=0;i<*M;i++) spf[i]=sp[i];
    spf[*M]=1.0; /* includes sp for fixed term */
  } 
  else {spf=sp;Mf = *M;} /* total number of components, including fixed one */

  /* Create a working copy of sqrtS, which can be modified  */
 
  rS = sqrtS; /* this routine modifies sqrtS */
  /* Explicitly form the Si (stored in a single block), so S_i is stored
     in Si + i * q * q (starting i from 0). As iteration progresses,
     blocks are shrunk -- always Q by Q */
  p = Si = (double *)CALLOC((size_t)*q * *q * Mf,sizeof(double));
  max_col = *q; /* need enough storage just in case square roots are over-sized */
  for (rSoff=i=0;i<Mf;p+= *q * *q,rSoff+=rSncol[i],i++) {
    bt=0;ct=1;mgcv_mmult(p,sqrtS+rSoff * *q,sqrtS+rSoff * *q,&bt,&ct,q,q,rSncol+i);
    if (rSncol[i]>max_col) max_col=rSncol[i];
  }

 
  /* Initialize the sub-dominant set gamma and the counters */
  K = 0; /* counter for coefs already deal with */
  Q = *q; /* How many coefs left to deal with */
  frob =  (double *)CALLOC((size_t)Mf,sizeof(double)); 
  gamma = (int *)CALLOC((size_t)Mf,sizeof(int));  /* terms remaining to deal with */
  gamma1 = (int *)CALLOC((size_t)Mf,sizeof(int)); /* new gamma */
  alpha = (int *)CALLOC((size_t)Mf,sizeof(int));  /* dominant terms */
  for (i=0;i<Mf;i++) gamma[i] = 1; /* no terms dealt with initially */
  
  /* Other storage... */

  U=Sb = (double *) CALLOC((size_t) Q * Q,sizeof(double)); /* summation storage */

  Sg = (double *) CALLOC((size_t) Q * Q,sizeof(double)); /* summation storage */
  ev = (double *) CALLOC((size_t) Q,sizeof(double));     /* eigenvalue storage */
  B = (double *) CALLOC((size_t) Q * max_col,sizeof(double)); /* Intermediate storage */
  C = (double *) CALLOC((size_t) Q * max_col,sizeof(double)); /* Intermediate storage */

  /* Start the main similarity transform loop */
  iter =0;
  while(1) {
    iter ++;

  /* Find the Frobenius norms of the Si in set gamma */
    max_frob=0.0;
    for (p=Si,i=0;i<Mf;i++,p += Q * Q) 
      if (gamma[i]) { /* don't bother if already dealt with */ 
        frob[i] = frobenius_norm(p,&Q,&Q);
        if (frob[i] *spf[i] >max_frob) max_frob=frob[i]  * spf[i];
    }
  /* Find sets alpha and gamma' */
    n_gamma1=0;
    for (i=0;i<Mf;i++) {
      if (gamma[i]) { /* term is still to be dealt with */
        if (frob[i]  * spf[i] > max_frob * *d_tol) { 
          alpha[i] = 1;gamma1[i] = 0; /* deal with it now */
        } else {
          alpha[i] = 0;gamma1[i] = 1; n_gamma1++; /* put it off */ 
        }
      } else { /* wasn't in gamma, so not in alpha or gamma1 */
        alpha[i] = gamma1[i] = 0;
      }
    }

  /* Form the scaled sum of the Si in alpha and eigen-decompose it to get its rank */
    if (n_gamma1) { /* stuff left in gamma1, so have to work out rank of contents of alpha */
      for (p=Sb,p1=p+Q*Q;p<p1;p++) *p=0.0; /* clear Sb */
      for (p=Si,i=0;i<Mf;i++,p += Q*Q) if (alpha[i]) { 
        x = frob[i];
        for (p1=p,p2=Sb,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 / x;
      } 
      mgcv_symeig(Sb,ev,&Q,&FALSE,&FALSE,&FALSE); /* get eigenvalues (ascending) of scaled sum over alpha */
      
      r=1;
      while(r<Q&&(ev[Q-r-1]>ev[Q-1] * *r_tol)) r++;
    } else { /* nothing left in gamma1, so... */
      r=Q;
    }
    /* ...  r is the rank of Sb, or any other positively weighted sum over alpha */


    /* If Q==r then terminate (form S first if it's the first iteration) */
    
    if (Q==r) { 
      if (iter==1 ) { /* form S and Qf*/
        for (p=Si,i=0;i<Mf;i++,p += Q*Q) { 
          x = spf[i];
          for (p1=p,p2=S,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
        }
        /* Qf is identity */
        for (p=Qf,p1=p+Q*Q;p<p1;p++) *p=0.0;
        for (p=Qf,i=0;i<Q;i++,p+=Q+1) *p=1.0;
        break; 
      } else break; /* just use current S */ 
    } /* end if (Q==r) */

  /* Form the dominant term and eigen-decompose it */
    for (p=Sb,p1=p+Q*Q;p<p1;p++) *p = 0.0; /* clear Sb */
    for (p=Si,i=0;i<Mf;i++,p += Q*Q) if (alpha[i]) { /* summing S[[i]]*sp[i] over i in alpha */
      x = spf[i];
      for (p1=p,p2=Sb,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
    } 

    mgcv_symeig(Sb,ev,&Q,&FALSE,&TRUE,&TRUE); /* get eigen decomposition of dominant term (ev descending) */
    
  /* .... U points to Sb, which now contains eigenvectors */
    if (iter==1) for (p=U,p1=Qf,p2 = p+Q*Q;p<p2;p++,p1++) *p1 = *p;
    else { 
      bt=0;ct=0;mgcv_mmult(B,Qf+K * *q,U,&bt,&ct,q,&Q,&Q);
      for (p=Qf + K * *q,p1=Qf + *q * *q,p2=B;p<p1;p++,p2++) *p = *p2;
    }
  /* Form the sum over the elements in gamma1, Sg */

    for (p=Sg,p1=p+Q*Q;p<p1;p++) *p=0.0; /* clear Sg */
    for (p=Si,i=0;i<Mf;i++,p += Q*Q) if (gamma1[i]) { /* summing S[[i]]*sp[i] over i in gamma1 */
      x = spf[i];
      for (p1=p,p2=Sg,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
    } 

  /* Form S' the similarity transformed S */
    if (K>0) { /* deal with upper right component B */
      /* first copy out K by Q matrix B  */ 
      for (j=0;j<Q;j++) for (i=0;i<K;i++) C[i + K * j] = S[i + *q * (j+K)];
      /* Now form BU (store in B)*/
    
      bt=0;ct=0;mgcv_mmult(B,C,U,&bt,&ct,&K,&Q,&Q);
      /* Replace B into S */
      for (j=0;j<Q;j++) for (i=0;i<K;i++) S[i + *q * (j+K)]= S[j + K + *q * i] = B[i + K * j];
    }

    /* Now deal with the lower right component, C */
    /* U'SgU  */
 
    bt=0;ct=0;mgcv_mmult(B,Sg,U,&bt,&ct,&Q,&Q,&Q); /* SgU is in B */
    bt=1;ct=0;mgcv_mmult(C,U,B,&bt,&ct,&Q,&Q,&Q);  /* U'SgU is in C */ 
    
    for (i=0;i<r;i++) C[i+i * Q] += ev[i];  /* Adding in the non-zero eigenvalues */
   
    /* Now copy U'SgU + D back into right part of S' */
    for (j=0;j<Q;j++) for (i=0;i<Q;i++) S[i + K + *q * (j+K)] = C[i + Q * j];
    
    /* Transform the square roots of Si in alpha and gamma1 (Can leave fixed term alone - not needed)*/
   
     for (p=rS,k=0;k<*M;p += rSncol[k] * *q,k++) if (alpha[k]) {  /* p points to the square root of S_i */    
        /* extract the part of rS_i to be modified */
        for (i=0;i<Q;i++) for (j=0;j<rSncol[k];j++) C[i + Q * j] = p[i + K + *q * j]; 
        bt=1;ct=0;mgcv_mmult(B,U,C,&bt,&ct,&r,rSncol+k,&Q); 
        for (i=0;i<r;i++) for (j=0;j<rSncol[k];j++) p[i + K + *q * j] = B[i + r * j];
        for (i=K+r;i<K+Q;i++) for (j=0;j<rSncol[k];j++) p[i + *q * j] = 0.0;
      } else if (gamma1[k]) { 
        for (i=0;i<Q;i++) for (j=0;j<rSncol[k];j++) C[i + Q * j] = p[i + K + *q * j];
       
        bt=1;ct=0;mgcv_mmult(B,U,C,&bt,&ct,&Q,rSncol+k,&Q);
        for (i=0;i<Q;i++) for (j=0;j<rSncol[k];j++) p[i + K + *q * j] = B[i + Q * j];
      }
   

  /* Transform the Si in gamma' */
    Qr = Q - r;Un = U + r * Q;
   
    for (p1=p=Si,i=0;i<Mf;i++,p += Q*Q,p1 +=Qr*Qr) if (gamma1[i]) { /* p points to old Si, and p1 to new */
      bt=1;ct=0;mgcv_mmult(B,Un,p,&bt,&ct,&Qr,&Q,&Q);
      bt=0;ct=0;mgcv_mmult(p1,B,Un,&bt,&ct,&Qr,&Qr,&Q); 
    }
  /* Update K, Q and gamma */   
    K = K + r; Q = Qr;
    for (i=0;i<Mf;i++) gamma[i] = gamma1[i];
  } /* end of Similarity Transfrom Loop */

  /* Now get the determinant and inverse of the transformed S (stored in C, which gets overwritten) 
     inverse of S returned in B */
  for (p=S,p1=S + *q * *q,p2=C;p<p1;p++,p2++) *p2 = *p; /* copy S to C */
  *det = qr_ldet_inv(C,q,B,deriv);
  /* finally, the dervivatives, based on transformed S inverse and transformed square roots */  
  
  if (*deriv) { /* get the first derivatives */
   for (p=rS,i=0;i<*M;p += *q *rSncol[i],i++) det1[i] = trBtAB(B,p,q,rSncol+i)*sp[i]; /* tr(S^{-1}S_i) */
  }
 
  if (*deriv==2) { /* get second derivatives, as well */
    for (p1=Si,p=rS,i=0;i<*M;p += *q *rSncol[i],i++,p1 += *q * *q) { /* loop through creating S^{-1} S_i and storing in Si*/
      bt=0;ct=0;mgcv_mmult(C,B,p,&bt,&ct,q,rSncol+i,q);
      bt=0;ct=1;mgcv_mmult(p1,C,p,&bt,&ct,q,q,rSncol+i);
    }
    for (i=0;i<*M;i++) for (j=i;j<*M;j++) 
      det2[i + *M * j] = det2[j + *M * i] = -sp[i]*sp[j]*trAB(Si + *q * *q *i,Si + *q * *q *j,q,q);
    for (i=0;i<*M;i++) det2[i + *M * i] += det1[i];
  }

 
  FREE(frob);
  FREE(gamma);
  FREE(gamma1);
  FREE(alpha);
  FREE(Sb);
  FREE(Sg);
  if (*fixed_penalty) {FREE(spf);}
  FREE(Si);
  FREE(ev);
  FREE(B);
  FREE(C);
}/* end get_stableS */


void get_ddetXWXpS(double *det1,double *det2,double *P,double *K,double *sp,
      double *rS,int *rSncol,double *Tk,double *Tkm,int *n,int *q,int *r,int *M,int *M0,
		   int *deriv,int nthreads)

/* obtains derivatives of |X'WX + S| wrt the log smoothing parameters, as required for REML. 
   The determinant itself has to be obtained during intial decompositions: see gdi().

   * P is q by r
   * K is n by r
  
   * sp is of length M, and there are M penalties, but derivatives are required w.r.t. M0 + M
     parameters [theta,sp] and Tk, Tkm contain derivatives w.r.t. all of these. The parameters
     theta are parameters of the likelihood, but are not involved in penalties.

   * this routine assumes that sp contains smoothing parameters, rather than log smoothing parameters.
 
   * Note that P and K are as in Wood (2008) JRSSB 70, 495-518.

   * uses nthreads via openMP - assumes nthreads already set and nthreads already reset to 1
     if openMP not present
*/

{ double *diagKKt,xx,*KtTK,*PtrSm,*PtSP,*trPtSP,*work,*pdKK,*p1,*pTkm;
    int m,k,bt,ct,j,one=1,km,mk,*rSoff,deriv2,max_col,Mtot;
  int tid;
  if (nthreads<1) nthreads = 1;
  
  Mtot = *M0 + *M; /* total length of sp */

  if (*deriv==2) deriv2=1; else deriv2=0;
  /* obtain diag(KK') */ 
  if (*deriv) {
    diagKKt = (double *)CALLOC((size_t)*n,sizeof(double));
    xx = diagABt(diagKKt,K,K,n,r); 
  } else { /* nothing to do */
      return;
  }
  /* set up work space */
  work =  (double *)CALLOC((size_t)*n * nthreads,sizeof(double));
  tid=0; /* thread identifier defaults to zero if openMP not available */
  /* now loop through the smoothing parameters to create K'TkK */
  if (deriv2) {
    KtTK = (double *)CALLOC((size_t)(*r * *r * Mtot),sizeof(double));
    #ifdef OPENMP_ON
    #pragma omp parallel private(k,j,tid) num_threads(nthreads)
    #endif
    { /* open parallel section */
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (k=0;k < Mtot;k++) {
	#ifdef OPENMP_ON
        tid = omp_get_thread_num(); /* thread running this bit */
	#endif    
        j = k * *r * *r;
        getXtWX(KtTK+ j,K,Tk + k * *n,n,r,work + *n * tid);
      }
    } /* end of parallel section */
  } else { KtTK=(double *)NULL;} /* keep compiler happy */

  /* start first derivative */ 
  bt=1;ct=0;mgcv_mmult(det1,Tk,diagKKt,&bt,&ct,&Mtot,&one,n); /* tr(TkKK') */ 

  /* Finish first derivative and create create P'SmP if second derivs needed */
  
  max_col = *q;
  for (j=0;j<*M;j++) if (max_col<rSncol[j]) max_col=rSncol[j]; /* under ML can have q < max(rSncol) */

  PtrSm = (double *)CALLOC((size_t)(*r * max_col * nthreads),sizeof(double)); /* storage for P' rSm */
  trPtSP = (double *)CALLOC((size_t) *M,sizeof(double));

  if (deriv2) {
    PtSP = (double *)CALLOC((size_t)(*M * *r * *r ),sizeof(double));
  } else { PtSP = (double *) NULL;}
  
  
  rSoff =  (int *)CALLOC((size_t)*M,sizeof(int));
  if (*M>0) {
    rSoff[0] = 0;for (m=0;m < *M-1;m++) rSoff[m+1] = rSoff[m] + rSncol[m];
  }
  tid = 0;
  #ifdef OPENMP_ON
  #pragma omp parallel private(m,bt,ct,tid) num_threads(nthreads)
  #endif
  { /* parallel section start */
    #ifdef OPENMP_ON
    #pragma omp for
    #endif
    for (m=0;m < *M;m++) { /* loop through penalty matrices */
      #ifdef OPENMP_ON
      tid = omp_get_thread_num(); /* thread running this bit */
      #endif    
      bt=1;ct=0;mgcv_mmult(PtrSm + tid * *r * max_col,P,rS+rSoff[m] * *q,&bt,&ct,r,rSncol+m,q);
      /*rSoff += rSncol[m];*/
      trPtSP[m] = sp[m] * diagABt(work + *n * tid,PtrSm + tid * *r * max_col,
                                 PtrSm + tid * *r * max_col,r,rSncol+m); /* sp[m]*tr(P'S_mP) */ 
      det1[m + *M0] += trPtSP[m]; /* completed first derivative */
      if (deriv2) { /* get P'S_mP */
        bt=0;ct=1;mgcv_mmult(PtSP+ m * *r * *r,PtrSm + tid * *r * max_col,
                            PtrSm+ tid * *r * max_col ,&bt,&ct,r,r,rSncol+m);
      }
    }
  } /* end of parallel section */
  FREE(rSoff);
  /* Now accumulate the second derivatives */

  //  #ifdef OPENMP_ON
  //#pragma omp parallel private(m,k,km,mk,xx,tid,pdKK,p1,pTkm) num_threads(nthreads)
  //#endif
  if (deriv2)
  { /* start of parallel section */ 
    //if (deriv2) 
    #ifdef OPENMP_ON
    #pragma omp parallel for private(m,k,km,mk,xx,tid,pdKK,p1,pTkm) num_threads(nthreads)
    #endif
    for (m=0;m < Mtot;m++) {
      #ifdef OPENMP_ON
      tid = omp_get_thread_num(); /* thread running this bit */
      #endif
      if (m==0) pTkm = Tkm; else pTkm = Tkm + (m * Mtot - (m*(m-1))/2) * *n;        
      for (k=m;k < Mtot;k++) {
        km=k * Mtot + m;mk=m * Mtot + k;
        /* tr(Tkm KK') */
        /*for (xx=0.0,pdKK=diagKKt,p1=pdKK + *n;pdKK<p1;pdKK++,Tkm++) xx += *Tkm * *pdKK;*/
        for (xx=0.0,pdKK=diagKKt,p1=pdKK + *n;pdKK<p1;pdKK++,pTkm++) xx += *pTkm * *pdKK;
        det2[km] = xx;

        /* - tr(KTkKK'TmK) */
        det2[km] -= diagABt(work + *n * tid,KtTK + k * *r * *r,KtTK+ m * *r * *r,r,r);

        /* sp[k]*tr(P'S_kP) */
        if (k >= *M0 && k==m) det2[km] += trPtSP[m - *M0];

        /* -sp[m]*tr(K'T_kKP'S_mP) */
        if (m >= *M0) det2[km] -= sp[m - *M0]*diagABt(work + *n * tid,KtTK + k * *r * *r,PtSP + (m - *M0) * *r * *r,r,r);
     
        /* -sp[k]*tr(K'T_mKP'S_kP) */
        if (k >= *M0) det2[km] -= sp[k - *M0]*diagABt(work + *n * tid,KtTK + m * *r * *r,PtSP + (k - *M0) * *r * *r,r,r);
 
        /* -sp[m]*sp[k]*tr(P'S_kPP'S_mP) */
        if (k >= *M0 && m >= *M0) det2[km] -= sp[m - *M0]*sp[k - *M0]*
                                  diagABt(work + *n * tid,PtSP + (k - *M0) * *r * *r,PtSP + (m - *M0) * *r * *r,r,r);

        det2[mk] = det2[km];
      }     
    }
  } /* end of parallel section */
 
  /* free up some memory */
  if (deriv2) {FREE(PtSP);FREE(KtTK);}
  FREE(diagKKt);FREE(work);
  FREE(PtrSm);FREE(trPtSP);

} /* end get_ddetXWXpS */


void get_trA2(double *trA,double *trA1,double *trA2,double *P,double *K,double *sp,
	      double *rS,int *rSncol,double *Tk,double *Tkm,double *w,int *n,int *q,
              int *r,int *M,int *deriv,int *nt)

/* obtains trA and its first two derivatives wrt the log smoothing parameters 
   * P is q by r
   * K is n by r
   * U1 is q by r
   * this routine assumes that sp contains smoothing parameters, rather than log smoothing parameters.

   * If deriv is 0 then only tr(A) is obtained here.
   * This version uses only K and P, and is for the case where expressions involve weights which
     are reciprocal variances, not the squares of weights which are reciprocal standard deviations.

   * Note that tr(A) = tr(KK') and it is tempting to view diag(K'K) as giving the edfs
     of the parameters, but this seems to be wrong. It gives the edfs for R \beta, where
     R is (pseudo) inverse of P. 

     * uses nt threads via openMP. Assumes thread number already set on entry and nt already reset to
       1 if no openMP support.
*/

{ double *diagKKt,*diagKKtKKt,xx,*KtTK,*KtTKKtK,*KKtK,*KtK,*work,*pTk,*pTm,*pdKKt,*pdKKtKKt,*p0,*p1,*p2,*p3,*pd,
    *PtrSm,*PtSP,*KPtrSm,*diagKPtSPKt,*diagKPtSPKtKKt,*PtSPKtK, *KtKPtrSm, *KKtKPtrSm,*Ip,*IpK/*,lowK,hiK*/;
    int i,m,k,bt,ct,j,one=1,km,mk,*rSoff,deriv2,neg_w=0,tid=0;
  #ifdef OMP_REPORT
    Rprintf("get_trA2 (d=%d)...",*deriv);
  #endif
  if (*deriv==2) deriv2=1; else deriv2=0;
  /* Get the sign array for negative w_i */
  Ip = (double *)CALLOC((size_t)*n,sizeof(double));
  for (p0=w,p1=p0+ *n,p2=Ip;p0<p1;p0++,p2++) if (*p0 < 0) {*p2 = -1.0;neg_w=1;} else *p2 = 1.0;

  /* obtain tr(A) and diag(A) = diag(KK'Ip) */ 
  diagKKt = (double *)CALLOC((size_t)*n,sizeof(double));
  *trA = diagABt(diagKKt,K,K,n,r); 
  if (neg_w) { /* correct trA */
    for (*trA=0.0,p0=diagKKt,p1=p0 + *n,p2=Ip;p0<p1;p0++,p2++) *trA += *p2 * *p0;
  }
  if (!*deriv) {
    FREE(Ip);FREE(diagKKt); 
    #ifdef OMP_REPORT
    Rprintf("done\n");
    #endif 
    return;
  }

  /* set up work space */
  work =  (double *)CALLOC((size_t)*n * *nt,sizeof(double));
  /* Get K'IpK and KK'IpK  */
  KtK = (double *)CALLOC((size_t)*r * *r,sizeof(double));
  if (neg_w) { 
    IpK = (double *)CALLOC((size_t) *r * *n,sizeof(double));
    for (p0=IpK,p3=K,i=0;i<*r;i++) 
      for (p1=Ip,p2=p1 + *n;p1<p2;p1++,p0++,p3++) *p0 = *p1 * *p3; 
  } else { 
    IpK = (double *)CALLOC((size_t) *r * *n,sizeof(double));
    for (p0=IpK,p1=K,p2=K+ *n * *r;p1<p2;p0++,p1++) *p0 = *p1; 
    /*IpK = K;*/
  }
  /*  lowK=hiK=*K;

  for (p1=K,i=0;i<*n;i++) for (j=0;j<*r;j++,p1++) {
      if (*p1>hiK) hiK= *p1; else if (*p1<lowK) lowK = *p1;
    }
    Rprintf("K range = %g - %g\n",lowK,hiK);*/
  bt=1;ct=0;mgcv_pmmult(KtK,K,IpK,&bt,&ct,r,r,n,nt);  
  if (neg_w) FREE(IpK); else FREE(IpK);
  KKtK = (double *)CALLOC((size_t)*n * *r,sizeof(double));
  bt=0;ct=0;mgcv_pmmult(KKtK,K,KtK,&bt,&ct,n,r,r,nt);  

  /* obtain diag(KK'KK') */
  diagKKtKKt = (double *)CALLOC((size_t)*n,sizeof(double));
  xx = diagABt(diagKKtKKt,KKtK,K,n,r);
 
  /* now loop through the smoothing parameters to create K'TkK and K'TkKK'K */
  if (deriv2) {
    KtTK = (double *)CALLOC((size_t)(*r * *r * *M),sizeof(double));
    KtTKKtK = (double *)CALLOC((size_t)(*r * *r * *M),sizeof(double));
    #ifdef OPENMP_ON
    #pragma omp parallel private(k,j,tid) num_threads(*nt)
    #endif
    { /* open parallel section */
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (k=0;k < *M;k++) {
        #ifdef OPENMP_ON
        tid = omp_get_thread_num(); /* thread running this bit */
        #endif      
        j = k * *r * *r;
        getXtWX(KtTK+ j,K,Tk + k * *n,n,r,work + tid * *n);
        bt=ct=0;mgcv_mmult(KtTKKtK + k * *r * *r ,KtTK + j,KtK,&bt,&ct,r,r,r);
      }
    } /* parallel section end */
  } else { KtTK=KtTKKtK=(double *)NULL;}
  
  /* evaluate first and last terms in first derivative of tr(F) */
  bt=1;ct=0;mgcv_mmult(trA1,Tk,diagKKt,&bt,&ct,M,&one,n); /* tr(KK'Tk) */ 
  bt=1;ct=0;mgcv_mmult(work,Tk,diagKKtKKt,&bt,&ct,M,&one,n); /* tr(KK'TkKK') */
  for (i=0;i<*M;i++) trA1[i] +=  - work[i];

  /* now evaluate terms in Hessian of tr(F) which depend on what's available so far */
  if (deriv2) for (m=0;m < *M;m++) for (k=m;k < *M;k++){
     km=k * *M + m;mk=m * *M + k;

     /* tr(KK'Tkm  - KK'TkKK') */
     for (xx=0.0,pdKKt=diagKKt,pdKKtKKt=diagKKtKKt,p1=pdKKt + *n;pdKKt<p1;pdKKt++,pdKKtKKt++,Tkm++) 
          xx += *Tkm * (*pdKKt - *pdKKtKKt);
     trA2[km] = xx;

     /* -2 tr(K'TkKK'TmK)*/
     trA2[km] -= 2*diagABt(work,KtTK + k * *r * *r,KtTK+ m * *r * *r,r,r);

     /* 2 tr(K'TkKK'TmKK'K) -- needs correction*/
     xx = 2*diagABt(work,KtTK+k * *r * *r,KtTKKtK+m * *r * *r,r,r);
    
     trA2[km] += xx;

     trA2[mk] = trA2[km];     
  }

  /* free up some memory */
  if (deriv2) {FREE(KtTKKtK);FREE(KtTK);} 

  FREE(diagKKtKKt);FREE(diagKKt);

  /* create KP'rSm, KK'KP'rSm and P'SmP */
  PtrSm = (double *)CALLOC((size_t)(*r * *q * *nt),sizeof(double)); /* transient storage for P' rSm */
  KPtrSm = (double *)CALLOC((size_t)(*n * *q * *nt),sizeof(double)); /* transient storage for K P' rSm */
  diagKPtSPKt = (double *)CALLOC((size_t)(*n * *M),sizeof(double));
  if (deriv2) {
    PtSP = (double *)CALLOC((size_t)(*M * *r * *r ),sizeof(double));
    PtSPKtK = (double *)CALLOC((size_t)(*M * *r * *r ),sizeof(double));
    KtKPtrSm = (double *)CALLOC((size_t)(*r * *q * *nt),sizeof(double));/* transient storage for K'K P'rSm */ 
    KKtKPtrSm = (double *)CALLOC((size_t)(*n * *q * *nt),sizeof(double));/* transient storage for K'K P'rSm */ 
    diagKPtSPKtKKt = (double *)CALLOC((size_t)(*n * *M),sizeof(double));
  } else {  KKtKPtrSm=PtSPKtK= PtSP=KtKPtrSm=diagKPtSPKtKKt=(double *)NULL; }
  
  rSoff =  (int *)CALLOC((size_t)*M,sizeof(int));
  rSoff[0] = 0;for (m=0;m < *M-1;m++) rSoff[m+1] = rSoff[m] + rSncol[m];
  tid = 0;
  #ifdef OPENMP_ON
  #pragma omp parallel private(m,bt,ct,tid,xx,p0,p1,p2) num_threads(*nt)
  #endif
  { /* open parallel section */
    #ifdef OPENMP_ON
    #pragma omp for
    #endif
    for (m=0;m < *M;m++) { 
      #ifdef OPENMP_ON
      tid = omp_get_thread_num(); /* thread running this bit */
      #endif
      bt=1;ct=0;mgcv_mmult(PtrSm + *r * *q * tid,P,rS+rSoff[m] * *q,&bt,&ct,r,rSncol+m,q);
      bt=0;ct=0;mgcv_mmult(KPtrSm + *n * *q * tid,K,PtrSm + *r * *q * tid ,&bt,&ct,n,rSncol+m,r); 
      if (deriv2) {
        bt=0;ct=0;mgcv_mmult(KtKPtrSm + *r * *q * tid,KtK,PtrSm + *r * *q * tid,&bt,&ct,r,rSncol+m,r); 
        bt=0;ct=1;mgcv_mmult(PtSP+ m * *r * *r,PtrSm + *r * *q * tid,PtrSm + *r * *q * tid,&bt,&ct,r,r,rSncol+m);
    
        bt=0;ct=0;mgcv_mmult(KKtKPtrSm + *n * *q * tid,KKtK,PtrSm + *r * *q * tid,&bt,&ct,n,rSncol+m,r);      
        bt=0;ct=1;mgcv_mmult(PtSPKtK + m * *r * *r,PtrSm + *r * *q * tid,KtKPtrSm+ *r * *q * tid,&bt,&ct,r,r,rSncol+m); 
        xx = diagABt(diagKPtSPKtKKt+ m * *n,KPtrSm + *n * *q * tid,KKtKPtrSm + *n * *q * tid,n,rSncol+m);
      }
      xx = sp[m] * diagABt(diagKPtSPKt+ m * *n,KPtrSm + *n * *q * tid,KPtrSm + *n * *q * tid,n,rSncol+m);
      if (neg_w) { /* have to correct xx for negative w_i */
        for (xx=0.0,p0=diagKPtSPKt+m * *n,p1=p0 + *n,p2=Ip;p0<p1;p0++,p2++) xx += *p0 * *p2;
        xx *= sp[m];
      }
      trA1[m] -= xx; /* finishing trA1 */
      if (deriv2) trA2[m * *M + m] -=xx; /* the extra diagonal term of trA2 */
    } 
  } /* end of parallel section */
  FREE(rSoff);
  

  if (!deriv2) { /* trA1 finished, so return */
    FREE(PtrSm);FREE(KPtrSm);FREE(diagKPtSPKt);
    FREE(work);FREE(KtK);FREE(KKtK); 
    #ifdef OMP_REPORT
    Rprintf("done\n");
    #endif 
    return;
  }
  /* now use these terms to finish off the Hessian of tr(F) */ 
   for (m=0;m < *M;m++) for (k=m;k < *M;k++){
     km=k * *M + m;mk=m * *M + k;

     /* 2 sp[m] tr(KK'TkKP'SmPK') */
     pTk = Tk + k * *n;
     for (xx=0.0,pd = diagKPtSPKtKKt + m * *n,p1=pd + *n;pd < p1;pd++,pTk++) xx += *pd * *pTk;
     trA2[km] += 2*sp[m] *xx;

     /* 2 sp[k] tr(KK'TmKP'SkPK') */
     pTm = Tk + m * *n;
     for (xx=0.0,pd = diagKPtSPKtKKt + k * *n,p1=pd + *n;pd < p1;pd++,pTm++) xx += *pd * *pTm;
     trA2[km] += 2*sp[k] *xx;
     
     /* - sp[m] tr(TkKP'SmPK') */
     pTk = Tk + k * *n;
     for (xx=0.0,pd = diagKPtSPKt + m * *n,p1=pd + *n;pd < p1;pd++,pTk++) xx += *pd * *pTk;
     trA2[km] -= sp[m] *xx;
     
     /* - sp[k] tr(TmKP'SkPK') */
     pTm = Tk + m * *n;
     for (xx=0.0,pd = diagKPtSPKt + k * *n,p1=pd + *n;pd < p1;pd++,pTm++) xx += *pd * *pTm;
     trA2[km] -= sp[k] *xx;

     /* 2 sp[m] sp[k] tr(KP'SkPP'SmPK') */
     trA2[km] += 2 * sp[k]*sp[m]*diagABt(work,PtSP + m * *r * *r,PtSPKtK + k * *r * *r,r,r);
      
     trA2[mk] =trA2[km];
   } 
   /* clear up */
   FREE(PtrSm);FREE(KPtrSm);FREE(PtSP);FREE(KtKPtrSm);FREE(diagKPtSPKt);
   FREE(diagKPtSPKtKKt);FREE(work);FREE(KtK);FREE(KKtK);FREE(PtSPKtK);FREE(KKtKPtrSm);
   FREE(Ip);  
   #ifdef OMP_REPORT
    Rprintf("done\n");
   #endif 
} /* end get_trA2 */






void rc_prod(double *y,double *z,double *x,int *xcol,int *n)
/* obtains element-wise product of z and each of the  xcol columns of x,
   returning the result in y. z and the columns of x are n-vectors. 
   (equivalent to y = diag(z)%*%x)
*/ 
{ int i;
  double *pz,*pz1;
  pz1 = z + *n;
  for (i=0;i < *xcol;i++) 
      for (pz=z;pz<pz1;pz++,y++,x++) *y = *pz * *x;
}


void Rinv(double *Ri,double *R,int *c,int *r, int *ri)
/* invert c by c upper triangular matrix R, actually stored in upper 
   part of r by c matrix. Result returned in top of  Ri (actually ri by c).
*/

{ int i,j,k,eye;
  double xx,*rc;
  rc=Ri;
  for (i=0;i<*c;i++) {
      for (eye=1,k=i;k>=0;k--) {
	  for (xx=0.0,j=k+1;j <=i;j++) xx += R[k + j * *r] * rc[j];
          rc[k]=(eye-xx)/R[k + k * *r];
          eye=0;
      }
      for (k=i+1;k<*c;k++) rc[k]=0.0;
      rc += *ri;
  }
}


void pearson2(double *P, double *P1, double *P2,
              double *y,double *mu,double *V, double *V1,double *V2,double *g1,double *g2,
              double *p_weights,double *eta1, double *eta2,int n,int M,int deriv, int deriv2)
/* Alternative calculation of the derivatives of the Pearson statistic, which avoids assuming that
   z and w are based on Fisher scoring */
{ double resid,xx,*Pe1,*Pe2,*pp,*p1,*p0,*v2,*Pi1,*Pi2;
  int i,k,m,n_2dCols=0,one=1;
  if (deriv) {
    Pe1 = (double *)CALLOC((size_t)n,sizeof(double)); /* for dP/deta */
    Pi1 = (double *)CALLOC((size_t) n * M,sizeof(double)); /* for dPi/drho */
    if (deriv2) { 
      n_2dCols = (M * (1 + M))/2;
      Pe2 = (double *)CALLOC((size_t)n,sizeof(double)); /* for d2P/deta2 */
      v2 = (double *)CALLOC((size_t)n,sizeof(double));
      Pi2 = (double *)CALLOC((size_t)n_2dCols*n,sizeof(double)); /* for d2P_i/drho */
    } else {Pe2=v2=Pi2=NULL;}
  } else {Pi1 = Pe2 = v2 = Pe1 = Pi2 = NULL;}
  *P=0.0;
  for (i=0; i < n;i++) {
    resid = y[i]-mu[i];
    xx = resid*p_weights[i]/V[i];
    *P += xx*resid;
    if (deriv) {
      Pe1[i] = - xx* (2 + resid*V1[i])/g1[i];
      if (deriv2) {
        Pe2[i] = - Pe1[i]*g2[i]/g1[i] + 
	  (2*p_weights[i]/V[i]+2*xx*V1[i] - Pe1[i]*V1[i]*g1[i] - xx*resid*(V2[i]-V1[i]*V1[i]))/(g1[i]*g1[i]);
      }
    }
  } /* derivs wrt eta completed */

  if (deriv) { /* transform to derivs wrt rho */
    rc_prod(Pi1,Pe1,eta1,&M,&n); /* Pi1 = dP_i/drho_k done */
      if (deriv2) {  
        rc_prod(Pi2,Pe1,eta2,&n_2dCols,&n);
        for (pp=Pi2,m=0;m < M;m++) for (k=m;k < M;k++) {
	    rc_prod(Pe1,eta1 + n *  m,eta1 + n * k,&one,&n);
            rc_prod(v2,Pe2,Pe1,&one,&n);
            p1=v2 + n;
            for (p0=v2;p0<p1;p0++,pp++) *pp += *p0;        
        } /* Pi2 update completed */
      }
  } /* derivatives of Pi wrt rho completed */

  /* now sum the derivatives over i */

  if (deriv) {
    pp = Pi1;
    for (k=0;k<M;k++) { xx=0.0; for (i=0;i<n;i++,pp++) xx += *pp;P1[k] = xx;}
    if (deriv2) {
        for (pp=Pi2,m=0;m < M;m++) for (k=m;k < M;k++) {
          xx=0.0;
          for (i=0;i<n;i++,pp++) xx += *pp;
          P2[k*M+m] = P2[m*M+k] = xx;
        } 
    }
  } /* end of derivative summation */
  
  /* clear up */
  if (deriv) {
    FREE(Pe1);FREE(Pi1);
    if (deriv2) {
      FREE(Pe2);FREE(Pi2);FREE(v2);
    }
  }

} /* end pearson2 */



void applyP(double *y,double *x,double *R,double *Vt,int neg_w,int nr,int r,int c,int right)
/* Forms y = Px. If neg_w==0 P = R^{-1} otherwise P = R^{-1}V  where V is 
   transpose of Vt (r by r). x is r by c if right==0 and c br otherwise. 
   R is in the r by r upper triangle of nr by r array R. */
{ double *x1;
  int bt,ct;
  if (neg_w) { /* apply V */
    x1 = (double *)CALLOC((size_t)r*c,sizeof(double));
    if (right) {
      mgcv_backsolve(R,&nr,&r,x,x1, &c,&right);         /* x1 = x R^{-1} */
      bt=0;ct=1;mgcv_mmult(y,x1,Vt,&bt,&ct,&c,&r,&r);   /* y = x R^{-1} V */
    } else {
      bt=1;ct=0;mgcv_mmult(x1,Vt,x,&bt,&ct,&r,&c,&r);   /* x1 = V x */    
      mgcv_backsolve(R,&nr,&r,x1,y, &c,&right);         /* y = R^{-1} V x */
    }
    FREE(x1);
  } else mgcv_backsolve(R,&nr,&r,x,y, &c,&right);            /* y = R^{-1} x */
} /* applyP */

void applyPt(double *y,double *x,double *R,double *Vt,int neg_w,int nr,int r,int c,int right)
/* Forms y = P'x. If neg_w==0 P = R^{-1} otherwise P = R^{-1}V  where V is 
   transpose of Vt (r by r). x is r by c if right==0, otherwise x is c by r. 
   R is in the r by r upper triangle of nr by r array R. */
{ double *x1;
  int bt,ct;
  if (neg_w) { /* apply V */
    x1 = (double *)CALLOC((size_t)r*c,sizeof(double));
    if (right) {
      bt=0;ct=0;mgcv_mmult(x1,x,Vt,&bt,&ct,&c,&r,&r);   /* x1 = x V' */
      mgcv_forwardsolve(R,&nr,&r,x1,y, &c,&right);       /* y = x V' R^{-T} */
    } else {
      mgcv_forwardsolve(R,&nr,&r,x,x1, &c,&right);      /* x1 = R^{-T} x */
      bt=0;ct=0;mgcv_mmult(y,Vt,x1,&bt,&ct,&r,&c,&r);   /* y = V'R^{-T} x */    
    }
    FREE(x1);
  } else mgcv_forwardsolve(R,&nr,&r,x,y, &c,&right);            /* y = R^{-T} x */
} /* applyPt */


void ift1(double *R,double *Vt,double *X,double *rS,double *beta,double *sp,double *w,
         double *dwdeta,double *b1, double *b2,double *eta1,double *eta2,
	  int *n,int *r, int *M,int *rSncol,int *deriv2,int *neg_w,int *nr)

/* Uses the implicit function theorem to get derivatives of beta wrt rho = log(sp) cheaply
   without iteration, when Newton or Fisher-canonical-link are used... 
   X is n by r
   P is r by r, but in this version is defined by R and Vt, and computed
     with using `applyP' and `applyPt' 
   there are M smoothing parameters (unlogged) in sp
   beta is a q vector
   b1 is r by M
   b2 is r by n_2dCols 
*/
{ int n_2dCols,i,j,k,one=1,bt,ct;
  double *work,*Skb,*pp,*p0,*p1,*work1;
  work = (double *) CALLOC((size_t)*n,sizeof(double));
  work1 = (double *) CALLOC((size_t)*n,sizeof(double));
  Skb = (double *) CALLOC((size_t)*r,sizeof(double));
  n_2dCols = (*M * (1 + *M))/2;
  for (i=0;i<*M;i++) { /* first derivative loop */
    multSk(Skb,beta,&one,i,rS,rSncol,r,work); /* get S_i \beta */
    for (j=0;j<*r;j++) Skb[j] *= -sp[i]; 
    applyPt(work,Skb,R,Vt,*neg_w,*nr,*r,1,0);
    applyP(b1 + i * *r,work,R,Vt,*neg_w,*nr,*r,1,0);   
  } /* first derivatives of beta finished */

  bt=0;ct=0;mgcv_mmult(eta1,X,b1,&bt,&ct,n,M,r); /* first deriv of eta */

  if (*deriv2) { /* then second derivatives needed */
    pp = b2;   
    for (i=0;i<*M;i++) for (k=i;k<*M;k++) { 
      p0 = eta1 + *n * i;p1 = eta1 + *n * k;
      for (j=0;j<*n;j++,p0++,p1++) work[j] = - *p0 * *p1 * dwdeta[j];
      bt=1;ct=0;mgcv_mmult(Skb,X,work,&bt,&ct,r,&one,n); /* X'f */
      multSk(work,b1+k* *r,&one,i,rS,rSncol,r,work1); /* get S_i dbeta/drho_k */
      for (j=0;j<*r;j++) Skb[j] += -sp[i]*work[j];
      multSk(work,b1+i* *r,&one,k,rS,rSncol,r,work1); /* get S_k dbeta/drho_i */
      for (j=0;j<*r;j++) Skb[j] += -sp[k]*work[j];
      applyPt(work,Skb,R,Vt,*neg_w,*nr,*r,1,0);
      applyP(pp,work,R,Vt,*neg_w,*nr,*r,1,0);
      if (i==k) for (j=0;j< *r;j++) pp[j] += b1[i * *r + j];
      pp += *r;
    }

    bt=0;ct=0;mgcv_mmult(eta2,X,b2,&bt,&ct,n,&n_2dCols,r); /* second derivatives of eta */
  }

  FREE(work);FREE(Skb);FREE(work1);
} /* end ift1 */

void ift2(double *R,double *Vt,double *X,double *rS,double *beta,double *sp,double *theta,
          double *Det_th,double *Det2_th,double *Det3,double *Det_th2,double *b1, double *b2,
          double *eta1,double *eta2,
	  int *n,int *r, int *M,int *n_theta,int *rSncol,int *deriv2,int *neg_w,int *nr)

/* Uses the implicit function theorem to get derivatives of beta wrt rho = log(sp) 
   and theta (extra parameters of likelihood), for extended GAMs. 

   X is n by r
   P is r by r, but in this version is defined by R and Vt, and computed
     with using `applyP' and `applyPt' 
   there are M smoothing parameters (unlogged) in sp
   there are n_theta theta parameters.

   beta is a q vector
   b1 is r by (M+n_theta)
   b2 is r by n_2dCols 
*/
{ int n_2dCols,i,j,k,one=1,bt,ct,ntot,kk;
  double *work,*Db_th,*pp,*p0,*p1,*work1;
  work = (double *) CALLOC((size_t)*n,sizeof(double));
  work1 = (double *) CALLOC((size_t)*n,sizeof(double));
  Db_th = (double *) CALLOC((size_t)*r,sizeof(double));
  ntot = *M + *n_theta;
  n_2dCols = (ntot  * (1 + ntot))/2;

  /* Db_th is d^2 pen dev/d beta d theta where theta is theta or log sp */
  for (i=0;i < ntot ;i++) { /* first derivative loop */
    if (i < *n_theta) { /* Db_th is from direct dependence of deviance on theta */
      bt=1;ct=0;mgcv_mmult(Db_th,X,Det_th + i * *n,&bt,&ct,r,&one,n);
      for (j=0;j<*r;j++) Db_th[j] *= -.5;
    } else { /* Db_th is from dependence of penalty on sp */
      multSk(Db_th,beta,&one,i - *n_theta,rS,rSncol,r,work); /* get S_i \beta */
      for (j=0;j<*r;j++) Db_th[j] *=  - sp[i - *n_theta];
    } 
    /* note that PPt = (X'WX+S)^{-1} i.e. *twice* the inverse Hessian,
       hence the factors of 0.5 introduced in Db_th above */ 
    applyPt(work,Db_th,R,Vt,*neg_w,*nr,*r,1,0);
    applyP(b1 + i * *r,work,R,Vt,*neg_w,*nr,*r,1,0);   
  } /* first derivatives of beta finished */

  bt=0;ct=0;mgcv_mmult(eta1,X,b1,&bt,&ct,n,&ntot,r); /* first deriv of eta */

  if (*deriv2) { /* then second derivatives needed */
    pp = b2; /* pointer to the second derivative of beta array */
    kk = 0; /* column counter for the second derivative of D wrt theta arrays */
    for (i=0;i<ntot;i++) for (k=i;k<ntot;k++) { 
      /* first term */
      p0 = eta1 + *n * i;p1 = eta1 + *n * k; 
      for (j=0;j<*n;j++,p0++,p1++) work[j] = - *p0 * *p1 * Det3[j];
      bt=1;ct=0;mgcv_mmult(Db_th,X,work,&bt,&ct,r,&one,n); /* D_bbb^lpq db_q/dtheta_j db_p/dtheta_k */
      /* second term */
      if (k < *n_theta) { 
        p0 = Det2_th + *n * k;
        p1 = eta1 + *n * i; 
        for (j=0;j<*n;j++,p0++,p1++) work[j] = *p0 * *p1 ;
        bt=1;ct=0;mgcv_mmult(work1,X,work,&bt,&ct,r,&one,n); 
      } else {
        multSk(work1,b1+ i * *r,&one,k - *n_theta,rS,rSncol,r,work);
        for (j=0;j<*r;j++) work1[j] *=  sp[k - *n_theta] * 2; 
      }
      for (j=0;j<*r;j++) Db_th[j] -=  work1[j];

      /* third term */
      if (i < *n_theta) { 
        p0 = Det2_th + *n * i;
        p1 = eta1 + *n * k; 
        for (j=0;j<*n;j++,p0++,p1++) work[j] = *p0 * *p1 ;
        bt=1;ct=0;mgcv_mmult(work1,X,work,&bt,&ct,r,&one,n); 
      } else {
        multSk(work1,b1+ k * *r,&one,i - *n_theta,rS,rSncol,r,work);
        for (j=0;j<*r;j++) work1[j] *=  sp[i - *n_theta] * 2; 
      }
      for (j=0;j<*r;j++) Db_th[j] -=  work1[j];
      /* final term */
      if (i< *n_theta && k < *n_theta) {
        p0 = Det_th2 + *n * kk;
        bt=1;ct=0;mgcv_mmult(work,X,p0,&bt,&ct,r,&one,n);
        for (j=0;j<*r;j++) Db_th[j] -=  work[j];
	kk++;
      } else if (i==k) {
        multSk(work1,beta,&one,i - *n_theta,rS,rSncol,r,work); /* get S_i \beta */
        for (j=0;j<*r;j++) Db_th[j] -= work1[j] * sp[i - *n_theta] *2 ;
      }
      for (j=0;j<*r;j++) Db_th[j] *= .5; /* since PPt is twice inv Hessian */
      applyPt(work,Db_th,R,Vt,*neg_w,*nr,*r,1,0);
      applyP(pp,work,R,Vt,*neg_w,*nr,*r,1,0);
      pp += *r;
    }
    bt=0;ct=0;mgcv_mmult(eta2,X,b2,&bt,&ct,n,&n_2dCols,r); /* second derivatives of eta */
  }

  FREE(work);FREE(Db_th);FREE(work1);
} /* end ift2 */




void drop_cols(double *X, int r, int c,int *drop, int n_drop) 
/* Routine to drop the columns in X indexed by drop (*ascending* order) 
   Result returned in X.
*/ 
{ int k,j,j0,j1;  
  double *p,*p1,*p2;
  if (n_drop<=0) return;
  if (n_drop) { /* drop the unidentifiable columns */
    for (k=0;k<n_drop;k++) {
      j = drop[k]-k; /* target start column */
      j0 = drop[k]+1; /* start of block to copy */
      if (k<n_drop-1) j1 = drop[k+1]; else j1 = c; /* end of block to copy */
      for (p=X + j * r,p1=X + j0 * r,p2=X + j1 * r;p1<p2;p++,p1++) *p = *p1;
    }      
  }
}

void drop_rows(double *X,int r, int c,int *drop,int n_drop)
/* Drops rows indexed by drop from X, returning result packed in 
   r-n_drop by c matrix X. `drop' *must* be in ascending order */
{ int i,j,k;
  double *Xs;
  if (n_drop<=0) return;
  Xs=X;
  for (j=0;j<c;j++) { /* work across columns */ 
    for (i=0;i<drop[0];i++,X++,Xs++) *X = *Xs;  
    Xs++;
    for (k=1;k<n_drop;k++) { 
      for (i=drop[k-1]+1;i<drop[k];i++,X++,Xs++) *X = *Xs;
      Xs++;
    }
    for (i=drop[n_drop-1]+1;i<r;i++,X++,Xs++) *X = *Xs;  
  }
}


void undrop_rows(double *X,int r,int c,int *drop,int n_drop)
/* Inserts extra zero rows in X in the rows indicated by drop,
   and shifts the others up accordingly. So, X ends up r by c, with 
   zero rows in the positions given in drop.

   The assumption is that X is densely packed as (r-n_drop) by c on 
   entry.
   
   `drop' *must* be in ascending order.
*/
{ double *Xs;
  int i,j,k;
  if (n_drop <= 0) return;
  Xs = X + (r-n_drop)*c - 1; /* position of the end of input X */
  X += r*c - 1;              /* end of final X */
  for (j=c-1;j>=0;j--) { /* back through columns */
    for (i=r-1;i>drop[n_drop-1];i--,X--,Xs--) *X = *Xs;
    *X = 0.0;X--;
    for (k=n_drop-1;k>0;k--) {
      for (i=drop[k]-1;i>drop[k-1];i--,X--,Xs--) *X = *Xs;
      *X = 0.0;X--;
    }
    for (i=drop[0]-1;i>=0;i--,X--,Xs--) *X = *Xs;
  }
} /* end undrop rows */




double MLpenalty1(double *det1,double *det2,double *Tk,double *Tkm,double *nulli, double *X,
double *R,double *Q, int *nind,double *sp,double *rS,int *rSncol,int *q,int *n,
		  int *Ms,int *M,int *M0,int *neg_w,double *rank_tol,int *deriv,
                  int *nthreads,int *type) {
/* Routine to obtain the version of log|X'WX+S| that applies to ML, rather than REML.
   This version assumes that we are working in an already truncated range-null separated space.

   * nulli is an array indicating whether a parameter (column) relates to the null 
     space (+ve) or range space (-ve) of the total penalty matrix. Because of pivoting
     they can be in any order.
   * Q, R are the QR factors of diag(abs(W))X augmenented by the square root of S
   * nind is the array indexing the locations of the `neg_w' -ve elements of W.
   * q is the number of model coefficients
   * Ms is the penalty null space dimension.
   * M is number of smoothing parameters, and M0 the number of theta parameters.
   * n is the number of rows in Q. 

   Basic task of the routine is to project Hessian of the penalized log likelihood 
   into the range space of the penalty, in order to obtain the correction term that 
   applies for ML.

   NOTE: rS is over-written by this. 
*/

  double *RU1,*tau,*work,*Ri,*Qb=NULL,*K,*P,*IQ,*IQQ,*Vt,*XU1=NULL,
         *d,*p0,*p1,*p2,*p3,ldetXWXS,ldetI2D=0.0;
  int ScS,bt,ct,qM,*pivot,i,j,k,left,tp,n_drop=0,*drop,FALSE=0; 

  drop = (int *)CALLOC((size_t)*Ms,sizeof(int));
  for (i=0;i < *q;i++) if (nulli[i]>0.0) { drop[n_drop] = i;n_drop++; }

  for (ScS=0.0,i=0;i<*M;i++) ScS += rSncol[i]; /* total columns of rS */

  qM = *q - n_drop;

  RU1 = (double *)CALLOC((size_t) *q * *q ,sizeof(double));
  for (p1=RU1,p2=R,p3=R+ *q * *q;p2 < p3;p1++,p2++) *p1 = *p2;
 
  drop_cols(RU1,*q,*q,drop,n_drop); /* drop the null space columns from R */ 

  /* A pivoted QR decomposition of RU1 is needed next */
  tau=(double *)CALLOC((size_t)qM,sizeof(double)); /* part of reflector storage */
  pivot=(int *)CALLOC((size_t)qM,sizeof(int));
  
  mgcv_qr(RU1,q,&qM,pivot,tau); /* RU1 and tau now contain the QR decomposition information */
  /* pivot[i] gives the unpivoted position of the ith pivoted parameter.*/
  
  /* Ri needed */

  Ri =  (double *)CALLOC((size_t) qM * qM,sizeof(double)); 
  Rinv(Ri,RU1,&qM,q,&qM); /* getting R^{-1} */
  
  if (*type==0||*neg_w) { /* new Q factor needed explicitly */
    Qb = (double *)CALLOC((size_t) *q * qM,sizeof(double)); 
    for (i=0;i< qM;i++) Qb[i * *q + i] = 1.0;
    left=1;tp=0;mgcv_qrqy(Qb,RU1,tau,q,&qM,&qM,&left,&tp); /* Q from the QR decomposition */
  } else { /* need X with null space cols dropped */ 
    XU1 = (double *)CALLOC((size_t) *n * *q,sizeof(double));
    for (p1=XU1,p2=X,p3=X + *n * *q;p2 < p3;p1++,p2++) *p1 = *p2;
    drop_cols(XU1,*n,*q,drop,n_drop); /* drop the null space columns from X */ 
  }
  FREE(tau);

  K = (double *)CALLOC((size_t) *n * qM,sizeof(double));
  P = (double *)CALLOC((size_t) qM * qM,sizeof(double));

  if (*neg_w) { /* need to deal with -ve weight correction */
    if (*neg_w < *q+1) k = *q+1; else k = *neg_w;
    IQ = (double *)CALLOC((size_t) k * *q,sizeof(double)); 
    for (i=0;i< *neg_w;i++) { /* Copy the rows of Q corresponding to -ve w_i into IQ */
      p0 = IQ + i;p1 = Q + nind[i];
      for (j=0;j<*q;j++,p0+=k,p1+= *n) *p0 = *p1;
    }
    /* Note that IQ may be zero padded, for convenience */
    IQQ = (double *)CALLOC((size_t) k * qM,sizeof(double)); 
    bt=0;ct=0;mgcv_mmult(IQQ,IQ,Qb,&bt,&ct,&k,&qM,q); /* I^-Q_1 \bar Q is k by rank */
    FREE(IQ);
     
    /* Get the SVD of IQQ */
    Vt = (double *)CALLOC((size_t) qM * qM,sizeof(double));
    d = (double *)CALLOC((size_t) qM,sizeof(double));
    mgcv_svd_full(IQQ,Vt,d,&k,&qM); /* SVD of IQ */
    FREE(IQQ);
    for (i=0;i<qM;i++) {
      d[i] = 1 - 2*d[i]*d[i];
      if (d[i]<=0) d[i]=0.0; 
      else {
        ldetI2D += log(d[i]); /* log|I-2D^2| */ 
        d[i] = 1/sqrt(d[i]);
      }
    } /* d now contains diagonal of diagonal matrix (I-2D^2)^{-1/2} (possibly pseudoinverse) */
    /* Now form (I-2D^2)^.5 Vt and store in Vt... */
    for (p0=Vt,i=0;i<qM;i++)
    for (p1=d,p2=d+qM;p1<p2;p1++,p0++) *p0 *= *p1;
    
    /* Form P */
    bt=0;ct=1;mgcv_mmult(P,Ri,Vt,&bt,&ct,&qM,&qM,&qM);

    /* Form K */
    if (*type==1) {
      bt=0;ct=0;mgcv_mmult(K,XU1,P,&bt,&ct,n,&qM,&qM); 
    } else {
      work = (double *)CALLOC((size_t) *q * qM,sizeof(double));
      bt=0;ct=1;mgcv_mmult(work,Qb,Vt,&bt,&ct,q,&qM,&qM); /* \bar Q V (I - 2D^2)^.5 */

      bt=0;ct=0;mgcv_mmult(K,Q,work,&bt,&ct,n,&qM,q);
      FREE(work);
    }
    FREE(d);FREE(Vt);   
    
  } else { /* no negative weights, so P and K can be obtained directly */
    ldetI2D = 0.0;
    /* Form P */
    for (p0=P,p1=Ri,p2=Ri+ qM * qM;p1<p2;p0++,p1++) *p0 = *p1; /* copy R^{-1} into P */
    /* Form K */
    if (*type==1) {
      bt=0;ct=0;mgcv_mmult(K,XU1,P,&bt,&ct,n,&qM,&qM);
    } else {
      bt=0;ct=0;mgcv_mmult(K,Q,Qb,&bt,&ct,n,&qM,q);
    }
  }

  FREE(Ri);

  /* Evaluate the required log determinant... */

  for (ldetXWXS=0.0,i=0;i<qM;i++) ldetXWXS += log(fabs(RU1[i + i * *q])); 
  ldetXWXS *= 2;
 
  ldetXWXS += ldetI2D; /* the negative weights correction */
  
  FREE(RU1);

  /* rS also needs to have null space parts dropped, and to be pivoted... */

  drop_rows(rS,*q,ScS,drop,n_drop);   /* rS now rank by ScS */ 
  pivoter(rS,&qM,&ScS,pivot,&FALSE,&FALSE); /* row pivot of rS */
  
  if (*type==0||*neg_w) FREE(Qb); else FREE(XU1);
  FREE(pivot);

  /* Now we have all the ingredients to obtain required derivatives of the log determinant... */
  
  if (*deriv)
    get_ddetXWXpS(det1,det2,P,K,sp,rS,rSncol,Tk,Tkm,n,&qM,&qM,M,M0,deriv,*nthreads);

  FREE(P);FREE(K);FREE(drop);
  return(ldetXWXS);
} /* end of MLpenalty1 */


int icompare (const void * a, const void * b)
/* integer comparison function for qsort */ 
{
  return ( *(int*)a - *(int*)b );
}



void gdiPK(double *work,double *X,double *E,double *Es,double *rS,double *U1,double *z,double *raw,double *R,
           double *Rh,double *nulli,double *dev_hess,double *P, double *K,double *Vt,double *PKtz,double *Q1,
           int *nind,int *pivot1,int *drop,
           int *n,int *q,int *Mp,int neg_w, int *nt,int *Enrow,int *rank,int *n_drop,int deriv2,int ScS, int *REML,
      double *rank_tol,double *ldetXWXS,int *type)
/* Does initial QR decomposition for gdi routines work must be at least 3*q 
   If type==1 uses wz instead of sqrt(|w|)z and returns K s.t. \hat \beta = PK'Wz. In this case
   z must contain Wz on entry.
*/
{ int i,j,k,*pivot,nt1,nr,left,tp,bt,ct,TRUE=1,FALSE=0,one=1;
  double *zz=NULL,*WX,*tau,*R1,Rnorm,Enorm,Rcond,*Q=NULL,*tau1,*Ri,ldetI2D,*IQ,*d,*p0,*p1,*p2,*p3,*p4,xx,norm1,norm2;
  nt1 = *nt;
  
  if (*type==0) {
    zz = (double *)CALLOC((size_t)*n,sizeof(double)); /* storage for z=[sqrt(|W|)z,0] */
    for (i=0;i< *n;i++) { 
      zz[i] = z[i]*raw[i]; /* form z itself */
    }
    for (i=0;i<neg_w;i++) { k=nind[i]; zz[k] = -zz[k];} 
  } /* otherwise z must actually contain wz */

  WX = (double *) CALLOC((size_t) ( (*n + *nt * *q) * *q),sizeof(double));
  for (j=0;j<*q;j++) 
  { for (i=0;i<*n;i++) /* form WX */
    { k = i + *n * j;
      WX[k]=raw[i]*X[k];
    }
  }
  /* get the QR decomposition of WX */
 
  tau=(double *)CALLOC((size_t) *q * (*nt + 1),sizeof(double)); /* part of reflector storage */

  pivot=(int *)CALLOC((size_t)*q,sizeof(int));
  
  mgcv_pqr(WX,n,q,pivot,tau,&nt1); /* WX and tau now contain the QR decomposition information */

  /* pivot[i] gives the unpivoted position of the ith pivoted parameter.*/
  
  /* copy out upper triangular factor R, and unpivot it */
  R1 = (double *)CALLOC((size_t)*q * *q,sizeof(double));
 
  getRpqr(R1,WX,n,q,q,&nt1);

  pivoter(R1,q,q,pivot,&TRUE,&TRUE); /* unpivoting the columns of R1 */
 
  /* Form a nicely scaled version of [R',Es']' for rank determination */ 
  Rnorm = frobenius_norm(R1,q,q);
  Enorm =  frobenius_norm(Es,Enrow,q);
  nr = *q + *Enrow;
  for (j=0;j<*q;j++) { 
    for (i=0;i< *q;i++) R[i + nr * j] = R1[i + *q * j]/Rnorm;
    for (i=0;i< *Enrow;i++) R[i + *q + nr * j] = Es[i + *Enrow * j]/Enorm;
  }
  
  /* ... and now use it to establish rank */
   
  tau1=(double *)CALLOC((size_t)*q,sizeof(double)); /* part of reflector storage */

  mgcv_qr(R,&nr,q,pivot1,tau1);
  
  /* now actually find the rank of R */
  *rank = *q;
  R_cond(R,&nr,rank,work,&Rcond);
  while (*rank_tol * Rcond > 1) { (*rank)--;R_cond(R,&nr,rank,work,&Rcond);}

  /* Now have to drop the unidentifiable columns from R1, E and the corresponding rows from rS
     The columns to drop are indexed by the elements of pivot1 from pivot1[rank] onwards.
     Before returning, zeros will need to be inserted in the parameter vector at these locations. 
  */


  for (i=0;i<*q - *Mp;i++) nulli[i] = -1.0;    /* parameter in penalty range space */
  for (i= *q - *Mp;i < *q;i++) nulli[i] = 1.0; /* parameter in penalty null space */ 
  *n_drop = *q - *rank;
  if (*n_drop) {
    for (i=0;i < *n_drop;i++) drop[i] = pivot1[*rank+i];
    qsort(drop,*n_drop,sizeof(int),icompare); /* key assumption of the drop/undrop routines is that `drop' is ascending */
    /* drop columns indexed in `drop'... */
    drop_cols(R1,*q,*q,drop,*n_drop);    /* R1 now q by rank */
    drop_cols(E,*Enrow,*q,drop,*n_drop); /* E now q by rank */ 
    drop_cols(X,*n,*q,drop,*n_drop);     /* X now n by rank */
    drop_rows(rS,*q,ScS,drop,*n_drop);   /* rS now rank by ScS */ 
    drop_rows(nulli,*q,1,drop,*n_drop);  /* keeps track of null space params */
  }

  /* At this stage the parameter space has been purged of terms that are
     theoretically unidentifiable, given WX and the penalties */


  /* Now augment R1 with the real square root penalty (not the nicely scaled version), result in R... */  
  for (j=0;j < *rank;j++) { 
    for (i=0;i< *q;i++) R[i + nr * j] = R1[i + *q * j];
      for (i=0;i< *Enrow;i++) R[i + *q + nr * j] = E[i + *Enrow * j];
  }
   
  mgcv_qr(R,&nr,rank,pivot1,tau1); /* The final QR decomposition */ 
  
  i=1;pivoter(nulli,rank,&i,pivot1,&FALSE,&FALSE); /* pivoting the rows of nulli */

  if (deriv2) { /* get first bit of X'WX (hessian of the deviance)*/
    pivoter(R1,q,rank,pivot1,&TRUE,&FALSE); /* pivot the columns of R1 */
    getXtX(dev_hess,R1,q,rank);    
  } 
    
  /* Form Q1 = Qf Qs[1:q,] where Qf and Qs are orthogonal factors from first and final QR decomps
     respectively ... */
  if (neg_w || *type==0) { /* Q1 needed if neg_w correction needed, and anyway for type==0 */
    Q = (double *)CALLOC((size_t) nr * *rank,sizeof(double)); 
    for (i=0;i < *rank;i++) Q[i * nr + i] = 1.0;
    left=1;tp=0;mgcv_qrqy(Q,R,tau1,&nr,rank,rank,&left,&tp); /* Q from the second QR decomposition */

    /* Q1 = Qb Q[1:q,]  where Qb from first QR decomposition... */
    for (i=0;i<*q;i++) for (j=0;j < *rank;j++) Q1[i + *q * j] = Q[i + nr * j];
    tp=0;mgcv_pqrqy(Q1,WX,tau,n,q,rank,&tp,&nt1);
    /* so, at this stage WX = Q1 R, dimension n by rank */
  }

  Ri =  (double *)CALLOC((size_t) *rank * *rank,sizeof(double)); 
  Rinv(Ri,R,rank,&nr,rank); /* getting R^{-1} */
  

  ldetI2D = 0.0; /* REML determinant correction */
  if (neg_w) { /* then the correction for the negative w_i has to be evaluated */
    if (neg_w < *rank + 1) k = *rank + 1; else k = neg_w;
    IQ = (double *)CALLOC((size_t) k * *rank,sizeof(double)); 
    for (i=0;i < neg_w;i++) { /* Copy the rows of Q corresponding to -ve w_i into IQ */
      p0 = IQ + i;p1 = Q1 + nind[i];
      for (j=0;j < *rank;j++,p0+=k,p1+= *n) *p0 = *p1;
    }
    /* Note that IQ may be zero padded, for convenience */
   
    d = (double *)CALLOC((size_t) *rank,sizeof(double));
    mgcv_svd_full(IQ,Vt,d,&k,rank); /* SVD of IQ */
    FREE(IQ);

    if (deriv2) { /* correct the Hessian of the deviance */
      /* put DV'R into P, temporarily */
      p1=P;
      for (j=0;j < *rank;j++,p1 += *rank) {
        p0 = R + j * nr; /* start of column j of R */
        for (p2=Vt,p3=p1,p4=p1 + *rank;p3<p4;p3++,p2++) *p3 = *p2 * *p0;
        p0++;
        for (k=1;k<=j;k++,p0++) 
	  for (p3=p1;p3<p4;p3++,p2++) *p3 += *p2 * *p0;
      } /* end of loop to form V'R */
      /* Now form DV'R */
      for (p0=P,j=0;j < *rank;j++) for (p1=d,p2=d + *rank;p1<p2;p1++,p0++) *p0 *= *p1; 
      /* Form K = R'VDDV'R --- the correction factor for X'WX */
      getXtX(K,P,rank,rank);
      for (p0=dev_hess,p1=p0+ *rank * *rank,p2=K;p0<p1;p0++,p2++) *p0 += -2 * *p2;
    }

    for (i=0;i < *rank;i++) {
      work[i] = d[i] = 1 - 2*d[i]*d[i];
      if (d[i]<=0) work[i] = d[i] = 0.0; 
      else {
        ldetI2D += log(d[i]); /* log|I-2D^2| */ 
        work[i] = sqrt(d[i]);
        d[i] = 1/work[i];
      }
    } /* d now contains diagonal of diagonal matrix (I-2D^2)^{-1/2} (possibly pseudoinverse) */
    
    /* copy uppertriangle of R (nr by rank) temporarily to P (rank by rank)... */
    for (p0=P,i=0;i < *rank;i++) {
      for (p1 = R + i*nr,p2 = p1 + i;p1 <= p2;p1++,p0++) *p0 = *p1;
      p2 = p0 + *rank - i - 1;
      for (;p0<p2;p0++) *p0 = 0.0;
    } 
    /* form Rh = Vt R... */
    bt=0;ct=0;mgcv_pmmult(Rh,Vt,P,&bt,&ct,rank,rank,rank,nt);
    /* finally Rh = (I-2D^2)^.5 Vt R ... */
    for (p0=Rh,i=0;i < *rank;i++)
    for (p1=work,p2=work + *rank;p1<p2;p1++,p0++) *p0 *= *p1;
    /* ... so now Rh'Rh = XWX+S */

    /* Now form (I-2D^2)^-.5 Vt and store in Vt... */
    for (p0=Vt,i=0;i < *rank;i++)
    for (p1=d,p2=d + *rank;p1<p2;p1++,p0++) *p0 *= *p1;

    /* Form P */
    bt=0;ct=1;mgcv_pmmult(P,Ri,Vt,&bt,&ct,rank,rank,rank,nt);
    
    /* Form K - or postpone if type!=0*/
    if (*type==0) mgcv_pmmult(K,Q1,Vt,&bt,&ct,n,rank,rank,nt);
  
    FREE(d);   
  } else { /* no negative weights so Rh, P and K much simpler */
    if (*type==0) { /* Form K now */
      for (p0=K,p1=Q1,j=0;j< *rank;j++,p1 += *n) /* copy just Q1 into K */
      for (p2 = p1,p3=p1 + *n;p2<p3;p0++,p2++) *p0 = *p2;
    } 
    /* Form P */
    for (p0=P,p1=Ri,j=0;j < *rank;j++,p0+= *rank) /* copy R^{-1} into P */
    for (p2=p0,p3=p0 + *rank;p2<p3;p1++,p2++) *p2 = *p1;
    /* form Rh s.t. Rh'Rh = XWX+S, by simply copying R... */  
    for (p0=Rh,i=0;i < *rank;i++) {
      for (p1 = R + i*nr,p2 = p1 + i;p1 <= p2;p1++,p0++) *p0 = *p1;
      p2 = p0 + *rank - i - 1;
      for (;p0<p2;p0++) *p0 = 0.0;
    }       
  }

  /* At this stage P is complete and K is complete for type==0 */

  if (*REML>0) {  
    for (*ldetXWXS=0.0,i=0;i < *rank;i++) *ldetXWXS += log(fabs(R[i + i * nr])); 
    *ldetXWXS *= 2;
    *ldetXWXS += ldetI2D; /* correction for negative weights */
  }


  /* Apply pivoting to the parameter space - this simply means reordering the cols of E and X and the 
     rows of the rS_i, and then unscrambling the parameter vector at the end (along with any covariance matrix)
     pivot1[i] gives the unpivoted position of the ith pivoted parameter.
  */

  pivoter(rS,rank,&ScS,pivot1,&FALSE,&FALSE); /* row pivot of rS */
  
  pivoter(E,Enrow,rank,pivot1,&TRUE,&FALSE);  /* column pivot of E */  

  pivoter(X,n,rank,pivot1,&TRUE,&FALSE);  /* column pivot of X */ 
 
  if (*type==1) { /* create K = XP... */
    applyP(K,X,R,Vt,neg_w,nr,*rank,*n,1);
  } else {
    /* start PK'z --- the pivoted coefficients...*/
    bt=1;ct=0;mgcv_mmult(work,K,zz,&bt,&ct,rank,&one,n); /* K'z */ 
    mgcv_mmult(work + *q *2,Q1,zz,&bt,&ct,rank,&one,n); /* Q1'z */
  }
  
  if (*type==1) { 
    bt=1;ct=0;mgcv_mmult(work,K,z,&bt,&ct,rank,&one,n); /* K'Wz */
    applyP(PKtz,work,R,Vt,neg_w,nr,*rank,1,0);
  } else {
    /* Create Wz (not sqrt(|W|)z)... */
    for (i=0;i<*n;i++) zz[i] = raw[i] * raw[i] * z[i]; 
    for (i=0;i<neg_w;i++) { k=nind[i]; zz[k] = -zz[k];} 
    bt=1;ct=0;mgcv_mmult(work + *q,X,zz,&bt,&ct,rank,&one,n); /* X'Wz */
    for (norm1=norm2=0,i=0;i < *rank;i++) {
      /* form R'Q1'z one row at a time */
      for (xx=0,j=0;j<=i;j++) xx += R[j + nr * i] * work[*q * 2 + j];
      xx -= work[i + *q]; /* difference to X'Wz */ 
      norm1 += xx*xx; norm2 += work[i + *q] * work[i + *q];   
    }
    if (norm1 > *rank_tol * norm2) {
      //Rprintf("gdi2 instability detected norm1= %g norm2 = %g\n",norm1,norm2);
      applyPt(zz,work + *q,R,Vt,neg_w,nr,*rank,1,0); /* P'X'Wz */
      applyP(PKtz,zz,R,Vt,neg_w,nr,*rank,1,0);
    } else applyP(PKtz,work,R,Vt,neg_w,nr,*rank,1,0);
  }

  FREE(WX);FREE(tau);FREE(Ri);FREE(R1); 
  FREE(tau1);
  if (neg_w || *type==0) FREE(Q); 
  FREE(pivot);
  if (*type==0) FREE(zz);
} /* gdiPK */



void gdi2(double *X,double *E,double *Es,double *rS,double *U1,
	  double *sp,double *theta,double *z,double *w,double *wz,double *wf,
          double *Dth,double *Det,double *Det2,double *Dth2,double *Det_th,
          double *Det2_th,double *Det3,double *Det_th2,
          double *Det4, double *Det3_th, double *Det2_th2,
          double *beta,double *b1,double *w1,
          double *D1,double *D2,double *P0,double *P1,double *P2,
          double *ldet, double *ldet1,double *ldet2,double *rV,
          double *rank_tol,int *rank_est,
	  int *n,int *q, int *M,int *n_theta, int *Mp,int *Enrow,int *rSncol,int *deriv,
	  int *fixed_penalty,int *nt,int *type,double *dVkk)     
/* Extended GAM derivative function, for independent data beyond exponential family.
   
   On entry *ldet < 0 indicates that ML ingredients should be computed, else REML 

   type == 0 is the original computation involving (|w|)^{-1}dw/drho whereas type == 1
   avoids the (|w|)^{-1} to avoid problems with zero weights (and/or badly scaled sqrt(|w|)z). 

   Identifiability truncation is based on the "well scaled" 
   penalty square root, Es, and is assuming that a stability enhancing 
   reparameterization and stable E are being employed.

   This version deals properly with negative weights, which can occur with Newton based 
   PIRLS. In consequence w's in this routine are proportional to reciprocal variances,
   not reciprocal standard deviations.
  
   The function is to be called at convergence of a P-IRLS scheme, estimating model coefficients 
   by P-IRLS.

   All names ending in 1,2 or 3 are derivatives of some sort, with the integer
   indicating the order of differentiation. 

   The arguments of this function point to the following:
   *i X is and n by q model matrix. On output this will contain K.
   *i E is a q by Enrow square root of the total penalty matrix, so E'E=S
   *i Es is the square root of a "well scaled" version of the total penalty,
     suitable for numerical determination of the theoretical rank of the problem.
   *i rS is a list of square roots of individual penalty matrices, packed
     in one array. The ith such matrix rSi, say, has dimension q by rSncol[i]
     and the ith penalty is [rSi][rSi]'.
   *i U1 is an (orthogonal) basis for the penalty range space (q by (q-Mp), where Mp
     is the null space dimension).
   *i sp is an M array of smoothing parameters (NOT log smoothing parameters)
   *i theta is the n_theta vector of extra parameters of the likelihood.
   *i z and w are pseudodata and iterative newton weights.
   *i wf are Fisher weights
   *i Dth, Dth2, Det, Det2, Det_th, Det2_th, Det3, Det_th2, Det4, Det3_th, Det2_th2
     give derivs of deviance wrt eta (linear predictor) and theta (extra params)
     in obvious notation, where e.g Detj_thk is derivative of deviance j times wrt
     eta and k times wrt theta. absence of a j or k implies they are 1. 
   *o beta - coefficients.
   *o b1 - first deriv of coefs w.r.t. sps (incl. theta)
   *o w1 - first deriv of weights w.r.t. sps (incl. theta)
   *o D1, D2 - first and second deriv of deviance w.r.t. sps (incl. theta)
   *o P0, P1, P2 - penalty b'Sb and its first and second derivs wrt sps. (incl theta)
   *o ldet, ldet1, ldet2, log|X'WX + S| & derivs wrt sp (incl theta)
   *o rV sqrt covariance matrix of coefs.
   *i rank_tol tol to use for rank estimation
   *o rank_est estiamted rank
   *i n , q, M, n_theta number of data, coefs, smoothing params and theta params.
   *i Mp penalty null space dimension
   *i Enrow rows of E.
   *i *rSncol array of number of cols in components of rS
   *i deriv order of deriv required (0, 1 or 2)
   *i fixed_penalty, non-zero indicates that S includes a fixed penalty component,
      the range space projected square root of which is in the final element of `UrS'.
      This information is used by get_detS2().
   *i nt number of threads to use, if supported.
   *i type 0 for computation using |w|^{-1} scaling, 1 to avoif this.    
   *o dVkk is M by M matrix containing curvature terms for objective w.r.t smoothing 
      params (is zero when second deriv is zero, but otherwise is not second deriv).

   The method has 4 main parts:

   1. The initial QR- decomposition and negative w correction SVD are performed, 
      and various quantities which are independent of derivatives are created

   2. IFT used to obtain derivatives of the coefficients wrt the log smoothing 
      parameters. 

   3. Evaluation of the derivatives of the deviance wrt the log smoothing parameters
      (i.e. calculation of D1 and D2)

   The method involves first and second derivatives of a number of k-vectors wrt
   log smoothing parameters (\rho), where k is q or n. Consider such a vector, v. 
   
   * v1 will contain dv/d\rho_0, dv/d\rho_1 etc. So, for example, dv_i/d\rho_j
     (indices starting at zero) is located in v1[q*j+i].
   
   * v2 will contain d^2v/d\rho_0d\rho_0, d^2v/d\rho_1d\rho_0,... but rows will not be
     stored if they duplicate an existing row (e.g. d^2v/d\rho_0d\rho_1 would not be 
     stored as it already exists and can be accessed by interchanging the sp indices).
     So to get d^2v_k/d\rho_id\rho_j: 
     i)   if i<j interchange the indices
     ii)  off = (j*m-(j+1)*j/2+i)*q 
     iii) v2[off+k] is the required derivative.       

*/
{ double *work,*p0,*p1,*p2,*p3,*p4,*p5,*p6,*p7,*K=NULL,
    *Vt,*b2=NULL,*P,xx=0.0,*eta1=NULL,*eta2=NULL,
    *PKtz,*wi=NULL,*w2=NULL,*Tk=NULL,*Tkm=NULL,
    *dev_hess=NULL,*R,*Rh,*raw,*Q1,*Q,*nulli,*WX,*tau,*R1;
  int i,j,k,*pivot1,ScS,*pi,rank,*pivot,
    ntot,n_2dCols=0,n_drop,*drop,tp,
    n_work,deriv2,neg_w=0,*nind,nr,TRUE=1,FALSE=0,ML=0; 
  
  #ifdef OPENMP_ON
  int m;
  m = omp_get_num_procs(); /* detected number of processors */
  if (*nt > m || *nt < 1) *nt = m; /* no point in more threads than m */
  omp_set_num_threads(*nt); /* set number of threads to use */
  #else
  *nt = 1;
  #endif

  if (*ldet<0) ML=1; /* require ML not REML */ 

  if (*deriv==2) deriv2=1; else deriv2=0;

  ScS=0;for (pi=rSncol;pi<rSncol + *M;pi++) ScS+= *pi;  /* total columns of input rS */

  /*d_tol = sqrt(*rank_tol * 100);*/
  /* first step is to obtain P and K */

  PKtz = (double *)CALLOC((size_t)*q,sizeof(double)); /* PK'z --- the pivoted coefficients*/
  nulli = (double *)CALLOC((size_t)*q,sizeof(double)); /* keep track of the params in null space */   
  drop = (int *)CALLOC((size_t)*q,sizeof(int)); /* original locations of dropped parameters */
  raw = (double *)CALLOC((size_t) *n,sizeof(double)); /* storage for sqrt(|w|) */
  n_work = (4 * *n + 2 * *q) * *M + 2 * *n;
  k = 5 * *q; if (n_work < k) n_work = k;
  k = (*M * (1 + *M))/2 * *n;
  if (n_work < k) n_work = k;
  work = (double *)CALLOC((size_t) n_work,sizeof(double)); /* work space for several routines*/
  nr = *q + *Enrow;
  R = (double *)CALLOC((size_t)*q * nr,sizeof(double));
  Rh = (double *)CALLOC((size_t)*q * *q,sizeof(double)); /* to hold Rh s.t. Rh'Rh = X'WX + S */
  pivot1=(int *)CALLOC((size_t)*q,sizeof(int));
  if (deriv2) dev_hess = (double *)CALLOC((size_t) *q * *q,sizeof(double)); else dev_hess=NULL;
  K = (double *)CALLOC((size_t) *n * *q,sizeof(double));
  P = (double *)CALLOC((size_t) *q * *q,sizeof(double));
  Q1 = (double *)CALLOC((size_t) *n * *q,sizeof(double)); 

  for (i=0;i< *n;i++) 
    if (w[i]<0) { neg_w++;raw[i] = sqrt(-w[i]);} 
    else if (w[i]>0) raw[i] = sqrt(w[i]);
    else { *type=1; } /* zero weights so we have to use type 1 method */

  if (neg_w) {  
    Vt = (double *)CALLOC((size_t) *q * *q,sizeof(double));
    nind = (int *)CALLOC((size_t)neg_w,sizeof(int)); /* index the negative w_i */
    k=0;for (i=0;i< *n;i++) if (w[i]<0) { nind[k]=i;k++;}
  } else { nind = (int *)NULL; Vt = (double *)NULL;}
  
  /* get  R,nulli,dev_hess,P,K,Vt,PKtz (== beta),Q1, nind,pivot1,drop,rank,n_drop,ldetXWXS */
  if (*type==1) z=wz; /* need to pass wz to gdiPK */ 

  gdiPK(work,X,E,Es,rS,U1,z,raw,
        R,Rh,nulli,dev_hess,P,K,Vt,PKtz,Q1,
        nind,pivot1,drop,
        n,q,Mp,neg_w,nt,Enrow,
        &rank,&n_drop,
        deriv2,ScS,&TRUE,
      rank_tol,ldet,type);
       
  FREE(raw);

  /* now call ift2 to get derivatives of coefs w.r.t. smoothing/theta parameters */
  ntot = *M + *n_theta;
  n_2dCols = (ntot  * (1 + ntot))/2;
  if (*deriv) {
    //b1 = (double *)CALLOC((size_t) rank * ntot,sizeof(double)); 
    eta1 = (double *)CALLOC((size_t) *n * ntot,sizeof(double)); 
    if (deriv2) {
      b2 = (double *)CALLOC((size_t) rank * n_2dCols,sizeof(double)); 
      eta2 = (double *)CALLOC((size_t) *n * n_2dCols,sizeof(double)); 
    }
    ift2(R,Vt,X,rS,PKtz,sp,theta,
          Det_th,Det2_th,Det3,Det_th2,
          b1,b2,eta1,eta2,
          n,&rank,M,n_theta,rSncol,&deriv2,&neg_w,&nr);
    
    if (*M>0) {
      i=0;mgcv_mmult(work,Rh,b1,&i,&i,&rank,M,&rank); /* Rh db/drho */
      /* Now obtain dVkk = db'/drho Rh' Rh db/drho ... */
      getXtX(dVkk,work,&rank,M);  
    }
    /* compute the grad of the deviance... */
    for (p4 = Dth,p0=D1,p1=eta1,i=0;i < *n_theta;i++,p0++) {
      for (*p0=0.0,p2 = Det,p3=Det + *n;p2<p3;p2++,p1++,p4++) *p0 += *p1 * *p2 + *p4;
    }
    for (;i<ntot;i++,p0++) {
      for (*p0=0.0,p2 = Det,p3=Det + *n;p2<p3;p2++,p1++) *p0 += *p1 * *p2;
    }

    if (deriv2) { /* get the Hessian of deviance as well */
       p0 = eta2;p6 = Dth2;
       for (i=0;i<ntot;i++) for (k=i;k<ntot;k++) {
           p1 = eta1 + i * *n; p2 = eta1 + k * *n;
           p3 = Det2; p4 = Det; p5 = Det + *n;
           for (xx=0.0;p4<p5;p0++,p1++,p2++,p3++,p4++) xx +=  *p3 * *p1 * *p2 + *p4 * *p0;
           if (k < *n_theta) {
	      for (j=0;j<*n;j++,p6++) xx += *p6; 
           }
           if (i < *n_theta) {
             p1 = Det_th + *n * i;p2 = p1 + *n;
             p3 = eta1 + *n *k;
             for (;p1<p2;p1++,p3++) xx += *p1 * *p3;
           }
           if (k < *n_theta) {
             p1 = Det_th + *n * k;p2 = p1 + *n;
             p3 = eta1 + *n *i;
             for (;p1<p2;p1++,p3++) xx += *p1 * *p3;
           }
           D2[i + ntot * k] = D2[k + ntot * i] = xx;
       }
    }
  } /* if (*deriv) */

  /* Obtain the derivatives of the weights w.r.t. theta, sp */
  
  if (*deriv) {
    /* first derivs... */
    p0 = w1; // = (double *)CALLOC((size_t) *n * ntot,sizeof(double)); 
    p3 = Det2_th;p4 = eta1;
    for (i=0;i<ntot;i++) {
      p1=Det3;p2 = p1 + *n;
      if (i < *n_theta) { 
	for (;p1<p2;p0++,p1++,p3++,p4++) *p0 = (*p1 * *p4 + *p3 ) * .5; 
      } else {
        for (;p1<p2;p0++,p1++,p4++) *p0 = *p1 * *p4 * .5;
      }
    }
    if (deriv2) { /* second derivs... */ 
      p0 = w2 = (double *)CALLOC((size_t)  *n * n_2dCols,sizeof(double)); 
      p1 = Det2_th2;p2 = eta2;
      for (i=0;i<ntot;i++) for (k=i;k<ntot;k++) {
      	p3 = Det3;p4 = Det4; 
        p5 = eta1 + i * *n;p6 = eta1 + k * *n;p7 = p3 + *n;
        for (;p3<p7;p0++,p2++,p3++,p4++,p5++,p6++) 
	   *p0 = *p4 * *p5 * *p6 + *p3 * *p2;
        if (i < *n_theta) {
           p0 -= *n;p3 = Det3_th + i * *n;
           p4 = eta1 + k * *n;p5 = p0 + *n;
           for (;p0 < p5;p0++,p3++,p4++) *p0 += *p3 * *p4;
        }
        if (k < *n_theta) {
           p0 -= *n;p3 = Det3_th + k * *n;
           p4 = eta1 + i * *n;p5 = p0 + *n;
           for (;p0 < p5;p0++,p3++,p4++) *p0 += *p3 * *p4;
           if (i < *n_theta) {
             p0 -= *n;
             for (;p0 < p5;p0++,p1++) *p0 += *p1;
           }
	}
        p0 -= *n;p5 = p0 + *n;
        for (;p0 < p5;p0++) *p0 *= .5;
      }
    } /* end of 2nd derivs */
    
    if (*type==0) { /* form |w|^{-1}dw/drho etc... */
      /* a useful array for Tk and Tkm */
      wi=(double *)CALLOC((size_t)*n,sizeof(double)); 
      for (i=0;i< *n;i++) { wi[i]=1/fabs(w[i]);}

      /* get Tk and Tkm */
      Tk = (double *)CALLOC((size_t) *n * ntot,sizeof(double)); 
      rc_prod(Tk,wi,w1,&ntot,n); 
      if (deriv2) { 
        Tkm = (double *)CALLOC((size_t)  *n * n_2dCols,sizeof(double));
        rc_prod(Tkm,wi,w2,&n_2dCols,n);
      }
      FREE(wi);
    } else { /* Tk and Tkm simply contain unscaled derivatives */
      Tk = w1;
      Tkm = w2;
    }
  } /* end of w derivs */

  if (ML) { /* ML is required (rather than REML), and need to save some stuff */
    /* save just R as rank by rank instead of nr by rank */
    for (p0=R,p1=R,j=0;j<rank;j++,p1+=nr) {
      for (p2=p1,i=0;i<=j;i++,p0++,p2++) *p0 = *p2;
      for (i=j+1;i<rank;i++,p0++) *p0 = 0.0; 
    }
  }
  /* the derivatives of b'S'b w.r.t. all parameters [theta,sp] */
  
  get_bSb(P0,P1,P2,sp,E,rS,rSncol,Enrow,&rank,M,n_theta,PKtz,b1,b2,deriv);

  /* Now get the derivatives of log|X'WX+S| w.r.t. all parameters [theta,sp] */ 
 
  if (ML) { 
    *ldet = MLpenalty1(ldet1,ldet2,Tk,Tkm,nulli,X,R,Q1,nind,sp,rS,rSncol,
		       &rank,n,Mp,M,n_theta,&neg_w,rank_tol,deriv,nt,type);

  } else get_ddetXWXpS(ldet1,ldet2,P,K,sp,rS,rSncol,Tk,Tkm,n,&rank,&rank,M,n_theta,deriv,*nt); 
  
  if (*deriv) { /* unpivot and zero pad b1 */
    
    for (j = ntot - 1;j>=0;j--) { 
      p0 = b1 + rank * j; /* start of source column */
      for (i=0;i< rank;i++) beta[pivot1[i]] = p0[i];
      undrop_rows(beta,*q,1,drop,n_drop); /* zero rows inserted */
      p1 = b1 + *q * j; /* start of target column */
      for (p0=beta,p2=p0 + *q;p0<p2;p0++,p1++) *p1 = *p0;
    }
  }
  /* PKtz into beta... */
  for (i=0;i< rank;i++) beta[pivot1[i]] = PKtz[i];
  undrop_rows(beta,*q,1,drop,n_drop); /* zero rows inserted */

  /* get rV and K using E(W)... these things are needed in order to compute 
     effective degrees of freedom safely, and for posterior inference, */
    
  /* form sqrt(wf)X augmented with E */
  nr = *n + *Enrow;
  WX = (double *) CALLOC((size_t) ( (nr + *nt * rank) * rank),sizeof(double));
  for (p0=w,p1=w + *n,p2=wf;p0<p1;p0++,p2++) *p0 = sqrt(*p2);
  for (p3=X,p0 = WX,i=0;i<rank;i++) {
    for (p1=w,p2=w+*n;p1<p2;p1++,p0++,p3++) *p0 = *p3 * *p1;
    for (j=0;j<*Enrow;j++,E++,p0++) *p0 = *E;
  }
  /* QR decompose it and hence get new P and K */
  pivot = (int *)CALLOC((size_t)rank,sizeof(int));
  tau = (double *)CALLOC((size_t)rank*(*nt+1),sizeof(double));
  mgcv_pqr(WX,&nr,&rank,pivot,tau,nt);
  R1 = (double *)CALLOC((size_t)rank*rank,sizeof(double));
  getRpqr(R1,WX,&nr,&rank,&rank,nt);
  Rinv(P,R1,&rank,&rank,&rank);
  FREE(R1); 
  Q = (double *)CALLOC((size_t) nr * rank,sizeof(double)); 
  for (i=0;i< rank;i++) Q[i * rank + i] = 1.0;
  tp=0;mgcv_pqrqy(Q,WX,tau,&nr,&rank,&rank,&tp,nt);
  for (p1=Q,p0=K,j=0;j<rank;j++,p1 += *Enrow) for (i=0;i<*n;i++,p1++,p0++) *p0 = *p1;
  FREE(Q);FREE(WX);FREE(tau);
  pivoter(P,&rank,&rank,pivot,&FALSE,&TRUE); /* unpivoting the rows of P */
  FREE(pivot);
  for (p1=P,i=0;i < rank; i++) for (j=0;j<rank;j++,p1++) rV[pivot1[j] + i * rank] = *p1;
  undrop_rows(rV,*q,rank,drop,n_drop); /* zero rows inserted */
  p0 = rV + *q * rank;p1 = rV + *q * *q;
  for (p2=p0;p2<p1;p2++) *p2 = 0.0; /* padding any trailing columns of rV with zeroes */

  /* Now unpack K into X -- useful for forming F = PK'W^.5X, diag of which is edf vector... */
  for (p0=X,p1=K,p2=K + rank * *n;p1<p2;p0++,p1++) *p0 = *p1;
  /* fill trailing columns with zero */ 
  for (p0 = X + rank * *n,p1 = X + *q * *n;p0<p1;p0++) *p0 = 0.0;

  /* ... now X contains K and rV%*%t(rV)*scale is the covariance matrix */
  
  /* and tidy up */

  *rank_est = rank;

  if (*deriv) { 
    FREE(eta1);
    if (*type==0) FREE(Tk);
    if (deriv2) {
      FREE(b2);FREE(eta2);FREE(w2);
      if (*type==0) FREE(Tkm);
      FREE(dev_hess);
    }
  }
  if (neg_w) {
    FREE(Vt);FREE(nind);
  }
  FREE(PKtz);FREE(nulli);FREE(drop);
  FREE(work);FREE(R);FREE(Rh);FREE(pivot1);FREE(K);
  FREE(P);FREE(Q1);

} /* gdi2 */



void gdi1(double *X,double *E,double *Es,double *rS,double *U1,
	  double *sp,double *z,double *w,double *wf,double *alpha,double *mu,double *eta, double *y,
	 double *p_weights,double *g1,double *g2,double *g3,double *g4,double *V0,
	  double *V1,double *V2,double *V3,double *beta,double *b1,double *w1,
         double *D1,double *D2,double *P0, double *P1,double *P2,double *trA,
    double *trA1,double *trA2,double *rV,double *rank_tol,double *conv_tol, int *rank_est,
	 int *n,int *q, int *M,int *Mp,int *Enrow,int *rSncol,int *deriv,
	  int *REML,int *fisher,int *fixed_penalty,int *nt,double *dVkk)     
/* 
   Version of gdi, based on derivative ratios and Implicit Function Theorem 
   calculation of the derivatives of beta. Assumption is that Fisher is only used 
   with canonical link, when it is equivalent to Newton anyway.

   This version does identifiability truncation on the basis of "well scaled" 
   penalty square root, Es, and is assuming that a stability enhancing 
   reparameterization and stable E are being employed.

   This version deals properly with negative weights, which can occur with Newton based 
   PIRLS. In consequence w's in this routine are proportional to reciprocal variances,
   not reciprocal standard deviations.
  
   The function is to be called at convergence of a P-IRLS scheme so that 
   z, w, mu and functions of these can be treated as fixed, and only the 
   derivatives need to be updated.

   All matrices are packed into arrays in column order (i.e. col1, col2,...)
   as in R. 

   All names ending in 1,2 or 3 are derivatives of some sort, with the integer
   indicating the order of differentiation. 

   The arguments of this function point to the following:
   * X is and n by q model matrix. On output this will contain K.
   * E is a q by Enrow square root of the total penalty matrix, so E'E=S
   * Es is the square root of a "well scaled" version of the total penalty,
     suitable for numerical determination of the theoretical rank of the problem.
   * rS is a list of square roots of individual penalty matrices, packed
     in one array. The ith such matrix rSi, say, has dimension q by rSncol[i]
     and the ith penalty is [rSi][rSi]'.
   * U1 is an (orthogonal) basis for the penalty range space (q by (q-Mp), where Mp
     is the null space dimension).
   * sp is an M array of smoothing parameters (NOT log smoothing parameters)
   * z, w and wf are n-vectors of the pseudodata iterative newton weights and iterative 
     fisher weights (only if `fisher' is zero) 
   * p_weights is an n-vector of prior weights (as opposed to the iterative weights in w)
   * mu and y are n-vectors of the fitted values and data.
   * g1,g2,g3,g4 are the n-vectors of the link derivatives 
     Note that g''(mu) g'''(mu) and g''''(mu) are *divided by* g'(mu)
   * V0, V1, V2, V3 are n-vectors of the variance function and first three derivatives,
     Note that V'(mu), V''(mu) & V'''(mu) are divided by V(mu)
   * beta is the vector for the returned coef vector and
     b1 is for returning the q by M array of derivs of coefs w.r.t. log sps 
   * D1 and D2 are an M-vector and M by M matrix for returning the first 
     and second derivatives of the deviance wrt the log smoothing parameters.
     if *REML is non zero then the derivs will be of the penalized deviance,
     and b'Sb will be returned in conv_tol  
   * trA1 and trA2 are an M-vector and M by M matrix for returning the first 
     and second derivatives of tr(A) wrt the log smoothing parameters.
     If *REML is non zero then the derivatives of the REML penalty, K, are 
     returned instead (with the REML penalty returned in `rank_tol', hack, hack).
   * P0,P1,P2 are for returning the Pearson statistic and its derivatives, or 
     the Pearson scale estimate and derivatives if *REML is non - zero. 
   * rank_est is for returning the estimated rank of the problem.
   * the remaining arguments are the dimensions already refered to except for:
   * deriv, which controls which derivatives are produced:
       deriv==0 for no derivatives: only trA, rV and beta returned
       deriv==1 for first derivatives only
       deriv==2 for gradient and Hessian
     -- on exit contains the number of iteration steps required.   

    * If REML is +ve non-zero, then the REML penalty is returned in rank_tol, with its 
      derivatives in trA1, trA2: it is to be added to the *deviance* to get D_r.
    * If REML is -ve non-zero, then the ML penalty is returned in place of the REML one.
    * non-zero `fisher' indicates that Fisher scoring, rather than full Newton,
      is the basis for iteration. 
    * non-zero `fixed_penalty' indicates that S includes a fixed penalty component,
      the range space projected square root of which is in the final element of `UrS'.
      This information is used by get_detS2().
    * nthreads tells routine how many threads to use for parallel code sections.

   The method has 4 main parts:

   1. The initial QR- decomposition and SVD are performed, various quantities which 
      are independent of derivatives are created

   2. IFT used to obtain derivatives of the coefficients wrt the log smoothing 
      parameters. 

   3. Evaluation of the derivatives of the deviance wrt the log smoothing parameters
      (i.e. calculation of D1 and D2)

   4. Evaluation of the derivatives of tr(A) (i.e. trA1 and trA2)

   The method involves first and second derivatives of a number of k-vectors wrt
   log smoothing parameters (\rho), where k is q or n. Consider such a vector, v. 
   
   * v1 will contain dv/d\rho_0, dv/d\rho_1 etc. So, for example, dv_i/d\rho_j
     (indices starting at zero) is located in v1[q*j+i].
   
   * v2 will contain d^2v/d\rho_0d\rho_0, d^2v/d\rho_1d\rho_0,... but rows will not be
     stored if they duplicate an existing row (e.g. d^2v/d\rho_0d\rho_1 would not be 
     stored as it already exists and can be accessed by interchanging the sp indices).
     So to get d^2v_k/d\rho_id\rho_j: 
     i)   if i<j interchange the indices
     ii)  off = (j*m-(j+1)*j/2+i)*q 
     iii) v2[off+k] is the required derivative.       

*/
{ double *WX,*tau,*work,*p0,*p1,*p2,*p3,*K=NULL,
    *R1,*Vt,xx,*b2,*P,*Q,
    *af1=NULL,*af2=NULL,*a1,*a2,*eta1=NULL,*eta2=NULL,
    *PKtz,*v1,*v2,*wi,*w2,*pw2,*Tk,*Tkm,*Tfk=NULL,*Tfkm=NULL,
         *pb2, *dev_grad,*dev_hess=NULL,
    ldetXWXS=0.0,reml_penalty=0.0,bSb=0.0,*R,*Rh,
    *alpha1,*alpha2,*raw,*Q1,*nulli;
  int i,j,k,*pivot=NULL,*pivot1,ScS,*pi,rank,tp,bt,ct,iter=0,m,one=1,
    n_2dCols=0,n_b2,n_drop,*drop,nt1,
      n_eta1=0,n_eta2=0,n_work,deriv2,neg_w=0,*nind,nr,TRUE=1,FALSE=0; 
  
  #ifdef OPENMP_ON
  m = omp_get_num_procs(); /* detected number of processors */
  if (*nt > m || *nt < 1) *nt = m; /* no point in more threads than m */
  omp_set_num_threads(*nt); /* set number of threads to use */
  #else
  *nt = 1;
  #endif
  nt1 = *nt; /* allows threading to be switched off for QR for debugging*/ 

  if (*deriv==2) deriv2=1; else deriv2=0;

  ScS=0;for (pi=rSncol;pi<rSncol + *M;pi++) ScS+= *pi;  /* total columns of input rS */

  /*d_tol = sqrt(*rank_tol * 100);*/
  /* first step is to obtain P and K */

  PKtz = (double *)CALLOC((size_t)*q,sizeof(double)); /* PK'z --- the pivoted coefficients*/
  nulli = (double *)CALLOC((size_t)*q,sizeof(double)); /* keep track of the params in null space */   
  drop = (int *)CALLOC((size_t)*q,sizeof(int)); /* original locations of dropped parameters */
  raw = (double *)CALLOC((size_t) *n,sizeof(double)); /* storage for sqrt(|w|) */
  n_work = (4 * *n + 2 * *q) * *M + 2 * *n; 
  k = 5 * *q; if (n_work < k) n_work = k;
  k = (*M * (1 + *M))/2 * *n; if (n_work < k) n_work = k;
 
  work = (double *)CALLOC((size_t) n_work,sizeof(double)); /* work space for several routines*/
  nr = *q + *Enrow;
  R = (double *)CALLOC((size_t)*q * nr,sizeof(double));
  Rh = (double *)CALLOC((size_t)*q * *q,sizeof(double)); /* to hold Rh s.t. Rh'Rh = X'WX + S */
  pivot1=(int *)CALLOC((size_t)*q,sizeof(int));
  if (deriv2) dev_hess = (double *)CALLOC((size_t) *q * *q,sizeof(double)); else dev_hess=NULL;
  K = (double *)CALLOC((size_t) *n * *q,sizeof(double));
  P = (double *)CALLOC((size_t) *q * *q,sizeof(double));
  Q1 = (double *)CALLOC((size_t) *n * *q,sizeof(double)); 

  for (i=0;i< *n;i++) 
    if (w[i]<0) { neg_w++;raw[i] = sqrt(-w[i]);} 
    else raw[i] = sqrt(w[i]);

  if (neg_w) {  
    Vt = (double *)CALLOC((size_t) *q * *q,sizeof(double));
    nind = (int *)CALLOC((size_t)neg_w,sizeof(int)); /* index the negative w_i */
    k=0;for (i=0;i< *n;i++) if (w[i]<0) { nind[k]=i;k++;}
  } else { nind = (int *)NULL; Vt = (double *)NULL;}
  
  /* get  R,nulli,dev_hess,P,K,Vt,PKtz,Q1, nind,pivot1,drop,rank,n_drop,ldetXWXS */
  gdiPK(work,X,E,Es,rS,U1,z,raw,
        R,Rh,nulli,dev_hess,P,K,Vt,PKtz,Q1,
        nind,pivot1,drop,
        n,q,Mp,neg_w,nt,Enrow,
        &rank,&n_drop,
        deriv2,ScS,REML,
        rank_tol,&ldetXWXS,&FALSE);
       
  /************************************************************************************/
  /* free some memory */                    
  /************************************************************************************/
  FREE(raw);
  /************************************************************************************/
  /* The coefficient derivative setup starts here */
  /************************************************************************************/
  /* set up some storage first */
  if (*deriv) {
    n_2dCols = (*M * (1 + *M))/2;
    n_b2 = rank * n_2dCols;
    b2 = (double *)CALLOC((size_t)n_b2,sizeof(double)); /* 2nd derivs of beta */
   
    //n_b1 = rank * *M;
    //b1 = (double *)CALLOC((size_t)n_b1,sizeof(double)); /* 1st derivs of beta */
   
    n_eta1 = *n * *M;
    eta1 = (double *)CALLOC((size_t)n_eta1,sizeof(double));
    Tk = (double *)CALLOC((size_t)n_eta1,sizeof(double));
   
    //w1 = (double *)CALLOC((size_t)n_eta1,sizeof(double));

    n_eta2 = *n * n_2dCols;
    eta2 = (double *)CALLOC((size_t)n_eta2,sizeof(double));
    Tkm = (double *)CALLOC((size_t)n_eta2,sizeof(double));
  
    w2 = (double *)CALLOC((size_t)n_eta2,sizeof(double));

 
    v1 = work;v2=work + *n * *M; /* a couple of working vectors */ 
   
    /* Set up constants involved updates (little work => leave readable!)*/
  
    a1=(double *)CALLOC((size_t)*n,sizeof(double));  
    a2=(double *)CALLOC((size_t)*n,sizeof(double));
    alpha1=alpha2 =(double *)NULL;
    if (*fisher) { /* Fisher scoring updates */
   
      /* set up constants involved in w updates */
      /* dw/deta = - w[i]*(V'/V+2g''/g')/g' */
      for (i=0;i< *n;i++) a1[i] = -  w[i] *(V1[i] + 2*g2[i])/g1[i];
     
      
      /* d2w/deta2 .... */
      for (i=0;i< *n;i++) 
        a2[i] = a1[i]*(a1[i]/w[i]-g2[i]/g1[i]) - w[i]*(V2[i]-V1[i]*V1[i] + 2*g3[i]-2*g2[i]*g2[i])/(g1[i]*g1[i]) ;

    } else { /* full Newton updates */
      
      alpha1 = (double *) CALLOC((size_t)*n,sizeof(double));
      alpha2 = (double *) CALLOC((size_t)*n,sizeof(double));
      for (i=0;i< *n;i++) {
        xx = V2[i]-V1[i]*V1[i]+g3[i]-g2[i]*g2[i]; /* temp. storage */
        alpha1[i] = (-(V1[i]+g2[i]) + (y[i]-mu[i])*xx)/alpha[i];
        alpha2[i] = (-2*xx + (y[i]-mu[i])*(V3[i]-3*V1[i]*V2[i]+2*V1[i]*V1[i]*V1[i]+g4[i]-3*g3[i]*g2[i]+2*g2[i]*g2[i]*g2[i]))/alpha[i];
      }
    
      /* end of preliminaries, now setup the multipliers that go forward */

   
      /* dw/deta ... */
      for (i=0;i<*n;i++) a1[i] = w[i]*(alpha1[i]-V1[i]-2*g2[i])/g1[i];
      /* d2w/deta2... */
      for (i=0;i<*n;i++) a2[i] = a1[i]*(a1[i]/w[i]-g2[i]/g1[i]) - 
                                 w[i]*(alpha1[i]*alpha1[i] - alpha2[i] + V2[i]-V1[i]*V1[i] + 2*g3[i]-2*g2[i]*g2[i])/(g1[i]*g1[i]) ;

      if (! *REML) { /* then Fisher versions of a1 and a2 also needed */
        af1=(double *)CALLOC((size_t)*n,sizeof(double));  
        af2=(double *)CALLOC((size_t)*n,sizeof(double));
        /* dwf/deta = - w[i]*(V'/V+2g''/g')/g' */
        for (i=0;i< *n;i++) af1[i] = -  wf[i] *(V1[i] + 2*g2[i])/g1[i];
        /* d2wf/deta2 .... */
        for (i=0;i< *n;i++) 
        af2[i] = af1[i]*(af1[i]/wf[i]-g2[i]/g1[i]) - wf[i]*(V2[i]-V1[i]*V1[i] + 2*g3[i]-2*g2[i]*g2[i])/(g1[i]*g1[i]) ;
      } 


      FREE(alpha1);FREE(alpha2);
      
    } /* end of full Newton setup */

    /* get gradient vector and Hessian of deviance wrt coefficients */
    for (i=0;i< *n ;i++) v1[i] = -2*p_weights[i]*(y[i]-mu[i])/(V0[i]*g1[i]);
    dev_grad=(double *)CALLOC((size_t) rank,sizeof(double));
    bt=1;ct=0;mgcv_mmult(dev_grad,X,v1,&bt,&ct,&rank,&one,n);
    
    if (deriv2) { /* get hessian of deviance w.r.t. beta */
       for (p0=dev_hess,p1=p0 + rank * rank;p0<p1;p0++) *p0 *= 2.0;      
    } 
  
  } /* end of if (*deriv) */ 
  else { /* keep compilers happy */
    v1=v2=b1=(double *)NULL;
    a1=a2=dev_grad=w2=b2=(double *)NULL;
    Tk=Tkm=(double *)NULL;
  }
  /************************************************************************************/
  /* End of the coefficient derivative preparation  */
  /************************************************************************************/


  /************************************************************************************/
  /* Implicit Function Theorem code */
  /************************************************************************************/
  if (*deriv) {
    /* obtain derivatives of beta (b1,b2) and eta (eta1,eta2) using the IFT (a1 = dw/deta) */

    /* Note that PKtz used as pivoted version of beta, but not clear that PKtz really essential if IFT used */

    ift1(R,Vt,X,rS,PKtz,sp,w,a1,b1,b2,eta1,eta2,n,&rank,M,rSncol,&deriv2,&neg_w,&nr);
    if (*M>0) {
      i=0;mgcv_mmult(work,Rh,b1,&i,&i,&rank,M,&rank); /* Rh db/drho */
      /* Now obtain dVkk = diag(db'/drho Rh' Rh db/drho) ... */
      getXtX(dVkk,work,&rank,M); 
    }
    /* Now use IFT based derivatives to obtain derivatives of W and hence the T_* terms */

    /* get derivatives of w */  
    rc_prod(w1,a1,eta1,M,n); /* w1 = dw/d\rho_k done */
    if (deriv2) {
      rc_prod(w2,a1,eta2,&n_2dCols,n); 
      for (pw2=w2,m=0;m < *M;m++) for (k=m;k < *M;k++) {
        rc_prod(v1,eta1 + *n * m,eta1 + *n * k,&one,n);
        rc_prod(v2,a2,v1,&one,n);
        p1=v2 + *n;
        for (p0=v2;p0<p1;p0++,pw2++) *pw2 += *p0;           
      } /* w2 completed */
    }


    /* a useful array for Tk and Tkm */
    wi=(double *)CALLOC((size_t)*n,sizeof(double)); 
    for (i=0;i< *n;i++) { wi[i]=1/fabs(w[i]);}

    /* get Tk and Tkm */
      
    rc_prod(Tk,wi,w1,M,n); 
    if (deriv2) rc_prod(Tkm,wi,w2,&n_2dCols,n);
    
    if (! *REML && ! *fisher) { /* then Fisher based versions of Tk and Tkm needed */ 
      rc_prod(w1,af1,eta1,M,n); /* w1 = dwf/d\rho_k done */
        Tfk = (double *)CALLOC((size_t)n_eta1,sizeof(double));
        Tfkm = (double *)CALLOC((size_t)n_eta2,sizeof(double));
      if (deriv2) {
        rc_prod(w2,af1,eta2,&n_2dCols,n); 
        for (pw2=w2,m=0;m < *M;m++) for (k=m;k < *M;k++) {
          rc_prod(v1,eta1 + *n * m,eta1 + *n * k,&one,n);
          rc_prod(v2,af2,v1,&one,n);
          p1=v2 + *n;
          for (p0=v2;p0<p1;p0++,pw2++) *pw2 += *p0;           
        } /* w2 completed */
      }
      for (i=0;i< *n;i++) { wi[i]=1/wf[i];}
      rc_prod(Tfk,wi,w1,M,n); 
      if (deriv2) rc_prod(Tfkm,wi,w2,&n_2dCols,n);
      FREE(af1);FREE(af2);
    } /* Fisher based Tk, Tkm, completed */
    else {Tfk = Tfkm = NULL;}

    /* evaluate gradient and Hessian of deviance */

    bt=1;ct=0;mgcv_mmult(D1,b1,dev_grad,&bt,&ct,M,&one,&rank); /* gradient of deviance is complete */
      
    if (deriv2) {       
      getXtMX(D2,b1,dev_hess,&rank,M,v1);
          
      for (pb2=b2,m=0;m < *M;m++) for (k=m;k < *M;k++) { /* double sp loop */
          p1 = dev_grad + rank;  
          for (xx=0.0,p0=dev_grad;p0<p1;p0++,pb2++) xx += *p0 * *pb2;
          D2[k * *M + m] += xx;
          D2[m * *M + k] = D2[k * *M + m];
      } /* Hessian of Deviance is complete !! */
    }

  } else {wi=NULL;} /* end of if (*deriv) */

  /* END of IFT */

  if (*REML<0) { /* ML is required (rather than REML), and need to save some stuff */
    /* save just R as rank by rank instead of nr by rank */
    for (p0=R,p1=R,j=0;j<rank;j++,p1+=nr) {
      for (p2=p1,i=0;i<=j;i++,p0++,p2++) *p0 = *p2;
      for (i=j+1;i<rank;i++,p0++) *p0 = 0.0; 
    }
  } else { FREE(R);FREE(Q1);FREE(nind); } /* needed later for ML calculation */


  /* REML NOTE: \beta'S\beta stuff has to be done here on pivoted versions.
     Store bSb in bSb, bSb1 in trA1 and bSb2 in trA2.
  */
  if (*REML) {
    get_bSb(&bSb,trA1,trA2,sp,E,rS,rSncol,Enrow,&rank,M,&FALSE,PKtz,b1,b2,deriv);
    if (*deriv) for (p2=D2,p1=trA2,i = 0; i< *M;i++) { /* penalized deviance derivs needed */
        D1[i] += trA1[i];
        if (deriv2) for (j=0;j<*M;j++,p1++,p2++) *p2 += *p1;   
    } 
  }

  pearson2(P0,P1,P2,y,mu,V0,V1,V2,g1,g2,p_weights,eta1,eta2,*n,*M,*deriv,deriv2);
  
  if (*REML) { /* really want scale estimate and derivatives in P0-P2, so rescale */
    j = *n - *Mp;
    i = 0; /* set to 1 to use penalized pearson statistic as basis for scale estimate */
    if (i) {
      *P0 += bSb;   /* penalized pearson statistic */
       if (*deriv) for (p0=trA1,p1 = P1,p2 = P1 + *M;p1<p2;p1++,p0++) *p1 += *p0;
       if (*deriv>1) for (p0=trA2,p1 = P2,p2 = P2 + *M * *M;p1<p2;p1++,p0++) *p1 += *p0;
    }
    *P0 /= j;     /* REML type scale estimate */
    /* derivatives also need penalized part added, before rescaling to get variance estimate...  */
    if (*deriv) for (p1 = P1,p2 = P1 + *M;p1<p2;p1++) *p1 /= j;
    if (*deriv>1) for (p1 = P2,p2 = P2 + *M * *M;p1<p2;p1++) *p1 /= j;
  }


  
  if (*deriv) { /* unpivot and zero pad b1 */
    
    for (j = *M-1;j>=0;j--) {
      p0 = b1 + rank * j; /* start of source column */
      for (i=0;i< rank;i++) beta[pivot1[i]] = p0[i];
      undrop_rows(beta,*q,1,drop,n_drop); /* zero rows inserted */
      p1 = b1 + *q * j; /* start of target column */
      for (p0=beta,p2=p0 + *q;p0<p2;p0++,p1++) *p1 = *p0;
    }
  }
  /* PKtz into beta... */
  for (i=0;i< rank;i++) beta[pivot1[i]] = PKtz[i];
  undrop_rows(beta,*q,1,drop,n_drop); /* zero rows inserted */


 
  /* Now get the REML penalty */

  if (*REML>0) { /* It's REML */
    /* Now deal with log|X'WX+S| */   
    reml_penalty = ldetXWXS;
    get_ddetXWXpS(trA1,trA2,P,K,sp,rS,rSncol,Tk,Tkm,n,&rank,&rank,M,&FALSE,deriv,*nt); /* trA1/2 really contain det derivs */
  } /* So trA1 and trA2 actually contain the derivatives for reml_penalty */

  if (*REML<0) { /* it's ML, and more complicated */
    
    /* get derivs of ML log det in trA1 and trA2... */

    reml_penalty =  MLpenalty1(trA1,trA2,Tk,Tkm,nulli,X,R,Q1,nind,sp,rS,rSncol,
			       &rank,n,Mp,M,&FALSE,&neg_w,rank_tol,deriv,nt,&FALSE);
    
    FREE(R);FREE(Q1);FREE(nind);
  } /* note that rS scrambled from here on... */


  /* clean up memory, except what's needed to get tr(A) and derivatives 
  */ 

  if (neg_w) FREE(Vt);   
  FREE(work);FREE(PKtz);
 
  if (*deriv) {
    //FREE(b1);
    FREE(eta1);
    FREE(eta2);
    FREE(a1);FREE(a2);FREE(wi);FREE(dev_grad);
    //FREE(w1);
    FREE(w2);FREE(b2);

    if (deriv2) { FREE(dev_hess);}
  }
  
  /* Note: the following gets only trA if REML is being used,
           so as not to overwrite the derivatives actually needed,
           which also means that it doesn't matter if MLpenalty
           has messed up rS */

  
  if (*fisher) { /* then all quantites are the ones required for EDF calculations */ 
    wf = w;Tfk=Tk;Tfkm=Tkm;
  } else { /* Need expected value versions of everything for EDF calculation */
    /* form sqrt(wf)X augmented with E */
    nr = *n + *Enrow;
    /* st WX = (double *)CALLOC((size_t)nr * rank,sizeof(double)); */
    WX = (double *) CALLOC((size_t) ( (nr + *nt * rank) * rank),sizeof(double));
    for (p0=w,p1=w + *n,p2=wf;p0<p1;p0++,p2++) *p0 = sqrt(*p2);
    for (p3=X,p0 = WX,i=0;i<rank;i++) {
      for (p1=w,p2=w+*n;p1<p2;p1++,p0++,p3++) *p0 = *p3 * *p1;
      for (j=0;j<*Enrow;j++,E++,p0++) *p0 = *E;
    }
    /* QR decompose it and hence get new P and K */
    pivot = (int *)CALLOC((size_t)rank,sizeof(int));
    /* st tau = (double *)CALLOC((size_t)rank,sizeof(double)); */
    tau = (double *)CALLOC((size_t)rank*(*nt+1),sizeof(double));
    /* st mgcv_qr(WX,&nr,&rank,pivot,tau); */
    mgcv_pqr(WX,&nr,&rank,pivot,tau,&nt1);

    /* st Rinv(P,WX,&rank,&nr,&rank); */ /* P= R^{-1} */
    R1 = (double *)CALLOC((size_t)rank*rank,sizeof(double));
    getRpqr(R1,WX,&nr,&rank,&rank,&nt1);

    Rinv(P,R1,&rank,&rank,&rank);
    FREE(R1); 

    /* there's something about the way you taste that makes me want to clear my throat, 
       there's a method to your madness, that really gets my goat */
    Q = (double *)CALLOC((size_t) nr * rank,sizeof(double)); 
    /* st for (i=0;i< rank;i++) Q[i * nr + i] = 1.0; */
    /* st left=1;tp=0;mgcv_qrqy(Q,WX,tau,&nr,&rank,&rank,&left,&tp); */ /* Q from the second QR decomposition */
    for (i=0;i< rank;i++) Q[i * rank + i] = 1.0;
    tp=0;mgcv_pqrqy(Q,WX,tau,&nr,&rank,&rank,&tp,&nt1);

    for (p1=Q,p0=K,j=0;j<rank;j++,p1 += *Enrow) for (i=0;i<*n;i++,p1++,p0++) *p0 = *p1;
    FREE(Q);FREE(WX);FREE(tau);
    if (*deriv)  pivoter(rS,&rank,&ScS,pivot,&FALSE,&FALSE); /* apply the latest pivoting to rows of rS */
    
  }

 
  if (*REML) i=0; else i = *deriv;
  get_trA2(trA,trA1,trA2,P,K,sp,rS,rSncol,Tfk,Tfkm,wf,n,&rank,&rank,M,&i,nt);


  /* unpivot P into rV.... */
  /* note that PP' and hence rV rV' are propto the cov matrix. */

  if (!*fisher) { /* first unpivot rows of P */
    pivoter(P,&rank,&rank,pivot,&FALSE,&TRUE); /* unpivoting the rows of P */
    FREE(pivot);
  }

  for (p1=P,i=0;i < rank; i++) for (j=0;j<rank;j++,p1++) rV[pivot1[j] + i * rank] = *p1;
  undrop_rows(rV,*q,rank,drop,n_drop); /* zero rows inserted */
  p0 = rV + *q * rank;p1 = rV + *q * *q;
  for (p2=p0;p2<p1;p2++) *p2 = 0.0; /* padding any trailing columns of rV with zeroes */

  /* Now unpack K into X -- useful for forming F = PK'W^.5X, diag of which is edf vector... */
  for (p0=X,p1=K,p2=K + rank * *n;p1<p2;p0++,p1++) *p0 = *p1;
  /* fill trailing columns with zero */ 
  for (p0 = X + rank * *n,p1 = X + *q * *n;p0<p1;p0++) *p0 = 0.0;

  *rank_est = rank;

  FREE(drop);
  FREE(nulli);
  FREE(pivot1);
  FREE(P);FREE(K);
  FREE(Rh);
  if (*deriv) { 
    FREE(Tk);FREE(Tkm);
    if (! *REML && ! *fisher) { FREE(Tfk);FREE(Tfkm);}
  }

  if (*REML) {*rank_tol = reml_penalty;*conv_tol = bSb;}

  *deriv = iter; /* the number of iteration steps taken */
} /* end of gdi1() */



void R_cond(double *R,int *r,int *c,double *work,double *Rcondition)
/* estimates the condition number of c by c upper triangular matrix 
   R stored in an r by c matrix, stored in column order. 
   work should be of minimum length c * 4. 
   Algorithm is due to Cline, Moler, Stewart and Wilkinson (1979) as reported in
   Golub and Van Loan (1996).
   n <- 100;m <- 10
   X <- matrix(runif(n*m),n,m)
   X[,1] <- X[,m]
   qrx <- qr(X,LAPACK=TRUE)
   R <- qr.R(qrx)
   .C("R_cond",R=as.double(R),r=as.integer(m),c=as.integer(m),
       work=as.double(rep(0,4*m)),Rcond=as.double(1),PACKAGE="mgcv")$Rcond
*/ 
{ double kappa,*pm,*pp,*y,*p,ym,yp,pm_norm,pp_norm,y_inf=0.0,R_inf=0.0;
  int i,j,k;
  /* allocate work */
  pp=work;work+= *c;pm=work;work+= *c;
  y=work;work+= *c;p=work;work+= *c;   
  for (i=0;i<*c;i++) p[i] = 0.0;
  for (k=*c-1;k>=0;k--) {
      yp = (1-p[k])/R[k + *r *k];
      ym = (-1-p[k])/R[k + *r *k];
      for (pp_norm=0.0,i=0;i<k;i++) { pp[i] = p[i] + R[i + *r * k] * yp;pp_norm += fabs(pp[i]);}
      for (pm_norm=0.0,i=0;i<k;i++) { pm[i] = p[i] + R[i + *r * k] * ym;pm_norm += fabs(pm[i]);}
      if (fabs(yp)+pp_norm >= fabs(ym)+pm_norm) {
	  y[k]=yp;
          for (i=0;i<k;i++) p[i] = pp[i];
      } else {
          y[k]=ym;
          for (i=0;i<k;i++) p[i] = pm[i];
      }
      kappa=fabs(y[k]);
      if (kappa>y_inf) y_inf=kappa;
  }
  for (i=0;i<*c;i++) { 
    for (kappa=0.0,j=i;j<*c;j++) kappa += fabs(R[i + *r * j]);  
    if (kappa>R_inf) R_inf = kappa;
  }
  kappa=R_inf*y_inf;
  *Rcondition=kappa;
} /* end R_cond */


void pls_fit1(double *y,double *X,double *w,double *wy,double *E,double *Es,int *n,int *q,int *rE,double *eta,
              double *penalty,double *rank_tol,int *nt,int *use_wy)
/* Fast but stable PLS fitter. Obtains linear predictor, eta, of weighted penalized linear model,
   without evaluating the coefficients, but also returns coefficients in case they are needed.

   Uses QR approach, but tests that X'Wz = R'Q_1'sqrt(\bar w)\bar z (in Wood 2011 notation), to
   ensure that rhs is stable, and uses R^{-T}X'Wy in plce of Q_1'sqrt(\bar w)\bar z if not. 
   The reason for this is that it is possible for sqrt(w)*z to be *very* badly scaled when 
   w*z is well scaled.... 

   Also has the option to not test, but simply use X'Wy directly, if *use_wy is non-zero. This is 
   useful in situations in which y is pseudodata involving a reciprocal w and some w_i is zero.  
   
   Note that here E'E = S, while Es'Es = `well scaled version of S'   

   In this version the w_i are the w_i in \sum_i w_i (y_i - X_i \beta)^2
   rather than being the square root of these. Some w_i may be negative (as
   may occur when using Newton, rather than Fisher updates on IRLS). Note that it is still 
   assumed that any zero weighted data will have been dropped before the call.

   If nt>1 and openMP is available then routine computes with the optimal number of threads up 
   to nt.

   On return:
   
   * if *n is -ve then X'WX+E'E was not +ve definite (which means that the routine should be
     called again with weights based on Fisher scoring).

   otherwise:

   * eta contains the linear predictor
   * penalty is the evaluated penalty
   * the first q elements of y are the coefficients.     

*/

{ int i,j,k,rank,one=1,*pivot,*pivot1,left,tp,neg_w=0,*nind,bt,ct,nr,n_drop=0,*drop,TRUE=1,FALSE=0,nz;
  double *z,*WX,*tau,Rcond,xx,zz,zz1,*work,*Q,*Q1,*IQ,*raw,*d,*Vt,*p0,*p1,
    *R1,*tau1,Rnorm,Enorm,*R,*Xp;
  #ifdef OPENMP_ON
  int m;
  m = omp_get_num_procs(); /* detected number of processors */
  if (*nt > m || *nt < 1) *nt = m; /* no point in more threads than m */
  omp_set_num_threads(*nt); /* set number of threads to use */
  #else
  *nt = 1; /* no openMP support - turn off threading */
  #endif

  nr = *q + *rE;
  nz = *n; if (nz<nr) nz=nr; /* possible for nr to be more than n */
  z = (double *)CALLOC((size_t) nz,sizeof(double)); /* storage for z=[sqrt(|W|)z,0] */
  
  raw = (double *)CALLOC((size_t) *n,sizeof(double)); /* storage for sqrt(|w|) */
  
  for (i=0;i< *n;i++) 
    if (w[i]<0) { neg_w++;raw[i] = sqrt(-w[i]);} 
    else raw[i] = sqrt(w[i]);

  if (neg_w) {
    nind = (int *)CALLOC((size_t)neg_w,sizeof(int)); /* index the negative w_i */
    k=0;for (i=0;i< *n;i++) if (w[i]<0) { nind[k]=i;k++;}
  } else { nind = (int *)NULL;}

  for (i=0;i< *n;i++) { 
    z[i] = y[i]*raw[i]; /* form z itself*/
  }

  for (i=0;i<neg_w;i++) {k=nind[i];z[k] = -z[k];} 

  /* st WX = (double *) CALLOC((size_t) ( *n * *q),sizeof(double)); */
  WX = (double *) CALLOC((size_t) ( (*n + *nt * *q) * *q),sizeof(double));
  for (Xp=X,p0=WX,j=0;j<*q;j++) { 
    for (p1=raw,i=0;i<*n;i++,p1++,p0++,Xp++) { /* form WX */
      *p0 = *Xp * *p1;
      /* k = i + nn * j;
         WX[k]=raw[i]*X[i + *n *j];*/
    }
  } 
  /* get the QR decomposition of WX */
  /* st tau=(double *)CALLOC((size_t)*q,sizeof(double)); */ /* part of reflector storage */
  tau=(double *)CALLOC((size_t) *q * (*nt + 1),sizeof(double)); 
  
  pivot=(int *)CALLOC((size_t)*q,sizeof(int));
  
  /* st mgcv_qr(WX,n,q,pivot,tau); */ /* WX and tau now contain the QR decomposition information */
  mgcv_pqr(WX,n,q,pivot,tau,nt);

  /* pivot[i] gives the unpivoted position of the ith pivoted parameter.*/
  
  /* copy out upper triangular factor R, and unpivot it */
  R1 = (double *)CALLOC((size_t)*q * *q,sizeof(double));
  /* st for (i=0;i<*q;i++) for (j=i;j<*q;j++) R1[i + *q * j] = WX[i + *n * j]; */ 
  getRpqr(R1,WX,n,q,q,nt);
  
  pivoter(R1,q,q,pivot,&TRUE,&TRUE); /* unpivoting the columns of R1 */
 
  /* Form a nicely scaled version of [R',Es']' for rank determination */ 
  Rnorm = frobenius_norm(R1,q,q);
  Enorm =  frobenius_norm(Es,rE,q);
 
  R = (double *)CALLOC((size_t)*q * nr,sizeof(double));
  for (j=0;j<*q;j++) { 
    for (i=0;i< *q;i++) R[i + nr * j] = R1[i + *q * j]/Rnorm;
    for (i=0;i< *rE;i++) R[i + *q + nr * j] = Es[i + *rE * j]/Enorm;
  }
  
  /* ... and now use it to establish rank */
   
  tau1=(double *)CALLOC((size_t)*q,sizeof(double)); /* part of reflector storage */
  pivot1=(int *)CALLOC((size_t)*q,sizeof(int));
  mgcv_qr(R,&nr,q,pivot1,tau1);
  
  /* now actually find the rank of R */
  work = (double *)CALLOC((size_t)(4 * *q),sizeof(double));
  rank = *q;
  R_cond(R,&nr,&rank,work,&Rcond);
  while (*rank_tol * Rcond > 1) { rank--;R_cond(R,&nr,&rank,work,&Rcond);}
  
  /* Now have to drop the unidentifiable columns from R1, E and the corresponding rows from rS
     The columns to drop are indexed by the elements of pivot1 from pivot1[rank] onwards.
     Before returning, zeros will need to be inserted in the parameter vector at these locations. 
  */

  n_drop = *q - rank;
  if (n_drop) {
    drop = (int *)CALLOC((size_t)n_drop,sizeof(int)); /* original locations of dropped parameters */
    for (i=0;i<n_drop;i++) drop[i] = pivot1[rank+i];
    qsort(drop,n_drop,sizeof(int),icompare); /* key assumption of the drop/undrop routines is that `drop' is ascending */
    /* drop columns indexed in `drop'... */
    drop_cols(R1,*q,*q,drop,n_drop);    /* R1 now q by rank */
    drop_cols(E,*rE,*q,drop,n_drop); /* E now q by rank */ 
  } else {drop=NULL;}

  /* At this stage the parameter space has been purged of terms that are
     theoretically unidentifiable, given WX and the penalties */

  /* Now augment R1 with the real square root penalty (not the nicely scaled version), result in R... */  
  for (j=0;j<rank;j++) { 
    for (i=0;i< *q;i++) R[i + nr * j] = R1[i + *q * j];
      for (i=0;i< *rE;i++) R[i + *q + nr * j] = E[i + *rE * j];
  }
  FREE(R1);
  mgcv_qr(R,&nr,&rank,pivot1,tau1); /* The final QR decomposition */ 
  

  if (neg_w) { /* then the correction for the negative w_i has to be evaluated */
    Q = (double *)CALLOC((size_t) nr * rank,sizeof(double)); 
    for (i=0;i< rank;i++) Q[i * nr + i] = 1.0;
    left=1;tp=0;mgcv_qrqy(Q,R,tau1,&nr,&rank,&rank,&left,&tp); /* Q from the second QR decomposition */

    Q1 = (double *)CALLOC((size_t) *n * rank,sizeof(double)); 
    /* st for (i=0;i<*q;i++) for (j=0;j<rank;j++) Q1[i + *n * j] = Q[i + nr * j]; */
    /* st left=1;tp=0;mgcv_qrqy(Q1,WX,tau,n,&rank,q,&left,&tp); */ /* Q1 = Qb Q[1:q,]  where Qb from first QR decomposition */   
    for (i=0;i<*q;i++) for (j=0;j<rank;j++) Q1[i + *q * j] = Q[i + nr * j];
    tp=0;mgcv_pqrqy(Q1,WX,tau,n,q,&rank,&tp,nt);/* Q1 = Qb Q[1:q,]  where Qb from first QR decomposition */   
    
    FREE(Q);

    if (neg_w < rank+1) k = rank+1; else k = neg_w;
    IQ = (double *)CALLOC((size_t) k * rank,sizeof(double)); 
    for (i=0;i<neg_w;i++) { /* Copy the rows of Q1 corresponding to -ve w_i into IQ */
      p0 = IQ + i;p1 = Q1 + nind[i];
      for (j=0;j<rank;j++,p0+=k,p1+= *n) *p0 = *p1;
    }
    FREE(Q1); 
    /* Note that IQ may be zero padded, for convenience */
    Vt = (double *)CALLOC((size_t) rank * rank,sizeof(double));
    d = (double *)CALLOC((size_t)  rank,sizeof(double));
    mgcv_svd_full(IQ,Vt,d,&k,&rank); /* SVD of IQ */
    FREE(IQ);
    for (i=0;i<rank;i++) {
      d[i] = 1 - 2*d[i]*d[i];
      if (d[i]< - *rank_tol) { /* X'WX not +ve definite, clean up and abort */
        *n = -1; 
        FREE(Vt);FREE(d);FREE(pivot);FREE(tau);
        FREE(nind);FREE(raw);FREE(z);FREE(work);FREE(WX);
        FREE(tau1);FREE(pivot1);FREE(R);if (n_drop) FREE(drop);
        return;
      }
      if (d[i]<=0) d[i]=0.0; else d[i] = 1/d[i];
    }
    /* d now contains diagonal of diagonal matrix (I-2D^2)^{-1} (possibly pseudoinverse) */
  } else {Vt = d = (double *)NULL; }
  /* The -ve w_i correction is now complete */

  if (!*use_wy) { /* Now get the fitted values X \beta, *without* finding \beta */
    /* st left=1;tp=1;mgcv_qrqy(z,WX,tau,n,&one,q,&left,&tp); */ /* z = Q'z */
    tp=1;mgcv_pqrqy(z,WX,tau,n,q,&one,&tp,nt);
    for (i=rank;i<nz;i++) z[i]=0.0;
    left=1,tp=1; 
    mgcv_qrqy(z,R,tau1,&nr,&one,&rank,&left,&tp); /* z = Q1'Q'z, where Q1 is the first rank rows of second orth factor */
 
    for (i=rank;i < nz;i++) z[i]=0.0;

    for (i=0;i<rank;i++) y[i] = work[i+*q] = z[i];        /* y = Q'z, or corrected version, for finding beta */ 
    
    if (neg_w) { /* apply the correction factor for negative w_i terms */
       bt=0;ct=0;mgcv_mmult(work,Vt,z,&bt,&ct,&rank,&one,&rank); /* V' Q_1' z */
       for (i=0;i<rank;i++) work[i] *= d[i];
       bt=1;ct=0;mgcv_mmult(z,Vt,work,&bt,&ct,&rank,&one,&rank); /* V (I-2D^2)^{-1} V' Q1' Q' z */
       for (i=0;i<rank;i++) y[i] = z[i];
    }
  
    left=1,tp=0; mgcv_qrqy(z,R,tau1,&nr,&one,&rank,&left,&tp); /* z = Q1 Q1'Q'z */

    for (*penalty=0.0,i=rank;i<nr;i++) *penalty += z[i]*z[i]; /* the penalty term */

    for (i=rank;i < *n;i++) z[i]=0.0;

    /* st left=1;tp=0;mgcv_qrqy(z,WX,tau,n,&one,q,&left,&tp); */ /* z = Q Q1 Q1'Q z */
    tp=0;mgcv_pqrqy(z,WX,tau,n,q,&one,&tp,nt);

    for (i=0;i<*n;i++) eta[i] = z[i]/raw[i]; /* the linear predictor */
  } /* if (!*use_wy) */
  /* form Wz (not sqrt(|W|)z)... */ 
  //  for (i=0;i< *n;i++) { 
  //  z[i] = y[i]*w[i]; /* form z itself*/
  //}

  /* form X'Wz, drop any entries in drop and pivot... */
  bt=1;ct=0;mgcv_mmult(work,X,wy,&bt,&ct,q,&one,n);     
  drop_rows(work,*q,1,drop,n_drop);
  pivoter(work,&rank,&one,pivot1,&FALSE,&FALSE);  
  
  if (!*use_wy) { /* Test R'Q'wz = X'Wz - not equal implies stability loss */ 
    for (zz=zz1=0,i=0;i<rank;i++) {
      for (xx=0,j=0;j<=i;j++) xx += R[j + nr * i] * work[*q + j];
      //Rprintf("XtWz = %g   RtQtwz = %g\n",work[i],xx);
      xx -= work[i];
      zz1 += xx*xx; zz += work[i]*work[i];   
    }
    if (zz1 > *rank_tol * zz) {
      *use_wy = 1;
      //Rprintf("instability detected zz1= %g zz = %g\n",zz1,zz);
    }
  }

  if (*use_wy) { /* then  R'Q'wz unstable or this computation signalled on entry */
    for (k=0;k<rank;k++) { /* forward sub loop (transposing R)*/
      for (xx=0.0,j=0;j < k;j++) xx += R[j + nr * k]*z[j];
      z[k] = (work[k] - xx)/R[k + nr * k];
    }
    if (neg_w) { /* apply the correction factor for negative w_i terms */
      bt=0;ct=0;mgcv_mmult(work,Vt,z,&bt,&ct,&rank,&one,&rank); /* V' Q_1' z */
      for (i=0;i<rank;i++) work[i] *= d[i];
      bt=1;ct=0;mgcv_mmult(z,Vt,work,&bt,&ct,&rank,&one,&rank); /* V (I-2D^2)^{-1} V' Q1' Q' z */
    } 
    for (i=0;i<rank;i++) y[i] = z[i];
  } /* end of alternative computation */
 
  /* now find  \hat \beta = R^{-1}Q'z, which are needed if P-IRLS starts to diverge
     in order to be able to evaluate penalty on step reduction */ 

  /* now back substitute to find \hat \beta = R^{-1}y */  
  for (k=rank-1;k>=0;k--) {
      for (xx=0.0,j=k+1;j < rank;j++) xx += R[k + nr * j]*z[j];
      z[k] = (y[k] - xx)/R[k + nr * k];
  }
  /* unpivot result (in z) into y */
  for (i=0;i< rank;i++) y[pivot1[i]] = z[i];
  /* insert zeroes for unidentifiables */
  undrop_rows(y,*q,1,drop,n_drop); 
  
  if (*use_wy) { /* re-compute other results from beta, as originals appear unstable or this method requested */
   
    bt=0;ct=0;mgcv_mmult(eta,X,y,&bt,&ct,n,&one,q);
    bt=0;ct=0;mgcv_mmult(work,E,y,&bt,&ct,rE,&one,q);
    for (*penalty=0.0,i=0;i < *rE;i++) *penalty += work[i]*work[i]; /* the penalty term */
  }

  FREE(z);
  FREE(WX);FREE(tau);FREE(pivot);FREE(raw);
  FREE(R);FREE(pivot1);FREE(tau1);  FREE(work);
  if (n_drop) FREE(drop);
  if (neg_w) { FREE(nind);FREE(d);FREE(Vt);}
} /* end pls_fit1 */



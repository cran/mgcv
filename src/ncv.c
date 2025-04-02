/* Copyright (C) 2022/3 Simon N. Wood  simon.wood@r-project.org

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

#ifdef _OPENMP // needs to precede R.h (and mgcv.h)
#include <omp.h>
#endif

#include "mgcv.h"
#include <R.h>



void minres0(double *R, double *u,double *b, double *x, int *p,int *m) {
/* Brute force alternative to minres for testing purposes */ 
  double *A,xx,zz,*work,workq;
  int p2,j,one=1,*ipiv,lwork=-1;
  char ntrans = 'N',trans='T',uplo='U',diag='N',side='L';
  p2 = *p * *p;
  A = (double *)CALLOC((size_t) p2,sizeof(double));
  ipiv = (int *)CALLOC((size_t) *p,sizeof(int));
  for (j=0;j<p2;j++) A[j] = R[j];
  xx=1.0;
  F77_CALL(dtrmm)(&side,&uplo,&trans,&diag,p,p,&xx,R,p,A,p FCONE FCONE FCONE FCONE); /* A = R'R */
  zz = -1.0;
  F77_CALL(dsyrk)(&uplo,&ntrans,p,m,&zz,u,p,&xx,A,p FCONE FCONE); /* A = R'R - uu' */
  for (j=0;j<*p;j++) x[j] = b[j];
  F77_CALL(dsysv)(&uplo,p,&one,A,p,ipiv,x,p,&workq,&lwork,&j FCONE);
  lwork=floor(workq);if (lwork<workq) lwork++;
  work = (double *)CALLOC((size_t) lwork,sizeof(double));
  F77_CALL(dsysv)(&uplo,p,&one,A,p,ipiv,x,p,work,&lwork,&j FCONE);
  FREE(A);FREE(ipiv);FREE(work);
}

void woodbury(double *R, double *u,double *b, double *x, int *p,int *m,double *work,int *iwork) {
/* Solves (R'R - uu')x=b where u is p by m and (R'R-uu') need not be positive definite, by
   direct use of Woodbury identity 
   (R'R - uu')^{-1} = R^{-1} (I_p - R^{-T}U(U^TR^{-1}R^{-T}U - I_m)^{-1}U^TR^{-1})R^{-T} )
   uses dgemv, dgemm, dtrsv, dtrsm, dsysv

   iwork is length m and nwork = iwork[0] on entry
   work is length p*m + m*m + m + nwork

   nwork should really be set by a preliminary call to ssysv. It is supposed to be m*nb 
   where nb is optimal block size. Ideally a one off call would be made setting nwork = -1 
     F77_CALL(dsysv)(&uplo,m,&one,B,m,iwork,v,m,work,nwork,&i FCONE FCONE);
   then *work is the optimal nwork on exit, while B and v are not accessed.
*/
  char side = 'L',trans = 'T',ntrans='N',uplo='U',diag='N';
  double done=1.0,*RitU,dzero=0.0,*B,*p0,*v,xx;
  int i,pm,one=1,nwork;
  pm = *p * *m;
  RitU = work;
  work += pm;
  B = work; work += *m * *m;
  v = work; work += *m;
  nwork = *iwork;
  for (i=0;i<pm;i++) RitU[i] = u[i];
  F77_CALL(dtrsm)(&side,&uplo,&trans,&diag,p,m,&done,R,p,RitU,p FCONE FCONE FCONE FCONE); /* form R^{-T}U */
  F77_CALL(dgemm)(&trans,&ntrans,m,m,p,&done,RitU,p,RitU,p,&dzero,B,m FCONE FCONE); /* form U^TR^{-1}R^{-T}U */
  for (p0=B,i=0;i < *m;i++,p0 += *m+1) *p0 += -1.0; /* B = U^TR^{-1}R^{-T}U - I */
  for (i=0;i<*p;i++) x[i] = b[i];
  F77_CALL(dtrsv)(&uplo,&trans,&diag,p,R,p,x,&one FCONE FCONE FCONE); /* x = R^{-T} b */
  F77_CALL(dgemv)(&trans,p,m,&done,RitU,p,x,&one,&dzero,v,&one FCONE); /* v = U^TR^{-1})R^{-T} b */
  F77_CALL(dsysv)(&uplo,m,&one,B,m,iwork,v,m,work,&nwork,&i FCONE); /* v = B^{-1}U^TR^{-1})R^{-T} b */
  xx = -1.0;
  F77_CALL(dgemv)(&ntrans,p,m,&xx,RitU,p,v,&one,&done,x,&one FCONE); /* x = (I-R^{-T}UB^{-1}U^TR^{-1})R^{-T} b */
  F77_CALL(dtrsv)(&uplo,&ntrans,&diag,p,R,p,x,&one FCONE FCONE FCONE); /* x = R^{-1}(I-R^{-T}UB^{-1}U^TR^{-1})R^{-T} b */
  
} /* woodbury */ 

void minres(double *R, double *u,double *b, double *x, int *p,int *m,double *work) {
/* solves (R'R - uu')x=b where u is p by m and (R'R-uu') need not be positive definite. 
   Use R as pre-conditioner. R'x0* = b, we solve (I - u*u*') x* = R^{-T}b
   where x = R^{-1}x* and u* = R^{-T}u. R is upper triangular.   

   work is p * (m+7) + m

   Note that m is over-written with number of iterations on output.
*/
  int one=1,i,j;
  double *v,*z,xx,zz,*u1,*dum,beta1,beta2,eta,epsilon,sig0,sig1,sig2,gamma0,gamma1,gamma2,
    *v1,*v2,alpha,delta,rho1,rho2,rho3,maxb,*w,*w1,*w2,*pp;
  char ntrans = 'N',trans='T',uplo='U',diag='N',side='L';
  u1 = work; work += *p * *m;
  v = work;work += *p;
  v1 = work;work += *p;
  v2 = work;work += *p;
  w = work;work += *p;
  w1 = work;work += *p;
  w2 = work;work += *p;
  z = work;work += *p;
  dum = work; work += *m;
  for (maxb=0.0,i=0;i<*p;i++) {
    xx = x[i] = b[i];maxb += xx*xx; 
  }
  maxb = sqrt(maxb);
  F77_CALL(dtrsv)(&uplo,&trans,&diag,p,R,p,x,&one FCONE FCONE FCONE); /* Solve R'x0* = b */
  xx = 1.0;
  for (i=0;i < *p * *m;i++) u1[i] = u[i];
  F77_CALL(dtrsm)(&side,&uplo,&trans,&diag,p,m,&xx,R,p,u1,p FCONE FCONE FCONE FCONE); /* Solve R'u1 = u */
  /* x currently contains R^{-T}b, form v = x - (I-u1u1')x = u1u1'x */
  zz = 0.0;
  F77_CALL(dgemv)(&trans,p,m,&xx,u1,p,x,&one,&zz,dum,&one FCONE); /* dum = u1'x */
  F77_CALL(dgemv)(&ntrans,p,m,&xx,u1,p,dum,&one,&zz,v1,&one FCONE); /* v1 = u1u1'x */
  for (beta1=0.0,i=0;i<*p;i++) { xx=v1[i]; beta1 += xx*xx;}
  epsilon = eta = beta1 = sqrt(beta1); /* beta1 = \\v1\\ */
  gamma0 = gamma1 = 1.0;sig0 = sig1 = 0.0;
  for (i=0;i < *p ;i++) v[i] = w1[i] = w[i] = 0.0;
  for (j=0;j<200;j++) {
    for (i=0;i<*p;i++) { v1[i] /= beta1;z[i]=v1[i];}
    xx = 1.0; zz = 0.0;
    F77_CALL(dgemv)(&trans,p,m,&xx,u1,p,v1,&one,&zz,dum,&one FCONE); /* dum = u1'v1 */
    zz = 1.0; xx = -1.0;
    F77_CALL(dgemv)(&ntrans,p,m,&xx,u1,p,dum,&one,&zz,z,&one FCONE); /* z = (I-u1u1')v1 */
    for (alpha=0.0,i=0;i<*p;i++) alpha += v1[i]*z[i];  /* alpha = v1'z = v1'(I-u1u1')v1 */
    for (beta2=0.0,i=0;i<*p;i++) {
      xx = v2[i] = z[i] - alpha * v1[i] - beta1 * v[i];
      beta2 += xx*xx;
    }
    delta = gamma1 * alpha - gamma0 * sig1 * beta1;
    rho1 = sqrt(delta*delta+beta2);
    beta2 = sqrt(beta2); /* beta_2 = ||v2|| */
    rho2 = sig1 * alpha + gamma0 * gamma1 * beta1;
    rho3 = sig0*beta1;
    gamma2 = delta/rho1;
    sig2 = beta2/rho1;
    xx = gamma2*eta;
    for (i=0;i<*p;i++) {
      w2[i] = (v1[i] - rho3*w[i] - rho2*w1[i])/rho1;
      x[i] += xx * w2[i];  
    }
    epsilon *= fabs(sig2);
    if (epsilon<maxb*1e-10) break;
    eta *= -sig2;
 
    pp = v; v = v1; v1 = v2; v2 = pp;
    pp = w; w = w1; w1 = w2; w2 = pp;
    sig0 = sig1; sig1 = sig2; beta1 = beta2;
    gamma0 = gamma1; gamma1 = gamma2;
  }
  F77_CALL(dtrsv)(&uplo,&ntrans,&diag,p,R,p,x,&one FCONE FCONE FCONE); /* Solve R'x = x* */
  *m=j;
} /* minres */

int CG(double *A,double *Mi,double *b, double *x,int n,double tol,double *cgwork) {
/* Basic pre-conditioned conjugate gradient solver for Ax = b where A is n by n 
   and Mi is the pre-conditioner. tol is the convergence tolerance. On exit x is 
   the approximate solution, and the return value is the number of iterations used.
   cgwork is an n*5 vector. 
*/
  int i,one=1,k;
  double *r,*z,*p,*r1,*z1,c1,c2,bmax=0.0,rmax,alpha,rz,r1z1,pAp,*dum,beta;
  char ntrans = 'N';
  p = cgwork;
  r = p + n;r1 = r + n;z=r1 + n;z1 = z + n;
  /* dgemv(char *trans,int *m, int *n,double a,double *A,int *lda,double *x,int *dx,
         double *b, double *y,int *dy)
       trans='T' to transpose A, 'N' not to. A is m by n. Forms y = a*A'x + b*y, or
       y = a*Ax + b*y. lda is number of actual rows in A (to allow for sub-matrices)
       dx and dy are increments of x and y indices.  */
  for (i=0;i<n;i++) {
    c1 = b[i]; r[i] = c1; /* copy b to r for gemv call below */
    c1 = fabs(c1);
    if (c1>bmax) bmax = c1; /* find max abs b for convergence testing */
  }
  c1 = -1.0;c2=1.0;
  F77_CALL(dgemv)(&ntrans,&n,&n,&c1,A,&n,x,&one,&c2,r,&one FCONE); /* r = b - Ax */
  c1 = 0.0;
  F77_CALL(dgemv)(&ntrans,&n,&n,&c2,Mi,&n,r,&one,&c1,z,&one FCONE); /* z = Mi r */
  for (i=0;i<n;i++) p[i] = z[i];
  c1=1.0;c2=0.0;
  for (k=0;k<200;k++) {  
    F77_CALL(dgemv)(&ntrans,&n,&n,&c1,A,&n,p,&one,&c2,z1,&one FCONE); /* z1 = Ap */
    for (rz=0.0,pAp=0.0,i=0;i<n;i++) {rz += r[i]*z[i];pAp += p[i] * z1[i];} /*r'z & p'Ap*/
    if (pAp==0.0) {k=-k;break;} /* failed, A possibly not +ve def */
    alpha = rz/pAp;
    for (rmax=0.0,i=0;i<n;i++) {
      x[i] += alpha * p[i];
      r1[i] = r[i] - alpha * z1[i];
      if (fabs(r1[i])>rmax) rmax = fabs(r1[i]);
    }
    if (rmax < tol*bmax) break;
    F77_CALL(dgemv)(&ntrans,&n,&n,&c1,Mi,&n,r1,&one,&c2,z1,&one FCONE); /* z1 = Mi r1 */
    for (r1z1=0.0,i=0;i<n;i++) r1z1 += r1[i] * z1[i];
    if (rz==0.0) {k=-k;break;} /* failed, A possibly not +ve def */
    beta = r1z1/rz;
    for (i=0;i<n;i++) p[i] = z1[i] + beta * p[i];
    dum = z1; z1 = z; z = dum;
    dum = r1; r1 = r; r = dum;
  }
  return(k); /* number of steps taken */
} /* CG */


SEXP nei_cov(SEXP v,SEXP d, SEXP d1, SEXP M, SEXP K) {
/* Computes a direct estimate of the parameter covariance matrix V, given neighbourhood structure 
   encoded in m and k, and parameters under leave one out perturbation in D, and D1. Otherwise 
   computing correction terms in similar way, but with D being the lett term and D1 the right.   
*/
  int *m,*k,n,p,i,j,i0,i1=0,ii,q;
  double *D,*D1,*V,*Ds;
  M = PROTECT(coerceVector(M,INTSXP));
  K = PROTECT(coerceVector(K,INTSXP)); /* otherwise R might be storing as double on entry */
  m = INTEGER(M); k = INTEGER(K);
  V = REAL(v);
  D = REAL(d);D1 = REAL(d1);
  n = length(M);
  p = ncols(d);
  for (ii=0;ii<p*p;ii++) V[ii]=0.0;
  Ds = (double *)CALLOC((size_t) p,sizeof(double));
  for (i=0;i<n;i++) { /* neibourhood loop */
    i0 = i1;i1 = m[i]; /* k[i0:(i1-1)] are neighbours of i */
    j=k[i0]; 
    for (q=0;q<p;q++) Ds[q] = D1[j+q*n]; /* first neighbour term to accumulate */  
    for (ii=i0+1;ii<i1;ii++) { /* remainder */ 
      j = k[ii];
      for (q=0;q<p;q++) Ds[q] += D1[j+q*n];
    }
    /* now add contribution to V */
    for (ii=0,j=0;j<p;j++) for (q=0;q<p;q++,ii++) V[ii] += D[i+j*n]*Ds[q]; 
  }
  FREE(Ds);
  UNPROTECT(2);
  return(R_NilValue); 
} /* nei_cov */  

SEXP ncv(SEXP x, SEXP hi, SEXP W1, SEXP W2, SEXP DB, SEXP DW, SEXP rS, SEXP IND, SEXP MI, SEXP M, SEXP K,SEXP BETA, SEXP SP, SEXP ETA, SEXP DETA,SEXP DLET,SEXP DERIV) {
/* Neighbourhood cross validation function. CG version - not optimal.
   Return: eta - eta[i] is linear predictor of y[ind[i]] when y[ind[i]] and its neighbours are ommited from fit
           deta - deta[i,j] is derivative of eta[ind[i]] w.r.t. log smoothing parameter j.
   Input: X - n by p model matrix. Hi inverse penalized Hessian. H penalized Hessian. w1 = w1[i] X[i,j] is dl_i/dbeta_j to within a scale parameter.
          w2 - -X'diag(w2)X is Hessian of log likelihood to within a scale parameter. db - db[i,j] is dbeta_i/d rho_j where rho_j is a log s.p. or 
          possibly other parameter. dw - dw[i,j] is dw2[i]/drho_j. rS[[i]] %*% t(rS[[i]]) is ith smoothing penalty matrix. 
          k[m[i-1]:(m[i]-1)] index the points in the ith neighbourhood. m[-1]=0 by convention. 
          Similarly ind[mi[i-1]:(mi[i]-1)] index the points whose linear predictors are to be predicted on dropping of the ith neighbourhood.
          beta - model coefficients (eta=X beta if nothing dropped). sp the smoothing parameters. deriv==0
          for no derivative calculations, deriv!=0 otherwise. 
         
   Basic idea: to approximate the linear predictor on omission of the neighbours of each point in turn, a single Newton step is taken from the full fit
               beta, using the gradient and Hessian implied by omitting the neighbours. To keep the cost at O(np^2) a pre-conditioned conjugate gradient 
               iteration is used to solve for the change in beta caused by the omission. 
               The gradient of this step w.r.t. to each smoothing parameter can also be obtained, again using CG to avoid O(p^3) cost for each obs.
               A point can be predicted several times with different omitted neighbourhoods. ind[i] is the point being predicted and eta[i] its prediction.
               LOOCV is recovered if ind = 0:(n-1) and each points neighbourhood is just itself.     
 */
  SEXP S,kr;
  int maxn,i,nsp,n,p,*m,*k,j,l,ii,i0,ki,q,p2,one=1,deriv,kk,error=0,jj,nm,*ind,nth,*mi,io,io0,no;
  double *X,*g,*g1,*gp,*p1,*Hp,*Hi,*Xi,xx,*xip,*xip0,z,*Hd,w1ki,w2ki,*wXi,*d,*w1,*w2,*eta,
    *deta,*beta,*dg,*dgp,*dwX,*wp,*wp1,*db=NULL,*dw=NULL,*rSj,*sp,*d1,*dbp,*dH=NULL,*xp,*wxp,*bp,*bp1,*dwXi,*cgwork,*dlet=NULL;
  char trans = 'T',ntrans = 'N';
  M = PROTECT(coerceVector(M,INTSXP));
  MI = PROTECT(coerceVector(MI,INTSXP));
  IND = PROTECT(coerceVector(IND,INTSXP));
  K = PROTECT(coerceVector(K,INTSXP)); /* otherwise R might be storing as double on entry */
  deriv = asInteger(DERIV);
  mi = INTEGER(MI);m = INTEGER(M); k = INTEGER(K);ind = INTEGER(IND);
  nsp = length(rS);
  nth = ncols(DETA)-nsp; /* how many non-sp parameters are there - first cols of db and dw relate to these */
  sp = REAL(SP);
  w1=REAL(W1);w2=REAL(W2);
  X = REAL(x);beta = REAL(BETA);
  Hi=REAL(hi);eta = REAL(ETA);deta = REAL(DETA);
  p = ncols(x); n = nrows(x);p2=p*p;
  no = length(IND); /* number of output lp values */
  nm = length(M); /* number of elements in cross validated eta - need not be n*/
  Hp = (double *)CALLOC((size_t) p2,sizeof(double));
  g = (double *)CALLOC((size_t) 3*p,sizeof(double));
  g1 = g + p;dg = g1 + p;
  d = (double *)CALLOC((size_t) 2*p,sizeof(double)); /* perturbation to beta on dropping y_i and its neighbours */
  cgwork = (double *)CALLOC((size_t) p*5,sizeof(double));
  d1 = d + p;
  /* need to know largest neighbourhood */
  maxn = ii = 0;
  for (j=0;j<nm;j++) {
    i = m[j]; if (i-ii>maxn) maxn = i-ii; ii = i;
  }  
  Xi = (double *)CALLOC((size_t) p*maxn,sizeof(double)); /* holds sub-matrix removed for this neighbourhood */
  wXi = (double *)CALLOC((size_t) p*maxn,sizeof(double)); /* equivalent pre-multiplied by diag(w2) */
  dwXi = (double *)CALLOC((size_t) p*maxn,sizeof(double)); /* equivalent pre-multiplied by d diag(w2)/d rho_j */
  Hd = (double *)CALLOC((size_t) p2,sizeof(double));
  dwX = (double *)CALLOC((size_t) p*n,sizeof(double));
  /* create Hessian X'diag(w2)X + S_lambda... */
  for (xip0 = X,xip=dwX,q=0;q<p;q++) for (wp=w2,wp1=wp+n;wp<wp1;wp++,xip++,xip0++) *xip = *xip0 * *wp;
  xx=1.0;z=0.0;
  F77_CALL(dgemm)(&trans,&ntrans,&p,&p,&n,&xx,X,&n,dwX,&n,&z,Hd,&p FCONE FCONE);
  for (j=0;j<nsp;j++) {
    S = VECTOR_ELT(rS, j);rSj = REAL(S);q = ncols(S);
    F77_CALL(dgemm)(&ntrans,&trans,&p,&p,&q,sp+j,rSj,&p,rSj,&p,&xx,Hd,&p FCONE FCONE);
  }  
  if (deriv>0) { /* derivarives of Hessian, dH/drho_j, needed */
    db=REAL(DB);dw=REAL(DW);dlet = REAL(DLET);
    dH = (double *)CALLOC((size_t) p2*(nsp+nth),sizeof(double));
    for (j=0;j<nsp+nth;j++) {
      for (xip0 = X,xip=dwX,q=0;q<p;q++) for (wp=dw+n*j,wp1=wp+n;wp<wp1;wp++,xip++,xip0++) *xip = *xip0 * *wp;    
      F77_CALL(dgemm)(&trans,&ntrans,&p,&p,&n,&xx,X,&n,dwX,&n,&z,dH+j*p2,&p FCONE FCONE); /* X'diag(dw[,j])X */
      if (j>=nth) { /* it's a smoothing parameter */
        S = VECTOR_ELT(rS, j-nth); /* Writing R Extensions 5.9.6 */
        rSj = REAL(S);q = ncols(S);
        F77_CALL(dgemm)(&ntrans,&trans,&p,&p,&q,sp+j-nth,rSj,&p,rSj,&p,&xx,dH+j*p2,&p FCONE FCONE); /* X'diag(dw[,j])X + lambda_j S_j */
      }
    } 
  } else if (deriv<0) dlet = REAL(DLET); /* storage for returning coeff changes per fold */
  FREE(dwX);
  for (io=ii=0,i=0;i<nm;i++) { /* loop over neighbourhoods, k[ii] is start of neighbourhood of i */
    p1 = g + p; /* fill accumulated g vector */
    ki = k[ii];w1ki = w1[ki];w2ki = w2[ki];
    i0=ii;io0=io; /* record of start needed in deriv calc */
    for (xip0=xip=Xi,gp=g,xp=X,wxp=wXi;gp<p1;gp++,xp += n,xip += maxn,wxp += maxn) {
      xx = xp[ki]; /* X[k[ii],j] */
      *gp = w1ki * xx; /* g gradient of log lik */
      *xip = xx; /* Xi matrix holding X[k[i],] */
      *wxp = xx*w2ki; /* wXi matrix holding w2[k[i]]*X[k[i],] */
    }
    q=1; /* count rows of Xi */
    for (xip0++,ii++;ii<m[i];ii++,xip0++,q++) { /* accumulate rest of g and Xi */ 
      ki = k[ii];w1ki = w1[ki];w2ki = w2[ki];
      for (xip=xip0,gp=g,xp=X,wxp=wXi+q;gp<p1;gp++,xp += n,xip += maxn,wxp += maxn) {
	xx = xp[ki];
	*gp += w1ki * xx;
	*xip = xx;
	*wxp = xx*w2ki; 
      }
    }  
    /* Now assemble the perturbed Hessian for this i */
    for (j=0;j<p2;j++) Hp[j] = Hd[j]; /* copy penalized Hessian */
    xx = -1.0;z=1.0; /* subtract part for neighbours of point i */ 
    /* dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *A,
         int *lda, double *B, int *ldb, double *beta,double *C,int *ldc) 
         transa/b = 'T' or 'N' for A/B transposed or not. C = alpha op(A) op(B) + beta C,
         where op() is transpose or not. C is m by n. k is cols of op(A). ldx is rows of X
         in calling routine (to allow use of sub-matrices) */
    F77_CALL(dgemm)(&trans,&ntrans,&p,&p,&q,&xx,Xi,&maxn,wXi,&maxn,&z,Hp,&p FCONE FCONE);
    /* dgemv(char *trans,int *m, int *n,double a,double *A,int *lda,double *x,int *dx,
         double *b, double *y,int *dy)
       trans='T' to transpose A, 'N' not to. A is m by n. Forms y = a*A'x + b*y, or
       y = a*Ax + b*y. lda is number of actual rows in A (to allow for sub-matrices)
       dx and dy are increments of x and y indices.  */
    xx = 0.0;
    F77_CALL(dgemv)(&ntrans,&p,&p,&z,Hi,&p,g,&one,&xx,d,&one FCONE); /* initial step Hi g */
    kk=CG(Hp,Hi,g,d,p,1e-13,cgwork); /* d is approx change in beta caused by dropping y_i and its neighbours */
    if (kk>error) error=kk;
    /* now create the linear predictors for target points */
    for (;io<mi[i];io++) {
      for (xx=0.0,xip=X+ind[io],j=0;j<p;j++,xip += n) xx += *xip * (beta[j]-d[j]);  
      eta[io] = xx; /* neighbourhood cross validated eta */
    }
    if (deriv<0) for (xip=dlet+p*i,j=0;j<p;j++) xip[j] = d[j]; /* dlet[,i] = d[] */
    /* now the derivatives */
    if (deriv>0) for (l=0;l<nsp+nth;l++) { /* loop over smoothing parameters */
	//	Rprintf(".");	
      /* compute sum_nei(i) dg/drho_l, start with first element of neighbourhood */
      for (xx=0.0,xip=Xi,bp=db+p*l,bp1=bp+p;bp<bp1;bp++,xip+=maxn) xx += *xip * *bp;
      for (dgp=dg,xip=wXi,p1=dg+p;dgp < p1;dgp++,xip+= maxn) *dgp = - *xip * xx;
      jj = i0;
      if (l<nth) for (xx=dlet[l*n+k[jj]],dgp=dg,xip=Xi,p1=dg+p;dgp < p1;dgp++,xip+= maxn) *dgp = - *xip * xx;
      for (jj++,j=1;j<q;j++,jj++) { /* loop over remaining neighbours */
        for (xx=0.0,xip=Xi+j,bp=db+p*l,bp1=bp+p;bp<bp1;bp++,xip+=maxn) xx += *xip * *bp;
        for (dgp=dg,xip=wXi+j,p1=dg+p;dgp < p1;dgp++,xip+= maxn) *dgp -= *xip * xx;
	if (l<nth) for (xx=dlet[l*n+k[jj]],dgp=dg,xip=Xi,p1=dg+p;dgp < p1;dgp++,xip+= maxn) *dgp = - *xip * xx;
      }
      /* Now subtract dH/drho_j d */
      /* First create diag(dw[,l])Xi */
      for (j=0;j<p;j++) for (xip0=Xi+j*maxn,xip=dwXi+j*maxn,wp=dw+l*n,jj=i0;jj<m[i];jj++,xip0++,xip++) *xip = *xip0 * wp[k[jj]]; 
      z=1.0;xx=0.0;
      F77_CALL(dgemv)(&ntrans,&q,&p,&z,dwXi,&maxn,d,&one,&xx,g,&one FCONE); /* g = diag(dw[,l])Xi d */
      F77_CALL(dgemv)(&trans,&q,&p,&z,Xi,&maxn,g,&one,&xx,g1,&one FCONE);  /* g1 = Xi'diag(dw[,l])Xi d */
      F77_CALL(dgemv)(&ntrans,&p,&p,&z,dH+l*p2,&p,d,&one,&xx,g,&one FCONE); /* g = dH_l d */
      for (j=0;j<p;j++) dg[j] += g1[j] - g[j]; /* sum_nei(i) dg/drho_l - dH/drho_l d */
      F77_CALL(dgemv)(&ntrans,&p,&p,&z,Hi,&p,dg,&one,&xx,d1,&one FCONE); /* initial step Hi dg */
      kk=CG(Hp,Hi,dg,d1,p,1e-13,cgwork); /* d1 is deriv wrt rho_l of approx change in beta caused by dropping the y_i and its neighbours */
      if (kk>error) error=kk;
      for (io=io0;io<mi[i];io++) {
        for (xx=0.0,xip=X+ind[io],dbp=db+p*l,j=0;j<p;j++,xip += n) xx += *xip * (dbp[j]-d1[j]);  
        deta[io+l*no] = xx;
      }	
    }
  }  
  FREE(Hp);FREE(Hd);
  FREE(g);FREE(d);FREE(cgwork);
  FREE(Xi);FREE(wXi);FREE(dwXi);
  if (deriv>0) FREE(dH);
  PROTECT(kr=allocVector(INTSXP,1));
  INTEGER(kr)[0] = error; /* max CG iterations used */
  UNPROTECT(5);
  return(kr);
} /* ncv */




SEXP Rncv(SEXP x, SEXP r, SEXP W1, SEXP W2, SEXP DB, SEXP DW, SEXP rS, SEXP IND, SEXP MI, SEXP M, SEXP K,SEXP BETA, SEXP SP, SEXP ETA,
	   SEXP DETA,SEXP DLET,SEXP DERIV,SEXP EPS,SEXP NT) {
/* Neighbourhood cross validation function, based on updating the Cholesky factor of the Hessian, rather than CG. 
   This is still O(np^2), but has the advantage of detecting any Hessian that is not positive definite. 
   
   OMP parallel version - scaling reasonable, as irreducibly level 2 dominated.

   Return: eta - eta[i] is linear predictor of y[ind[i]] when y[ind[i]] and its neighbours are omitted from fit
           deta - deta[i,j] is derivative of eta[ind[i]] w.r.t. log smoothing parameter j.
   Input: X - n by p model matrix. R chol factor of penalized Hessian. w1 = w1[i] X[i,j] is dl_i/dbeta_j to within a scale parameter.
          w2 - -X'diag(w2)X is Hessian of log likelihood to within a scale parameter. db - db[i,j] is dbeta_i/d rho_j where rho_j is a log s.p. or 
          possibly other parameter. dw - dw[i,j] is dw2[i]/drho_j. rS[[i]] %*% t(rS[[i]]) is ith smoothing penalty matrix. 
          k[m[i-1]:(m[i])] index the points in the ith neighbourhood. m[-1]=0 by convention. 
          Similarly ind[mi[i-1]:mi[i]] index the points whose linear predictors are to be predicted on dropping of the ith neighbourhood.
          beta - model coefficients (eta=X beta if nothing dropped). sp the smoothing parameters. deriv==0
          for no derivative calculations, deriv>0 to obtain first derivatives.
 
          deriv < 0 to compute NCV score without derivatives and return perturbations of beta in columns of DLET.  
         
   Basic idea: to approximate the linear predictor on omission of the neighbours of each point in turn, a single Newton step is taken from the full fit
               beta, using the gradient and Hessian implied by omitting the neighbours. To keep the cost at O(np^2) an O(p^2) update of the Cholesky
               factor is made to obtain Cholesky factor with dropped observations, enabling O(p^2) solution for the updated parameter.

               The gradient of this step w.r.t. to each smoothing parameter can also be obtained, again using the updated Cholesky factor to avoid 
               O(p^3) cost for each obs.

               If the updated Hessian is not positive definite the Cholesky update will detect this, and the routing falls back on Woodbury or minres.
               
               A point can be predicted several times with different omitted neighbourhoods. ind[i] is the point being predicted and eta[i] its prediction.
               LOOCV is recovered if ind = 0:(n-1) and each points neighbourhood is just itself.     
 */
  SEXP S,kr;
  int maxn,i,nsp,n,p,pg,*m,*k,j,l,ii,i0,ki,q,p2,one=1,deriv,*error,jj,nm,*ind,nth,*mi,io,io0,no,pdef,nddbuf,nwork = 0,
    nt,tid=0,pmaxn,*iwork=NULL,use_minres=1,niwork=0;
  double *X,*g,*g1,*gp,*p1,*R0,*R,*Xi,xx,*xip,*xip0,z,w1ki,w2ki,*wXi,*d,*w1,*w2,*eta,*p0,*p3,*ddbuf,*Rb,*work=NULL,
    *deta,*beta,*dg,*dgp,*dwX=NULL,*wp,*wp1,*db=NULL,*dw=NULL,*rSj,*sp,*d1,*dbp,*dH=NULL,*xp,*wxp,*bp,*bp1,*dwXi,*dlet=NULL,*dp,eps,alpha;
  char trans = 'T',ntrans = 'N',uplo='U',diag='N';
  M = PROTECT(coerceVector(M,INTSXP));
  MI = PROTECT(coerceVector(MI,INTSXP));
  IND = PROTECT(coerceVector(IND,INTSXP));
  K = PROTECT(coerceVector(K,INTSXP)); /* otherwise R might be storing as double on entry */
  deriv = asInteger(DERIV);
  nt = asInteger(NT);
  mi = INTEGER(MI);m = INTEGER(M); k = INTEGER(K);ind = INTEGER(IND);
  nsp = length(rS);
  nth = ncols(DETA)-nsp; /* how many non-sp parameters are there - first cols of db and dw relate to these */
  sp = REAL(SP);
  w1=REAL(W1);w2=REAL(W2);
  X = REAL(x);beta = REAL(BETA);
  R=REAL(r);eta = REAL(ETA);deta = REAL(DETA);
  eps = asReal(EPS);
  p = ncols(x); n = nrows(x);p2=p*p;
  no = length(IND); /* number of output lp values */
  nm = length(M); /* number of elements in cross validated eta - need not be n*/
  #ifndef _OPENMP
  nt = 1;
  #endif
 
  d = (double *)CALLOC((size_t) 2*p*nt,sizeof(double)); /* perturbation to beta on dropping y_i and its neighbours */
  d1 = d + p*nt;
  /* need to know largest neighbourhood */
  maxn = ii = 0;
  for (j=0;j<nm;j++) {
    i = m[j]; if (i-ii>maxn) maxn = i-ii; ii = i;
  }
  pmaxn = p*maxn;
  pg = p; if (pg<maxn) pg=maxn;
  g = (double *)CALLOC((size_t) (2*p+pg)*nt,sizeof(double));
  g1 = g + pg*nt;dg = g1 + p*nt;
  
  Xi = (double *)CALLOC((size_t) pmaxn*nt,sizeof(double)); /* holds sub-matrix removed for this neighbourhood */
  wXi = (double *)CALLOC((size_t) pmaxn*nt,sizeof(double)); /* equivalent pre-multiplied by diag(w2) */
  dwXi = (double *)CALLOC((size_t) pmaxn*nt,sizeof(double)); /* equivalent pre-multiplied by d diag(w2)/d rho_j */
  R0 = (double *)CALLOC((size_t) p2*nt,sizeof(double));Rb = (double *)CALLOC((size_t) p2*nt,sizeof(double));
  ddbuf = (double *)CALLOC((size_t) pmaxn*nt,sizeof(double)); /* buffer for downdates that spoil +ve def */ 
  if (use_minres) {
    nwork =  p*(maxn+7)+maxn;
  } else {
    niwork = maxn;
    iwork = (int *)CALLOC((size_t)niwork*nt,sizeof(int));
    /* workspace query (d not referenced) ... */
    nwork = -1;z = 1.0;i=1;
    F77_CALL(dsysv)(&uplo,&niwork,&i,d,&niwork,iwork,d,&niwork,&xx,&nwork,&j FCONE);
    nwork = (int) ceil(xx);
    nwork += maxn*(p+1+maxn);
  }  
  work = (double *)CALLOC((size_t) nwork*nt,sizeof(double));
  error = (int *)CALLOC((size_t) nt,sizeof(int));   
  xx=1.0;z=0.0;

  if (deriv>0) { /* derivarives of Hessian, dH/drho_j, needed */
    dwX = (double *)CALLOC((size_t) p*n,sizeof(double));
    db=REAL(DB);dw=REAL(DW);dlet = REAL(DLET);
    dH = (double *)CALLOC((size_t) p2*(nsp+nth),sizeof(double));
    for (j=0;j<nsp+nth;j++) {
      for (xip0 = X,xip=dwX,q=0;q<p;q++) for (wp=dw+n*j,wp1=wp+n;wp<wp1;wp++,xip++,xip0++) *xip = *xip0 * *wp;    
      F77_CALL(dgemm)(&trans,&ntrans,&p,&p,&n,&xx,X,&n,dwX,&n,&z,dH+j*p2,&p FCONE FCONE); /* X'diag(dw[,j])X */
      if (j>=nth) { /* it's a smoothing parameter */
        S = VECTOR_ELT(rS, j-nth); /* Writing R Extensions 5.9.6 */
        rSj = REAL(S);q = ncols(S);
        F77_CALL(dgemm)(&ntrans,&trans,&p,&p,&q,sp+j-nth,rSj,&p,rSj,&p,&xx,dH+j*p2,&p FCONE FCONE); /* X'diag(dw[,j])X + lambda_j S_j */
      }
    }
    FREE(dwX);
  } else if (deriv<0) dlet = REAL(DLET); /* for storing the coefficient changes for each fold */
  
 
  #ifdef _OPENMP
  /* schedule: static, dynamic or guided - seems hard to do better than guided */
#pragma omp parallel for schedule(guided) private(i,ii,io,i0,io0,j,jj,l,ki,q,alpha,nddbuf,pdef,xx,z,wxp,xip,xp,xip0,p0,p1,p3,dp,gp,dgp,bp,bp1,wp,dbp,w1ki,w2ki,tid) num_threads(nt)
  #endif
  for (i=0;i<nm;i++) { /* loop over neighbourhoods, k[ii] is start of neighbourhood of i */
    if (!i) io=ii=0;else {io=mi[i-1];ii=m[i-1];} 
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #endif
    p1 = g + tid*pg + p; /* fill accumulated g vector */
    ki = k[ii];w1ki = w1[ki];w2ki = w2[ki];
    i0=ii;io0=io; /* record of start needed in deriv calc */
    for (p0=R0+tid*p2,p3=R,j=0;j<p;j++,p0+=p,p3+=p) for (q=0;q<=j;q++) p0[q] = p3[q]; /* copy Cholesky factor*/
    alpha = sqrt(fabs(w2ki));
    nddbuf=0; /* counter for number of updates to store as they cause loss of definiteness */
    for (dp=d+tid*p,xip0=xip=Xi+pmaxn*tid,gp=g+tid*pg,xp=X,wxp=wXi+pmaxn*tid;gp<p1;dp++,gp++,xp += n,xip += maxn,wxp += maxn) { /* first element of neighbourhood */
      xx = xp[ki]; /* X[k[ii],j] */
      *gp = w1ki * xx; /* g gradient of log lik */
      *xip = xx; /* Xi matrix holding X[k[i],] */
      *wxp = xx*w2ki; /* wXi matrix holding w2[k[i]]*X[k[i],] */
      *dp = xx*alpha;
    }
    if (alpha<0) j=1; else j=0; /* update or downdate ? */
    chol_up(R0+p2*tid,d+tid*p,&p,&j,&eps);
    if (*(R0+p2*tid+1) < -0.5) { /* is update positive definite? */
      pdef=0;*(R0+p2*tid+1)=0.0; 
      for (p0=R0+p2*tid,p3=R,j=0;j<p;j++,p0+=p,p3+=p) for (q=0;q<=j;q++) p0[q] = p3[q]; /* restore factor to state before update attempt */
      for (p0=ddbuf+p*nddbuf+tid*pmaxn,p3=d+tid*p,j=0;j<p;j++) p0[j] = p3[j]; /* store the skipped update */
      nddbuf++;
    } else pdef=1;
    q=1; /* count rows of Xi */
    for (xip0++,ii++;ii<m[i];ii++,xip0++,q++) { /* accumulate rest of g and Xi for rest of neighbourhood*/ 
      ki = k[ii];w1ki = w1[ki];w2ki = w2[ki];alpha = sqrt(fabs(w2ki));
      for (dp=d+tid*p,xip=xip0,gp=g+pg*tid,xp=X,wxp=wXi+q+pmaxn*tid;gp<p1;dp++,gp++,xp += n,xip += maxn,wxp += maxn) {
	xx = xp[ki];
	*gp += w1ki * xx;
	*xip = xx;
	*wxp = xx*w2ki;
	*dp = xx*alpha;
      }
      if (alpha<0) {
	j=1; /* update */
      } else {/* downdate */
        for (p0=Rb+p2*tid,p3=R0+p2*tid,j=0;j<p;j++,p0+=p,p3+=p) for (l=0;l<=j;l++) p0[l] = p3[l]; /* backup state of R0 before attempting downdate */
        j=0;
      }
      chol_up(R0+p2*tid,d+tid*p,&p,&j,&eps);
      if (*(R0+p2*tid+1)< -0.5) { /* is update positive definite? */
	pdef=0;*(R0+p2*tid+1) = 0.0;
	for (p0=R0+tid*p2,p3=Rb+tid*p2,j=0;j<p;j++,p0+=p,p3+=p) for (l=0;l<=j;l++) p0[l] = p3[l]; /* restore factor to state before update attempt */
        for (p0=ddbuf+p*nddbuf+pmaxn*tid,p3=d+p*tid,j=0;j<p;j++) p0[j] = p3[j]; /* store the skipped update */
        nddbuf++;
      } 	
    }  
    
    if (pdef) { /* solve for R0'R0 d = g - the change in beta caused by dropping neighbourhood i */
      for (p0=d+tid*p,p3=g+tid*pg,j=0;j<p;j++) p0[j] = p3[j]; /* copy g to d */
      F77_CALL(dtrsv)(&uplo,&trans,&diag,&p,R0+tid*p2,&p,d+tid*p,&one FCONE FCONE FCONE);
      F77_CALL(dtrsv)(&uplo,&ntrans,&diag,&p,R0+tid*p2,&p,d+tid*p,&one FCONE FCONE FCONE);
    } else {  /* fallback solve (R0'R0 - uu') d= g via minres iteration, u, the skipped downdates are in ddbuf*/
     
      j = nddbuf; /* modified to number of iterations on exit */     
      if (use_minres) minres(R0+p2*tid,ddbuf+tid*pmaxn,g+pg*tid,d+p*tid,&p,&j,work+tid*nwork); else {
        *(iwork+niwork*tid) = nwork;
	woodbury(R0+p2*tid,ddbuf+tid*pmaxn,g+pg*tid,d+p*tid,&p,&j,work+tid*nwork,iwork+niwork*tid);
      }	
      error[tid]++; /* count the number of non +ve def cases */
    }
    for (;io<mi[i];io++) {
      for (p0=d+tid*p,xx=0.0,xip=X+ind[io],j=0;j<p;j++,xip += n) xx += *xip * (beta[j]-p0[j]);  
      eta[io] = xx; /* neighbourhood cross validated eta */
    }
    if (deriv<0) { /* return the coeff change */
      for (p0=d+tid*p,p1=dlet+p*i,j=0;j<p;j++) p1[j] = p0[j]; /* dlet[,i] = d[] */
    }  
    /* now the derivatives */
    if (deriv>0) for (l=0;l<nsp+nth;l++) { /* loop over smoothing parameters */
      /* compute sum_nei(i) dg/drho_l, start with first element of neighbourhood */
      for (xx=0.0,xip=Xi+tid*pmaxn,bp=db+p*l,bp1=bp+p;bp<bp1;bp++,xip+=maxn) xx += *xip * *bp;
      for (dgp=dg+tid*p,xip=wXi+tid*pmaxn,p1=dgp+p;dgp < p1;dgp++,xip+= maxn) *dgp = - *xip * xx;
      jj = i0;
      if (l<nth) for (xx=dlet[l*n+k[jj]],dgp=dg+p*tid,xip=Xi+tid*pmaxn,p1=dgp+p;dgp < p1;dgp++,xip+= maxn) *dgp += - *xip * xx;
      for (jj++,j=1;j<q;j++,jj++) { /* loop over remaining neighbours */
        for (xx=0.0,xip=Xi+j+pmaxn*tid,bp=db+p*l,bp1=bp+p;bp<bp1;bp++,xip+=maxn) xx += *xip * *bp;
        for (dgp=dg+tid*p,xip=wXi+j+tid*pmaxn,p1=dgp+p;dgp < p1;dgp++,xip+= maxn) *dgp -= *xip * xx;
	if (l<nth) for (xx=dlet[l*n+k[jj]],dgp=dg+tid*p,xip=Xi+j+tid*pmaxn,p1=dgp+p;dgp < p1;dgp++,xip+= maxn) *dgp += - *xip * xx;
      }
      /* Now subtract dH/drho_j d */
      /* First create diag(dw[,l])Xi */
      for (j=0;j<p;j++) for (xip0=Xi+j*maxn+pmaxn*tid,xip=dwXi+j*maxn+tid*pmaxn,wp=dw+l*n,jj=i0;jj<m[i];jj++,xip0++,xip++) *xip = *xip0 * wp[k[jj]]; 
      z=1.0;xx=0.0;
      F77_CALL(dgemv)(&ntrans,&q,&p,&z,dwXi+pmaxn*tid,&maxn,d+tid*p,&one,&xx,g+tid*pg,&one FCONE); /* g = diag(dw[,l])Xi d */
      F77_CALL(dgemv)(&trans,&q,&p,&z,Xi+pmaxn*tid,&maxn,g+tid*pg,&one,&xx,g1+tid*p,&one FCONE);  /* g1 = Xi'diag(dw[,l])Xi d */
      F77_CALL(dgemv)(&ntrans,&p,&p,&z,dH+l*p2,&p,d+tid*p,&one,&xx,g+tid*pg,&one FCONE); /* g = dH_l d */
      for (dgp = dg+tid*p,p0=g1+tid*p,p3=g+tid*pg,j=0;j<p;j++) dgp[j] += p0[j] - p3[j]; /* sum_nei(i) dg/drho_l - dH/drho_l d */
      if (pdef) { /* solve for R0'R0 d = g - the change in beta caused by dropping neighbourhood i */
	for (p0=d1+tid*p,p3=dg+tid*p,j=0;j<p;j++) p0[j] = p3[j];
        F77_CALL(dtrsv)(&uplo,&trans,&diag,&p,R0+tid*p2,&p,d1+tid*p,&one FCONE FCONE FCONE);
        F77_CALL(dtrsv)(&uplo,&ntrans,&diag,&p,R0+tid*p2,&p,d1+tid*p,&one FCONE FCONE FCONE);
      } else {  /* fallback solve (R0'R0 + uu') d= g where u are skipped downdates in ddbuf */
	j = nddbuf;
	if (use_minres) minres(R0+tid*p2,ddbuf+tid*pmaxn,dg+tid*p,d1+tid*p,&p,&j,work+nwork*tid); else {
          *(iwork+niwork*tid) = nwork;
	  woodbury(R0+tid*p2,ddbuf+tid*pmaxn,dg+tid*p,d1+tid*p,&p,&j,work+nwork*tid,iwork+niwork*tid);
	}  
      }
      for (io=io0;io<mi[i];io++) {
        for (p0=d1+tid*p,xx=0.0,xip=X+ind[io],dbp=db+p*l,j=0;j<p;j++,xip += n) xx += *xip * (dbp[j]-p0[j]);  
        deta[io+l*no] = xx;
      }	
    }
  }  
  FREE(R0);FREE(ddbuf);FREE(Rb);
  FREE(g);FREE(d);
  FREE(Xi);FREE(wXi);FREE(dwXi);
  if (deriv>0) FREE(dH);
  FREE(work);if (!use_minres) FREE(iwork);
  for (j=0,i=0;i<nt;i++) j += error[i];
  FREE(error);
  PROTECT(kr=allocVector(INTSXP,1));
  INTEGER(kr)[0] = j; /* max CG iterations used */
  UNPROTECT(5);
  return(kr);
} /* Rncv */



static inline int i2f(int i,int j,int K) { /* note *static* inline to avoid external image and potential link failure */
  /* see i3f - same idea for 2 indices */  
  int ii;
  if (i>j) { ii=j;j=i;i=ii;}
  ii = (i*(2*K-i+1))/2 + j-i;
  return(ii);
} /* i2f */

static inline int i3f(int i,int j,int k,int K) {
/* Suppose we fill an array... 
   for (m=i=0;i<K;i++) for (j=i;j<K;j++) for (k=j;k<K;k++,m++) a[m] = whatever
   ... the idea being that the same 'whatever' is stored for all permutations 
   (i,j,k), (j,k,i) etc. e.g. indices (0,2,3), (3,2,0), (2,0,3) etc all 
   result in the same stored quantity. We might then want to compute the 
   m for any such permutation of indices. This does it. 
   Note that all indices start at 0 here.
*/
  int ii;
  while (!(k>=j&&j>=i)) {
    if (i>j) {ii=j;j=i;i=ii;}
    if (j>k) {ii=j;j=k;k=ii;}
  }  
  ii = (i*(3*K*(K+1)+(i-1)*(i-3*K-2)))/6 + ((j-i)*(2*K+1-i-j))/2+k-j;
  return(ii);
} /* i3f */  

SEXP ncvls(SEXP x,SEXP JJ,SEXP h,SEXP hi,SEXP dH,SEXP L1, SEXP L2,SEXP L3,SEXP IND, SEXP MI, SEXP M, SEXP K,SEXP BETA,
	   SEXP ETACV,SEXP DETACV,SEXP DETA,SEXP DB,SEXP DERIV) {
/* This computes the NCV for GAMLSS families, using a Conjugate Gradient iteration. 
   Only intended as a fall back: Rncvls is better.

   Inputs:
   X[,jj[[i]]] is the model matrix for the ith linear predictor.
   H is the penalized Hessian, and Hi its inverse (or an approximation to it since its used as a pre-conditioner).
   dH[[i]] is the derivarive of H w.r.t. log(sp[i]);   
   lj contains jth derivatives of the likelihood w.r.t. the lp for each datum.
   k[m[i-1]:(m[i])] indexes the neighbours dropped for the ith neighbourhood.
   ind[mi[i-1]:m[i]] indexes the points predicted for the ith neighbourhood.
   m[-1]=mi[-1]=0 by convention. 
   beta - model coefficients. 
   deta and dbeta are matrices with derivaives of the lp's and coefs in their cols. deta has the lps stacked in each 
   column. 
   deriv=0 for no derivs, 1 to compute derivs.   
   
   Outputs:
   The perturbed etas will be returned in etacv: if nm is the length of ind,then eta[q*nm+i] is the ith element of 
   qth perturbed linear predictor.
   If deriv>0 then the the derivatives of the perturbed linear predictors are in deta: detacv[q*nm+i + l*(np*nlp)]] is the ith element of 
   deriv of qth lp w.r.t. lth log sp. If deriv<0 the detacv contains the perturbtations of the coefs for each neighbourhood in its rows.
 
 
   BUG? Offset handling!!
*/
  double *X,*H,*Hi,*l1,*l2,*l3=NULL,*beta,*g,*Hp,xx,z,*d,*d1,*cgwork,*eta,*deta,v,*db=NULL,*dbp,*detacv=NULL,*dh;
  int **jj,*jjl,*jjq,*ind,*m,*k,n,p,nm,nlp,*plp,ii,i,j,i0,i1,l,ln,ki,p2,q,r,l2i,one=1,kk,nsp,iter1=0,iter=0,deriv,*mi,io,io0,no;
  SEXP JJp,kr,DH;
  char ntrans = 'N';
  p = length(BETA);p2 = p*p;
  n = nrows(x);deriv = asInteger(DERIV);
  M = PROTECT(coerceVector(M,INTSXP));
  MI = PROTECT(coerceVector(MI,INTSXP));
  IND = PROTECT(coerceVector(IND,INTSXP));
  K = PROTECT(coerceVector(K,INTSXP)); /* otherwise R might be storing as double on entry */
  mi = INTEGER(MI); m = INTEGER(M); k = INTEGER(K);ind = INTEGER(IND);
  l1 = REAL(L1);l2=REAL(L2);
  nm = length(M);nlp = length(JJ);nsp = length(dH);no=length(IND);
  eta = REAL(ETACV);H=REAL(h);Hi = REAL(hi);
  beta = REAL(BETA);

  if (deriv>0) {
    l3=REAL(L3);deta=REAL(DETA);detacv=REAL(DETACV);db = REAL(DB);
  }  
  /* unpack the jj indices to here in order to avoid repeated list lookups withn loop */
  jj = (int **)CALLOC((size_t) nlp,sizeof(int *));
  plp = (int *)CALLOC((size_t) nlp,sizeof(int));
  for (l=0;l<nlp;l++) {
    JJp = VECTOR_ELT(JJ, l);JJp = PROTECT(coerceVector(JJp,INTSXP)); /* see R extensions 5.9.1 Handling the effects of garbage collection */
    plp[l] = length(JJp); jj[l] = INTEGER(JJp); /* jj[l][1:plp[l]] indexes cols of X for this lp */
    jjl=jj[l];for (i=0;i<plp[l];i++) jjl[i]--; /* 5.9.3 Details of R types. In fact coerceVector creates new vector only if type needs to change. */
  }  
  X = REAL(x);
  g = (double *)CALLOC((size_t) 3*p,sizeof(double)); /* gradient change */
  d = g+p; /* change in beta */
  d1 = d + p; /* deriv of above */
  cgwork = (double *)CALLOC((size_t) p*5,sizeof(double));
  Hp = (double *)CALLOC((size_t) p2,sizeof(double)); /* perturbed Hessian */
  for (io=ii=0,i=0;i<nm;i++) { /* loop over obs, k[ii] is start of neighbourhood of i */
    i0=ii;io0=io; /* record start of neigbourhood record */
    /* start with the change in gradient term */
    for (l=0;l<p;l++) g[l] = 0.0; /* have to clear first as multiple lps may be added */
    for (l=0;l<p2;l++) Hp[l] = H[l]; /* Hessian before dropping of neighbours */
    for (;ii<m[i];ii++) { /* loop over neighbours */
      ki=k[ii]; /* neighbour index */
      for (l2i=l=0;l<nlp;l++) { /* loop over linear predictors */
        jjl = jj[l];ln=l*n;
        for (j=0;j<plp[l];j++) g[jjl[j]] += X[n*jjl[j]+ki] * l1[ln+ki]; /* accumulate grad change */
	for (q=l;q<nlp;q++,l2i++) { /* Hessian block loop X[,jj[l]]' %*% l2*X[,jj[q]] */
	  jjq = jj[q];ln = l2i*n;
	  for (r=0;r<plp[q];r++) { /* NOTE: substantial optimization possible here once checked */
	    for (j=0;j<plp[l];j++) {
	      /* H is penalized *negative* Hessian, so dropped part needs to be *added* */
	      xx = X[jjl[j]*n+ki]*X[jjq[r]*n+ki]*l2[ln+ki]; 
              Hp[jjl[j]+jjq[r]*p] += xx;
	      if (q>l) Hp[jjl[j]*p+jjq[r]] += xx;
	    }
	  }  
	}  
      } /* lp loop */
    } /* neighbour loop */
    xx = 0.0;z=1.0;
    F77_CALL(dgemv)(&ntrans,&p,&p,&z,Hi,&p,g,&one,&xx,d,&one FCONE); /* initial step Hi g */
    kk=CG(Hp,Hi,g,d,p,1e-13,cgwork); /* d is approx change in beta caused by dropping y_i and its neighbours */
    if (kk<0) {
      Rprintf("npd! ");kk = -kk;
    }   
	
    if (iter < kk) iter=kk;
    /* now create the linear predictors for the ith point */
    for (;io<mi[i];io++)
    for (l=0;l<nlp;l++) {
      ln = no*l;jjl = jj[l];
      for (xx=0.0,j=0;j<plp[l];j++) {
	  q = jjl[j];
	  xx += X[n*q+ind[io]] * (beta[q]-d[q]); 
      }
      eta[ln+io] = xx;	
    }
    if (deriv>0) { /* get the derivatives of the linear predictors */
      for (l=0;l<nsp;l++) {
	ln = l*nlp;
	DH = VECTOR_ELT(dH, l);dh = REAL(DH);
	xx=0.0;z=1.0;
        F77_CALL(dgemv)(&ntrans,&p,&p,&z,dh,&p,d,&one,&xx,g,&one FCONE); /* g = dH_l d */
	//for (j=0;j<p;j++) dg[j] = 0.0; /* clear deriv of g before accumulation - NOT needed - straight into g! */
        for (i1=i0;i1<m[i];i1++) { /* neighbour loop */
	  ki=k[i1]; /* current neighbour of interest */
	  for (j=0;j<nlp;j++) for (q=0;q<nlp;q++) { /* have to do both triangles to avoid full matrix */
	    v = l3[i3f(j,q,0,nlp)*n+ki]*deta[ln*n+ki];
	    for (r=1;r<nlp;r++) v += l3[i3f(j,q,r,nlp)*n+ki]*deta[(ln+r)*n+ki];
	    jjl=jj[j];jjq=jj[q];
	    for (xx=0.0,kk=0;kk<plp[q];kk++) xx += X[ki+n*jjq[kk]]*d[jjq[kk]]; 
            xx *= v;
	    for (kk=0;kk<plp[j];kk++) g[jjl[kk]] -= X[ki+n*jjl[kk]]*xx;
	    /* now work on deriv of grad */
	    for (xx=0.0,kk=0;kk<plp[q];kk++) xx += X[ki+n*jjq[kk]]*db[jjq[kk]+p*l];
	    xx *= l2[ki+i2f(j,q,nlp)*n];
	    for (kk=0;kk<plp[j];kk++) g[jjl[kk]] += X[ki+n*jjl[kk]]*xx;
	  }  
	} /* neighbour loop */
	/* at this point g contains deriv of perturned Hessian w.r.t. log sp l multiplied by
           the change in beta vector, d. Now the derivative of grad vect w.r.t. log(sp[l])
           is also required */
	xx=0.0;z=1.0;
	F77_CALL(dgemv)(&ntrans,&p,&p,&z,Hi,&p,g,&one,&xx,d1,&one FCONE); /* initial step Hi d */
        kk=CG(Hp,Hi,g,d1,p,1e-13,cgwork); /* d1 is deriv wrt rho_l of approx change in beta caused by dropping the y_i and its neighbours */
	
	if (kk<0) {
          Rprintf("npdg! ");kk = -kk;
	}   
	
	if (iter1 < kk) iter1 = kk;
	for (io=io0;io<mi[i];io++)
        for (q=0;q<nlp;q++) {
	  jjq=jj[q];dbp = db + p*l;
	  for (xx=0.0,j=0;j<plp[q];j++) {
	      kk = jjq[j];
	      xx += X[ind[io]+kk*n] * (dbp[kk]-d1[kk]);
	  }  
          detacv[io + q*no +l*(no*nlp)] = xx;  
        }	  
      } /* l loop - smoothing parameters */	
    } else if (deriv<0) { /* the output beta perturbations to rows of detacv */
      for (j=0;j<p;j++) detacv[i+nm*j] = d[j];
    }  
  } /* main obs loop */
  for (l=0;l<nlp;l++) {
    /* iff coerceVector did not have to create a new vector then subtracting 1 from index will have changed original object in R, so need to 
       guard against that by adding one back on again */
    jjl=jj[l];for (i=0;i<plp[l];i++) jjl[i]++; 
  }  
  FREE(jj);FREE(plp);FREE(g);FREE(Hp);FREE(cgwork);
  PROTECT(kr=allocVector(INTSXP,2));
  INTEGER(kr)[0] = iter; /* max CG iterations used */
  INTEGER(kr)[1] = iter1; /* max CG iterations used for derivs*/
  UNPROTECT(5+nlp);
  return(kr);
} /* ncvls */  







SEXP Rncvls(SEXP x,SEXP JJ,SEXP R1,SEXP dH,SEXP L1, SEXP L2,SEXP L3,SEXP IND, SEXP MI, SEXP M, SEXP K,SEXP BETA,
	    SEXP ETACV,SEXP DETACV,SEXP DETA,SEXP DB,SEXP DERIV,SEXP EPS,SEXP NT) {
/* This computes the NCV for GAMLSS families, using the Cholesky factor updating approach. 

   l.p.s can share coefs in this version

   OMP parallelization
   
   Inputs:
   X[,jj[[i]]] is the model matrix for the ith linear predictor.
   R1 is the Cholesky factor of the penalized Hessian.
   dH[[i]] is the derivarive of H w.r.t. log(sp[i]); 
   lj contains jth derivatives of the likelihood w.r.t. the lp for each datum.
   w.r.t. the linear predictors. 
   k[m[i-1]:m[i]] index the neighbours omitted for the ith neighbourhood.
   ind[mi[i-1]:mi[i]] index the points to predict when the ith neighbourhood is omitted. 
   mi[-1] = m[-1]=0 by convention.
   beta - model coefficients. 
   deta and dbeta are matrices with derivatives of the lp's and coefs in their cols. deta has the lps stacked in each 
   column. 
   deriv = 0 for no derivative computations, 1 to compute derivatives (first order).
   eps is the machine precision.
   nt is the number of threads to use. 

   Outputs:
   The perturbed etas will be returned in etacv: if nm is the length of ind,then eta[q*nm+i] is the ith element of 
   qth perturbed linear predictor.
   If deriv>0 then the derivatives of the perturbed linear predictors are in detacv: detacv[q*nm+i + l*(np*nlp)]] is the ith element of 
   deriv of qth lp w.r.t. lth log sp.
   If deriv<0 then detacv is used to return the perturbed coefficint vectors for each neighbourhood/fold (in its rows)
*/
  double *X,*R,*l1,*l2,*l3,*beta,*g,*R0,xx,z,*d,*d1,*eta,*deta,v,*db,*dbp,*detacv,*dh,*b,alpha,alpha0,eps,*Rb,*ddbuf,*p0,*p3,*work;
  int **jj,*jjl,*jjq,*ind,*m,*k,n,p,nm,nlp,*plp,ii,i,j,i0,i1,l,ln,ki,p2,q,r,l2i,one=1,kk,nsp,*error,deriv,nddbuf,
    *mi,io,io0,no,pdef=1,maxn,buffer_size=0,nwork=0,nt=1,tid=0,*iwork=NULL,use_minres=1,niwork=0;
  SEXP JJp,kr,DH;
  char trans = 'T',ntrans = 'N',uplo='U',diag='N';
  p = length(BETA);p2 = p*p;
  n = nrows(x);deriv = asInteger(DERIV);
  eps = asReal(EPS);
  M = PROTECT(coerceVector(M,INTSXP));
  MI = PROTECT(coerceVector(MI,INTSXP));
  IND = PROTECT(coerceVector(IND,INTSXP));
  K = PROTECT(coerceVector(K,INTSXP)); /* otherwise R might be storing as double on entry */
  mi = INTEGER(MI); m = INTEGER(M); k = INTEGER(K);ind = INTEGER(IND);
  l1 = REAL(L1);l2=REAL(L2);
  nm = length(M);nlp = length(JJ);nsp = length(dH);no=length(IND);
  eta = REAL(ETACV);beta = REAL(BETA);R=REAL(R1);

  if (deriv>0) {
    l3=REAL(L3);deta=REAL(DETA);db = REAL(DB);
  } else l3=deta=db=NULL;
  if (deriv!=0) detacv=REAL(DETACV); else detacv = NULL;
  /* unpack the jj indices to here in order to avoid repeated list lookups withn loop */
  jj = (int **)CALLOC((size_t) nlp,sizeof(int *));
  plp = (int *)CALLOC((size_t) nlp,sizeof(int));
  for (l=0;l<nlp;l++) {
    JJp = VECTOR_ELT(JJ, l);JJp = PROTECT(coerceVector(JJp,INTSXP)); /* see R extensions 5.9.1 Handling the effects of garbage collection */
    plp[l] = length(JJp); jj[l] = INTEGER(JJp); /* jj[l][1:plp[l]] indexes cols of X for this lp */
    jjl=jj[l];for (i=0;i<plp[l];i++) jjl[i]--; /* 5.9.3 Details of R types. In fact coerceVector creates new vector only if type needs to change. */
  }
  /* need to know largest neighbourhood to create storage buffer for problematic updates*/
  maxn = ii = 0;
  for (j=0;j<nm;j++) {
    i = m[j]; if (i-ii>maxn) maxn = i-ii; ii = i;
  }
  #ifdef _OPENMP
  nt = asInteger(NT);
  #endif
  
  X = REAL(x);
  g = (double *)CALLOC((size_t) (3*p+nlp)*nt,sizeof(double)); /* gradient change */
  d = g + p*nt; /* change in beta */
  d1 = d + p*nt; /* deriv of above */
  b = d1 + p*nt; /* multipliers on spurious leading diagonal blocks - length nlp*nt */
  R0 = (double *)CALLOC((size_t) p2*nt,sizeof(double)); /* chol factor of perturbed Hessian */
  Rb = (double *)CALLOC((size_t) p2*nt,sizeof(double)); /* back up of chol factor of perturbed Hessian, in case of downdate failures */
  niwork = maxn*nlp*(nlp+1)/2; /* max rank of low rank downdate */
  if (use_minres) {
    nwork =  p*(niwork+7)+niwork;   
  } else {
    iwork = (int *)CALLOC((size_t)niwork*nt,sizeof(int));
    /* workspace query (d not referenced) ... */
    nwork = -1;z = 1.0;i=1;
    F77_CALL(dsysv)(&uplo,&niwork,&i,d,&niwork,iwork,d,&niwork,&xx,&nwork,&j FCONE);
    nwork = (int) ceil(xx);
    nwork += niwork*(p+1+niwork);
  }
  work = (double *)CALLOC((size_t) nwork*nt,sizeof(double));
  /* need to create buffer for skipped down-dates */
  buffer_size = p*niwork; 
  ddbuf = (double *)CALLOC((size_t) buffer_size*nt,sizeof(double));
  error = (int *)CALLOC((size_t)nt,sizeof(int));	
  #ifdef _OPENMP
  /* schedule: static, dynamic or guided - seems hard to do better than guided */
#pragma omp parallel for schedule(guided) private(i,ii,i0,io,io0,tid,l,pdef,nddbuf,ki,j,l2i,jjl,ln,alpha0,jjq,alpha,xx,p0,p3,q,i1,DH,dh,z,v,r,kk,dbp) num_threads(nt)
  #endif
  for (i=0;i<nm;i++) { /* loop over folds, k[ii] is start of neighbourhood of i */
    if (!i) io=ii=0;else {io=mi[i-1];ii=m[i-1];} 
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #endif
    i0=ii;io0=io; /* record start of neigbourhood record */
    /* start with the change in gradient term */
    for (p0=g+p*tid,l=0;l<p;l++) p0[l] = 0.0; /* have to clear g first as multiple lps may be added */
    for (p0=R0+p2*tid,l=0;l<p2;l++) p0[l] = R[l]; /* R0=R, Hessian Chol Factor before dropping of neighbours */
    pdef=1; /* positive definiteness status of R0  */
    nddbuf = 0; /* number of skipped down-dates */
    for (;ii<m[i];ii++) { /* loop over neighbours */
      ki=k[ii]; /* neighbour index */
      for (p0=b+nlp*tid,j=0;j<nlp;j++) p0[j] = 0.0; /* clear b */
      for (kk=0;kk<2;kk++) /* do updates when kk==0 and downdates when kk==1 */
      for (l2i=l=0;l<nlp;l++) { /* loop over linear predictors */
        jjl = jj[l];ln=l*n;
        if (!kk) for (p0=g+p*tid,j=0;j<plp[l];j++) p0[jjl[j]] += X[n*jjl[j]+ki] * l1[ln+ki]; /* accumulate grad change into g */
        alpha0 = l2[l2i*n+ki];l2i++;
	for (q=l+1;q<nlp;q++,l2i++) { /* Hessian block loop X[,jj[l]]' %*% l2*X[,jj[q]] */
          /* H is penalized *negative* Hessian, so dropped part needs to be *added* */
          jjq = jj[q];alpha = l2[l2i*n+ki];
	  /* Now record multipliers on spuriously generated leading diagonal blocks... */
	  if (!kk) { *(b+l+nlp*tid) += alpha;*(b+q+nlp*tid) += alpha;}
          if ((!kk&&alpha>0)||(kk&&alpha<0)) { /* updates only on first pass, downdates only on second */
	    xx = sqrt(fabs(alpha));
	    for (p0=d+p*tid,j=0;j<p;j++) p0[j] = 0.0; /* clear d */
	    for (j=0;j<plp[l];j++) p0[jjl[j]] =  X[jjl[j]*n+ki]*xx; /* fill d */
	    for (j=0;j<plp[q];j++) p0[jjq[j]] +=  X[jjq[j]*n+ki]*xx; /* fill d */
            if (alpha>0) j=1; else {
	      for (p0=Rb+p2*tid,p3=R0+tid*p2,j=0;j<p;j++,p0+=p,p3+=p) for (q=0;q<=j;q++) p0[q] = p3[q]; /* backup state of R0 before attempting downdate */
              j=0; /* downdate indicator*/
	    }  
	    chol_up(R0+p2*tid,d+tid*p,&p,&j,&eps);/* add/subtract correction */
	    if (!j && *(R0+1+p2*tid) < -0.5) { /* downdate not pdef - restore R0 and store downdate */
              pdef=0;*(R0+1+p2*tid) = 0.0;
	      for (p0=Rb+tid*p2,p3=R0+tid*p2,j=0;j<p;j++,p0+=p,p3+=p) for (q=0;q<=j;q++) p3[q] = p0[q]; /* restore to state before update attempt */
	      for (p0=ddbuf+p*nddbuf+tid*buffer_size,p3=d+tid*p,j=0;j<p;j++) p0[j] = p3[j]; /* store the skipped update */
              nddbuf++;
	    }  
	  }  
	}  /* lp loop q */
	/* now update the lth leading diagonal block, removing the nuisance block already formed */
	alpha = alpha0 - *(b+l+tid*nlp);
	if ((!kk&&alpha>0)||(kk&&alpha<0)) { /* updates only on first pass, downdates only on second */
	  xx = sqrt(fabs(alpha));for (p0=d+tid*p,j=0;j<p;j++) p0[j] = 0.0;
	  for (j=0;j<plp[l];j++) p0[jjl[j]] =  xx*X[jjl[j]*n+ki];
	
	  if (alpha>0) { /* add or subtract correction? */
	    j=1; /* add */
	  } else { /* subtract */
            for (p0=Rb+p2*tid,p3=R0+tid*p2,j=0;j<p;j++,p0+=p,p3+=p) for (q=0;q<=j;q++) p0[q] = p3[q]; /* backup state of R0 before attempting downdate */
            j=0; 
	  }
	  chol_up(R0+tid*p2,d+tid*p,&p,&j,&eps);
	  if (!j && *(R0+1+p2*tid) < -0.5) { /* did update fail to be positive definite? */
	    pdef=0;*(R0+1+p2*tid) = 0.0;
	 
	    for (p0=Rb+tid*p2,p3=R0+tid*p2,j=0;j<p;j++,p0+=p,p3+=p) for (q=0;q<=j;q++) p3[q] = p0[q]; /* restore to state before update attempt */
	    for (p0=ddbuf+p*nddbuf+tid*buffer_size,p3=d+tid*p,j=0;j<p;j++) p0[j] = p3[j]; /* store the skipped update */
            nddbuf++;
	  }  
        } 
      } /* lp loop l */
    
    } /* neighbour loop */
    
    if (pdef) { /* solve for R0'R0 d = g - the change in beta caused by dropping neighbourhood i */
      for (p0=d+tid*p,p3=g+tid*p,j=0;j<p;j++) p0[j] = p3[j];
      F77_CALL(dtrsv)(&uplo,&trans,&diag,&p,R0+tid*p2,&p,d+tid*p,&one FCONE FCONE FCONE);
      F77_CALL(dtrsv)(&uplo,&ntrans,&diag,&p,R0+tid*p2,&p,d+tid*p,&one FCONE FCONE FCONE);
    } else {  /* fallback solve (R0'R0 - uu') d= g where u is stored in ddbuf */
      error[tid]++; /* count the number of non +ve def cases */ 
    
      j = nddbuf; /* modified to number of iterations on exit */
      if (use_minres) minres(R0+tid*p2,ddbuf+buffer_size*tid,g+tid*p,d+tid*p,&p,&j,work+nwork*tid); else {
        *(iwork+niwork*tid) = nwork;
	woodbury(R0+tid*p2,ddbuf+buffer_size*tid,g+tid*p,d+tid*p,&p,&j,work+nwork*tid,iwork+niwork*tid);
      }	
      //Rprintf(" fold=%d iter=%d updates=%d\n",i,j,nddbuf);
    }
	
    /* now create the linear predictors for the ith fold. d contains the perturbation to beta */
    for (;io<mi[i];io++)
    for (l=0;l<nlp;l++) {
      ln = no*l;jjl = jj[l];i1=ind[io];
      for (p0=d+tid*p,xx=0.0,j=0;j<plp[l];j++) {
	  q = jjl[j];
	  xx += X[n*q+i1] * (beta[q]-p0[q]); 
      }
      eta[ln+io] = xx;	
    }
    if (deriv>0) { /* get the derivatives of the linear predictors */
      for (l=0;l<nsp;l++) {
	ln = l*nlp;
	DH = VECTOR_ELT(dH, l);dh = REAL(DH);
	xx=0.0;z=1.0;
        F77_CALL(dgemv)(&ntrans,&p,&p,&z,dh,&p,d+tid*p,&one,&xx,g+tid*p,&one FCONE); /* g = dH_l d */
        for (i1=i0;i1<m[i];i1++) { /* neighbour loop */
	  ki=k[i1]; /* current neighbour of interest */
	  for (j=0;j<nlp;j++) for (q=0;q<nlp;q++) { /* have to do both triangles to avoid needing full matrix */
	    v = l3[i3f(j,q,0,nlp)*n+ki]*deta[ln*n+ki];
	    for (r=1;r<nlp;r++) v += l3[i3f(j,q,r,nlp)*n+ki]*deta[(ln+r)*n+ki];
	    jjl=jj[j];jjq=jj[q];
	    for (p0=d+tid*p,xx=0.0,kk=0;kk<plp[q];kk++) xx += X[ki+n*jjq[kk]]*p0[jjq[kk]]; 
            xx *= v;
	    for (p0=g+tid*p,kk=0;kk<plp[j];kk++) p0[jjl[kk]] -= X[ki+n*jjl[kk]]*xx;
	    /* now work on deriv of grad */
	    for (xx=0.0,kk=0;kk<plp[q];kk++) xx += X[ki+n*jjq[kk]]*db[jjq[kk]+p*l];
	    xx *= l2[ki+i2f(j,q,nlp)*n];
	    for (kk=0;kk<plp[j];kk++) p0[jjl[kk]] += X[ki+n*jjl[kk]]*xx;
	  }  
	} /* neighbour loop */
	/* at this point g contains deriv of perturbed Hessian w.r.t. log sp l multiplied by
           the change in beta vector, d. Now the derivative of grad vect w.r.t. log(sp[l])
           is also required */
        if (pdef) { /* solve for R0'R0 d1 = g1 */
	  for (p0=d1+tid*p,p3=g+tid*p,j=0;j<p;j++) p0[j] = p3[j];
          F77_CALL(dtrsv)(&uplo,&trans,&diag,&p,R0+tid*p2,&p,d1+tid*p,&one FCONE FCONE FCONE);
          F77_CALL(dtrsv)(&uplo,&ntrans,&diag,&p,R0+tid*p2,&p,d1+tid*p,&one FCONE FCONE FCONE);
        } else {  /* fallback solve (R0'R0 -uu')d1 = g1 where u stored in ddbuf*/
          j = nddbuf; 
	  if (use_minres) minres(R0+tid*p2,ddbuf+tid*buffer_size,g+tid*p,d1+tid*p,&p,&j,work+nwork*tid); else {
            *(iwork+niwork*tid) = nwork;
	    woodbury(R0+tid*p2,ddbuf+tid*buffer_size,g+tid*p,d1+tid*p,&p,&j,work+nwork*tid,iwork+niwork*tid);
	  }  
	  //Rprintf(" %d",j);
        }
	
	for (io=io0;io<mi[i];io++)
        for (q=0;q<nlp;q++) {
	  jjq=jj[q];dbp = db + p*l;
	  for (p0=d1+tid*p,xx=0.0,j=0;j<plp[q];j++) {
	      kk = jjq[j];
	      xx += X[ind[io]+kk*n] * (dbp[kk]-p0[kk]);
	  }  
          detacv[io + q*no +l*(no*nlp)] = xx;  
        }	  
      } /* l loop - smoothing parameters */	
    } else if (deriv<0) { /* write the perturbations of the coefs out to deta.cv (nm by p) */
      for (p0=d+tid*p,j=0;j<p;j++) detacv[i+nm*j] = p0[j]; 
    }  
  } /* main obs loop */
  for (l=0;l<nlp;l++) {
    /* iff coerceVector did not have to create a new vector then subtracting 1 from index will have changed original object in R, so need to 
       guard against that by adding one back on again */
    jjl=jj[l];for (i=0;i<plp[l];i++) jjl[i]++; 
  }  
  FREE(jj);FREE(plp);FREE(g);FREE(R0);
  FREE(Rb);FREE(ddbuf);FREE(work);if (!use_minres) FREE(iwork);
  for (j=i=0;i<nt;i++) j+=error[i];
  FREE(error);
  PROTECT(kr=allocVector(INTSXP,1));
  INTEGER(kr)[0] = j; /* number of Hessians not positive def */
  UNPROTECT(5+nlp);
  return(kr);
} /* Rncvls */  

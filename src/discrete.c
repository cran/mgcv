/* (c) Simon N. Wood (2015) Released under GPL2 */

/* Routines to work with discretized covariate models 
   Data structures:
   * X contains sub matrices making up each term 
   * nx is number of sub-matrices.
   * jth matrix is m[j] by p[j]
   * k contains the index vectors for each matrix in Xb.
   * ts[i] is starting matrix of ith term    
   * dt[i] is number of matrices making up ith term.
   * n_terms is the number of terms. 
   * Q contains reparameterization matrices for each term
   * qs[i] is starting address within Q for ith repara matrix
   R CMD SHLIB discrete2.c NOTE: *Must* put CFLAGS = -fopenmp -O2 in .R/Makevars 

   Here is an interesting bit of code for converting an index kk
   running over the upper triangle of an nt by nt matrix to rows 
   and columns (first row corresponds to kk = 0, 1, 2 ,...)
   i=kk;r=0;while (i >= *nt-r) { i -= *nt - r; r++;} c = r + i; 
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <Rconfig.h>
#include "mgcv.h"
#ifdef OPENMP_ON
#include <omp.h>
#endif



/* basic extraction operations */ 

void singleXj(double *Xj,double *X, int *m, int *k, int *n,int *j) {
/* Extract a column j of matrix stored in compact form in X, k into Xj. 
   X has m rows. k is of length n. ith row of result is Xj = X[k(i),j]
   (an n vector). This function is O(n). Thread safe.
*/
  double *pe; 
  X += *m * *j; /* shift to start of jth column */ 
  for (pe = Xj + *n;Xj < pe;Xj++,k++) *Xj = X[*k]; 
} /* singleXj */

void tensorXj(double *Xj, double *X, int *m, int *p,int *dt, 
              int *k, int *n, int *j, int *kstart,int *koff) {
/* Extract a column j of tensor product term matrix stored in compact 
   form in X, k into Xj. There are dt sub matrices in Xj. The ith is 
   m[i] by p[i]. There are dt index n - vectors stacked end on end in 
   k. The ith component (starting at 0) has index vector at column 
   kstart[i] + *koff of k. This function is O(n*dt)

   This routine performs pure extraction only if Xj is a vector of 1s on 
   entry. Otherwise the jth column is multiplied element wise by the 
   contents of Xj on entry. Thread safe.
*/ 
  int q=1,l,i,jp,*kp;
  double *p0,*p1,*M;
  p1 = Xj + *n; /* pointer for end of Xj */
  for (i = 0;i < *dt;i++) q *= p[i];
  jp = *j; 
  for (i = 0;i < *dt; i++) {
    q /= p[i]; /* update q */
    l = jp/q; /* column of current marginal */
    jp = jp%q;
    M = X + m[i] * l; /* M now points to start of col l of ith marginal model matrix */
    kp = k + (kstart[i] + *koff) * (ptrdiff_t) *n;
    for (p0=Xj;p0<p1;p0++,kp++) *p0 *= M[*kp];
    //k += *kjump * (ptrdiff_t) *n; /*k cols relating to different marginals may be separated by kjump cols*/
    X += m[i] * p[i]; /* move to the next marginal matrix */
  }
} /* tensorXj */

void singleXty(double *Xy,double *temp,double *y,double *X, int *m,int *p, int *k, int *n, int *add) {
/* forms X'y for a matrix stored in packed form in X, k, with ith row 
   X[k[i],]. X as supplied is m by p, while k is an n vector.
   Xy and temp are respectively p and m vectors and do not need to be cleared on entry.
   Thread safe.
*/
  double *p0,*p1,done=1.0,dzero=0.0; 
  char trans = 'T';
  int one=1;
  for (p0=temp,p1 = p0 + *m;p0<p1;p0++) *p0 = 0.0;
  for (p1=y + *n;y<p1;y++,k++) temp[*k] += *y;
  if (*add) dzero = 1.0;
  F77_CALL(dgemv)(&trans, m, p,&done,X,m,temp,&one,&dzero,Xy,&one);
} /*singleXty*/

void tensorXty(double *Xy,double *work,double *work1, double *y,double *X, 
               int *m, int *p,int *dt,int *k, int *n, int *add, int *kstart,int *koff) {
/* forms X'y for a matrix stored as a compact tensor product in X, k.
   There are dt maginal matrices packed in X, the ith being m[i] by p[i].
   y and work are n vectors. work does not need to be cleared on entry.
   work1 is an m vector where m is the dimension of the final marginal.
   Note: constraint not dealt with here. Thread safe.
   k[,kstart[i] + *koff] is the index for the ith term (starting at i=0)
*/
  int pb=1,i,j,pd;
  double *p1,*yn,*p0,*M; 
  yn = y + *n;
  M = X;
  for (i=0;i<*dt-1;i++) { 
    pb *= p[i]; 
    M += m[i] * p[i]; /* shift M to final marginal */ 
  }
  pd = p[*dt-1];
  for (i=0;i<pb;i++) {
    /* extract the ith column of M_0 * M_1 * ... * M_{d-2} * y, where 
       '*' is the row tensor product here */
    for (p0=y,p1=work;p0<yn;p0++,p1++) *p1 = *p0; /* copy y to work */ 
    j = *dt - 1; 
    tensorXj(work,X,m,p,&j,k,n,&i,kstart,koff);
    /* now form M'work */
    singleXty(Xy+i*pd,work1,work,M,m + *dt-1,&pd, k + (ptrdiff_t) *n * (kstart[j] + *koff) ,n,add);
  }
} /* tensorXty */

void singleXb(double *f,double *work,double *X,double *beta,int *k,int *m, int *p,int *n,
              int *kstart, int *kstop) {
/* Forms X beta, where beta is a p - vector, and X is stored in compact form in 
   X, k, with ith row sum_j X[k[i,j],]. X is stored as an m by p matrix. k is the matrix
   containing the indices. indices for this term start at column kstart and end 
   at column kstop-1 (of k). often there is only one j. 
   work is an m vector. Thread safe.
*/
  char trans='N';
  double done=1.0,dzero=0.0,*p1,*fp;
  int one=1,j;
  F77_CALL(dgemv)(&trans, m, p,&done,X,m,beta,&one,&dzero,work,&one);
  p1 = f + *n;
  fp = f;
  k += *kstart * (ptrdiff_t) *n;
  for (;fp < p1;fp++,k++) *fp = work[*k];
  for (j=1;j < *kstop - *kstart;j++) for (fp=f;fp < p1;fp++,k++) *fp += work[*k];
} /* singleXb */

void tensorXb(double *f,double *X, double *C,double *work, double *beta,
              int *m, int *p,int *dt,int *k, int *n,double *v,int *qc,
              int *kstart, int *kstop) {
/* for X* beta where X* is a tensor product term with nt marginal model matrices
   stored in X. ith such is m[i] by p[i]. 
   work is an n vector. C is an m[d] by pb working matrix.  
   v is a vector such that if Q=I-vv' and Z is Q with the first column dropped then 
   Z is a null space basis for the identifiability constraint. If *qc <= 0 then
   no constraint is applied. Thread safe.
   k[,kstart[i]...kstop[i]-1] are the indices for the ith component of this term 
   (starting at i=0).
*/ 
  char trans='N';
  int pb=1,md,*kp,*kd,pd,i,j,q;
  double *M,done=1.0,dzero=0.0,*p0,*p1,*p2,*p3,*pf,*pc,x;
  M = X;
  for (i=0;i<*dt-1;i++) {
    pb *= p[i];
    M += m[i] * (ptrdiff_t) p[i]; /* final marginal */ 
  }
  md = m[*dt - 1];
  pd = p[*dt -1];
  //kd = k + (*dt-1) * (ptrdiff_t) *n; /* index vector for final term */
  kd = k + kstart[*dt-1] * (ptrdiff_t) *n;
  /* form work = M B, where vec(B) = beta */
  if (*qc<=0) { /* no constraint supplied */
    F77_CALL(dgemm)(&trans,&trans,&md,&pb, &pd, &done,
		    M,&md,beta,&pd,&dzero,C,&md); 
  } else { /* there is a constraint matrix */
    /* first map supplied beta to unconstrained parameterization */ 
    j = pb * pd; /* total number of coeffs - length of unconstrained beta */
    
    *work = 0.0;x=0.0;
    for (p0=work+1,p1=p0+j-1,p2=beta,p3=v+1;p0<p1;p0++,p2++,p3++) { 
      *p0 = *p2; 
      x += *p0 * *p3; /* v'beta where beta padded with extra zero at start */ 
    }
    for (p0=work,p1=p0+j,p2=v;p0<p1;p0++,p2++) *p0 -= *p2 * x; /* (I-vv')(0,beta')' */

    /*F77_CALL(dgemv)(&trans, &j, qc,&done,Q,&j,beta,&one,&dzero,work,&one); old when Q full matrix */
    F77_CALL(dgemm)(&trans,&trans,&md,&pb, &pd, &done,
		    M,&md,work,&pd,&dzero,C,&md);
  }
  p1 = work + *n;
  for (pf=f,p0=f + *n;pf<p0;pf++) *pf = 0.0;
  for (q = 0;q < *kstop - *kstart;q++) /* loop over possible multiple indices */
  for (i=0;i<pb;i++) {
    /* extract ith col of truncated tensor product (final marginal omitted) */
    for (p0=work;p0<p1;p0++) *p0 = 1.0; /* set work to 1 */ 
    j = *dt - 1; 
    tensorXj(work,X,m,p,&j,k,n,&i,kstart,&q);
    pc = C + i * md; /* ith col of C */
    for (kp=kd + q  *(ptrdiff_t) *n,pf=f,p0=work;p0<p1;p0++,pf++,kp++) {
      *pf += pc[*kp] * *p0;
    }
  }
} /* tensorXb */

void Xbd(double *f,double *beta,double *X,int *k,int *ks, int *m,int *p, int *n, 
	 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *bc) {
/* Forms f = X beta for X stored in the packed form described in function XWX
   bc is number of cols of beta and f... 
   LIMITED thread safety. Allocates/frees memory using R_chk_calloc/R_chk_free (usually),
   protected within critical sections. Not safe to use with other routines allocating memory using these routines
   within a parallel section. Safe to use as the only allocator of memory within a parallel section. 
*/
  ptrdiff_t *off,*voff;
  int i,j,q,*pt,*tps,dC=0,c1,first;
  double *f0,*pf,*p0,*p1,*p2,*C=NULL,*work,maxp=0,maxrow=0;
  /* obtain various indices */
#pragma omp critical (xbdcalloc)
  { pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
    off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
    voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvector starts */
    tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts in param vector or XWy */
  }
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) { /* work through components of each term */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (maxrow<m[q]) maxrow=m[q];
      if (j>0 && j==dt[i]-1) {
        c1 = pt[i] * (ptrdiff_t) m[q]; 
        if (c1>dC) dC = c1; /* dimension of working matrix C */
      }
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
    } 
    if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else voff[i+1] = voff[i]; /* start of ith v matrix */
    if (maxp < pt[i]) maxp = pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
  }
  /* now form the product term by term... */ 
  i = *n; if (i<maxp) i=maxp; if (i<maxrow) i=maxrow;
#pragma omp critical (xbdcalloc)
  { pf=f0 = (double *)CALLOC((size_t)*n,sizeof(double));
    work = (double *)CALLOC((size_t)i,sizeof(double));
    if (dC) C = (double *)CALLOC((size_t)dC,sizeof(double));
  }
  for (j=0;j < *bc;j++) { /* loop over columns of beta */
    first = 1;
    for (i=0;i < *nt;i++) { /* work through terms */  
        if (first) f0 = f; /* result written straight to f for i==0 */
        if (dt[i]==1) singleXb(f0,work,X+off[ts[i]],beta+tps[i],k /*+ (ptrdiff_t) *n *ts[i]*/,
			       m+ts[i], p+ts[i],n,ks + ts[i], ks + ts[i] + *nx); 
        else tensorXb(f0,X+off[ts[i]],C,work, beta+tps[i],m+ts[i], p+ts[i],dt+i,k 
		      /*+ (ptrdiff_t) *n * ts[i]*/,n,v+voff[i], qc+i,ks + ts[i], ks + ts[i] + *nx);
        if (!first) {
          for (p0=f,p1=f + *n,p2=f0;p0<p1;p0++,p2++) *p0 += *p2; /* f <- f + f0 */     
        } else { 
          f0=pf; /* restore f0 */
          first=0;
        }
    } /* term loop */
    f += *n;beta += tps[*nt]; /* move on to next column */
  } /* col beta loop */
#pragma omp critical (xbdcalloc)
  { if (dC) FREE(C);
    FREE(work);FREE(f0);
    FREE(pt);FREE(off);FREE(voff);FREE(tps);
  }
} /* Xb */

void diagXVXt(double *diag,double *V,double *X,int *k,int *ks,int *m,int *p, int *n, 
	      int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pv,int *nthreads) {
/* Forms diag(XVX') where X is stored in the compact form described in XWXd.
   V is a pv by pv matrix. 
   Parallelization is by splitting the columns of V into nthreads subsets.
   Currently inefficient. Could be speeded up by a factor of 2, by supplying a
   square root of V in place of V.
*/
  double *xv,*dc,*p0,*p1,*p2,*p3,*ei,*xi;
  ptrdiff_t bsj,bs,bsf,i,j,kk;
  int one=1;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  if (*nthreads > *pv) *nthreads = *pv;
  xv = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for cols of XV */
  xi = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for cols of X */
  ei = (double *) CALLOC((size_t) *nthreads * *pv,sizeof(double)); /* storage for identity matrix cols */
  dc = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for components of diag */
  if (*nthreads>1) {
    bs = *pv / *nthreads;
    while (bs * *nthreads < *pv) bs++;
    while (bs * *nthreads - bs >= *pv) (*nthreads)--;
    bsf = *pv - (bs * *nthreads - bs); 
  } else {
    bsf = bs = *pv;
  }
  #ifdef OPENMP_ON
  #pragma omp parallel for private(j,bsj,i,kk,p0,p1,p2,p3) num_threads(*nthreads)
  #endif  
  for (j=0;j < *nthreads;j++) {
    if (j == *nthreads - 1) bsj = bsf; else bsj = bs;
    for (i=0;i<bsj;i++) { /* work through this block's columns */
      kk = j * bs + i;
      ei[j * *pv + kk] = 1;if (i>0) ei[j * *pv + kk - 1] = 0;
      /* Note thread safety of XBd means this must be only memory allocator in this section*/
      Xbd(xv + j * *n,V + kk * *pv,X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,&one); /* XV[:,kk] */
      Xbd(xi + j * *n,ei + j * *pv,X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,&one); /* X[:,kk] inefficient, but deals with constraint*/
      p0 = xi + j * *n;p1=xv + j * *n;p2 = dc + j * *n;p3 = p2 + *n;
      for (;p2<p3;p0++,p1++,p2++) *p2 += *p0 * *p1; /* elementwise product of XV[:,kk] X[:,kk] */
    } 
  } /* parallel loop end */
  /* sum the contributions from the different threads into diag... */
  for (p0=diag,p1=p0+ *n,p2=dc;p0<p1;p0++,p2++) *p0 = *p2;
  for (i=1;i< *nthreads;i++) for (p0=diag,p1=p0+ *n;p0<p1;p0++,p2++) *p0 += *p2;
  FREE(xv);FREE(dc);FREE(xi);FREE(ei);
} /* diagXVXt */


void XWyd(double *XWy,double *y,double *X,double *w,int *k,int *ks, int *m,int *p, int *n, 
         int *nx, int *ts, int *dt, int *nt,double *v,int *qc,
         int *ar_stop,int *ar_row,double *ar_weights) {
/* NOT thread safe */
  double *Wy,*p0,*p1,*p2,*p3,*Xy0,*work,*work1,x;
  ptrdiff_t i,j,*off,*voff;
  int *tps,maxm=0,maxp=0,one=1,zero=0,*pt,add,q;
  if (*ar_stop>=0) { /* model has AR component, requiring sqrt(weights) */
    for (p0 = w,p1 = w + *n;p0<p1;p0++) *p0 = sqrt(*p0);
  }
  /* obtain various indices */
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvector starts */
  tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts in param vector or XWy */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) { /* work through components of each term */
      off[q+1] = off[q] + p[q]*m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
      if (maxm<m[q]) maxm=m[q];
    } 
    if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else voff[i+1] = voff[i];  /* start of ith Q matrix */
    if (maxp < pt[i]) maxp=pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
  }
  Xy0 =  (double *) CALLOC((size_t)maxp,sizeof(double));
  work =  (double *) CALLOC((size_t)*n,sizeof(double));
  work1 = (double *) CALLOC((size_t)maxm,sizeof(double));
  /* apply W to y */
  Wy = (double *) CALLOC((size_t)*n,sizeof(double)); /* Wy */
  for (p0=Wy,p1=Wy + *n,p2=w;p0<p1;p0++,y++,p2++) *p0 = *y * *p2; 
  if (*ar_stop>=0) { /* AR components present (weights are sqrt, therefore) */
    rwMatrix(ar_stop,ar_row,ar_weights,Wy,n,&one,&zero,work);
    rwMatrix(ar_stop,ar_row,ar_weights,Wy,n,&one,&one,work); /* transpose of transform applied */
    for (p0=w,p1=w + *n,p2=Wy;p0<p1;p0++,p2++) *p2 *= *p0; /* sqrt weights again */
  }
  /* now loop through terms applying the components of X'...*/
  for (i=0;i<*nt;i++) { /* term loop */ 
   add=0; 
   if (dt[i]>1) { /* it's a tensor */
      //tensorXty(Xy0,work,work1,Wy,X+off[ts[i]],m+ts[i],p+ts[i],dt+i,k+ts[i] * (ptrdiff_t) *n,n);
      for (q=0;q<ks[ts[i] + *nx]-ks[ts[i]];q++) {  /* loop through index columns */
        tensorXty(Xy0,work,work1,Wy,X+off[ts[i]],m+ts[i],p+ts[i],dt+i,k,n,&add,ks+ts[i],&q);
        add=1;
      }
      if (qc[i]>0) { /* there is a constraint to apply Z'Xy0: form Q'Xy0 and discard first row... */
        /* Q' = I - vv' */
        for (x=0.0,p0=Xy0,p1=p0 + pt[i],p2=v+voff[i];p0<p1;p0++,p2++) x += *p0 * *p2; /* x = v'Xy0 */
        p0=XWy + tps[i];p1 = p0 + pt[i]-1;p2 = v+voff[i] + 1;p3=Xy0+1;
        for (;p0<p1;p0++,p2++,p3++) *p0 = *p3 - x * *p2; /* (I-vv')Xy0 less first element */
      } else { /* straight copy */
        for (p0=Xy0,p1=p0+pt[i],p2=XWy+tps[i];p0<p1;p0++,p2++) *p2 = *p0;
      }
    } else { /* it's a singleton */
      //singleXty(XWy+tps[i],work1,Wy,X+off[ts[i]], m+ts[i],p+ts[i], k+ts[i] * (ptrdiff_t) *n,n);
      for (q=ks[ts[i]];q<ks[ts[i] + *nx];q++) { /* loop through index columns */  
        singleXty(XWy+tps[i],work1,Wy,X+off[ts[i]], m+ts[i],p+ts[i], k + q * (ptrdiff_t) *n,n,&add);
        add=1;
      }
    }
  } /* term loop */
  FREE(Wy); FREE(Xy0); FREE(work); FREE(work1); 
  FREE(pt); FREE(off); FREE(voff); FREE(tps);
} /* XWy */

void XWXd(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, int *n, int *nx, int *ts, 
          int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,int *ar_row,double *ar_weights) {
/* This version uses open MP by column, within sub-block.

   An alternative in which one thread did each sub-block scaled poorly, largely 
   because the work is very uneven between sub-blocks. 

   Forms Xt'WXt when Xt is divided into blocks of columns, each stored in compact form
   using arguments X and k. 'X' contains 'nx' blocks, the ith is an m[i] by p[i] matrix. 
   k contains nx index vectors. If kj and Xj represent the jth sub-matrix and index vector 
   then the ith row of the corresponding full sub-matrix is Xj[kj[i],]. The blocks of Xt 
   need not be directly given by terms like Xj[kj[i],], some are actually given by row 
   tensors of several such matrices. Conceptually Xt has nt <= nx blocks (terms) and the ith 
   block starts at sub-matrix  ts[i] of X, and is made up of the row tensor product of 
   dt[i] consecutive sub-matrices. 
    
   Tensor product terms may have constraint matrices Z, which post multiply the tensor product 
   (typically imposing approximate sum-to-zero constraints). Actually Z is Q with the first column 
   dropped where Q =  I - vv'. qc[i]==0 for singleton terms.  

   AR models are handled via the 3 ar_* arrays. ar_stop[0] < 0 signals no AR. 
  
*/  
  int r,c,i,j,q,*pt,*pd,a,b,*tps,ptot,maxp=0,maxm=0,pa,pb,kk,dk,*start,one=1,zero=0,add; 
  ptrdiff_t *off,*voff;
  double *p0,*p1,*p2, *Xi, *Xj, *temp,*tempn,*xwx,*xwx0,
    *XiB,*XjB,*tempB,*tempnB,*x0,*x1,x;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  start = (int *) CALLOC((size_t)(*nthreads+1),sizeof(int));
  XiB = (double *) CALLOC((size_t)(*n * *nthreads),sizeof(double)); /* column of X */ 
  XjB = (double *) CALLOC((size_t)(*n * *nthreads),sizeof(double)); /* column of X */
  tempnB = (double *) CALLOC((size_t)(*n * *nthreads),sizeof(double));
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  pd = (int *) CALLOC((size_t)*nt,sizeof(int)); /* storage for last marginal size */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvectors starts */
  tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts */
  if (*ar_stop>=0) { /* model has AR component, requiring sqrt(weights) */
    for (p0 = w,p1 = w + *n;p0<p1;p0++) *p0 = sqrt(*p0);
  }
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) {
      if (j==dt[i]-1) pd[i] = p[q]; /* the relevant dimension for deciding product ordering */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
      if (maxm<m[q]) maxm=m[q];
    } 
    if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else voff[i+1] = voff[i]; /* start of ith v vector */
    if (maxp<pt[i]) maxp=pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
  }
  tempB = (double *) CALLOC((size_t) maxm * *nthreads,sizeof(double));
  xwx = (double *) CALLOC((size_t) maxp * maxp,sizeof(double)); /* working cross product storage */
  xwx0 = (double *) CALLOC((size_t) maxp * maxp,sizeof(double)); /* working cross product storage */
  ptot = tps[*nt]; /* total number of parameters */
  for (r=0;r < *nt;r++) for (c=r;c< *nt;c++) { /* the block loop */
    if (pd[r]>pd[c]) { /* Form Xr'WXc */
      a=r;b=c;
    } else { /* Form Xc'WXr */
      a=c;b=r; 
    }
    /* split cols between threads... */  
    dk = pt[b] / *nthreads; //rk = pt[b] % *nthreads;
    if (dk * *nthreads < pt[b]) dk++;
    start[0]=0; 
    for (i=0;i<*nthreads;i++) { 
      start[i+1] = start[i] + dk;
      if (start[i+1]>pt[b]) start[i+1]=pt[b];
    }
    #ifdef OPENMP_ON
    #pragma omp parallel private(Xi,Xj,i,q,add,temp,tempn,p0,p1,p2) num_threads(*nthreads)
    #endif 
    { /* begin parallel section */
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (kk=0;kk<*nthreads;kk++) { 
        /* allocate thread specific storage... */
        temp = tempB + kk * (ptrdiff_t) maxm;
        Xi = XiB + kk * (ptrdiff_t) *n; Xj = XjB + kk * (ptrdiff_t) *n;
        tempn = tempnB + kk * (ptrdiff_t) *n;
        for (i=start[kk];i<start[kk+1];i++) { /* loop over columns of Xb */ 
          /* extract Xb[:,i]... */
          if (ks[ts[b]]==ks[ts[b] + *nx]-1) { /* single index column */
            if (dt[b]>1) { /* tensor */
              for (p0=Xi,p1=p0+*n;p0<p1;p0++) *p0 = 1.0; /* set Xi to 1 */
	      tensorXj(Xi,X+off[ts[b]], m + ts[b] , p + ts[b], dt+b, 
	      		    k, n,&i,ks + ts[b],&zero);
            } else { /* singleton */
              singleXj(Xi,X+off[ts[b]], m + ts[b], k + (ptrdiff_t)*n * ks[ts[b]] /* ts[b]*/, n,&i);
            }
          } else { /* need to sum over multiple cols of k */
            for (q = 0;q<ks[ts[b] + *nx]-ks[ts[b]];q++) { /* k column loop */
              if (dt[b]>1) { /* tensor */
                for (p0=Xj,p1=p0+*n;p0<p1;p0++) *p0 = 1.0; /* set Xj to 1 */
                tensorXj(Xj,X+off[ts[b]], m + ts[b] , p + ts[b], dt+b, 
			    k, n,&i,ks+ts[b],&q);
              } else { /* singleton */
                singleXj(Xj,X+off[ts[b]], m + ts[b], k + (ptrdiff_t)*n * (q + ks[ts[b]]) /* ts[b]*/, n,&i);
              }
              if (q == 0) for (p0=Xj,p1=p0+*n,p2=Xi;p0<p1;p0++,p2++) *p2 = *p0;
              else for (p0=Xj,p1=p0+*n,p2=Xi;p0<p1;p0++,p2++) *p2 += *p0;
            } /* k column loop */
          }
          /* apply W to Xi - for AR models this is sqrt(W) */
          for (p0=w,p1=w + *n,p2=Xi;p0<p1;p0++,p2++) *p2 *= *p0; 
          if (*ar_stop>=0) { /* AR components present (weights are sqrt, therefore) */
            rwMatrix(ar_stop,ar_row,ar_weights,Xi,n,&one,&zero,tempn);
            rwMatrix(ar_stop,ar_row,ar_weights,Xi,n,&one,&one,tempn); /* transpose of transform applied */
            for (p0=w,p1=w + *n,p2=Xi;p0<p1;p0++,p2++) *p2 *= *p0; /* sqrt weights again */
          }
          /* now form Xa'WXb[:,i]... */ 
          add=0; 
          for (q=0;q<ks[ts[a] + *nx]-ks[ts[a]];q++) {
            if (dt[a]>1) { /* tensor */
              tensorXty(xwx + i * pt[a],tempn,temp,Xi,X+off[ts[a]],m+ts[a],p+ts[a],
	      			    dt+a,k, n,&add,ks+ts[a],&q);
            } else { /* singleton */
              singleXty(xwx + i * pt[a],temp,Xi,X+off[ts[a]],m+ts[a],p+ts[a],k + (ptrdiff_t)*n * (q + ks[ts[a]]),n,&add);
            }
            add = 1; /* for q>0 accumulate result */
          }  
        } /* loop over columns of Xb */
      }
      /* so now xwx contains pt[a] by pt[b] matrix Xa'WXb */
    } /* end parallel section */

    /* if Xb is tensor, may need to apply constraint */
    if (dt[a]>1&&qc[a]>0) { /* first term is a tensor with a constraint */
      x0=x1=xwx; /* pointers to columns of xwx */ 
      /* col by col form (I-vv')xwx, dropping first row... */
      for (j=0;j<pt[b];j++) {
        for (x=0.0,p0=x0,p1=x0+pt[a],p2=v+voff[a];p0<p1;p0++,p2++) x += *p0 * *p2;
        for (p2=v+voff[a]+1,p1=x0+pt[a],x0++;x0<p1;x1++,x0++,p2++) *x1 = *x0 - *p2 *x;
      }
      pa = pt[a] -1;
    } else pa = pt[a];
    if (dt[b]>1&&qc[b]>0) { /* second term is a tensor with a constraint */
      /* copy xwx to xwx0 */
      for (p0=xwx,p1=p0 + pt[b] * (ptrdiff_t) pa,p2=xwx0;p0<p1;p0++,p2++) *p2 = *p0;
      /* row by row form xwx(I-vv') dropping first col... */
      for (j=0;j<pa;j++) {
        for (x=0.0,p0=xwx0+j,p1=v+voff[b],p2=p1+pt[b];p1<p2;p0 += pa,p1++) x += *p0 * *p1; 
        x1 = xwx + j;
        for (p0=xwx0+j+pa,p1=v+voff[b],p2=p1+pt[b],p1++;p1<p2;x1+=pa,p1++,p0+=pa) *x1 = *p0 - *p1 *x; 
      }
      pb=pt[b]-1;
    } else pb=pt[b];
    /* copy result into overall XWX*/  
  
    if (pd[r]>pd[c]) { /* xwx = Xr'WXc */
      for (i=0;i<pa;i++) for (j=0;j<pb;j++) 
      XWX[i+tps[r]+(j+tps[c]) * (ptrdiff_t)ptot] = XWX[j+tps[c]+(i+tps[r]) * (ptrdiff_t)ptot] 
                                                 = xwx[i + pa * (ptrdiff_t)j];
    } else { /* Form xwx = Xc'WXr */
      for (i=0;i<pa;i++) for (j=0;j<pb;j++) 
      XWX[i+tps[c]+(j+tps[r])*(ptrdiff_t)ptot] = XWX[j+tps[r]+(i+tps[c])*(ptrdiff_t)ptot] 
                                               = xwx[i + pa * (ptrdiff_t)j];
    }
  } /* end of block loop */
  FREE(start);
  FREE(XiB); FREE(XjB); FREE(tempnB); FREE(pt); FREE(pd); FREE(off);
  FREE(voff); FREE(tps); FREE(tempB); FREE(xwx); FREE(xwx0);
} /* XWX */


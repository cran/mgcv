/* (c) Simon N. Wood (2019) Released under GPL2 
   Sparse matrix utility routines. 

   In particular provides sparse discrete method routines:
   sXWXd, sdiagXWXt, sXbd & sXyd as well as sparse row tensor product
   routine stmm.
   
   Note that there is little attempt to hand optimize the code at the basic
   avoiding-look-up level. Also sdiagXWXt is not optimized and is single threaded 
   at present.
  
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mgcv.h"
#include <R.h>

/* A sparse matrix structure useful for storing discrete sparse marginals */
typedef struct {
  int m,c, // matrix is m by c
    n, // full matrix has n rows
    nk, // number of index vectors 
    *p, // column j starts at element p[j] of i and x. p[c] is total NZ count.
    *i, // record of non-zero rows
    *k, // index n-vectors ith full row in row k[i] of compact
    *r, // reverse index n-vectors  
    *off, // m-vectors of offsets in reverse n-vector
    nzmax; // entries allocated in i an x (could be more than p[c])
  double *x; // record of non-zero entries
} spMat;

typedef struct {
  int rb,cb,r,c,r_start,c_start,r_size,c_size;
} XWXblock;



SEXP getListEl(SEXP list, const char *str) {
/* get the list element named str, or return NULL 
   - From `Writing R extensions' 5.9.6 */
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    for (int i = 0; i < length(list); i++)
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
        }
    return elmt;
} /* getListElement */

/****************************************************************************** 
 tedious infrastructure routines for handling basic sparse matrix operations... 
*******************************************************************************/

void tri_to_cs(int *Ti,int *Tj,double *Tx,int *Cp,int *Ci,double *Cx,int *w, int nz, int c) {
/* Convert a sparse c column matrix represented in triplet form by Ti, Tj, Tx to 
   compressed column format. nz is the number of non-zero entries.
   Ti, Tj, Tx, Ci and Cx are nz vectors. 
   Cp is a c+1 vectors. w is a c vector it should be zero on entry and will be zeroed on exit.
   Recall that Cx[Cp[j]:(Cp[j+1]-1)] contains the non-zero elements of column j, which
   are rows Ci[Cp[j]:(Cp[j+1]-1)] of the full matrix.
   Algorithm is essentially the one given in Davis (2006) 2.4.  
*/
  int i,j;
  
  for (i=0;i<nz;i++) w[Tj[i]]++; /* count the entries per column */
  /* form cumsum of w to get the indices of the column starts in Ci and Cx... */ 
  for (j=0,i=0;i<c;i++) {
    Cp[i] = j;
    j+=w[i];
    w[i]=Cp[i];
  }
  Cp[c] = j;
  /* now w can be re-used to index where each column's next row to accumulate should go */
  for (i=0;i<nz;i++) {
    j = w[Tj[i]]++; /* position of column Tj[i]'s current entry in Cx/Ci */ 
    Ci[j] = Ti[i];
    Cx[j] = Tx[i];
  }
  for (i=0;i<c;i++) w[i] = 0;  
} /* tri_to_cs */

void cs_trans(int *Ap,int *Ai,double *Ax,int *Cp,int *Ci, double *Cx,int *w,int r,int c) {
/* computes C=A^T where A and C are compressed column sparse matrices. A is r by c.
   w is an r-vector. 
   Note that the columns of C are sorted, irrespective of the ordering of the cols of A.
   Recall that Cx[Cp[j]:(Cp[j+1]-1)] contains the non-zero elements of column j, which
   are rows Ci[Cp[j]:(Cp[j+1]-1)] of the full matrix.
   Algorithm is essentially the one given in Davis (2006) 2.5.   */
  int i,j,q;
  for (i=0;i<r;i++) w[i] = 0;
  for (i = 0 ; i<Ap[c];i++) w[Ai[i]]++; /* counting the entries in A's rows */
  /* form cumsum of w to get the indices of the column starts in Ci and Cx... */ 
  for (j=0,i=0;i<r;i++) {Cp[i] = j;j+=w[i];w[i]=Cp[i];} Cp[r] = j;
  /* now w[j] used to keep track of entries added to col j of C, basically w is an array
     of indices w[j] indicating where the next entry of col j is stored in Ci/Cx.  */
  for (j=0;j<c;j++) { /* work across the columns of A */
    for (i=Ap[j];i<Ap[j+1];i++) { /* work through NZ rows of A[,j] */
      q = w[Ai[i]]++; /* row j of col i of C goes here in Ci, Cx */
      Ci[q] = j;
      Cx[q] = Ax[i];
    }  
  }  
} /* cs_trans */  


void sprealloc(spMat *A, int nzmax) {
/* Alter the storage in A to accomodate at most nzmax non-zero entries. 
   Note that if using openMP then allocation of storage needs to be locked.
*/
#pragma omp critical
  { A->i = (int *)REALLOC(A->i,sizeof(int) * nzmax);
    A->x = (double *)REALLOC(A->x,sizeof(double) * nzmax);
    A->nzmax = nzmax;
  }  
} /* sprealloc */  

void spalloc(spMat *A,int c,int nzmax) {
  A->p = (int *)CALLOC((size_t)c+1,sizeof(int));
  A->c = c;
  A->nzmax = nzmax;
  A->i = (int *)CALLOC((size_t)nzmax,sizeof(int));
  A->x = (double *)CALLOC((size_t)nzmax,sizeof(double));
} /* spalloc */

void spfree(spMat *A,int m) {
/* clear m sparse matrices in A */
  spMat *a;
  for (a=A;a<A+m;a++) {
    FREE(a->p);FREE(a->i);FREE(a->x);
  }  
} /* spfree */  


void cs_mult(spMat *A,spMat *B,spMat *C,int *w, double *x,int init) {
/* Computes C=AB where A,B and C are compressed column sparse matrices using the 
   algorithm of section 2.8 of Davis (2006). Note that this approach involves *no*
   redundant multiplications by (or checking of) zeroes, which appear unvoidable 
   in a direct computation of A'B. For this reason A'B should be performed by first
   explicitly transposing A.
   If init!=0 then this routine will 
   allocate or reallocate storage in C->i and C->x.  C should have p pre-allocated 
   (to a B->c+1 vector), and i and x set to NULL with nzmax = 0, or i, x and nzmax initialized.  
   This routine only accesses the m,c,i,p and x elements of its sparse matrix arguments, and 
   also the nzmax argument of C.
   w and x are A->m vectors.   
 */
  int i,j,k,r,c,nz=0,*Cp,*Ci,*Bp,*Bi,*Ap,*Ai,ii,jj;
  double *Cx,*Bx,*Ax,b; 
  C->c = c = B->c; /* columns of B and result C */
  C->m = r = A->m; /* rows of A and C */
  Cp = C->p;Ci = C->i;Cx = C->x;
  Bp = B->p;Bi = B->i;Bx = B->x;
  Ap = A->p;Ai = A->i;Ax = A->x;
  /* w[i] is set to j the first time that there is a row i entry in column j.
     If w[i] <j there is not yet an entry at row i.
  */
  for (i=0;i<r;i++) w[i] = -1; 
  for (j=0;j<c;j++) { /* work through columns of B and C */
    if (init && nz+r > C->nzmax) {
      sprealloc(C,C->nzmax*2+r); /* expand storage */
      Ci = C->i;Cx = C->x;
    }
    Cp[j]=nz; /* location of jth column start in Ci/Cx */
    for (i=Bp[j];i<Bp[j+1];i++) { /* work through B[:,j]'s NZ elements */
      b = Bx[i]; /* B[jj,j] to multiply A[:,jj] */
      jj = Bi[i]; /* next non-zero row in B[:,j] */
      for (k=Ap[jj];k<Ap[jj+1];k++) { /* work through A[:,jj]'s NZ elements */
        ii = Ai[k]; /* A[ii,jj] is next NZ */
	if (w[ii]<j) { /* first component of row ii in column C[:,j] */
	  w[ii] = j; /* signal that column C[:,j] has an entry at row ii */ 
	  Ci[nz] = ii;nz++; /* record row in C[:,j] sparse storage */
	  x[ii] = b * Ax[k]; /* C[:,j] accumulated to dense (but indexed) storage in x */
	} else x[ii] += b * Ax[k]; 
      }  
    } /* done with B[:,j] */
    for (i=Cp[j];i<nz;i++) Cx[i] = x[Ci[i]]; /* read C[:,j] from dense storage to sparse */
  }
  Cp[c] = nz;
  /* shrink back to exact size needed if init==1, but not otherwise */
  if (init==1 && nz!=C->nzmax) { if (!nz) nz=1; sprealloc(C,nz); C->nzmax=nz;}
} /* cs_mult */


int sum_dup(int *Cp,int *Ci,double *Cx,int *w, int r, int c) {
/* Sum the duplicate entries in a compressed column matrix.
   Cx and Ci are Cp[c+1] vectors. Cp is a c+1 vector. w is an r vector where r is the 
   number of rows of C (it will be zeroed on exit)
   Recall that Cx[Cp[j]:(Cp[j+1]-1)] contains the non-zero elements of column j, which
   are rows Ci[Cp[j]:(Cp[j+1]-1)] of the full matrix.
   Algorithm is essentially the one given in Davis (2006) 2.6.
*/
  int i,j,k,nz=0,q,Cpj=0;
  /* w is used to signal whether a row has been seen yet in this column. If row i
     has been seen already then w[i] records where the original entry is in the 
     output Cx/Ci. If it has not been seen already then w[i] will be before the 
     start index for the current output column. Next line initializes... */ 
  for (i=0;i<r;i++) w[i] = -1;
  for (j=0;j<c;j++) { /* work over the columns */
    q=nz; /* start of current column in output Ci/Cx */
    for (k=Cpj;k<Cp[j+1];k++) { /* work down NZ entries of input col j */
      i = Ci[k];     /* entry at i,j is NZ */
      if (w[i]>=q) { /* entry is a duplicate */
        Cx[w[i]] += Cx[k]; /* so just add it to original */ 
      }	else { /* new entry */
        w[i] = nz; /* record where row i of col j is located in Ci/Cx */
	Ci[nz] = i;
	Cx[nz] = Cx[k];
	nz++;
      }	
    } /* NZ of col j loop */
    Cpj = Cp[j+1];
    Cp[j+1] = nz; /* update the output column record */
  } /* column loop */
  for (i=0;i<r;i++) w[i] = 0;
  return(nz);  /* retun the total number of non-zeroes */
}  /* sum_dup */

void cs_accumulate(spMat *A,spMat *B,int *iwork) {
/* computes A += B by first appending elements of B to A, and 
   then summing the duplicate entries in A. 
   iwork is an A->m vector. */
  int nzmax,i,j,k,*Ap,*Bp,*Ai,*Bi;
  double *Bx,*Ax;
  Ap = A->p;Bp = B->p;
  Ax = A->x;Bx = B->x;
  Ai = A->i;Bi = B->i;
  nzmax = Ap[A->c] + Bp[B->c]; /* maximum storage required */
  if (A->nzmax < nzmax) sprealloc(A,nzmax); /* need to expand storage */
  k = nzmax-1;
  for (j=A->c;j>0;j--) { /* work down columns */
    for (i=Bp[j]-1;i >= Bp[j-1];i--) {
      Ax[k] = Bx[i];
      Ai[k] = Bi[i];
      k--;
    }
    for (i=Ap[j]-1;i >= Ap[j-1];i--) {
      Ax[k] = Ax[i];
      Ai[k] = Ai[i];
      k--;
    }
    Ap[j] = nzmax; /* new start of column j */
    nzmax = k + 1;
  }
  /* now sum the duplicates */
  nzmax = sum_dup(Ap,Ai,Ax,iwork,A->m,A->c);
} /* cs_accumulate */

/**********************************************************************
  Sparse-dense product routines...
***********************************************************************/

void spMtv(spMat *M,double *v,double *u,int add) {
/* Forms u = M'v where M is a sparse matrix and v a dense M.m vector */
  int i,j,Mc,*Mp,*Mi;double *Mx;
  Mc = M->c;Mi = M->i;Mp = M->p;Mx = M->x;
  if (!add) for (j=0;j<Mc;j++) u[j] = 0.0;
  for (j=0;j<Mc;j++,u++) { /* cols of M, rows of u */
    for (i=Mp[j];i<Mp[j+1];i++) *u += v[Mi[i]] * Mx[i];
  }  
} /* spMtv */ 

void spMtA(spMat *M,double *A,double *B,int c,int add) {
/* Forms B = M'A where M is a sparse matrix and A a dense M.m by c matrix 
   B is obviously M.c by c. 
*/
  int i,j,k,Mm,Mc,*Mi,*Mp;
  double *Mx;
  Mm = M->m;Mc = M->c;Mi = M->i;Mp = M->p;Mx = M->x;
  if (!add) for (j=0;j < Mm * c ;j++) B[j] = 0.0;
  for (j=0;j < Mc;j++) { /* cols of M */
    for (i=Mp[j];i<Mp[j+1];i++) for (k=0;k<c;k++) B[j + k* Mc] += A[Mi[i]+Mm*k] * Mx[i];
  }  
} /* spMv */ 

void spMv(spMat *M,double *v,double *u) {
/* Forms u = Mv where M is a sparse matrix and v a dense M.c vector */
  int i,j,Mm,Mc,*Mi,*Mp;
  double *Mx;
  Mm = M->m;Mc = M->c;Mi = M->i;Mp = M->p;Mx = M->x;
  for (j=0;j<Mm;j++) u[j] = 0.0;
  for (j=0;j<Mc;j++) { /* cols of M */
    for (i=Mp[j];i<Mp[j+1];i++) u[Mi[i]] += v[j] * Mx[i];
  }  
} /* spMv */ 

void spMA(spMat *M,double *A,double *B,int c) {
/* Forms B = MA where M is a sparse matrix and A a dense M.c by c matrix 
   B is obviously M.m by c. 
*/
  int i,j,k,Mm,Mc,*Mi,*Mp;
  double *Mx;
  Mm = M->m;Mc = M->c;Mi = M->i;Mp = M->p;Mx = M->x;
  for (j=0;j < Mm * c ;j++) B[j] = 0.0;
  for (j=0;j < Mc;j++) { /* cols of M */
    for (i=Mp[j];i<Mp[j+1];i++) for (k=0;k<c;k++) B[Mi[i] + k* Mm] += A[j+Mc*k] * Mx[i];
  }  
} /* spMv */ 


/*****************************************************
  routines to deal with constraints, which involves 
    sparse -> dense -> transform -> sparse 
******************************************************/

void sp_to_dense_insitu(spMat *A,int r){
/* Sparse matrix A has had its x component realloced to a accomodate at least an 
   r by A->c dense matrix. This routine moves the existing elements of A->x into 
   the position required in a dense representation, assuming that the original A
   is in the first block within the re-sized A.
   Assumes that columns are sorted into row order.
*/
  int i,j,Ac,*Ap,*Ai;
  double x,*Ax;
  Ac = A->c;Ai=A->i;Ap=A->p;Ax=A->x;
  for (j = Ac-1;j >= 0;j--) /* down the columns */
    for (i = Ap[j+1]-1;i>=Ap[j];i--) { /* and from last to first nz row of each column */ 
      /* 3 steps (or other check) needed to avoid scrubbing source location coinciding with destination! */   
      x = Ax[i]; /* copy */
      Ax[i] = 0.0; /* scrub from old location */
      Ax[Ai[i] + j * r] = x; /* move to new location */
  }
  Ap[0] = -1; /* signal a dense matrix packed in A->x */
} /* sp_to_dense_insitu */

void sp_to_dense(spMat *A,double *x,int r,int c,int nr) {
/* copy sparse matrix A to a block of dense matrix stored in x starting at r,c. 
   x has nr rows in total.
*/
  int i,j,Ac,*Ap,*Ai,jc;double *Ax;
  Ac = A->c;Ai=A->i;Ap=A->p;Ax=A->x;
  for (j=0;j<Ac;j++) {
    jc = j+c;
    for (i=Ap[j];i<Ap[j+1];i++) x[r + Ai[i] + jc * nr] = Ax[i];
  }
} /* sp_to_dense */

void dense_to_sp(spMat *A) {
/* Represent a dense A->m by A->c matrix stored in A->x in sparse form, by filling in 
   A->i and A->p appropriately. 
*/
  int i,j,*Ai,Am;
  Am = A->m;
  A->i = REALLOC(A->i,sizeof(int)*Am * A->c);
  A->p = REALLOC(A->p,sizeof(int)*(A->c+1));
  Ai = A->i;
  for (j=0;j < A->c;j++) {
    A->p[j] = Am*j;
    for (i=0;i<Am;i++,Ai++) *Ai = i;
  }
  A->p[A->c] = Am * A->c;
} /* dense_to_sp */

void left_con_vec(double *x,double *v,double *y,int n,int rev) {
/* Forms y = (I-vv')x and drops first row of y if rev=0.
   Otherwise forms y=(I-vv')[0,x']'  
   n is vector dimension
*/
  double vx=0.0;
  int i;
  if (rev) rev=1;
  for (i=rev;i<n;i++) vx += v[i]*x[i-rev];
  if (rev) {
    y[0] = -v[0]*vx;
    for (i=1;i<n;i++) y[i] = x[i-1] - v[i]*vx; 
  } else for (i=1;i<n;i++) y[i-1] = x[i] - v[i]*vx;
} /* left_con_vec */

void left_con(spMat *A,double *v,double *work) {
/* Forms (I-vv')A and drops first row.
   A is an A->m by A->c dense matrix stored in A->x
   work is length A->c and v is length A->m 
*/
  char trans='T';
  int one=1,Am,i,j,k;
  double done=1.0,dzero=0.0,*Ax,x,*Ax1;
  Ax = A->x;Am = A->m;
  /* compute v'A... */
  F77_CALL(dgemv)(&trans, &(A->m), &(A->c),&done,Ax,&Am,v,&one,&dzero,work,&one FCONE);
  /* compute A - vv'A */
  for (k=j=0;j<A->c;j++) for (x=work[j],i=0;i<Am;i++,k++) Ax[k] -= v[i] * x;
  /* drop first row of A */
  for (Ax1=Ax,j=0;j<A->c;j++) {
    Ax++; /* move pointer to second row of column */
    for (i=1;i<Am;i++,Ax++,Ax1++) *Ax1 = *Ax;
  }
  A->m--;
} /* left_con */

void right_con(spMat *A,double *v,double *work) {
/* Forms A(I-vv') and drops first col.
   A is an A->m by A->c dense matrix stored in A->x
   work is length A-> m and v is length A->c 
*/
  char trans='N';
  int one=1,Am,i,j,k;
  double done=1.0,dzero=0.0,*Ax,x,*Ax1,*Ax2;
  Ax = A->x;Am = A->m;
  /* compute Av... */
  F77_CALL(dgemv)(&trans, &(A->m), &(A->c),&done,Ax,&Am,v,&one,&dzero,work,&one FCONE);
  /* compute A - Avv' */
  for (k=j=0;j<A->c;j++) for (x=v[j],i=0;i<Am;i++,k++) Ax[k] -= work[i] * x;
  /* drop the first column of A */
  Ax1 = Ax; Ax += Am; Ax2 = Ax1 + (A->c-1)*A->m;
  for (;Ax1<Ax2;Ax1++,Ax++) *Ax1 = *Ax;
  A->c--;
} /* right_con */

/********************************************************
 The cross product routines. By far the most involved 
 algorithms. sXWXdij does the main algorithmic work, and
 sXWXd is in some sense a driver - but it still has to 
 handle constraints, a large quantity of book keeping and
 parallelization
 ********************************************************/

void sXWXdij(double *w, double *d,spMat *Xs,spMat *Xt,int rb,int cb,int r,int c,int *dn,int *ts,
	     int *dt,spMat *W,spMat *V,spMat *XWX,int *init,int *iwork,double *xwork) {
/*
   Forms the r, c sub-block of the rb, cb term block of X'WX, where W = diag(w), and the 
   sparse marginal matrices defining the terms are in sparse matrix array Xs.
   The ith term starts at Xs[ts[i]] and has dt[i] component marginal model matrices.   
   dn and d are zero on entry and exit and are both n vectors. Xt stores the transposes 
   of matrices in Xs, to avoid redundant transposition here (this matters, since many 
   blocks can be empty, in which case transposition cost completely dominate all other costs).  

   W and V should have i and x initialized to at least n vectors on entry and nzmax set 
   accordingly. W->p and V->p should be m- vectors where m is the maximum of n, p and the 
   largest nz count in any element of Xs.

   iwork and xwork should be max(n,p)-vectors. 

   if *init!=0 then W and XWX will have their x and i components expanded as needed, and 
   on exit *init with contain the largest amount of NZ storage required by either returned.
   When used with *init==0 then the storage in W and XWX should be at least as large as the 
   largest such number required.  

   NOTE: *init=0 is not fully tested, as currently unused. 
*/
  double x,*Wx,*Xsx;
  int q,qr,qc,qt,l,i,k,ri,ci,ii,jj,k1,j,j1,first_mat,
    first_col=0,first_s=0,*kr,*kc,n,Wnz,*Wi,*Wj,*Xsi,*Xsr,*Xsoff,s,t,r0,qr0,qc0,c0;
  
  n = Xs[0].n;
  ri = ts[rb]; /* element of Xs at which term rb starts */
  if (dt[rb]==1) qr=0; else 
  for (qr=1,i = ri;i < ri+dt[rb]-1;i++) qr *= Xs[i].c; /* pre-final column count for term rb */
  r0=r;qr0=qr;
  ci = ts[cb]; /* element of Xs at which term cb starts */
  if (dt[cb]==1) qc=0; else 
  for (qc=1,i = ci;i < ci+dt[cb]-1;i++) qc *= Xs[i].c; /* pre-final column count for term cb */
  c0=c;qc0=qc;
  for (s=0;s<Xs[ri].nk;s++) for (t=0;t<Xs[ci].nk;t++) {
    kr = Xs[ri+dt[rb]-1].k + n*s; /* index vector for final marginal of term rb */
    kc = Xs[ci+dt[cb]-1].k + n*t; /* index vector for final marginal of term cb */
    r=r0;qr=qr0;c=c0;qc=qc0;Wnz=0;first_mat = -1;
    qt = dt[cb]+dt[rb] - 2; /* number of columns contributing to the diagonal */
    q=0; /* total pre-column count and counter */
    /* accumulate the diagonal for the rb term */
    Wi = W->i;Wx = W->x;Wj = W->p;
    for (i = ri;i < ri+dt[rb]-1;i++) { /* loop over columns before final */
      qr /= Xs[i].c; /* update q */
      l = r/qr; /* column of current (ith) marginal */
      r = r%qr;
      Xsi = Xs[i].i;Xsx = Xs[i].x;Xsr = Xs[i].r+n*s;Xsoff = Xs[i].off+Xs[i].m*s+s;
      if (first_mat<0) {first_mat=i;first_col=l;first_s=s;} /* save for later clearing of dn and d */
      k1 = Xs[i].p[l+1];
      for (k=Xs[i].p[l];k<k1;k++) { /* loop over NZ elements of this column */
        j = Xsi[k]; /* non-zero row of col l of current marginal */
        x = Xsx[k]; /* the corresponding element */
        j1 = Xsoff[j+1];
        for (jj=Xsoff[j];jj<j1;jj++) { /* loop over full rows corresponding to row j */
          ii = Xsr[jj]; /* current full row */
	  if (q == dn[ii]) { /* condition for product to be non-zero */ 
	    if (q==0) d[ii] = x * w[ii]; else d[ii] *= x;
	    dn[ii]++;
	    if (q==qt-1) { /* no more columns to contribute to d accumulate entry to \bar W*/
              Wi[Wnz] = kr[ii];Wj[Wnz] = kc[ii];Wx[Wnz]=d[ii];Wnz++;
	    }  /* end of sparse \bar W accumulation */
	  } /* if (i == dn[ii]) */  
        }	 
      }
      q++;
    } /* rb term pre-final column loop */
    /* accumulate the diagonal for the cb term */
    for (i = ci;i < ci+ dt[cb]-1;i++) { /* loop over columns before final */
      qc /= Xs[i].c; /* update q */
      l = c/qc; /* column of current (ith) marginal */
      c = c%qc;
      Xsi = Xs[i].i;Xsx = Xs[i].x;Xsr = Xs[i].r+n*t;Xsoff = Xs[i].off+Xs[i].m*t+t;
      if (first_mat<0) {first_mat=i;first_col=l;first_s=t;}/* save for later clearing of dn and d */
      k1 = Xs[i].p[l+1];
      for (k=Xs[i].p[l];k<k1;k++) { /* loop over NZ elements of this column */
        j = Xsi[k]; /* non-zero row of col l of current marginal */
        x = Xsx[k]; /* the corresponding element */
        j1 = Xsoff[j+1];
        for (jj=Xsoff[j];jj<j1;jj++) { /* loop over full rows corresponding to row j */
          ii = Xsr[jj]; /* current full row */
	  if (q == dn[ii]) { /* condition for product to be non-zero */ 
	    if (q==0) d[ii] = x * w[ii]; else d[ii] *= x;
	    dn[ii]++;
	    if (q==qt-1) { /* no more columns to contribute to d - accumulate to \bar W */
	      Wi[Wnz] = kr[ii];
	      Wj[Wnz] = kc[ii];
	      Wx[Wnz]=d[ii];
	      Wnz++;
	    }  
	  } /* if (i == dn[ii]) */  
        }	 
      }
      q++; /* update total column counter */
    } /* cb term pre-final column loop */
  
    /* now clear dn and d back to zero */
    if (qt) {
      Xsi = Xs[first_mat].i;Xsr = Xs[first_mat].r+first_s*n;
      Xsoff = Xs[first_mat].off+first_s*(Xs[first_mat].m+1);
      k1 = Xs[first_mat].p[first_col+1];
      for (k=Xs[first_mat].p[first_col];k<k1;k++) { /* loop over NZ elements of this column */
        j = Xsi[k]; /* non-zero row of col l of current marginal */
        j1 = Xsoff[j+1];
        for (jj=Xsoff[j];jj<j1;jj++) { /* loop over full rows of discrete row j */
          ii = Xsr[jj]; /* current full row */
	  d[ii] = 0.0;dn[ii] = 0;
        }	 
      } 
    } else { /* both terms are non-tensor, and so \bar W has to be directly accumulated */
      for (ii=0;ii<n;ii++) {
        Wi[Wnz] = kr[ii];Wj[Wnz] = kc[ii];Wx[Wnz]=w[ii];Wnz++;
      }   
    }
    /* at this stage there are duplicate entries in the sparse matrix stored in triplet form in 
       Wi, Wj, Wx (each with Wnz elements). The duplicates need to be summed. This can be efficiently 
       accomplished by first converting to compressed column format, and then stripping the duplicates
       from each column. See Davis (2016) sections 2.4 & 2.6 for basic algorithms for this. */
    c = Xs[ts[cb]+dt[cb]-1].m; r = Xs[ts[rb]+dt[rb]-1].m; /* \bar W is r by c */ 
    tri_to_cs(Wi,Wj,Wx,V->p,V->i,V->x,dn,Wnz,c); /* note requirement that dn zero on entry */
    Wnz = sum_dup(V->p,V->i,V->x,dn,r,c); /* sum the duplicates */
    V->m = r;V->c=c;
    /* V now contains \bar W in column compressed format.
     The actual crossproduct can now be formed. 
     A substantial problem for parallelization is that the NZ size of the crossproduct is not 
     known in advance (it would be possible to find the NZ count for \bar W B without forming 
     sparse matrix storage, but not A'\bar W B). One way to handle this is to have a switch 
     that allows the routine to expand memory as needed (reporting the memory requirement), but 
     protect the expansion using (in the obvious way)... 
     omp_lock_t lck;
     omp_init_lock(&lck); 
     omp_set_lock(&lck); 
     omp_unset_lock(&lck); 
     - Actually it's simpler to use the approach now in sprealloc, above.
     When the switch is 0 then the memory supplied is assumed to be adequate. Hence memory can 
     be allocated at first call, and assumed correct at subsequent calls.  
    */
   
    W->m = r;W->c= Xs[ts[cb]+dt[cb]-1].c;
    /* iwork and xwork are V->m vectors ... */
    j  = *init; if (j==1) j = 2; // signal that W is not to be shrunk, only expanded (matters above!)...
    cs_mult(V,Xs+ts[cb]+dt[cb]-1,W,iwork, xwork,j); /* W=\bar W B where B is final marginal of term cb */
   
    /* V->i and V->x must be at least as large as A->i and A->x, but will be if it meets requirements given
       at start... */
  
    if (!s && !t) cs_mult(Xt+ts[rb]+dt[rb]-1,W,XWX,iwork, xwork,*init); else {
      V->m = Xs[ts[rb]+dt[rb]-1].c;V->c = W->c;
      cs_mult(Xt+ts[rb]+dt[rb]-1,W,V,iwork, xwork,j);
      cs_accumulate(XWX,V,iwork);
    }  
    /* return the largest non-zero memory allocation required in init */
    if (*init) { if (XWX->nzmax > W-> nzmax) *init = XWX->nzmax; else *init = W->nzmax;}
  } /* s,t summation loop */
  /* there is no guarantee that XWX has column indices sorted. Double transpose of (XWX) 
     will sort this out... */
  k = XWX->p[XWX->c]; /* NZ elements in XWX */
  if (k>1) { /* more than one entry so maybe something to do */
    if (V->nzmax<k) sprealloc(V,k);
    cs_trans(XWX->p,XWX->i,XWX->x,V->p,V->i,V->x,iwork,XWX->m,XWX->c);
    V->m = XWX->c;V->c = XWX->m;
    cs_trans(V->p,V->i,V->x,XWX->p,XWX->i,XWX->x,iwork,V->m,V->c); 
  }
  /* so XWX now contains the required result */
} /* sXWXdij */


SEXP sXWXd(SEXP X,SEXP W,SEXP LT, SEXP RT,SEXP NT) {
/* X is a list defining the model matrix. Its elements are
   Xd, kd, ks, v, ts, dt, qc, off and r.
   See ?XWXd for definitions of these and other arguments except:
   elements of list Xd are sparse matrices, 
   r[off[[i]][j]:(off[[i]][j+1]-1),i] indexes the full row of ith submatrix 
   corresponding to the jth row stored in Xd[[i]]   

   Note that kd, ks, off and r should be zero based on entry and 
   off should be unlisted before entry.

   W is the diagonal vector in X'WX

   LT and RT index which smooth terms to include from left and right.

   NT is the number of threads to use.

   NOTE: it is assumed that only tensor product terms are constrained using 
   the qc/v mechanism since it is always more efficient to constrain singletons
   before calling this routine.

*/
  spMat *Xs,*Xt,*xwx,*V1,*W1,*Mp;
  int mx,i,j,k,b,i1,j1,*dim,n,*kd,*ks,*r,*off_start,*off,nt,*ts,*dt,*qc,*lt,*rt,nlt,nrt,nb,nr,nc,nrc,ncc,
    brs0,brs1,bcs0,bcs1=0,is,js,ic,jc,is0,ic0,js0,jc0,init=1,p,ii,jj,*str,*stc,symmetric,*sub_blocks,*block_size,
    rcum=0,ccum=0,is1,rcum0,ccum0,*dn,*iwork,*XWXi,*XWXp,nzmax=0,rb,cb,*ip,n_threads,tid,*B,nprot=0;
  SEXP Xd,M, i_sym,x_sym,dim_sym,p_sym,ul_sym,KD,R,KS,OFF,TS,DT,QC,V,XWX,OFFS;
  double **v,*w,*d,*xwork,*XWXx,*xp,*xp1,*cost;
  XWXblock *block,*blp;
  /* register the names of the slots in the sparse matrices */
  p_sym = install("p");
  dim_sym = install("Dim");
  i_sym = install("i");
  x_sym = install("x");
  ul_sym = install("uplo");
  Xd = getListEl(X,"Xd"); /* the sparse matrix list */
  KD = getListEl(X,"kd"); /* the matrix of index vectors */
  n = nrows(KD); /* number of data */
  /* get the row indices (forward then reverse) */
  KD = PROTECT(coerceVector(KD,INTSXP));nprot++;
  kd = INTEGER(KD);
  R = getListEl(X,"r");
  R = PROTECT(coerceVector(R,INTSXP));nprot++;
  r = INTEGER(R);
  OFF =  getListEl(X,"off"); /* the offset list for the reverse indices */
  OFF = PROTECT(coerceVector(OFF,INTSXP));nprot++;
  off = INTEGER(OFF);
  OFFS =  getListEl(X,"offstart"); /* the start points in the offset array */
  OFFS = PROTECT(coerceVector(OFFS,INTSXP));nprot++;
  off_start = INTEGER(OFFS);
  /* get the matrix defining the range of k vectors for each matrix */
  KS = getListEl(X,"ks");
  KS = PROTECT(coerceVector(KS,INTSXP));nprot++;
  ks = INTEGER(KS);
  
  mx = length(Xd); /* list length */
  /* read sparse matrix list into structure suitable for passing on to 
     other routines */
  Xs = (spMat *) CALLOC((size_t)mx,sizeof(spMat)); // sparse matrix list
  Xt = (spMat *) CALLOC((size_t)mx,sizeof(spMat)); // chached sparse matrix transpose list
  iwork = (int *)CALLOC((size_t)n,sizeof(int));
  for (i=0;i<mx;i++) { // work through the sparse matrices
    M = VECTOR_ELT(Xd, i); // the ith sparse matrix
    Xs[i].x = REAL(R_do_slot(M,x_sym)); // the elements of ith matrix
    Xs[i].p = INTEGER(R_do_slot(M,p_sym)); // .p[j] is where jth col starts in .i and .x
    Xs[i].i = INTEGER(R_do_slot(M,i_sym)); // row index for each non-zero element
    dim = INTEGER(R_do_slot(M,dim_sym));
    Xs[i].m = dim[0];Xs[i].c = dim[1]; // matrix is .m by .c
    Xs[i].k = kd + n * ks[i]; // the .nk index n vectors (end to end)
    Xs[i].r = r + n * ks[i];  // ... equivalent reverse indices
    Xs[i].n = n; // rows in full matrix
    Xs[i].nk = ks[mx+i] - ks[i]; // number of index vectors for this term
    Xs[i].off = off + off_start[ks[i]]; // the .nk offset .m+1 vectors for .r 
    j = Xs[i].p[Xs[i].c]; /* number of NZs in this matrix */
    if (j>nzmax) nzmax=j; /* largest amount of NZ storage required by any element of Xs */
    /* now Cache its transpose... */
    spalloc(Xt+i,Xs[i].m,j);
    cs_trans(Xs[i].p,Xs[i].i,Xs[i].x,Xt[i].p,Xt[i].i,Xt[i].x,iwork,Xs[i].m,Xs[i].c);
    Xt[i].m = Xs[i].c; Xt[i].c = Xs[i].m;
  } /* sparse matrix array filling loop */    
  FREE(iwork);
  /* now deal with the smooth term information... */
  TS = getListEl(X,"ts");
  nt = length(TS); /* number of smooth terms */
  TS = PROTECT(coerceVector(TS,INTSXP));nprot++;
  ts = INTEGER(TS); /* term starts in matrix array */
  DT = getListEl(X,"dt");
  DT = PROTECT(coerceVector(DT,INTSXP));nprot++;
  dt= INTEGER(DT); /* number of marginal matrices for this term. */
  QC = getListEl(X,"qc");
  QC = PROTECT(coerceVector(QC,INTSXP));nprot++;
  qc= INTEGER(QC); /* the constraint indicator nt-vector */
  V = getListEl(X,"v");
  v = (double **) CALLOC((size_t)nt,sizeof(double *));
  for (i=0;i<nt;i++) if (qc[i]) {
      v[i] = REAL(VECTOR_ELT(V, i));
  }    
  nlt = length(LT);lt = INTEGER(LT);
  nrt = length(RT);rt = INTEGER(RT);
  n_threads = asInteger(NT);
  if (n_threads<1) n_threads = 1;
  w = REAL(W);
  /* is result symmetric? */
  if (nlt==nrt) { for (symmetric=1,i=0;i<nlt;i++) if (lt[i]!=rt[i]) { symmetric=0;break;}} else symmetric = 0; 
  /* compute the number of sub-blocks, and their sizes for each smooth term
     If term i is involved in a transposed model matrix applied from the 
     left, then it contributes sub_blocks[i] row blocks each with block_size[i] rows.
     Applied untransposed from the right it contributes sub_blocks[i] col blocks 
     each of block_size[i] cols. 
  */
  sub_blocks = (int *)CALLOC((size_t)nt,sizeof(int)); /* number of sub-blocks for this term */
  block_size = (int *)CALLOC((size_t)nt,sizeof(int)); 
  
  for (p=0,i=0;i<nt;i++) { /* term loop (all, not just lt,rt selected) */
    sub_blocks[i] = 1;
    if (dt[i]>1) for (j=ts[i];j<ts[i]+dt[i]-1;j++) sub_blocks[i] *= Xs[j].c;
    block_size[i] = Xs[ts[i]+dt[i]-1].c; // in use may require adjust by -1 if there is a constraint 
    p += block_size[i] * sub_blocks[i]; // total parameter count
  }
  /* adjust maximum NZ storage upwards if needed... */
  if (nzmax<n) nzmax = n;
  if (nzmax<p) nzmax=p;
  /* now work out the total number of actual blocks */
  for (j=0,i=0;i<nlt;i++) j += sub_blocks[lt[i]];
  if (symmetric) { k=j;nb = j*(j+1)/2;} else { 
    for (k=0,i=0;i<nrt;i++) k += sub_blocks[rt[i]];
    nb = j*k;
  }
  /* allocate storage for the structure matrix. struct[i,j] stores the 
     location in block of the ij'th sub-block of the total requested 
     cross-product */
  nr = j;nc=k; 
  str = (int *)CALLOC((size_t)nr*nc,sizeof(int)); 
  
  /* next create an array storing the sub-block information. i.e.
     rb,cb - the term block row and col.
     r,c - the sub block within the current term block.
     r_size, c_size - the number of rows and cols in block.
     r_start, c_start - the row and column at which block starts in full result
     ... this is the information required to form the block, and part of the information 
     required to then to place it correctly in the final result (the remaining information
     has to wait until the NZP for all blocks is known).  
  */ 
  
  block = (XWXblock *)CALLOC((size_t)nb,sizeof(XWXblock)); // NOTE: free this
  
  brs1=rcum0=0; 
  for (ii=0,blp = block,i=0;i<nlt;i++) { // term block-rows
    if (symmetric) i1 = i+1; else i1 = nrt;
    brs0=brs1; /* start row of current term block */
    for (ccum0=0,bcs0=j=0;j<i1;j++) { // term block-cols
      for (rcum=rcum0,brs1=brs0,is=0;is<sub_blocks[lt[i]];is++,rcum++) { // sub block-rows
        if (symmetric&&i==j) is1 = is+1; else is1 = sub_blocks[rt[j]];
	for (ccum=ccum0,bcs1=bcs0,js=0;js<is1;js++,blp++,ccum++) { // sub block-cols
          blp->rb = lt[i];blp->cb = rt[j];
	  blp->r = is;blp->c = js; /* r,c sub-block within rb,cb block */
	  blp->r_size = block_size[lt[i]]; /* number of rows in sub-block */
	  blp->c_size = block_size[rt[j]]; /* number of cols in sub-block */
	  blp->r_start = brs1; /* at which row of full result does this start? */
	  blp->c_start = bcs1; /* at which col of full result does this start? */
	  bcs1 += block_size[rt[j]];
	  str[rcum + nr * ccum] = ii++; 
	} // sub block-cols
	brs1 += block_size[lt[i]]; /* update start row of current sub block */
      }	// sub block rows
      ccum0 = ccum; // update starting col-block
      bcs0 = bcs1; /* start of next term block - final sub block loop runs through all sub-blocks*/
    } // term block columns
    rcum0 = rcum; // update starting block-row
  }  // block rows
  /* reset rcum and ccum from counting (sub) blocks of output to counting total rows and columns 
     of output (because of lt and rt, these need not be total parameter count)... */
  blp--;
  rcum = blp->r_start+blp->r_size;
  ccum = blp->c_start+blp->c_size;
  /* given the block array, it is easy to re-order the block processing order, by re-ordering the 
     indices 0:(nb-1) */

  /* now process the blocks - the major decision is how to deal with storage allocation
     and de-allocation
  */
  xwx = (spMat *)CALLOC((size_t)nb,sizeof(spMat)); /* sparse sub matrix storage */
  for (b=0;b<nb;b++) spalloc(xwx+b,block[b].c_size,250); // NOTE: is a less arbitrary starting size possible?
  d = (double *)CALLOC((size_t)n*n_threads,sizeof(double));
  dn = (int *)CALLOC((size_t)n*n_threads,sizeof(int));
  if (n>p) j=n; else j=p;
  iwork = (int *)CALLOC((size_t)j*n_threads,sizeof(int));
  xwork = (double *)CALLOC((size_t)j*n_threads,sizeof(double));
  V1 = (spMat *)CALLOC((size_t)n_threads,sizeof(spMat));
  W1 = (spMat *)CALLOC((size_t)n_threads,sizeof(spMat));
  for (i=0;i<n_threads;i++) { spalloc(V1+i,n,nzmax);spalloc(W1+i,n,nzmax);}
  /* the main work of the routine, by a substantial margin, is the following loop... */		
  B = (int *) CALLOC((size_t)nb,sizeof(int));
  cost = (double *)CALLOC((size_t)nb,sizeof(double));
  for (b=0;b<nb;b++) {
    B[b] = b;
    cost[b] = block[b].r_size*block[b].c_size;
  }
  revsort(cost,B,nb); /* R reverse sort on cost, to re-order B - see R.h*/
  /* idea of re-ordering is that expensive blocks should be processed early, 
     which helps balance - ordering only approximate since true cost depends 
     on NZ pattern of components, not just result size */ 
  tid=0; /* single thread default */
  #ifdef _OPENMP
  /* schedule: static, dynamic or guided - seems hard to do better than guided
     static leads to quite substantial load imbalance.*/
  #pragma omp parallel for schedule(guided) private(tid,b,i) num_threads(n_threads)
  #endif
  for (i=0;i<nb;i++) {
    b=B[i];
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #endif
    sXWXdij(w, d+tid*n,Xs,Xt,block[b].rb,block[b].cb,block[b].r,block[b].c,
	    dn+tid*n,ts,dt,W1+tid,V1+tid,xwx+b,&init,iwork+tid*j,xwork+tid*j);

  }
  FREE(B);FREE(cost);
  /* now deal with any constraints. This involves merging the sub-blocks of constrained tensor 
     products. Singletons should not be constrained this way. */
  for (nrc=b=i=0;i<nlt;i++) {
    rb=lt[i];k = 1; // number of sub-blocks in this block 
    if (!qc[rb]) k *= sub_blocks[rb];
    nrc += k; /* counts number of row blocks in constrained structure matrix */
    for (j=0;j<nrt;j++) {
      cb=rt[j];
      if (!qc[cb]) b += k*sub_blocks[cb]; else b += k;
    }
  }
  ncc = b/nrc; // columns of constraint matrix
  
  if (nr*nc != b) { /* then there is block merging to do */
    stc = (int *)CALLOC((size_t)b,sizeof(int)); /* constrained structure matrix */
    for (i=0;i<b;i++) stc[i] = -1; /* unfilled signal for checking */
    /* now work through the elements of stc, creating the merged blocks,
       ic,jc index row col in stc matrix
    */
    for (is0=ic0=i=0;i<nlt;i++) { // term row block loop
      if (symmetric) j1 = i+1; else j1 = nrt;
      rb = lt[i];
      for (js0=jc0=j=0;j<j1;j++) { // term col block loop
	/* if either model matrix block involved in this crossprod block is 
           constrained, then all its sub-blocks have to be combined, and the 
           constraint applied. Furthermore, blocks below or right of this 
           will require r_start and/or c_start decremented.
	*/
	cb = rt[j];//Rprintf("rb=%d cb=%d\n",rb,cb);
	if (!qc[rb]&&!qc[cb]) // simply copy entries from str to stc
	  for (ic = ic0,is = is0; is < is0 + sub_blocks[rb];is++,ic++)
	    for (jc = jc0, js = js0; js < js0 + sub_blocks[cb];js++,jc++) stc[ic + nrc * jc] = str[is + nr * js];
	else if (!qc[cb]){ /* qc[rb] only, so consolidate rows */
	  for (jc = jc0, js = js0; js < js0 + sub_blocks[cb];js++,jc++) { /* loop over sub block columns */
             k = sub_blocks[rb]*block_size[rb]*block_size[cb]; // total combined block size (dense)
	     b = str[is0+nr*js]; // starting sub-block
	     stc[ic0 + nrc * jc] = b;
	     block[b].r_size = sub_blocks[rb]*block_size[rb];
	     xwx[b].x = REALLOC(xwx[b].x,sizeof(double)*k);
	     // need to clear new part of xwx[b].x to zero...
	     for (xp=xwx[b].x+xwx[b].p[xwx[b].c],xp1=xwx[b].x+k;xp<xp1;xp++) *xp=0.0;
	     sp_to_dense_insitu(xwx+b,block[b].r_size); /* unpack the xwx[b] to dense storage in situ */ 
             xwx[b].m =  block[b].r_size;
	     for (is = is0+1; is < is0 + sub_blocks[rb];is++) { // copy in remaining blocks...
	       sp_to_dense(xwx + str[is + nr*js],xwx[b].x,(is-is0)*block_size[rb],0,xwx[b].m);
	     }
	     left_con(xwx + b,v[rb],xwork); // apply the constraint
	     block[b].r_size--;
	     dense_to_sp(xwx+b); /* convert back to a sparse representation for easier copying into final result */
	  }   
	} else if (!qc[rb]) { /* qc[cb] only so consolidate cols */
	  for (ic = ic0,is = is0; is < is0 + sub_blocks[rb];is++,ic++) { /* loop over sub block rows */
             k = sub_blocks[cb]*block_size[rb]*block_size[cb]; // total combined block size (dense)
	     b = str[is+nr*js0]; // starting sub-block
	     stc[ic + nrc * jc0] = b;
	     block[b].c_size = sub_blocks[cb]*block_size[cb];
	     xwx[b].x = REALLOC(xwx[b].x,sizeof(double)*k);
	     // need to clear new part of xwx[b].x to zero...
	     for (xp=xwx[b].x+xwx[b].p[xwx[b].c],xp1=xwx[b].x+k;xp<xp1;xp++) *xp=0.0;
	     sp_to_dense_insitu(xwx+b,block[b].r_size); /* unpack the xwx[b] to dense storage in situ */ 
	     xwx[b].c =  block[b].c_size;
             for (js = js0+1; js < js0 + sub_blocks[cb];js++) {
               sp_to_dense(xwx + str[is + nr*js],xwx[b].x,0,(js-js0)*block_size[cb],xwx[b].m);
	     }
	     right_con(xwx + b,v[cb],xwork);
	     block[b].c_size--;
	     dense_to_sp(xwx+b); /* convert back to a sparse representation for easier copying into final result */
	  }  
        } else {
	  /*  the sub-blocks must be consolidated and constrained...
          1. Write all the sub-blocks to a dense matrix.
          2. Apply the constraints to the dense matrix.
          3. Put the dense matrix into the xwx/block pointed to by the
             first str[is,js] and have stc[ic,jc] point to this.
          4. update the trailing block starts.   	  
	  */
	  k = sub_blocks[rb]*block_size[rb]*sub_blocks[cb]*block_size[cb]; // total combined block size (dense)
	  b = str[is0+nr*js0]; // starting sub-block
          stc[ic0 + nrc * jc0] = b;
	  block[b].r_size = sub_blocks[rb]*block_size[rb];
	  block[b].c_size = sub_blocks[cb]*block_size[cb];  
	  // now expand the storage for xwx[b].x to allow full dense storage
	  xwx[b].x = REALLOC(xwx[b].x,sizeof(double)*k);
	  // need to clear new part of xwx[b].x to zero...
	  for (xp=xwx[b].x+xwx[b].p[xwx[b].c],xp1=xwx[b].x+k;xp<xp1;xp++) *xp=0.0;
	  sp_to_dense_insitu(xwx+b,block[b].r_size); /* unpack the xwx[b] to dense storage in situ */ 
          xwx[b].m =  block[b].r_size;
	  xwx[b].c =  block[b].c_size;
	  // copy in remaining blocks...
	  for (is = is0; is < is0 + sub_blocks[rb];is++) {
	    for (js = js0; js < js0 + sub_blocks[cb];js++) {
	      if ((is==is0&&js==js0) ||(symmetric&&j==i&&js>is)) { // nothing to do
	      }	else {
		/* copy the current sub-block into xwx[b].x at row is*block_size[r]
                   js*block_size[c]  */
                sp_to_dense(xwx + str[is + nr*js],xwx[b].x,(is-is0)*block_size[rb],(js-js0)*block_size[cb],xwx[b].m);
              }
	    } // js/jc
	  } // is/ic
	  if (symmetric && i==j) { /* fill in upper triangle of block */
	    k = xwx[b].m;
            for (is=0;is<k;is++) for (js=0;js<is;js++) xwx[b].x[js + k * is] = xwx[b].x[is + k * js]; 
	  }
	  /* now apply the constraints and adjust the row/col starts */
	  if (qc[rb]) {
	    left_con(xwx + b,v[rb],xwork);
	    block[b].r_size--;
	  }
	  if (qc[cb]) {
	    right_con(xwx + b,v[cb],xwork);
	    block[b].c_size--;
	  }
	  dense_to_sp(xwx+b); /* convert back to a sparse representation for easier copying into final result */
	} /* sub block consolidation */
	js0 += sub_blocks[cb]; // str col start update
	if (qc[cb]) jc0++; else jc0 += sub_blocks[cb]; // stc col start update
      } // j loop (term block cols)
      is0 += sub_blocks[rb];  // str row start update
      if (qc[rb]) ic0++; else ic0 += sub_blocks[rb];// stc row start update
      
    } // i loop  (term block rows)
    /* now run through stc enforcing consistency of block starts */
    k = stc[0];
    for (i=1;i<nrc;i++) {
      /* set r_start to correct value for first block on each row (first row always ok)... */
      b = stc[i-1];k = stc[i];
      block[k].r_start = block[b].r_start + block[b].r_size;
      if (symmetric) j1=i; else j1=ncc-1;
      for (j=0;j<j1;j++) { /* work along cols of row */
        b = stc[i + j * nrc];k = stc[i + (j+1) * nrc];
	block[k].c_start = block[b].c_start + block[b].c_size;
	block[k].r_start = block[b].r_start; /* copy r_start from previous block on row */
      }	
    }
    /* get total rows and cols of result... */
    rcum = block[k].r_start + block[k].r_size;
    ccum = block[k].c_start + block[k].c_size;
  } else { /* structure matrix is unchanged */
    stc = str;
  }  
  /* finally copy the sub blocks into the main result matrix */ 

  for (k=0,i=0;i<nrc;i++) { /* obtain the storage requirement (slightly oversized in symmetric case) */
    if (symmetric) j1 = i+1; else j1 = ncc;
    for (j=0;j<j1;j++) {
      b = stc[i + nrc * j];
      if (xwx[b].p[0]<0) k += xwx[b].m * xwx[b].c; /* dense matrix */
      else k += xwx[b].p[xwx[b].c]; /* sparse matrix */
    }  
  }

  /* write out results to a sparse matrix for return */
  
  if (symmetric) {
    XWX = PROTECT(R_do_new_object(PROTECT(R_getClassDef("dsCMatrix"))));nprot++;nprot++; // create symmetric sparse matrix
    SET_STRING_ELT(PROTECT(R_do_slot(XWX,ul_sym)),0,PROTECT(mkChar("L")));nprot++;nprot++; // set to lower triangle storage - note charecter setting special
  } else { XWX = PROTECT(R_do_new_object(PROTECT(R_getClassDef("dgCMatrix"))));nprot++;nprot++;}
  dim = INTEGER(R_do_slot(XWX,dim_sym));
  dim[0] = rcum;dim[1] = ccum;
  if (symmetric) { /* don't know exact size yet - make temporary over-sized storage */
    XWXi = (int *) CALLOC((size_t)k,sizeof(int));
    XWXx = (double *) CALLOC((size_t)k,sizeof(double));
  } else { /* not symmetric - storage size already correct */
    R_do_slot_assign(XWX,x_sym,allocVector(REALSXP,k));
    R_do_slot_assign(XWX,i_sym,allocVector(INTSXP,k));
    XWXi = INTEGER(R_do_slot(XWX,i_sym));
    XWXx = REAL(R_do_slot(XWX,x_sym));
  }  
  R_do_slot_assign(XWX,p_sym,allocVector(INTSXP,dim[1]+1));
  XWXp = INTEGER(R_do_slot(XWX,p_sym));
  
  /* Need to work through column by column, employing the nr by nc (constrained) structure matrix to get 
     the right components.... */
  XWXp[0]= 0;
  for (ii=0,jj=0,j=0;j<ncc;j++) { /* cols of (constrained) structure matrix */
    for (js=0;js<xwx[stc[j*nrc+j*symmetric]].c;js++,jj++) { /* loop over columns of this xwx block, jj is output column */
      for (i=j*symmetric;i<nrc;i++) { /* the rows of the structure matrix */
	k = stc[i+j*nrc];
	if (k<0) Rprintf("\n WARNING: incomplete constrained structure matrix\n ");
	Mp = xwx + k; /* the current xwx matrix */
	blp = block + k; /* the current block information */
	is = Mp->p[js]; /* start of current column of xwx */
	if (i==j&&symmetric) { /* main diagonal block - find first entry below diagonal*/
	  while (Mp->i[is]<js&&is<Mp->p[js+1]) is++;
	}
	//Rprintf("Mp->c = %d js = %d j = %d i = %d\n",Mp->c,js,j,i);
	while (is < Mp->p[js+1]) { 
	  XWXi[ii] = Mp->i[is] + blp->r_start;
	  XWXx[ii] = Mp->x[is];ii++;is++;
        }	   
      } /* rows of structure matrix */
      XWXp[jj+1] = ii;
    } /* cols of this block */
  } /* cols structure matrix */
  if (symmetric) { /* initial allocation oversized - copy to exact sized return object */
    R_do_slot_assign(XWX,x_sym,allocVector(REALSXP,ii));
    R_do_slot_assign(XWX,i_sym,allocVector(INTSXP,ii));
    ip = INTEGER(R_do_slot(XWX,i_sym));
    xp = REAL(R_do_slot(XWX,x_sym));
    for (i=0;i<ii;i++) {xp[i]=XWXx[i];ip[i]=XWXi[i];}
    FREE(XWXi);FREE(XWXx);
  }  
  if (ncc*nrc!=nr*nc) FREE(stc);
  FREE(block);//xwx,
  for (i=0;i<n_threads;i++) {
    spfree(V1+i,1);spfree(W1+i,1);
  }
  FREE(W1);FREE(V1);
  spfree(Xt,mx);FREE(Xt);
  spfree(xwx,nb);FREE(xwx);
  FREE(iwork); FREE(xwork);FREE(dn);FREE(d); FREE(str);FREE(sub_blocks);FREE(block_size);
  FREE(Xs); // free the sparse matrix array
  FREE(v);
  UNPROTECT(nprot);
  return(XWX);
} /* sXWXd */  


/******************************************************************************
  X'y is relatively straightforward............
*******************************************************************************/

SEXP sXyd(SEXP X,SEXP Y,SEXP LT) {
/* Forms X'y where X is defined by sparse discrete storage blocks and y is a vector.
   lt selects the blocks to include - it relates to whole smoothers. 
*/
  spMat *Xs,*Xj;
  int mx,i,j,k,q,out_start,n,bb,nzmax=0,*dim,*kd,*ks,*r,*off_start,*off,nt,
    *ts,*dt,*qc,*lt,nlt,*dn,*ki,*c,*p,P,s,l,ii,b,c0,os,qq,nc=0,no=0,*tps,cy,yj;
  SEXP Xd,M,XY,
    i_sym,x_sym,dim_sym,p_sym,
    KD,R,KS,OFF,TS,DT,QC,V,OFFS;
  double **v,*y,*Xy,*yb,*My,*Xyo,*yp;
  p_sym = install("p");
  dim_sym = install("Dim");
  i_sym = install("i");
  x_sym = install("x");
  Xd = getListEl(X,"Xd"); /* the sparse matrix list */
  KD = getListEl(X,"kd"); /* the matrix of index vectors */
  n = nrows(KD); /* number of data */
  /* get the row indices (forward then reverse) */
  KD = PROTECT(coerceVector(KD,INTSXP));
  kd = INTEGER(KD);
  R = getListEl(X,"r");
  R = PROTECT(coerceVector(R,INTSXP));
  r = INTEGER(R);
  OFF =  getListEl(X,"off"); /* the offset list for the reverse indices */
  OFF = PROTECT(coerceVector(OFF,INTSXP));
  off = INTEGER(OFF);
  OFFS =  getListEl(X,"offstart"); /* the start points in the offset array */
  OFFS = PROTECT(coerceVector(OFFS,INTSXP));
  off_start = INTEGER(OFFS);
  /* get the matrix defining the range of k vectors for each matrix */
  KS = getListEl(X,"ks");
  KS = PROTECT(coerceVector(KS,INTSXP));
  ks = INTEGER(KS);
  
  mx = length(Xd); /* list length */
  /* read sparse matrix list into structure suitable for passing on to 
     other routines */
  Xs = (spMat *) CALLOC((size_t)mx,sizeof(spMat)); // sparse matrix list
  for (i=0;i<mx;i++) { // work through the sparse matrices
    M = VECTOR_ELT(Xd, i); // the ith sparse matrix
    Xs[i].x = REAL(R_do_slot(M,x_sym)); // the elements of ith matrix
    Xs[i].p = INTEGER(R_do_slot(M,p_sym)); // .p[j] is where jth col starts in .i and .x
    Xs[i].i = INTEGER(R_do_slot(M,i_sym)); // row index for each non-zero element
    dim = INTEGER(R_do_slot(M,dim_sym));
    Xs[i].m = dim[0];Xs[i].c = dim[1]; // matrix is .m by .c
    Xs[i].k = kd + n * ks[i]; // the .nk index n vectors (end to end)
    Xs[i].r = r + n * ks[i];  // ... equivalent reverse indices
    Xs[i].n = n; // rows in full matrix
    Xs[i].nk = ks[mx+i] - ks[i]; // number of index vectors for this term
    Xs[i].off = off + off_start[ks[i]]; // the offsets for index vectors into r 
    j = Xs[i].p[Xs[i].c]; /* number of NZs in this matrix */
    if (j>nzmax) nzmax=j; /* largest amount of NZ storage required by any eleemnt of Xs */
  } /* sparse matrix array filling loop */    

  /* now deal with the smooth term information... */
  TS = getListEl(X,"ts");
  nt = length(TS); // number of smooth terms
  TS = PROTECT(coerceVector(TS,INTSXP));
  ts = INTEGER(TS); // term starts in matrix array
  DT = getListEl(X,"dt");
  DT = PROTECT(coerceVector(DT,INTSXP));
  dt= INTEGER(DT); // number of marginal matrices for this term.
  QC = getListEl(X,"qc");
  QC = PROTECT(coerceVector(QC,INTSXP));
  qc= INTEGER(QC); // the constraint indicator nt-vector
  V = getListEl(X,"v");
  v = (double **) CALLOC((size_t)nt,sizeof(double *));
  for (i=0;i<nt;i++) if (qc[i]) {
    v[i] = REAL(VECTOR_ELT(V, i));
  }    
  nlt = length(LT);lt = INTEGER(LT);
  for (i=0;i<nlt;i++) if (qc[lt[i]]) nc++; /* count actual constraints */
  y = REAL(Y);
  cy = length(Y)/n; // cols of y
  out_start = 0;
  /* figure out length of return vector, Xy, and allocate it... */
  for (ii=bb=no=i=0;i<nlt;i++) {
    b = lt[i]; /* selected smooth block */
    for (q=1,j=0;j<dt[b];j++) q *= Xs[ts[b]+j].c; /* dim of this term */
    no += q; /* add to total count */
    if (Xs[ts[b]+dt[b]-1].m > bb) bb = Xs[ts[b]+dt[b]-1].m; /* max number of rows in a final marginal */ 
    if (ii<dt[b]) ii=dt[b]; /* maximum number of marginals */
  }  
 
  Xy = (double *) CALLOC((size_t) no*cy,sizeof(double));
  yb = (double *) CALLOC((size_t) bb,sizeof(double));
  c = (int *) CALLOC((size_t)ii,sizeof(int)); 
  p = (int *) CALLOC((size_t)ii,sizeof(int)); /* storage for marginal numbers of columns */ 
  My = (double *) CALLOC((size_t) n * ii,sizeof(double)); /* storage for partial column products */
  dn = (int *) CALLOC((size_t) n,sizeof(int));
  tps = (int *) CALLOC((size_t) nlt+1,sizeof(int)); /* smooth term parameter starts */
  for (bb=0;bb<nlt;bb++) {
    b = lt[bb]; /* block to proccess next */
    tps[bb]=out_start;
    if (dt[b] == 1) { /* singleton */
      for (yj=0;yj<cy;yj++) { /* loop over cols of y */
	yp = y + yj * n; // current col of y
        for (i=0;i<Xs[ts[b]].m;i++) yb[i] = 0.0;
        for (s=0;s<Xs[ts[b]].nk;s++) { /* loop for any summation convention */
          ki = Xs[ts[b]].k+s*n;
          for (i=0;i<n;i++) yb[ki[i]] += yp[i]; /* accumulate to yb using index for compressed matrix */ 
        }
        /* now form the singleton product */
        spMtv(Xs + ts[b],yb,Xy + yj*no + out_start,0);
      }	
      out_start += Xs[ts[b]].c;
    } else { /* tensor */
      for (P=1,i=0;i<dt[b]-1;i++) {
	p[i] = Xs[ts[b]+i].c;
	P *= p[i]; /* total cols divided by number of cols of final marginal */
        c[i] = 0; /* c[j] indexes current column of jth marginal */ 
      }
      p[i] = Xs[ts[b]+i].c;
      /* for tensor products the summation handling should be external to the product 
         formation, otherwise storage of partial products can't work without potentially
         very large storage costs */
      os = out_start;
      for (yp=y,yj=0;yj<cy;yj++,yp += n)  /* loop over cols of y */
      for (s=0;s<Xs[ts[b]].nk;s++) { /* loop for any summation convention */
        k = 0;out_start=os;
	for (i=0;i<dt[b];i++) c[i] = 0; /* clear the column index counter */
	if (s) for (i=0;i<n;i++) dn[i] = 0;
        ki = Xs[ts[b]+dt[b]-1].k + n *s; /* index for final marginal */
        for (l=0;l<P;l++) { /* loop over columns of row kronecker product of marginals excluding final */
	  /* need to clear yb */
	  for (i=0;i<Xs[ts[b]+dt[b]-1].m;i++) yb[i] = 0.0;
	  for (j=k;j<dt[b]-1;j++) {
	    Xj = Xs + ts[b] + j;
	    for (q=Xj->p[c[j]];q<Xj->p[c[j]+1];q++)
	      for (qq = Xj->i[q],ii=Xj->off[qq + (Xj->m+1)*s];ii<Xj->off[qq+1+ (Xj->m+1)*s];ii++) { 
              i = Xj->r[ii + s*n]; /* full matrix row */
	      if (dn[i]==j) { /* not already zeroed*/
	        dn[i]++;
	        /* get required partial column products with y... */
	        if (j>0) My[i+n*j] = My[i+n*(j-1)]*Xj->x[q];
	        else My[i] = yp[i] * Xj->x[q];
	        if (j==dt[b]-2) yb[ki[i]] += My[i+j*n]; /* accumulate according to final marg index */ 
	      }  
	    }
	  }  
	
	  spMtv(Xs + ts[b]+dt[b]-1,yb,Xy + yj*no + out_start,s); /* clear at s=0, add thereafter */
	  out_start += Xs[ts[b]+dt[b]-1].c; 
	  /* now update the column counter and figure out which partial products
             have to be updated.. */
	  k = dt[b]-2;c[k]++;
	  while (c[k]==p[k])  {
	    c[k]=0; if (k) { k--;c[k]++;}
	  }
	  /* roll back dn to state before previous column starts */
	  for (j=dt[b]-2;j>=k;j--) {
	    if (c[j]>0) c0 = c[j] - 1; else c0 = p[j] - 1;
	    Xj = Xs + ts[b] + j;
	    for (q=Xj->p[c0];q<Xj->p[c0+1];q++)
	    for (qq = Xj->i[q],ii=Xj->off[qq];ii<Xj->off[qq+1];ii++) {
              i = Xj->r[ii + s*n]; /* full matrix row */
	      if (dn[i]==j+1) dn[i]--; /* role back counter to previous state */
	    }	
	  }  
	} /* col loop (l) */
      }	/* summation loop (s) */
    } /* tensor */ 
  } /* smooth terms */
  tps[bb] = out_start;
  /* Now copy result to an output vector, dealing with any constraints at the same time */
  PROTECT(XY=allocVector(REALSXP,cy*(no-nc))); /* vector of dimension of Xy, less number of constraints */
  Xyo=REAL(XY);
  for (yp=Xy,yj=0;yj<cy;yj++,Xyo += no-nc,yp += no)
  for (i=bb=0;bb<nlt;bb++) {
    b = lt[bb]; /* current block */
    if (qc[b]) { /* apply constraint */
      k = tps[bb+1]-tps[bb]; /* term dimension, pre-constraint */
      left_con_vec(yp+tps[bb],v[b],Xyo+i,k,0);
      i += k-1; /* update output vector index */
    } else { /* no constraint - just copy */
      for (j=tps[bb];j<tps[bb+1];j++,i++) Xyo[i] = yp[j];
    }  
  }  
  
  FREE(Xy);FREE(Xs);FREE(v);FREE(yb);FREE(c);FREE(p);FREE(My);FREE(tps);FREE(dn);
  UNPROTECT(9);
  return(XY);
} /* sXyd */

/************************************************************************
 Functions for Xbeta and diag(XVX'). R interfaces are sXbd and sdiagXVXd,
 which are drivers for sXbdwork and sXbsdwork that handle the actual 
 multiplication when 'beta' is dense or sparse, respectively.
************************************************************************/

int spac(int *i,int i0,int j0,int j1,int r,int c,int *ii,int *q) {
/* subsets a sparse vector, and packs it into an r by c sparse matrix, with 
   column start index vector q and row index ii. q is length c+1. 
   i gives the non-zero rows of the input matrix. i0 is the starting row. 
   ii is usually
   same length as i to guarantee adequate length. j0 is the sugested starting 
   entry in i, but may be below the actual starting value. j1 is the length of i.
   Return value is actual j0, while returned j0 + q[c] - 1 gives final entry. 
*/
  int i1,ik,j00,k;
  i1 = i0 + r*c; // one beyond final entry in the input
  while (j0>0&&i[j0]>i0) j0--;
  while (i[j0]<i0&&j0<j1) j0++;
  j00 = j0;
  k = 0; /* column counter */
  q[0] = 0;ik = 0; 
  while (i[j0]<i1&&j0<j1) {
    while (i[j0]-i0 < r+k*r && j0<j1) ii[ik++] = i[j0++] - i0 - k*r;
    if (j0<j1) while (k<c && i[j0]-i0 >= r+k*r) {k++;q[k]=ik;}
  }
  while (k<c) {k++;q[k]=ik;}
  return(j00);
} /* spac */  

void sXbsdwork(double *Xb,double *a,spMat beta0,int bp,spMat *Xs,double **v,int *qc,int nt,
	       int *ts,int *dt,int *lt,int nlt,int n,double *work,int *worki,int unit_a) {
/* Sparse beta0/b version of routine that actually does the work of forming Xb for sXbd and sdiagXVX. 
   Note that Xb is not cleared to zero in this function, and the result is added to its initial value.
   This version uses a partial product approach to reduce workload on tensor products.
   * Xb is a n-vector to contain the result
   * a is the n-vector of weights to use in Xb accumulation.
   * beta0 is the vector to multiply, stored as a sparse matrix
   * Xs is the array of sparse marginal discretized matrices.
   * v[i] defines the ith constraint if qc[i]!=0
   * nt is the number of model terms
   * ts is the term starts in Xs 
   * dt is how many Xs matrices are involved in each term
   * lt is the nlt-array of terms to include 
    * work and worki are workspace of size n*maxd+2*(bp+nc) + nm and nm + n + 2*(nt+bp+nc+1) + maxdim + 3*maxd respectively
     maxd is largest number of marginals in a term, maxdim the largest number of coefs in a term, nc is number of constraints
     nm is largest number of rows in a marginal matrix. 
     -- see sXbd for how to initialize.    
   * unit_a=0 to use a and 1 not to (treat as 1).
*/
  spMat *Xj,beta,bsub,C; 
  int nc,i,j,k,k0,q,bb,*dim,maxd,maxdim,//*ki,
    *dn,*c,*c0,*tps,*p,P,s,l,ii,b,qq,m,/*mx,*/ j0,j1,*ri;
  double *d;
  for (nc=i=0;i<nt;i++) if (qc[i]) nc++; // count constraints
  beta.x = work;work += bp+nc;
  beta.i = worki;worki += bp + nc;
  beta.p = worki;worki += 2;
  dim = worki;worki += nt; // smooth term dimensions
  tps = worki;worki += nt; // smooth term start in param vector
  dn = worki;worki += n;
  /*for (maxdim=maxd=0,i=0;i<nt;i++) { // compute smooth term sizes 
    b = ts[i];
    for (dim[i]=1,j=0;j<dt[i];j++) dim[i] *= Xs[b+j].c;
    if (i) tps[i] = tps[i-1] + dim[i-1];
    if (dt[i] > maxd) maxd = dt[i];
    if (dim[i] > maxdim) maxdim = dim[i];
    }*/
  for (maxdim=maxd=k=0,bb=0;bb<nlt;bb++) { // compute smooth term sizes 
    i = lt[bb];
    b = ts[i];
    for (dim[i]=1,j=0;j<dt[i];j++) dim[i] *= Xs[b+j].c;
    tps[i] = k; k += dim[i];
    if (dt[i] > maxd) maxd = dt[i];
    if (dim[i] > maxdim) maxdim = dim[i];
  }
  /* constraint handling q is start row for beta0, k start row for beta, j0 is current start element of 
     sparse beta0, j1 start element in beta... */
  /*for (j0=j1=q=k=i=0;i<nt;i++) { 
    if (qc[i]) { // constraints make that subsection of beta dense  
      // copy relevant sections of beta0 to a dense vector 
      for (j=0;j<dim[i];j++) work[j] = 0.0; // This is re-cycled below - sizing ok as require bp sized allocation later 
      while (j0<beta0.p[1] && beta0.i[j0]<q+dim[i]-1) {
	work[beta0.i[j0]-q] = beta0.x[j0];j0++; 
      }
      left_con_vec(work,v[i],beta.x+j1,dim[i],1); // undo constraint  
      for (j=k+dim[i];k<j;k++,j1++) { beta.i[j1] = k;}
      q += dim[i]-1;// note k updated in above loop
    } else { // no constraint
      while (j0<beta0.p[1] && beta0.i[j0]<q+dim[i]) {
        beta.i[j1] = beta0.i[j0] + k - q; 
	beta.x[j1] = beta0.x[j0];
	j1++;j0++;
      }
      q += dim[i];k += dim[i];
    }
    } */
  for (j0=j1=q=k=bb=0;bb<nlt;bb++) {
    i = lt[bb];
    if (qc[i]) { // constraints make that subsection of beta dense  
      // copy relevant sections of beta0 to a dense vector 
      for (j=0;j<dim[i];j++) work[j] = 0.0; // This is re-cycled below - sizing ok as require bp sized allocation later 
      while (j0<beta0.p[1] && beta0.i[j0]<q+dim[i]-1) {
	work[beta0.i[j0]-q] = beta0.x[j0];j0++; 
      }
      left_con_vec(work,v[i],beta.x+j1,dim[i],1); // undo constraint  
      for (j=k+dim[i];k<j;k++,j1++) { beta.i[j1] = k;}
      q += dim[i]-1;// note k updated in above loop
    } else { // no constraint
      while (j0<beta0.p[1] && beta0.i[j0]<q+dim[i]) {
        beta.i[j1] = beta0.i[j0] + k - q; 
	beta.x[j1] = beta0.x[j0];
	j1++;j0++;
      }
      q += dim[i];k += dim[i];
    }
  }
  
  beta.p[1] = j1;beta.p[0]=0;  
 
  d = work; work += n*maxd; /* partial product array */
  
  p = worki;worki += maxd;
  c = worki;worki += maxd;
  c0 = worki;worki += maxd; 
  bsub.i = worki; worki += bp + nc; 
  bsub.x = work; work += bp + nc; 
  bsub.p = worki;worki += maxdim + 1; 
  spalloc(&C,maxdim,1000);
  j0 = 0; /* starting entry in beta */
  for (bb=0;bb<nlt;bb++) { // now actually do the multiplying.
    b = lt[bb];
    if (dt[b]==1) { // singleton
      /* call spac to extract sparse subvector equivalent to dense vector beta+tps[b] */
      Xj = Xs + ts[b];
      j0 = spac(beta.i,tps[b],j0,beta.p[1],Xj->c,1,bsub.i,bsub.p);
      for (i=0;i<bsub.p[1];i++) bsub.x[i] = beta.x[i+j0];
      j0 += bsub.p[1];bsub.m = Xj->c;bsub.c=1;
      cs_mult(Xj,&bsub,&C,worki,work,1); // NOTE: work, iwork Xs[ts[b]].m used

      /* needs proper unpacking here. that is work through non-zeros of C
         putting them in correct place using the reverse index of Xs[ts[b]] */
      for (s=0;s<Xj->nk;s++) { // summation over multiple indices
        ri = Xj->r + s * n;
	for (q=0;q<C.p[1];q++)
	  for (j=Xj->off[q+s*Xj->m];j<Xj->off[q+1+s*Xj->m];j++)
	    if (unit_a) Xb[ri[j]] += C.x[q]; else
	      Xb[ri[j]] += a[ri[j]] * C.x[q];
      }		
    } else { // tensor product
      for (j=0;j<dt[b];j++) p[j]=Xs[ts[b]+j].c;
      for (P=1,j=0;j<dt[b]-1;j++) P *= p[j];
      /* get product of final marginal with B, w and x need to be A.m vectors where A
         is the marginal (first argument) */
      // First pack beta + tps[b] into a P column matrix, bsub
      Xj = Xs + ts[b] + dt[b] - 1;
      j0 = spac(beta.i,tps[b],j0,beta.p[1],Xj->c,P,bsub.i,bsub.p);
      for (i=0;i<bsub.p[P];i++) bsub.x[i] = beta.x[i+j0];
      j0 += bsub.p[P]; bsub.m = Xj->c;bsub.c = P;
      cs_mult(Xj,&bsub,&C,worki,work,1); 
      for (s=0;s<Xs[ts[b]].nk;s++) { // summation loop
	for (i=0;i<dt[b];i++) c0[i]=c[i] = 0; /* clear the column index counters */
	k0=0;
        for (l=0;l<P;l++) {
	  if (C.p[l]<C.p[l+1]) { /* then there is something to do as C[,l] non-zero */
            k = k0; /* marginal that we need to start from */
	    k0 = dt[b]-2; /* now update to this marginal */
	    for (m=k;m<dt[b]-1;m++) { /* loop over marginals */
	      c0[m] = c[m];
	      Xj = Xs + ts[b] + m;
	      for (q = Xj->p[c[m]];q<Xj->p[c[m]+1];q++) {
                for (qq = Xj->i[q],ii=Xj->off[qq + (Xj->m+1)*s];ii<Xj->off[qq+1+ (Xj->m+1)*s];ii++) { 
                  i = Xj->r[ii + s*n]; /* full matrix row */
	          if (dn[i]==m) { /* not already zeroed */
	            dn[i]++;
	            if (m==0) d[i] = Xj->x[q]; else d[i+m*n] = d[i+m*n-n] * Xj->x[q];
		  }  
	        }	
	      }  
	    } /* m loop */
	    /* now multiply d by correct C element to update Xb */
	    m = dt[b] - 2; // just to make this explicit!
	    Xj = Xs + ts[b] + m + 1;
	    for (q = C.p[l];q<C.p[l+1];q++) {
	      for (qq = C.i[q],ii=Xj->off[qq + (Xj->m+1)*s];ii<Xj->off[qq+1+(Xj->m+1)*s];ii++) {
                i = Xj->r[ii + s*n]; /* full matrix row */
	        if (dn[i]==dt[b]-1) {
                  if (unit_a) Xb[i] += d[i+m*n] * C.x[q]; else Xb[i] += a[i] * d[i+m*n] * C.x[q];
	        }	
	      }    
	    }
	  } /* non zero C[:,l] code */
	  /* update the marginal column indices and downdate dn accordingly */
	  k = dt[b]-2;c[k]++;
	  while (c[k]==p[k]) {
            c[k]=0;
	    if (k>0) {
              k--;c[k]++;
	    }  
	  }
	  if (k<k0) k0=k; /* keep track of first marginal in need of update */
	  if (l==P-1 || C.p[l+1]<C.p[l+2]) {
	    for (m=dt[b]-2;m>=k0;m--) { 
	      Xj = Xs + ts[b] + m;
	      for (q = Xj->p[c0[m]];q<Xj->p[c0[m]+1];q++)
	      for (qq = Xj->i[q],ii=Xj->off[qq + (Xj->m+1)*s];ii<Xj->off[qq+1+ (Xj->m+1)*s];ii++) {
                i = Xj->r[ii + s*n];
	        if (dn[i]==m+1) dn[i]--;
	      }
	    }	 
	  }
        } /* l loop */
      } /* s loop */
    } /* tensor product */
  } /* term (bb) loop */
  spfree(&C,1);
} /* sXsbdwork */  

void sXbdwork(double *Xb,double *a,double *beta0,int bp,spMat *Xs,double **v,int *qc,int nt,
	      int *ts,int *dt,int *lt,int nlt,int n,double *work,int *worki,int unit_a) {
/* Routine that actually does the work of forming Xb for sXbd and sdiagXVX 
   Note that Xb is not cleared to zero in this function.
   This version uses a partial product approach to reduce workload on tensor products.
   * Xb is a n-vector to contain the result
   * a is the n-vector of weights to use in Xb accumulation.
   * beta0 is the bp-vector to multiply
   * Xs is the array of sparse marginal discretized matrices.
   * v[i] defines the ith constraint if qc[i]!=0
   * nt is the number of model terms
   * ts is the term starts in Xs 
   * dt is how many Xs matrices are involved in each term
   * lt is the nlt-array of terms to include    
   * work and worki are workspace of size n*maxd+bp+nc and n + 2*(nt+maxd) respectively
     maxd is largest number of marginals in a term, nc is number of constraints 
     -- see sXbd for how to initialize.
   * unit_a=0 to use a and 1 not to (treat as 1).
*/
  spMat *Xj;
  int nc,i,j,k,q,bb,*dim,maxd,ok,*tps,
    *dn,*ki,*c,*p,P,s,l,ii,b,qq,m,mx;
  double *beta,*C,*d;
  for (nc=i=0;i<nt;i++) if (qc[i]) nc++; // count constraints
  beta = work;work += bp+nc;//(double *)CALLOC((size_t)bp+nc,sizeof(double));
  dim = worki;worki += nt;//(int *)CALLOC((size_t)nt,sizeof(int)); // smooth term dimensions
  tps = worki;worki += nt;//(int *)CALLOC((size_t)nt,sizeof(int)); // smooth term start in param vector
  dn = worki;worki += n;//(int *)CALLOC((size_t)n,sizeof(int));
  /*  for (maxd=0,i=0;i<nt;i++) { // compute smooth term sizes 
    b = ts[i];
    for (dim[i]=1,j=0;j<dt[i];j++) dim[i] *= Xs[b+j].c;
    if (i) tps[i] = tps[i-1] + dim[i-1];
    if (dt[i] > maxd) maxd = dt[i];
  }
  for (q=k=i=0;i<nt;i++) { // constraint handling 
    if (qc[i]) {
      left_con_vec(beta0+q,v[i],beta+k,dim[i],1); // undo constraint 
      q += dim[i]-1;k += dim[i];
    } else { // no constraint
      for (j=0;j<dim[i];j++,k++,q++) beta[k] = beta0[q];
    }  
  } */ 
  for (k=maxd=0,bb=0;bb<nlt;bb++) { // compute smooth term sizes 
    i = lt[bb]; // selected term
    b = ts[i];
    for (dim[i]=1,j=0;j<dt[i];j++) dim[i] *= Xs[b+j].c;
    tps[i] = k; k += dim[i];
    if (dt[i] > maxd) maxd = dt[i];
  }
  for (q=k=bb=0;bb<nlt;bb++) { // constraint handling
    i = lt[bb];
    if (qc[i]) {
      left_con_vec(beta0+q,v[i],beta+k,dim[i],1); // undo constraint 
      q += dim[i]-1;k += dim[i];
    } else { // no constraint
      for (j=0;j<dim[i];j++,k++,q++) beta[k] = beta0[q];
    }  
  }
  d = work;//(double *)CALLOC((size_t) n*maxd, sizeof(double)); /* partial product array */
  
  p = worki;worki += maxd;//(int *)CALLOC((size_t)maxd,sizeof(int));
  c = worki;//(int *)CALLOC((size_t)maxd,sizeof(int));
  
  for (bb=0;bb<nlt;bb++) { // now actually do the multiplying.
    b = lt[bb];
    if (dt[b]==1) { // singleton
      spMv(Xs+ts[b],beta+tps[b],d);
     for (s=0;s<Xs[ts[b]].nk;s++) { // summation over multiple indices
        ki = Xs[ts[b]].k + s * n;
	if (unit_a) for (j=0;j<n;j++) Xb[j] += d[ki[j]]; // unpack to full vector.
	else for (j=0;j<n;j++) Xb[j] += a[j] * d[ki[j]];
      }	
    } else { // tensor product
      for (j=0;j<dt[b];j++) p[j]=Xs[ts[b]+j].c;
      for (P=1,j=0;j<dt[b]-1;j++) P *= p[j];
      /* get product of final marginal with B, w and x need to be A.m vectors where A
         is the marginal (first argument) */
      mx = Xs[ts[b]+dt[b]-1].m;
      C = (double *)CALLOC((size_t)P * mx,sizeof(double));
      spMA(Xs+ts[b]+dt[b]-1,beta+tps[b],C,P); /* C = Xs[ts[b]+dt[b]-1] B, vec(B) = beta+tps[b] */
      for (s=0;s<Xs[ts[b]].nk;s++) { // summation loop
        ki = Xs[ts[b]+dt[b]-1].k + s * n;
	for (i=0;i<dt[b];i++) c[i] = 0; /* clear the column index counter */
	k=0;
        for (l=0;l<P;l++) {
          for (m=k;m<dt[b]-1;m++) {
	    Xj = Xs + ts[b] + m;
	    for (q = Xj->p[c[m]];q<Xj->p[c[m]+1];q++) {
              for (qq = Xj->i[q],ii=Xj->off[qq + (Xj->m+1)*s];ii<Xj->off[qq+1+ (Xj->m+1)*s];ii++) { 
                i = Xj->r[ii + s*n]; /* full matrix row */
	        if (dn[i]==m) { /* not already zeroed*/
	          dn[i]++;
	          if (m==0) d[i] = Xj->x[q]; else d[i+m*n] = d[i+m*n-n] * Xj->x[q];
	          if (m==dt[b]-2) {
		    if (unit_a) Xb[i] += d[i+m*n] * C[ki[i] + l*mx];
		    else Xb[i] += a[i] * d[i+m*n] * C[ki[i] + l*mx];
		  }  
	        }	
	      }  
	    }	
          } /* m loop */
	 
	  /* update the marginal column indices and downdate dn accordingly */
	  k = dt[b]-2;ok = 0;
	  while (!ok) {
	    Xj = Xs + ts[b] + k;
	    for (q = Xj->p[c[k]];q<Xj->p[c[k]+1];q++)
	    for (qq = Xj->i[q],ii=Xj->off[qq + (Xj->m+1)*s];ii<Xj->off[qq+1+ (Xj->m+1)*s];ii++) {
              i = Xj->r[ii + s*n];
	      if (dn[i]==k+1) dn[i]--;
	    }	
	    c[k]++;ok=1;
	    if (c[k]==p[k]) {
	      c[k]=0; if (k>0) {k--;ok=0;}
	    }  
	  }
        } /* l loop */
      } /* s loop */
      FREE(C);
    } /* tensor product */
  } /* term (bb) loop */
} /* sXbdwork */  


SEXP sXbd(SEXP X,SEXP BETA,SEXP LT) {
/* Form X beta (beta dense) */
  spMat *Xs/*,sbeta0*/;
  int bp,mx,i,j,k,n,nzmax=0,*dim,*kd,*ks,*r,*off_start,*off,nt,
    *ts,*dt,*qc,*lt,nlt,maxd,nc,*worki,maxm=0,bc;
  SEXP Xd,M,
    i_sym,x_sym,dim_sym,p_sym,
    KD,R,KS,OFF,TS,DT,QC,V,XB,OFFS;
  double **v,*beta0,*Xb,*work,a=1;
  p_sym = install("p");
  dim_sym = install("Dim");
  i_sym = install("i");
  x_sym = install("x");
  Xd = getListEl(X,"Xd"); /* the sparse matrix list */
  KD = getListEl(X,"kd"); /* the matrix of index vectors */
  n = nrows(KD); /* number of data */
  /* get the row indices (forward then reverse) */
  KD = PROTECT(coerceVector(KD,INTSXP));
  kd = INTEGER(KD);
  R = getListEl(X,"r");
  R = PROTECT(coerceVector(R,INTSXP));
  r = INTEGER(R);
  OFF =  getListEl(X,"off"); /* the offset list for the reverse indices */
  OFF = PROTECT(coerceVector(OFF,INTSXP));
  off = INTEGER(OFF);
  OFFS =  getListEl(X,"offstart"); /* the start points in the offset array */
  OFFS = PROTECT(coerceVector(OFFS,INTSXP));
  off_start = INTEGER(OFFS);
 
  /* get the matrix defining the range of k vectors for each matrix */
  KS = getListEl(X,"ks");
  KS = PROTECT(coerceVector(KS,INTSXP));
  ks = INTEGER(KS);
  
  mx = length(Xd); /* list length */
  /* read sparse matrix list into structure suitable for passing on to 
     other routines */
  Xs = (spMat *) CALLOC((size_t)mx,sizeof(spMat)); // sparse matrix list

  for (i=0;i<mx;i++) { // work through the sparse matrices
    M = VECTOR_ELT(Xd, i); // the ith sparse matrix
    Xs[i].x = REAL(R_do_slot(M,x_sym)); // the elements of ith matrix
    Xs[i].p = INTEGER(R_do_slot(M,p_sym)); // .p[j] is where jth col starts in .i and .x
    Xs[i].i = INTEGER(R_do_slot(M,i_sym)); // row index for each non-zero element
    dim = INTEGER(R_do_slot(M,dim_sym));
    Xs[i].m = dim[0];Xs[i].c = dim[1]; // matrix is .m by .c
    if (maxm<dim[0]) maxm=dim[0]; // largest number of rows in a marginal matrix
    Xs[i].k = kd + n * ks[i]; // the .nk index n vectors (end to end)
    Xs[i].r = r + n * ks[i];  // ... equivalent reverse indices
    Xs[i].n = n; // rows in full matrix
    Xs[i].nk = ks[mx+i] - ks[i]; // number of index vectors for this term
    Xs[i].off = off + off_start[ks[i]]; // the .nk offset .m+1 vectors for .r 
    j = Xs[i].p[Xs[i].c]; /* number of NZs in this matrix */
    if (j>nzmax) nzmax=j; /* largest amount of NZ storage required by any eleemnt of Xs */
  } /* sparse matrix array filling loop */    

  /* now deal with the smooth term information... */
  TS = getListEl(X,"ts");
  nt = length(TS); // number of smooth terms
  TS = PROTECT(coerceVector(TS,INTSXP));
  ts = INTEGER(TS); // term starts in matrix array
  DT = getListEl(X,"dt");
  DT = PROTECT(coerceVector(DT,INTSXP));
  dt= INTEGER(DT); // number of marginal matrices for this term.
  QC = getListEl(X,"qc");
  QC = PROTECT(coerceVector(QC,INTSXP));
  qc= INTEGER(QC); // the constraint indicator nt-vector
  V = getListEl(X,"v");
  v = (double **) CALLOC((size_t)nt,sizeof(double *));
  for (nc=i=0;i<nt;i++) if (qc[i]) {
    v[i] = REAL(VECTOR_ELT(V, i));
    nc++;
  }
  for (maxd=i=0;i<nt;i++) if (maxd<dt[i]) maxd=dt[i];
  nlt = length(LT);lt = INTEGER(LT);
  beta0 = REAL(BETA); // raw beta
  //bp = length(BETA);
  bp = nrows(BETA);bc = ncols(BETA);
  PROTECT(XB=allocVector(REALSXP,n*bc)); /* vector for X beta */
  Xb=REAL(XB);
  for (i=0;i<n*bc;i++) Xb[i] = 0.0;
 
  //k = n*maxd+2*(bp+nc) + maxm; // sparse version
  k = n*maxd+bp+nc; // dense version 
  work = (double *)CALLOC((size_t) k,sizeof(double));
  //k = 3*maxd+bp+2*(bp+nc+nt+1) + maxm + n; // sparse version
  k = n+2*(nt+maxd); // dense version
  worki = (int *)CALLOC((size_t) k,sizeof(int));
  /* // following is for sparse code testing only... 
  for (k=0,i=0;i<bp;i++) if (beta0[i]!=0.0) k++;
  spalloc(&sbeta0,2,k);
  for (k=0,i=0;i<bp;i++) if (beta0[i]!=0.0) {
    sbeta0.x[k] = beta0[i];sbeta0.i[k]=i;k++;
  }
  sbeta0.p[1] = k;
  sXbsdwork(Xb,&a,sbeta0,bp,Xs,v,qc,nt,ts,dt,tps,lt,nlt,n,work,worki,1);
  spfree(&sbeta0,1);
  */
  for (j=0;j<bc;j++,beta0 += bp,Xb += n) 
    sXbdwork(Xb,&a,beta0,bp,Xs,v,qc,nt,ts,dt,lt,nlt,n,work,worki,1);
  
  FREE(worki);FREE(work);FREE(v);FREE(Xs);
  UNPROTECT(9);
  return(XB);
} /* sXbd */ 

SEXP sdiagXVXt(SEXP X, SEXP V, SEXP LT, SEXP RT) {
/* Computes diag(XVX') when V is sparse and the model matrix is stored in discrete sparse form.
   The approach is based on using sXsbdwork to extract X[:,i] and XV[:,i] and accumulate 
   diag(XVX') = sum_i XV[:,i]*X[:,i]. LT and RT are indices of model terms to include in X and X'
   respectively. The cost is O(np) even under very high sparsity, because of fully clearing 'a' for
   each i. An alternative would record the non-zero entries in a, and wipe only those. 
*/
  spMat Vs,Ii,Vi,*Xs;
  int i,j,k,p,n,mx,*r,*kd,*off,*ks,*off_start,nt,*ts,*dt,*qc,nc,maxd,maxm=0,nlt,*lt,nrt,*rt,*dim,*worki;
  SEXP Xd,i_sym,x_sym,dim_sym,p_sym,KD,R,OFF,OFFS,KS,M,TS,DT,QC,Vc,DXVX;
  double **v,*dxvx,*a,*work;
  p_sym = install("p");
  dim_sym = install("Dim");
  i_sym = install("i");
  x_sym = install("x");
  Xd = getListEl(X,"Xd"); /* the sparse matrix list */
  KD = getListEl(X,"kd"); /* the matrix of index vectors */
  n = nrows(KD); /* number of data */
  KD = PROTECT(coerceVector(KD,INTSXP));
  kd = INTEGER(KD);
  R = getListEl(X,"r");
  R = PROTECT(coerceVector(R,INTSXP));
  r = INTEGER(R);
  OFF =  getListEl(X,"off"); /* the offset list for the reverse indices */
  OFF = PROTECT(coerceVector(OFF,INTSXP));
  off = INTEGER(OFF);
  OFFS =  getListEl(X,"offstart"); /* the start points in the offset array */
  OFFS = PROTECT(coerceVector(OFFS,INTSXP));
  off_start = INTEGER(OFFS);
  /* get the matrix defining the range of k vectors for each matrix */
  KS = getListEl(X,"ks");
  KS = PROTECT(coerceVector(KS,INTSXP));
  ks = INTEGER(KS);
  
  mx = length(Xd); /* list length */
  /* read sparse matrix list into structure suitable for passing on to 
     other routines */
  Xs = (spMat *) CALLOC((size_t)mx,sizeof(spMat)); // sparse matrix list
  for (i=0;i<mx;i++) { // work through the sparse matrices
    M = VECTOR_ELT(Xd, i); // the ith sparse matrix
    Xs[i].x = REAL(R_do_slot(M,x_sym)); // the elements of ith matrix
    Xs[i].p = INTEGER(R_do_slot(M,p_sym)); // .p[j] is where jth col starts in .i and .x
    Xs[i].i = INTEGER(R_do_slot(M,i_sym)); // row index for each non-zero element
    dim = INTEGER(R_do_slot(M,dim_sym));
    Xs[i].m = dim[0];Xs[i].c = dim[1]; // matrix is .m by .c
    if (maxm<dim[0]) maxm=dim[0]; // largest number of rows in a marginal matrix
    Xs[i].k = kd + n * ks[i]; // the .nk index n vectors (end to end)
    Xs[i].r = r + n * ks[i];  // ... equivalent reverse indices
    Xs[i].n = n; // rows in full matrix
    Xs[i].nk = ks[mx+i] - ks[i]; // number of index vectors for this term
    Xs[i].off = off + off_start[ks[i]]; // the .nk offset .m+1 vectors for .r 
    //j = Xs[i].p[Xs[i].c]; /* number of NZs in this matrix */
    //if (j>nzmax) nzmax=j; /* largest amount of NZ storage required by any eleemnt of Xs */
  } /* sparse matrix array filling loop */
  /* read V into a sparse matrix stucture */
  Vs.x =  REAL(R_do_slot(V,x_sym));Vs.p = INTEGER(R_do_slot(V,p_sym));
  Vs.i = INTEGER(R_do_slot(V,i_sym)); // row index for each non-zero element
  dim = INTEGER(R_do_slot(V,dim_sym));
  p = Vs.m = dim[0];Vs.c = dim[1]; // matrix is .m by .c
  
  /* now deal with the smooth term information... */
  TS = getListEl(X,"ts");
  nt = length(TS); // number of smooth terms
  TS = PROTECT(coerceVector(TS,INTSXP));
  ts = INTEGER(TS); // term starts in matrix array
  DT = getListEl(X,"dt");
  DT = PROTECT(coerceVector(DT,INTSXP));
  dt= INTEGER(DT); // number of marginal matrices for this term.
  QC = getListEl(X,"qc");
  QC = PROTECT(coerceVector(QC,INTSXP));
  qc= INTEGER(QC); // the constraint indicator nt-vector
  Vc = getListEl(X,"v");
  v = (double **) CALLOC((size_t)nt,sizeof(double *));
  for (nc=i=0;i<nt;i++) if (qc[i]) {
    v[i] = REAL(VECTOR_ELT(Vc, i));
    nc++;
  }
  for (maxd=i=0;i<nt;i++) if (maxd<dt[i]) maxd=dt[i];
  nlt = length(LT);lt = INTEGER(LT);
  nrt = length(RT);rt = INTEGER(RT);
  PROTECT(DXVX=allocVector(REALSXP,n)); /* vector for diag(XVX') */
  dxvx=REAL(DXVX);
  a = (double *)CALLOC((size_t)n,sizeof(double));
  for (i=0;i<n;i++) a[i] = dxvx[i] = 0.0;
  /* Vi holds ith col of V, and Ii the ith col of identity */
  Vi.p = (int *)CALLOC((size_t)2,sizeof(int));
  Ii.m = Vi.m = p; Ii.c = Vi.c = 1;
  spalloc(&Ii,1,1);Ii.x[0]=1.0;Ii.p[0]=0;Ii.p[1]=1;
  k = n*maxd+2*(p+nc) + maxm;
  work = (double *)CALLOC((size_t) k,sizeof(double));
  k = 3*maxd+p+1+2*(p+nc+nt+1) + maxm + n +1; /* p+1 is upper bound on max coefs per term (single tp with constraint could be p+1) */ 
  worki = (int *)CALLOC((size_t) k,sizeof(int));
  for (i=0;i<p;i++) { /* work through the columns */
    /* extract XV[:,i] ... */
    Vi.i = Vs.i + Vs.p[i]; Vi.x = Vs.x + Vs.p[i];
    Vi.p[1] = Vs.p[i+1] - Vs.p[i];
    for (j=0;j<n;j++) a[j]=0.0; // NOTE: inefficient
    sXbsdwork(a,a,Vi,p,Xs,v,qc,nt,ts,dt,lt,nlt,n,work,worki,1);
    /* form element wise product of X[:,i] and XV[:,i] and accumulate to result */
    Ii.i[0] = i; /* put the single entry here! */
    sXbsdwork(dxvx,a,Ii,p,Xs,v,qc,nt,ts,dt,rt,nrt,n,work,worki,0);

    //sXbsdwork(double *Xb,double *a,spMat beta0,int bp,spMat *Xs,double **v,int *qc,int nt,
    //	      int *ts,int *dt,int *lt,int nlt,int n,double *work,int *worki,int unit_a);
  }
  spfree(&Ii,1);FREE(Vi.p);FREE(work);FREE(worki);FREE(a);FREE(v);FREE(Xs);
  UNPROTECT(9);
  return(DXVX);
} /* sdiagXVXt */


/***************************************************************************
 Sparse row tensor product routine...
***************************************************************************/


SEXP stmm(SEXP X) {
/* X is a list of dgCMatrix sparse matrices whose row tensor product is required.
   Basic idea is... 
   1. to maintain partial products of columns, so that we only need to
      update part of product that changes for each column.
   2. to count columns processed in forming a column product until we encounter 
      a zero at each row. If a zero already encountered on some row, then 
      no furhter processing of row is needed. 

   set env R_HOME /usr/local/lib/R
   break stmm

   dyn.load("spardisc.so")
   library(Matrix);library(mgcv)
   set.seed(1)
   B <- A <- list()
   n <- 10; p <- c(2,3,2)
   A[[1]] <- matrix(sample(c(0,0,0,0,0,1:6),n*p[1],replace=TRUE),n,p[1])
   A[[2]] <- matrix(sample(c(0,0,0,0,0,1:6),n*p[2],replace=TRUE),n,p[2])
   A[[3]] <- matrix(sample(c(0,0,0,0,0,1:6),n*p[3],replace=TRUE),n,p[3])
   for (i in 1:3) B[[i]] <- as(A[[i]],"dgCMatrix")
   R <- .Call("stmm",B)
   Rf <- tensor.prod.model.matrix(A)
   Rf;R;Rf-R
*/
  spMat *Xs;
  int mx,i,ii,j,k,l,*c,n,p,*Rp,*Ri=NULL,rj,*dn,*dim,op,*ip;
  SEXP M,R,i_sym,x_sym,dim_sym,p_sym;
  double *Rx=NULL, *pp,*pp0,*pp1;
  /* register the names of the slots in the sparse matrices */
  p_sym = install("p");
  dim_sym = install("Dim");
  i_sym = install("i");
  x_sym = install("x");
  mx = length(X); /* list length */
  if (mx==1) return(VECTOR_ELT(X, 0));
  /* read sparse matrix list into structure suitable for passing on to 
     other routines */
  Xs = (spMat *) CALLOC((size_t)mx,sizeof(spMat)); // sparse matrix list
  //off_off=0;
  for (p=1,i=0;i<mx;i++) { // work through the sparse matrices
    M = VECTOR_ELT(X, i); // the ith sparse matrix
    Xs[i].x = REAL(R_do_slot(M,x_sym)); // the elements of ith matrix
    Xs[i].p = INTEGER(R_do_slot(M,p_sym)); // .p[j] is where jth col starts in .i and .x
    Xs[i].i = INTEGER(R_do_slot(M,i_sym)); // row index for each non-zero element
    dim = INTEGER(R_do_slot(M,dim_sym));
    Xs[i].m = dim[0];Xs[i].c = dim[1]; // matrix is .m by .c
    p *= dim[1]; // getting total column count
  } /* sparse matrix array filling loop */
  n = Xs[mx-1].m;
  c = (int *)CALLOC((size_t)mx,sizeof(int)); /* c[i] indexes which col of Xs[i] is corrently being processed */
  pp = (double *)CALLOC((size_t)mx*n,sizeof(double)); /* partial column products */
  dn = (int *)CALLOC((size_t)n,sizeof(int)); /* the non-zero tracker */
  R =  PROTECT(R_do_new_object(PROTECT(R_getClassDef("dgCMatrix")))); /* the result sparse matrix */
  dim = INTEGER(R_do_slot(R,dim_sym));dim[0] = n;dim[1] = p; /* set dimensions to n by p */
  R_do_slot_assign(R,p_sym,allocVector(INTSXP,p+1)); /* column starts vector */
  Rp = INTEGER(R_do_slot(R,p_sym));
  /* The total entries in a sparse row Kronecker product is not known in advance (unlike a Kronecker product) 
     The following therefore does two passes, the first to compute size and the second to compute the product.
     It might be faster to usa a realloc and copy approach instead (there does not seem to be an equivalent to 
     realloc for allocVector).
  */
  for (op=0;op<2;op++) { /* op==0 to count, 1 to output */
    k = 0; /* first matrix in list to have experienced a column change */
    for (l=0;l<mx;l++) c[l] = 0;
    if (op) { /* assign storage counted on previous (op==0) pass */
      R_do_slot_assign(R,x_sym,allocVector(REALSXP,rj));
      R_do_slot_assign(R,i_sym,allocVector(INTSXP,rj));
      Ri = INTEGER(R_do_slot(R,i_sym));
      Rx = REAL(R_do_slot(R,x_sym));
    }
    rj = 0; /* location in output arrays */ 
    for (l=0;l<p;l++) { /* loop over cols of result */
      Rp[l] = rj; /* starting entry of col l in output sparse matrix R */
      for (i=k;i<mx;i++) { /* marginal loop */
	ip = Xs[i].p + c[i] + 1; /* points to next col start for Xs[i] */
	pp1 = pp + i * n; /* points to ith partial product array */
	if (i) pp0 = pp + (i-1)*n; /* points to (i-1)th partial product array */
        for (j = *(ip-1);j<*ip;j++) { /* loop over NZ of c[i] of Xs[i] */
	  ii = Xs[i].i[j]; /* row of product */
          if (dn[ii]==i) { /* not yet zero */
            dn[ii]++; /* and still not zero */
            if (op) {
	      if (!i) pp[ii] = Xs[i].x[j]; else if (i<mx-1) pp1[ii] = pp0[ii]*Xs[i].x[j]; else { /* write directly into result */
                Rx[rj] = pp0[ii]*Xs[i].x[j];
		Ri[rj] = ii;
		rj++;
              }
	    } else if (i==mx-1) rj++; /* just counting */
	  } /* column not zero */ 
        } /* loop over NZ of c[i] of Xs[i] */	
      } /* marginal loop */
      /* now update the column counters - when a counter hits its limit it is 
         reset to zero and the preceding column counter is updated. Each counter 
         update requires a corresponding down data of dn, so that it ends up in 
         the state it was in before processing of the column we are moving on from */ 
      k = mx-1;
      /* clear previous dn additions for this col */
      ip = Xs[k].p + c[k] + 1; /* points to next col start for Xs[k] */
      for (j= *(ip-1);j<*ip;j++) { ii = Xs[k].i[j]; if (dn[ii]==k+1) dn[ii]--;}
      c[k]++;
      while (c[k]==Xs[k].c) { /* update the column counter */
        c[k]=0;
        if (k>0) {
	  k--;
	  ip = Xs[k].p + c[k] + 1; /* points to next col start for Xs[k] */
	  for (j= *(ip-1);j<*ip;j++) { ii = Xs[k].i[j]; if (dn[ii]==k+1) dn[ii]--;}
	  c[k]++;
        } 
      }
      /* adjust dn */
    } /* result column l loop */
    Rp[l] = rj;
  } /* op loop */
  FREE(Xs);FREE(pp);FREE(dn);FREE(c);
  UNPROTECT(2);
  return(R);
} /* stmm */


SEXP AddBVB(SEXP A,SEXP bt, SEXP vbt) {
/* A is a column compressed sparse matrix (symmetric or otherwise)
   Bt and VBt are p by n dense matrices. Forms
   A + t(Bt) %*% VBt on NZP(A). Modifies A@x only.
*/
  SEXP i_sym,x_sym,dim_sym,p_sym;
  int *dim,*Ai,*Ap,n,p,i,j,k;
  double *Ax,*B,*VB,x,*b,*b1,*vb;
  p_sym = install("p");
  dim_sym = install("Dim");
  i_sym = install("i");
  x_sym = install("x");
  dim = INTEGER(R_do_slot(A,dim_sym));
  n = dim[0];
  Ap = INTEGER(R_do_slot(A,p_sym));
  Ai = INTEGER(R_do_slot(A,i_sym));
  Ax = REAL(R_do_slot(A,x_sym));
  B = REAL(bt);
  p = nrows(bt);
  VB = REAL(vbt);
  for (j=0;j<n;j++) { /* loop over cols of A */
    for (k=Ap[j];k<Ap[j+1];k++) {
      i = Ai[k]; /* current row to A */
      x = 0.0;
      for (b=B+i*p,vb=VB+j*p,b1=b+p;b<b1;b++,vb++) x += *b * *vb;
      Ax[k] += x;
    }  
  }
  return(R_NilValue);
} /* AddBVB */



/*** inverse subset algorithm ***/

static inline int kij(int *Ap,int *Ai,int i, int j) {
/* find location of A[i,j] in Ax */ 
  int k0,k1,kt;
  k0 = Ap[j];k1 = Ap[j+1]-1;
  if (Ai[k0]==i) return(k0);
  if (Ai[k1]==i) return(k1);
  while(1) {
    //global_ops++;
    kt = (k0+k1)/2;
    if (Ai[kt] == i) return(kt);
    if (Ai[kt] < i) k0 = kt; else k1 = kt;
  }
} /* kij */


SEXP isa1p(SEXP L,SEXP S,SEXP NT) {
/* Parallel version. Fine scale single column parallelization.

   Inverse subset algorithm adapted from Rue (2005). Idea is to compute the 
   elements of S = (LL')^{-1} on the non-zero Pattern on L+L'. On entry 
   L and S are of class "dgCMatrix" and S has the non-zero pattern of L+L'.
   On exit S contains the required elements of the inverse.
   The rate limiting step is the search for the matching non-zero elements in 
   the intermost summation loop. This version uses an algorithm that 
   simultaneously finds search brackets in S[,j] for each element in L[,i] 
   (always the smaller) in the summation. Bisection is then used within the 
   search bracket. It's only about 10% faster than simple bisection in reality!
*/
  
  SEXP i_sym,x_sym,dim_sym,p_sym,kr;
  int *Lp,*Li,*Sp,*Si,i,j,k,q,*dim,s,k0,k1,l0,l1,n,mm,
    *li0,*li1,s0,s1,s2,m,*ul,*ll,*ul0,*ll0,kk,*ulq,*llq,*llq1,*liq,nt,tid;
  //long long int ops=0; 
  double *Lx,*Sx,x=0.0,Lii,*lxp;
  //global_ops=0;
  /* register the names of the slots in X and XWX */
  p_sym = install("p");
  dim_sym = install("Dim");
  i_sym = install("i");
  x_sym = install("x");
  nt = asInteger(NT); 
  /* Get pointers to the relevant slots in L*/ 
  Lp = INTEGER(R_do_slot(L,p_sym));
  dim = INTEGER(R_do_slot(L,dim_sym));
  n=dim[1];
  Li = INTEGER(R_do_slot(L,i_sym));
  Lx = REAL(R_do_slot(L,x_sym));

  /* Now get pointers to slots in S */
  Sp = INTEGER(R_do_slot(S,p_sym)); /* col j's elements lie between Sp[j] and Sp[j+1]-1 in x*/
  Si = INTEGER(R_do_slot(S,i_sym)); /* row index corresponding to elements in x*/
  Sx = REAL(R_do_slot(S,x_sym)); /* non-zero elements */

  /* we need the maximum column length in L... */
  for (mm=0,i=0;i<n;i++) {
    j = Lp[i+1] - Lp[i]; if (j>mm) mm=j;
  }  
  /* allocate storage for the search interval limits... */
  ll0 = ll = (int *)CALLOC((size_t)mm*nt,sizeof(int));
  ul0 = ul = (int *)CALLOC((size_t)mm*nt,sizeof(int));
  
  for (i=n-1;i>=0;i--) { /* work down columns */
    Lii = Lx[Lp[i]]; /* Lii is first element in ith col of L */
    l0 = Lp[i]+1;l1 = Lp[i+1]; /* limits of L[,i] over which to sum */
    li0 = Li + l0; li1 = Li + l1; /* pointers to row indices */
    k1 = Sp[i+1]-1; // k0 = Sp[i];
    k0 = kij(Sp,Si,i,i); /* get index for S[i,i] as we do not want to go further than this */
    /* in fact elements j/k can be processed in any order, so openMP
       parallel section could be used here. However calc is not 
       block oriented */
    tid = 0;
    #ifdef _OPENMP
#pragma omp parallel private(k,tid,ul,ll,j,m,q,s,s1,kk,ulq,llq,liq,llq1,s2,s0,x,lxp) num_threads(nt)
    #endif 
    { /* there seems to be little in it between blocking and by column threading */
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (k=k1;k>k0;k--) { /* work up rows of col i */
      #ifdef _OPENMP
      tid = omp_get_thread_num();
      #endif
      ul = ul0 + tid*mm;
      ll = ll0 + tid*mm;
      j = Si[k]; /* the row corresponding to stored element k */ 
      /* Now compute S[i,j] which is S[j,i] which is in Sx[k] */
      /* Loop over the ith column of L */
      m = l1-l0; /* number of non zero elements in L[,i] */
      if (m>0) {
	s = kij(Sp,Si,*li0,j); /* location of Sp[Li[l0],j] */
        s1 = kij(Sp,Si,li1[-1],j); /* location of Sp[Li[l1-1],j] */
      }	else s1=s=Sp[j];
      for (q=0;q<m;q++) { /* fill out initial search bracket ends */
        ul[q] = s1;ll[q] = s;
      }
     	
      kk = 0; /* interval we are working on */

      while (kk<m-1) { /* iterate for bracketing intervals */
        s = (ll[kk] + ul[kk])/2; /* current interval mid-point */
	s1 = Si[s]; /* row number at s */

	for (q=kk;q<m;q++) { /* work through remaining intervals */
          if (s1 > li0[q]) { /* is new point above point q? */
            if (s<ul[q]) ul[q] = s; /* is new point closer than previous upper limit ? */
          } else { /* new point is below or equal to point q */
            if (s>ll[q]) ll[q] = s; else break; /* either it is a higher lower limit, or we can stop updating */
          }   
	}
	/* now test whether interval kk is complete and we can move on... */
	if (ul[kk] <= ll[kk+1]||ul[kk] == ll[kk]+1) kk++; 
      }	/* bracketing interval loop */			       
      /* at this stage we have non overlapping intervals within which to search for the
         elements of S matching the non-zeores of L[,i], can now finish the matching 
         by simple bisection within each interval. */
      ulq=ul;llq=ll;liq=li0;llq1=ll+m;
      for (x=0.0,lxp=Lx+l0;llq<llq1;llq++,ulq++,liq++,lxp++) { 
        s = *llq;s2 = *ulq;kk = *liq;
	while (Si[s] != kk) { /* search for S element corresponding to L element */
          s0 = (s+s2+1)/2;
	  if (Si[s0] > kk) s2=s0; else s=s0;
	}
	x -=  *lxp * Sx[s];
      }	

      x /= Lii;//ops++;
      Sx[k] = x; /* S[j,i] */
      s = kij(Sp,Si,i,j);
      Sx[s] = x;/* S[i,j] */
    } /***** k loop end *****/
    } /* parallel section end */ 
    /* now do k0, to fill in S[i,i] */
    ul=ul0;ll=ll0;
      
    /* Now compute S[i,j] which is S[j,i] which is in Sx[k] */
    /* Loop over the ith column of L */
    m = l1-l0; /* number of non zero elements in L[,i] */
    if (m>0) {
      s = kij(Sp,Si,*li0,i); /* location of Sp[Li[l0],j] */
        s1 = kij(Sp,Si,li1[-1],i); /* location of Sp[Li[l1-1],j] */
    }  else s1=s=Sp[i];
    for (q=0;q<m;q++) { /* fill out initial search bracket ends */
      ul[q] = s1;ll[q] = s;
    }
     	
    kk = 0; /* interval we are working on */

    while (kk<m-1) { /* iterate for bracketing intervals */
      s = (ll[kk] + ul[kk])/2; /* current interval mid-point */
      s1 = Si[s]; /* row number at s */

      for (q=kk;q<m;q++) { /* work through remaining intervals */
        if (s1 > li0[q]) { /* is new point above point q? */
          if (s<ul[q]) ul[q] = s; /* is new point closer than previous upper limit ? */
        } else { /* new point is below or equal to point q */
          if (s>ll[q]) ll[q] = s; else break; /* either it is a higher lower limit, or we can stop updating */
        }   
      }
      /* now test whether interval kk is complete and we can move on... */
      if (ul[kk] <= ll[kk+1]||ul[kk] == ll[kk]+1) kk++; 
    }	/* bracketing interval loop */			       
    /* at this stage we have non overlapping intervals within which to search for the
       elements of S matching the non-zeores of L[,i], can now finish the matching 
       by simple bisection within each interval. */
    ulq=ul;llq=ll;liq=li0;llq1=ll+m;
    for (x=0.0,lxp=Lx+l0;llq<llq1;llq++,ulq++,liq++,lxp++) { 
      s = *llq;s2 = *ulq;kk = *liq;
      while (Si[s] != kk) { /* search for S element corresponding to L element */
        s0 = (s+s2+1)/2;
	if (Si[s0] > kk) s2=s0; else s=s0;
      }
      x -=  *lxp * Sx[s];
    }	
    x += 1/Lii;
    x /= Lii;//ops++;
    Sx[k0] = x; /* S[j,i] */
    
  }
  FREE(ul0);FREE(ll0);
  PROTECT(kr=allocVector(REALSXP,1));
  REAL(kr)[0] = 0.0;//ops/(double) global_ops;
  UNPROTECT(1);
  return(kr);
} /* isa1p */

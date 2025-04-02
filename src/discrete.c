/* (c) Simon N. Wood (2015-2024) Released under GPL2 */

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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mgcv.h"
#include <R.h>

struct SM_el { /* stores element of sparse matrix */
  int i,j; /* row and column */
  double w; /* entry */
  struct SM_el *next; /* pointer to next record */
};  

typedef struct SM_el SM;

void SMinihash(unsigned long long *ht) {
/* initialize hash generator ht (length 256),
   used in converting matrix row column pairs to hash keys.    
*/ 
  unsigned long long h;
  int j,k;
  h = 0x987564BACF987454LL;
  for (j=0;j<256;j++) {
    for (k=0;k<31;k++) {
        h = (h>>7)^h;
        h = (h<<11)^h;
        h = (h>>10)^h;
    }
    ht[j] = h;
  }
} /* SMinihash */ 

void indReduce(int *ka,int *kb,double *w,int tri,int *n,
	       unsigned long long *ht,SM **sm,SM * SMstack,double *C,double *A,
	       int rc,int cc,int ra,int trans,int *worki,int buffer) {
/* on input ka, kb and w are n-vectors. Let W be a matrix initialized to zero. w[i] 
   is to be added to element ka[i], kb[i] of W. If tri!=0 then ws[i] is added to element 
   ka[i], kb[i+1] and wl[i] to ka[i+1], kb[i]. ws is stored at w + n and wl at ws + n;

   C is rc by cc and A is ra by cc. If trans!=0 then form C+=W'A otherwise C+=WA 

   On entry SMstack should be an n-vector if tri==0 and a 3n-vector otherwise. 
   sm is an n vector.

   This routine accumulates W in a sparse, i,j,w structure constructed using a hash table. 
   After accumulation  the hash table contains n_u <= n unique matrix entries, which can then 
   be used directly to form the matrix product.       

   Accumulation cost is O(n). Access cost is then O(n_u) while matrix product cost is O(n_u * cc). 
   In comparison direct accumulation costs would be O(n*cc). 

   If buffer!=0 then the routine will access worki (dimension 6 * n) and will modify w. It does 
   this because it massively improves data locality and cache performance to read the sparse matrix 
   out of the hash table structure into 3 arrays, before using it for multiplication.
*/
  SM **sm1, **sm2,*smk;
  int bpp,nstack,nkey,ij[2],k,l,i,j,t,*kao,*kbo;
  char *key,*keyend;
  unsigned long long h;
  double Wij,*Cq,*Aq,*Cq1,*ws,*wl;
  if (tri) { ws = w + *n;wl = ws + *n;} else {ws=wl=w; /* avoid compiler warning */}
  bpp = sizeof(int)/sizeof(char);
  nkey = bpp * 2; /* length of key in char */ 
  if (tri) nstack = 3 * *n-1; else nstack = *n-1; /* top of the SM element stack */
  /* clear the hash table */
  for (sm1=sm,sm2=sm + *n;sm1<sm2;sm1++) *sm1 = (SM *)NULL;
  if (tri) tri=3; else tri=1;
  for (l=0;l<*n;l++) {
    for (t=0;t<tri;t++) {
      if (!t) { /* leading diagonal */
	Wij = w[l];i=ka[l];j=kb[l];
      } else if (t==1) { /* super */
        i=ka[l];j=kb[l+1];Wij = ws[l];
      } else {
	i=ka[l+1];j=kb[l];Wij = wl[l];
        if (l == *n-2) tri=1; /* no sub-diagonals beyond here */
      }	
      /* add Wij to row i col j */
      /* generate hash key */
      ij[0]=i;ij[1]=j;
      key = (char *) ij;
      h = 0x99EFB145DAA48450LL;
      for (keyend = key + nkey;key<keyend;key++)
      h = (h * 7664345821815920749LL)^ht[(unsigned char)(*key)];
      k = (int) (h % *n); /* get location in hash table */
      //Rprintf("k=%d, i=%d, j=%d, w=%g\n",k,i,j,Wij);
      if (sm[k]) { /* already something in this slot */
        smk = sm[k];
        /* traverse linked list looking for a match */
        while (smk) { /* not at end of linked list */
          if (smk->i == i && smk->j==j) {
            smk->w += Wij;break;
          }
          smk = smk->next; /* next element in list */
        }
        if (!smk) { /* no match found - add to linked list */
          smk = SMstack + nstack;nstack--; /* remove an entry from the stack */
          smk->next = sm[k];sm[k]=smk;smk->i = i;smk->j = j;smk->w=Wij;
        }  
      } else { /* slot was empty */
        smk = sm[k] = SMstack + nstack;nstack--; /* remove an entry from the stack */
        smk->i = i;smk->j = j;smk->w=Wij;smk->next = (SM *)NULL; 
      }
    } /* t loop */
  } /* l loop */
  /* Now form C = \bar W A -- There is a problem here - the access to C and  A are cache inefficient ---
     The only solutions to this are either to read out the elements of the linked list to temporary arrays,
     which allows working over cols to be outer to working though the list, or transposing C and A so that 
     they are in row-major order. 
  */
  if (buffer) { /* read data structure out to arrays */
    ka = kao = worki; kb = kbo= worki + 3 * *n;ws=w;
    for (k=0,sm1=sm,sm2=sm + *n;sm1<sm2;sm1++) if (*sm1) { /* there is something in this slot */
      smk = *sm1; /* base of linked list at this slot */
      while (smk) { /* traverse the linked list writing the compressed ka,kb,w vectors coding \bar W*/
        //kao[k] = smk->i;kbo[k] = smk->j;w[k] = smk->w;k++;smk = smk->next;
	*ka = smk->i;*kb = smk->j;*ws = smk->w;k++;ka++;kb++;ws++;smk = smk->next;
      }  
    }
    if (trans) for (Aq=A,Cq = C,Cq1=C+rc*cc;Cq<Cq1;Cq+=rc,Aq+=ra)
      for (wl=w,ws=w+k,ka=kao,kb=kbo;wl<ws;wl++,ka++,kb++) Cq[*kb] += Aq[*ka] * *wl;
    else for (Aq=A,Cq = C,Cq1=C+rc*cc;Cq<Cq1;Cq+=rc,Aq+=ra)
      for (wl=w,ws=w+k,ka=kao,kb=kbo;wl<ws;wl++,ka++,kb++) Cq[*ka] += Aq[*kb] * *wl;

  } else /* work directly from hash table and linked list */
    for (k=0,sm1=sm,sm2=sm + *n;sm1<sm2;sm1++) if (*sm1) {  /* there is something in this slot */
    smk = *sm1; /* base of linked list at this slot */
     while (smk) { /* traverse the linked list extracting matrix elements */
       if (trans) { j = smk->i;i = smk->j;} else { i = smk->i;j = smk->j;}
       Wij = smk->w;smk = smk->next;
       for (Cq = C + i,Aq = A + j,Cq1 = C + rc*cc;Cq<Cq1;Cq += rc,Aq += ra) *Cq += Wij * *Aq;
     }
  }   
} /* indReduce */

void Cdgemv(char *trans, int *m, int *n, double *alpha, double *a, int *lda,
	    double *x, int *incx, double *beta, double *y, int *incy) {
/* miserably vanilla implementation of dgemv to plug in in place of BLAS call for
   debugging purposes (written to examine openblas thread safety problem) */
  int i,j,q;
  double *yp,*ap,*xp;
  if (*trans == 'T') q = *n; else q =*m; /* length of y */
  if (*alpha==0) { /* matrix part does not contribute */
    for (i=0;i<q;i++,y += *incy) *y *= *beta;
    return;
  } else *beta /= *alpha;
  if (*trans == 'N') { /* y <- alpha a x + beta y */
    for (yp=y,i=0;i<*m;i++,a++,yp += *incy) *yp = *beta * *yp + *a * *x;
    x += *incx;
    for (j=1;j<*n;j++,x += *incx) 
      for (ap=a + *lda * j,yp=y,i=0;i<*m;i++,ap++,yp += *incy) *yp += *ap * *x;
  } else { /* y <- alpha a' x + beta y */
    for (yp=y,i=0;i<*n;i++,yp++) {
      *yp *= *beta;
      for (ap=a + *lda * i,xp=x,j=0;j<*m;j++,ap++,xp += *incx) *yp += *ap * *xp;
    }  
  }
  for (i=0;i<q;i++,y += *incy) *y *= *alpha;
} /* Cdgemv */ 


/* constraint application helper functions */

void Zb(double *b1,double *b0,double *v,int *qc, int *p,double *w) {
/* Form b1 = Z b0 where constraint matrix Z has more rows than columns. 
   b1 and b0 must be separate storage.
   p is dim(b1)
   qc > 0 for simple HH sum to zero constraint.
   qc < 0 for Kronecker product of sum to zero contrasts.
   v is either the householder vector, or a vector giving: 
   [the number of contrasts, the leading dimension of each, the total number of 
    linear constraints implied]
   w is work space at least 2p long in the qc<0 case, but unused otherwise.
      
*/
  double x,*p0,*p1,*p2,*p3,*w0,*w1,z;
  int M,k0,pp,i,j,k,q,mk,p0m;
  if (*qc>0) { 
    *b1 = 0.0;x=0.0;
    for (p0 = b1+1,p1=b1+ *p,p2 = b0,p3=v+1;p0<p1;p0++,p2++,p3++) {
      *p0 = *p2; x += *p0 * *p3; /* v'(0,b0')' */
    }  
    for (p0=b1,p2=v;p0<p1;p0++,p2++) *p0 -= *p2 * x; /* (I-vv')(0,b0')' */
  } else if (*qc < 0){ /* Z is sequential Kronecker product of M sum-to-zero contrasts and an identity matrix */ 
    M = (int)round(v[0]); /* number of sum to zero contrasts making up contrast */
    for (pp = *p,k0=1,i=0; i<M;i++) {
      mk = (int) round(v[i+1]); /* stz contrast dimension */
      k0 *= mk - 1;
      pp /= mk; /* dimension of final I */
    }
    k0 *= pp; /* dim b0 */
    w1 = w + *p;w0=w;
    for (k=0;k<=M;k++) {
      if (k<M) mk = (int) round(v[k+1])-1; else {
	mk = pp;
	w1=b1; /* final w1 is b1 */ 
      }	
      p0m = k0/mk;q = 0;
      for (i=0;i<p0m;i++) {
        x=0.0;
	for (j=0;j<mk;j++) {
          z = b0[j*p0m+i];
	  x += z;
	  w1[q] = z; q++;
	}
	if (k<M) {
	  w1[q] = -x; q++;
	}  
      }
      if (k<M) k0 += p0m;
      b0 = w1;w1=w0;w0=b0;
    }
    //for (i=0;i<p;i++) b1[i] = b0[i]; //- redundant?
  } /* Kronecker stz end */
} /* Zb */  

void Ztb(double *b1,double *b0,double *v,int *qc,int *di, int *p,double *w) {
/* Apply constraint to vector p-vector b0. b1 = Z'b0. b0 and b1 may be spaced 
   within arrays b0 and b1. b0[0:(p-1)*di] and b1[0:(p-k)*di]  are the elements 
   to be operated on (where k-1 is the number of constraints). 
   This facilitates computation of b0'Z where b0' is the row of a matrix. 
   If qc==1 then v is the p-vector defining a Housholder rotation (I-vv') and 
   Z'b0 is (I-vv')b0 with the first element dropped.
   if (qc == -1) then the constraint is the Kronecker product of a sequence of 
   sum-to-zero contrasts and a final identity matrix. In that case v[0] is the 
   number of sum to zero contrasts, the ith of which is v[1+i] by (v[1+i]-1).
   NOTE: b1 and b0 may point to same array - care needed in implementing new 
         constraints.
   w is workspace and should be double the length of b0. It is not required 
   and not referenced if qc==1.
*/
  double z,x,*p0,*p1,*p2,*w0,*w1;
  int k1,M,i,j,q,k,p1m,mk,pp,mk1;
  if (*qc>0) { /* signals a rank 1 HH constraint */
    for (x=0.0,p0=b0,p1=v,p2=v + *p;p1<p2;p0 += *di,p1++) x += *p0 * *p1; /* b0'v */
    for (p1=v+1,b0 += *di;p1<p2;b0 += *di,b1 += *di,p1++) *b1 = *b0 - *p1 * x;
  } else if (*qc<0) { /* signals a constraint based on sum to zero contrasts */
    /* Z is a Kronecker product of a sequence of sum-to-zero contrasts and an identity matrix */
    for (p0=w,p1=w + *p,p2=b0;p0<p1;p0++,p2 += *di) *p0 = *p2; /* copy b0 to w */
    M = (int) round(v[0]); /* number of sum to zero contrasts */
    for (pp = *p,i=0; i<M;i++) {
      mk =  (int) round(v[i+1]);
      pp /= mk;
    } 
    /* pp is dimension of final identity matrix in constraint */
    k1 = *p;w0 = w;w1 = w + *p;
    for (k=0;k<=M;k++) {
      if (k<M) { mk = (int) round(v[k+1]);mk1 = mk-1;} else mk1=mk = pp;
      p1m = k1/mk;
      q = 0;
      for (i=0;i<p1m;i++) {
        if (k<M) z = w0[i+(mk-1)*p1m]; else z = 0.0;
	for (j=0;j<mk1;j++) {
	  w1[q] = w0[i+j*p1m] - z; q++;
	}  
      }
      if (k<M) k1 -= p1m;
      p0 = w0;w0 = w1;w1 = p0;
    }
    for (p0=w0,p1=w0+k1,p2=b1;p0<p1;p0++,p2 += *di) *p2 = *p0; /* copy w0 to b1 */
  } /* sum-to-zero Kronecker constraints */
} /* Ztb */  


/* basic extraction operations */ 

void singleXj(double *Xj,double *X, int *m, int *k, ptrdiff_t  *n,int *j) {
/* Extract a column j of matrix stored in compact form in X, k into Xj. 
   X has m rows. k is of length n. ith row of result is Xj = X[k(i),j]
   (an n vector). This function is O(n). Thread safe.
*/
  double *pe; 
  X += *m * *j; /* shift to start of jth column */ 
  for (pe = Xj + *n;Xj < pe;Xj++,k++) *Xj = X[*k]; 
} /* singleXj */

void tensorXj(double *Xj, double *X, int *m, int *p,int *dt, 
              int *k, ptrdiff_t *n, int *j, int *kstart,int *koff) {
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
    kp = k + (kstart[i] + *koff) * *n;
    for (p0=Xj;p0<p1;p0++,kp++) *p0 *= M[*kp];
    X += m[i] * p[i]; /* move to the next marginal matrix */
  }
} /* tensorXj */

void singleXty(double *Xy,double *temp,double *y,double *X, int *m,int *p, int *k, ptrdiff_t  *n, int *add) {
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
  F77_CALL(dgemv)(&trans, m, p,&done,X,m,temp,&one,&dzero,Xy,&one FCONE);
  /* dgemm call equivalent to dgemv call above */ 
  //F77_CALL(dgemm)(&trans,&ntrans,p,&one,m,&done,X,m,temp,m,&dzero,Xy,m); 
} /*singleXty*/

void tensorXty(double *Xy,double *work,double *work1, double *y,double *X, 
               int *m, int *p,int *dt,int *k, ptrdiff_t *n, int *add, int *kstart,int *koff) {
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
    singleXty(Xy+i*pd,work1,work,M,m + *dt-1,&pd, k +  *n * (kstart[j] + *koff) ,n,add);
  }
} /* tensorXty */

void singleXb(double *f,double *work,double *X,double *beta,int *k,int *m, int *p,ptrdiff_t *n,
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
  F77_CALL(dgemv)(&trans, m, p,&done,X,m,beta,&one,&dzero,work,&one FCONE);
  /* dgemm call equivalent to dgemv call above */ 
  //F77_CALL(dgemm)(&trans,&trans,m,&one,p,&done,X,m,beta,p,&dzero,work,m); 
  p1 = f + *n;
  fp = f;
  k += *kstart * *n;
  for (;fp < p1;fp++,k++) *fp = work[*k];
  for (j=1;j < *kstop - *kstart;j++) for (fp=f;fp < p1;fp++,k++) *fp += work[*k];
} /* singleXb */

void tensorXb(double *f,double *X, double *C,double *work, double *beta,
              int *m, int *p,int *dt,int *k, ptrdiff_t *n,double *v,int *qc,
              int *kstart, int *kstop) {
/* for X* beta where X* is a tensor product term with nt marginal model matrices
   stored in X. ith such is m[i] by p[i]. 
   work is an n vector. C is an m[d] by pb working matrix.  
   v is a vector such that if Q=I-vv' and Z is Q with the first column dropped then 
   Z is a null space basis for the identifiability constraint. If *qc == 0 then
   no constraint is applied. Thread safe.
   k[,kstart[i]...kstop[i]-1] are the indices for the ith component of this term 
   (starting at i=0).
*/
  char trans='N';
  int pb=1,md,*kp,*kd,pd,i,j,q;
  double *M,done=1.0,dzero=0.0,*p0,*p1,*pf,*pc;
  M = X;
  for (i=0;i<*dt-1;i++) {
    pb *= p[i];
    M += m[i] * (ptrdiff_t) p[i]; /* final marginal */ 

  }
  md = m[*dt - 1];
  pd = p[*dt -1];
  kd = k + kstart[*dt-1] *  *n;
  /* form work = M B, where vec(B) = beta */
  if (*qc==0) { /* no constraint supplied */
    F77_CALL(dgemm)(&trans,&trans,&md,&pb, &pd, &done,
		    M,&md,beta,&pd,&dzero,C,&md FCONE FCONE); 
  } else { /* there is a constraint matrix */
    /* first map supplied beta to unconstrained parameterization */ 
    j = pb * pd; /* total number of coeffs - length of unconstrained beta */
    Zb(work,beta,v,qc,&j,work+j);
    F77_CALL(dgemm)(&trans,&trans,&md,&pb, &pd, &done,
		    M,&md,work,&pd,&dzero,C,&md FCONE FCONE);
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
    for (kp=kd + q  * *n,pf=f,p0=work;p0<p1;p0++,pf++,kp++) {
      *pf += pc[*kp] * *p0;
    }
  }
} /* tensorXb */



SEXP CXbd(SEXP fr, SEXP betar, SEXP Xr, SEXP kr, SEXP ksr, SEXP mr, SEXP pr,
	  SEXP tsr, SEXP dtr,SEXP vr,SEXP qcr,SEXP bcr,SEXP csr) {
/* .Call wrapper for Xbd allowing R long vector storage for k. Note that 
  this does not allow more than maxint data - that would require re-writting R code 
  to avoid storing k in a matrix (which is only allowed maxint rows).

  n is length of diag, nx is length of m or p, nt is the length of ts or dt, 
  ncs is length of cs.
*/
  double *f,*X,*v,*beta,*dwork;
  int *k,*ks,*m,*p,*ts,*dt,*qc,*cs,nx,nt,ncs,*bc,*iwork;
  ptrdiff_t n,space[3],*pwork;
  n = (ptrdiff_t) nrows(kr); f = REAL(fr);
  beta = REAL(betar);
  X = REAL(Xr);
  k = INTEGER(kr); ks = INTEGER(ksr);
  m = INTEGER(mr); nx = length(mr);
  p = INTEGER(pr);
  ts = INTEGER(tsr); dt = INTEGER(dtr); nt = length(tsr);
  v = REAL(vr);qc = INTEGER(qcr);
  bc = INTEGER(bcr);
  cs = INTEGER(csr);ncs = length(csr);
  Xbdspace(space,m,p,&n, &nx,dt,&nt);
  iwork = (int *) CALLOC((size_t) space[0],sizeof(int));
  pwork = (ptrdiff_t *) CALLOC((size_t)space[1],sizeof(ptrdiff_t));
  dwork = (double *)CALLOC((size_t)space[2],sizeof(double));
  Xbd(f,beta,X,k,ks,m,p,&n,&nx,ts,dt,&nt,v,qc,bc,cs,&ncs,iwork,pwork,dwork);
  FREE(iwork);FREE(pwork);FREE(dwork);
  return(R_NilValue);
} /*  CXbd */

void Xbdspace(ptrdiff_t *space,int *m,int *p, ptrdiff_t *n, int *nx, int *dt, int *nt) {
  /* computes space needed by Xbd returning in 3 vector space [int,ptrdiff_t,double] */
  int i,q,j,c1,dC=0,maxrow=0,maxp=0,pti=0;
  ptrdiff_t ii;
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) { /* work through components of each term */
      if (maxrow<m[q]) maxrow=m[q];
      if (j>0 && j==dt[i]-1) {
        c1 = pti * (ptrdiff_t) m[q]; 
        if (c1>dC) dC = c1; /* dimension of working matrix C */
      }
      if (j==0) pti = p[q]; else pti *= p[q]; /* term dimension */
    } 
    if (maxp < pti) maxp = pti;
  }
  space[0] = *nt *2 + 1; /* int */
  space[1] = *nx + *nt + 2; /* ptrdiff_t */
  ii = *n; if (ii<3*maxp) ii= (ptrdiff_t) 3*maxp; if (ii<maxrow) ii=(ptrdiff_t)maxrow;
  if (!dC) dC=1; /* avoid pointing to unallocated */
  space[2] = ii + *n + dC; /* double */
} /* Xbdspace */


void Xbd(double *f,double *beta,double *X,int *k,int *ks, int *m,int *p, ptrdiff_t *n, 
	 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *bc,int *cs,int *ncs,
	 int *iwork, ptrdiff_t *pwork,double *dwork) {
/* Forms f = X beta for X stored in the packed form described in function XWX
   bc is number of cols of beta and f... 
   This version allows selection of terms via length ncs array cs. If ncs<=0 then all terms
   are returned (but the physical storage for cs must then be of length nt).
   BUT: the selection mechanism is clunky as it assumes that beta has been re-ordered and 
        truncated to match cs, so that the working through the selected terms in order 
        we also only need to move through the supplied beta in order.
   No memory allocation, but Xbdspace must be called first and iwork,pwork and dwork 
   allocated accordingly.

   Note that if ncs <=0 then ncs and cs are reset to nt and 0:(nt-1). Requires care in parallel!
*/
  ptrdiff_t *off,*voff,ii;
  int i,j,q,*pt,*tps,first,kk;
  double *f0,*pf,*p0,*p1,*p2,*C=NULL,*work,maxp=0,maxrow=0;
  pt = iwork; /* the term dimensions */
  off = pwork; /* offsets for X submatrix starts */
  pf=f0 = dwork;
  tps = pt + *nt; /* the term starts in param vector */
  voff = off + *nx + 1;/* offsets for v subvector starts */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) { /* work through components of each term */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (maxrow<m[q]) maxrow=m[q];
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
    } 
    if (qc[i]==0) voff[i+1] = voff[i]; else if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else {  /* start of ith v matrix */
      kk = (int) round(v[voff[i]]); /* number of contrasts in this KP contrast */
      voff[i+1] = voff[i] + kk + 2; 
    }  
    if (maxp < pt[i]) maxp = pt[i];
  }
  if (*ncs<=0) { /* return everything */
    for (j=0;j<*nt;j++) cs[j] = j;
    *ncs = *nt;
  }  
  for (kk=j=0;j<*ncs;j++) { /* get the offsets for the returned terms in the output */
    i = cs[j];tps[i] = kk;
    if (qc[i]==0) kk += pt[i]; /* where cth terms starts in param vector */ 
    else if (qc[i]>0) kk += pt[i] - 1; else { /* there is a tensor constraint to apply - reducing param count*/
      /* Kronecker product of sum-to-zero contrasts */
      q = (int) round(v[voff[i]]); /* number of contrasts */
      kk += pt[i] - (int) round(v[voff[i]+q+1]); /* subtracting number of constraints */
    }
  }
  tps[*nt] = kk;
  /* now form the product term by term... */ 
  ii = *n; if (ii<3*maxp) ii=(ptrdiff_t)3*maxp; if (ii<maxrow) ii=(ptrdiff_t)maxrow;
  work = pf + *n;C = work + ii;
  for (j=0;j < *bc;j++) { /* loop over columns of beta */
    first = 1;
    for (kk=0;kk<*ncs;kk++) {
      i = cs[kk]; /* term to use */ 
        if (first) f0 = f; /* result written straight to f for i==0 */
        if (dt[i]==1) singleXb(f0,work,X+off[ts[i]],beta+tps[i],k /*+ (ptrdiff_t) *n *ts[i]*/,
			       m+ts[i], p+ts[i], n,ks + ts[i], ks + ts[i] + *nx); 
        else tensorXb(f0,X+off[ts[i]],C,work, beta+tps[i],m+ts[i], p+ts[i],dt+i,k 
		      /*+ (ptrdiff_t) *n * ts[i]*/, n,v+voff[i], qc+i,ks + ts[i], ks + ts[i] + *nx);
        if (!first) {
          for (p0=f,p1=f + *n,p2=f0;p0<p1;p0++,p2++) *p0 += *p2; /* f <- f + f0 */     
        } else { 
          f0=pf; /* restore f0 */
          first=0;
        }
    } /* term loop */
    f += *n;beta += tps[*nt]; /* move on to next column */
  } /* col beta loop */
} /* Xb */


SEXP CdiagXVXt(SEXP DIAG, SEXP Vp, SEXP x, SEXP K, SEXP KS, SEXP M, SEXP P, SEXP TS, SEXP DT,
	       SEXP vp,SEXP QC, SEXP NTHREADS, SEXP CS, SEXP RS) {
/* .Call wrapper for diagXVXt allowing R long vector storage for k. Note that 
  this does not allow more than maxint data - that would require re-writing R code 
  to avoid storing k in a matrix (which is only allowed maxint rows).

  n is length of diag, nx is length of m or p, nt is the length of ts or dt, pv is the
  number of rows of V, cv the number of cols. ncs and nrs are length of cs/rs.
*/
  double *diag,*V,*X,*v;
  int *k,*ks,*m,*p,*ts,*dt,*qc,*nthreads,*cs,*rs,nx,nt,pv,cv,ncs,nrs;
  ptrdiff_t n;
  n = (ptrdiff_t) length(DIAG); diag = REAL(DIAG);
  V = REAL(Vp);pv = nrows(Vp);cv = ncols(Vp);
  X = REAL(x);
  k = INTEGER(K); ks = INTEGER(KS);
  m = INTEGER(M); nx = length(M);
  p = INTEGER(P);
  ts = INTEGER(TS); dt = INTEGER(DT); nt = length(TS);
  v = REAL(vp);qc = INTEGER(QC);
  nthreads = INTEGER(NTHREADS);
  cs = INTEGER(CS); rs = INTEGER(RS);
  nrs = length(RS); ncs = length(CS);
  diagXVXt(diag,V,X,k,k,ks,m,p,&n,&nx,ts,dt,&nt,v,qc,&pv,&cv,nthreads,cs,&ncs,rs,&nrs);
  return(R_NilValue);
} /*  CdiagXVXt */

SEXP CijXVXt(SEXP DIAG, SEXP Vp, SEXP x, SEXP K,SEXP K1, SEXP KS, SEXP M, SEXP P, SEXP TS, SEXP DT,
	       SEXP vp,SEXP QC, SEXP NTHREADS, SEXP CS, SEXP RS) {
/* .Call wrapper for diagXVXt called for computing scattered elements XWX[i,j] allowing R long 
  vector storage for k and k1. Note that this does not allow more than maxint data - that would 
  require re-writing R code to avoid storing k in a matrix (which is only allowed maxint rows).

  n is length of diag, nx is length of m or p, nt is the length of ts or dt, pv is the
  number of rows of V, cv the number of cols. ncs and nrs are length of cs/rs.
*/
  double *diag,*V,*X,*v;
  int *k,*k1,*ks,*m,*p,*ts,*dt,*qc,*nthreads,*cs,*rs,nx,nt,pv,cv,ncs,nrs;
  ptrdiff_t n;
  n = (ptrdiff_t) length(DIAG); diag = REAL(DIAG);
  V = REAL(Vp);pv = nrows(Vp);cv = ncols(Vp);
  X = REAL(x);
  k = INTEGER(K); k1 = INTEGER(K1); ks = INTEGER(KS);
  m = INTEGER(M); nx = length(M);
  p = INTEGER(P);
  ts = INTEGER(TS); dt = INTEGER(DT); nt = length(TS);
  v = REAL(vp);qc = INTEGER(QC);
  nthreads = INTEGER(NTHREADS);
  cs = INTEGER(CS); rs = INTEGER(RS);
  nrs = length(RS); ncs = length(CS);
  diagXVXt(diag,V,X,k,k1,ks,m,p,&n,&nx,ts,dt,&nt,v,qc,&pv,&cv,nthreads,cs,&ncs,rs,&nrs);
  return(R_NilValue);
} /*  CdiagXVXt */

void diagXVXt(double *diag,double *V,double *X,int *k1,int *k2,int *ks,int *m,int *p, ptrdiff_t *n, 
	      int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pv,int *cv,int *nthreads,
	      int *cs,int *ncs,int *rs,int *nrs) {
/* If k1 and k2 are identical and k1[i] is storage row for ith obs, forms diag(XVX') where 
   X is stored in the compact form described in XWXd. Otherwise used to form selected elements
   of XVX'. e.g. k1[l] might be storage location of row/obs i and k2[l] the storage location 
   for row/obs j. Then lth element computed is XVX'[i,j].  

   V is a pv by pv matrix, if all terms are selected, otherwise pv by cv and V has to have been 
   manipulated to contain only elements relevant to retained terms. 
   
   Parallelization is by splitting the columns of V into nthreads subsets.
   Currently inefficient. Could be speeded up by a factor of 2, by supplying a
   square root of V in place of V.
   
   cs and rs are ncs and nrs vectors specifying which terms should be included 
   in the cross product. zero or negative ncs or nrs signals to include all, but then
   cs and rs must both be of length nt (Xbd requires this).

   Basic algorithm is to compute XV and then the row sums of XV.X. 
   In practice this is done by computing one column of XV and X at a time.
   When only some terms are selected then only those columns of XV corresponding 
   to selected columns of X need to be computed. The terms selected in rs determine
   which columns of X are required. 

   Called using diagXVXd in R
   library(mgcv)
   X <- list(matrix(runif(4*3),4,3),matrix(runif(8*4),8,4),matrix(runif(5*2),5,2))
   k <- cbind(sample(1:4,100,replace=TRUE),sample(1:8,100,replace=TRUE),sample(1:5,100,replace=TRUE))
   ks <- cbind(1:3,2:4); ts <- 1:3;dt <- rep(1,3)
   attr(X,"lpip") <- list(1:3,4:7,8:9)
   Xf <- cbind(X[[1]][k[,1],],X[[2]][k[,2],],X[[3]][k[,3],])
   V <- crossprod(matrix(runif(81)-.5,9,9))
   diagXVXd(X,V,k,ks,ts,dt,v=NULL,qc=rep(-1,3))
   rowSums((Xf%*%V)*Xf)

   lt <- c(1,3);rt <- c(2,3) 
   mt <- as.numeric(1:3 %in% lt)
   X1 <- cbind(X[[1]][k[,1],]*mt[1],X[[2]][k[,2],]*mt[2],X[[3]][k[,3],]*mt[3])
   mt <- as.numeric(1:3 %in% rt)
   X2 <- cbind(X[[1]][k[,1],]*mt[1],X[[2]][k[,2],]*mt[2],X[[3]][k[,3],]*mt[3])
   
   d1 <- diagXVXd(X,V,k,ks,ts,dt,v=NULL,qc=rep(-1,3),lt=lt,rt=rt)
   d2 <- rowSums((X1%*%V)*X2)
   range(d1-d2)

   (X1%*%V)%*%t(X2) -> XVX
   i <- c(2,4,8,25);j <- c(1,7,4,56)
   XVX[i+(j-1)*100]
   mgcv:::ijXVXd(i,j,X,V,k,ks,ts,dt,v=NULL,qc=rep(0,3),lt=lt,rt=rt)

   i <- sample(1:100,150,replace=TRUE);j <- sample(1:100,150,replace=TRUE)
   d1 <- XVX[i+(j-1)*100]
   d2 <- mgcv:::ijXVXd(i,j,X,V,k,ks,ts,dt,v=NULL,qc=rep(0,3),lt=lt,rt=rt)
   range(d1-d2)

   ## check Xbd
   beta <- runif(9)
   Xbd(X,beta,k,ks,ts,dt,v=NULL,qc=rep(-1,3),drop=NULL,lt=NULL)
   drop(Xf%*%beta)
  
   Xbd(X,beta,k,ks,ts,dt,v=NULL,qc=rep(-1,3),drop=NULL,lt=lt)
   drop(X1%*%beta)

   ## check out XWXd
   w <- runif(100)-.1
   XWX <- XWXd(X,w,k,ks,ts,dt,v=NULL,qc=rep(0,3))
   XWXf <- t(Xf)%*%(w*Xf)
   range(XWXf-XWX)
   
   ## XWXd with selection
   X1 <- matrix(0,100,0) 
   for (i in 1:length(X)) if (mt[i]>0) X1 <- cbind(X1,X[[i]][k[,i],]) 
   XWX <- XWXd(X,w,k,ks,ts,dt,v=NULL,qc=rep(0,3),lt=rt,rt=rt)
   XWXf <- t(X1)%*%(w*X1)
   range(XWXf-XWX)
*/
  double *xv,*dc,*p0,*p1,*p2,*p3,*ei,*xi,*dwork;
  ptrdiff_t bsj,bs,bsf,i,j,kk,space[3],*pwork;
  int one=1,ii,*iwork;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif

  if (*nthreads<1) *nthreads = 1;
  if (*nthreads > *cv) *nthreads = *cv;
  Xbdspace(space,m,p,n,nx,dt,nt);
  iwork = (int *) CALLOC((size_t) space[0] * *nthreads,sizeof(int));
  pwork = (ptrdiff_t *) CALLOC((size_t)space[1] * *nthreads,sizeof(ptrdiff_t));
  dwork = (double *)CALLOC((size_t)space[2] * *nthreads,sizeof(double));

  xv = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for cols of XV */
  xi = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for cols of X */
  ei = (double *) CALLOC((size_t) *nthreads * *cv,sizeof(double)); /* storage for identity matrix cols */
  dc = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for components of diag */
  if (*nthreads>1) {
    bs = *cv / *nthreads;
    while (bs * *nthreads < *cv) bs++;
    while (bs * *nthreads - bs >= *cv) (*nthreads)--;
    bsf = *cv - (bs * *nthreads - bs); 
  } else {
    bsf = bs = *cv;
  }
  if (*ncs <=0 ) { *ncs = *nt; for (ii=0;ii<*nt;ii++) cs[ii] = ii;} /* set here to avoid reset in parallel sections by Xbd */
  if (*nrs <=0 ) { *nrs = *nt; for (ii=0;ii<*nt;ii++) rs[ii] = ii;}
  #ifdef OPENMP_ON
  #pragma omp parallel for private(j,bsj,i,kk,p0,p1,p2,p3) num_threads(*nthreads)
  #endif  
  for (j=0;j < *nthreads;j++) {
    if (j == *nthreads - 1) bsj = bsf; else bsj = bs;
    for (i=0;i<bsj;i++) { /* work through this block's columns */
      kk = j * bs + i; /* column being worked on */
      ei[j * *pv + kk] = 1;if (i>0) ei[j * *pv + kk - 1] = 0;
      /* Note thread safety of XBd means this must be only memory allocator in this section*/
      Xbd(xv + j * *n,V + kk * *pv,X,k1,ks,m,p,n,nx,ts,dt,nt,v,qc,&one,cs,ncs,
	  iwork+j*space[0],pwork+j*space[1],dwork+j*space[2]); /* XV[:,kk] */
      /* Get X[:,kk] - inefficient when k1==k2, but deals with constraint*/
      Xbd(xi + j * *n,ei + j * *pv,X,k2,ks,m,p,n,nx,ts,dt,nt,v,qc,&one,rs,nrs,
	  iwork+j*space[0],pwork+j*space[1],dwork+j*space[2]); 
      p0 = xi + j * *n;p1=xv + j * *n;p2 = dc + j * *n;p3 = p2 + *n;
      for (;p2<p3;p0++,p1++,p2++) *p2 += *p0 * *p1; /* element-wise product of XV[:,kk] X[:,kk] */
    } 
  } /* parallel loop end */
  /* sum the contributions from the different threads into diag... */
  for (p0=diag,p1=p0+ *n,p2=dc;p0<p1;p0++,p2++) *p0 = *p2;
  for (i=1;i< *nthreads;i++) for (p0=diag,p1=p0+ *n;p0<p1;p0++,p2++) *p0 += *p2;
  FREE(xv);FREE(dc);FREE(xi);FREE(ei); FREE(iwork);FREE(pwork);FREE(dwork);
} /* diagXVXt */


static int kucomp(const void *a, const void *b) {
  /* comparison function for qsort as used in kunique. a and b are pointers to integer pointers. Each points in particular to an
     integer pair to sort, on the basis of the first element (the second element is an index tracker).
  */   
  int k;
  k = *(*(int **)a) - *(*(int **)b); /* difference between first elements pointed to */
  return(k);
} /* kucomp */  

void kunique(int *k, int *ind,ptrdiff_t *n) {
/* finds the unique integers in k returning the number of these in *n (which is array length on entry).
   The unique entries are returned in the early elements of k and ind[i] gives the elements of these unique entries
   giving the ith input entry. i.e. k[i] on input is located in k[ind[i]] on output.
   - could be optimized once tested.
*/
  int **a,*ap,*a0,i,j;
  a = (int **)CALLOC(*n,sizeof(int *));
  a0 = ap = (int *)CALLOC(*n * 2,sizeof(int));
  for (i=0;i<*n;i++,ap += 2) { /* a[i] points to an integer pair (k[i],i) */
    a[i] = ap;ap[0] = k[i];
    ap[1] = i;
  }  
  qsort(a,*n,sizeof(int *),kucomp); /* sort the pair pointers to order k */
  /* read sorted back into k and ind */
  for (i=0;i<*n;i++,ap += 3) {
    ap = a[i];
    k[i] = *ap;ind[i] = ap[1];
  }
  for (i=0;i<*n;i++) a0[i] = ind[i]; /* need a copy for re-arranging into ind */
  /* now compress unique into k, while recording index into ind */
  j = 0;
  ind[*a0] = j;
  for (i=1;i<*n;i++) {
    if (k[i]!=k[i-1]) { /* k change from previous */
      j++; k[j] = k[i];
    }
    ind[a0[i]] = j; /* record that original index maps to unique row j */
  }
  *n = j+1; /* unique element count */
  FREE(a0);FREE(a);
} /* kunique */  

void idiagXLLtXt(double *diag,double *L,double *X,int *k,int *ks,int *m,int *p, int *n, 
		 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pl,int *cl,
		 int *ri,int *ci,int *nrc,int *nthreads) {
  ptrdiff_t ni,nrci; ni = (ptrdiff_t) *n;nrci = (ptrdiff_t) *nrc;
  diagXLLtXt(diag,L,X,k,ks,m,p,&ni,nx,ts,dt,nt,v,qc,pl,cl,ri,ci,&nrci,nthreads);
}
  
void diagXLLtXt(double *diag,double *L,double *X,int *k,int *ks,int *m,int *p, ptrdiff_t *n, 
	      int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pl,int *cl,
	      int *ri,int *ci,ptrdiff_t *nrc,int *nthreads) {
/* Extracts elements indexed by ri and ci of XLL'X'. So on exit diag[i] = (XLL'X')[ri[i],ci[i]];

   L is a pl by cl matrix. ri and ci are nrc vectors. X is packed in standard disrete form, other arguments 
   relating to this.
   
   Parallelization is by splitting the columns of L into nthreads subsets.

   Basic algorithm is to compute XL and then the row sums of XL.XL. 
 
   In practice this is done by computing one column of XL at a time, and only the required rows 
   are computed.

   To test from R it is easiest to create a calling wrapper to deal with the ptrdiff_t *n and defined 
   in mgcv.h and init.c.

   Called using diagXVXd in R
   library(mgcv)
   X <- list(matrix(runif(4*3),4,3),matrix(runif(8*4),8,4),matrix(runif(5*2),5,2))
   k <- cbind(sample(1:4,100,replace=TRUE),sample(1:8,100,replace=TRUE),sample(1:5,100,replace=TRUE))
   ks <- cbind(1:3,2:4); ts <- 1:3;dt <- rep(1,3)
   attr(X,"lpip") <- list(1:3,4:7,8:9)
   Xf <- cbind(X[[1]][k[,1],],X[[2]][k[,2],],X[[3]][k[,3],])
   V <- tcrossprod(L <- matrix(runif(45)-.5,9,5))
   diagXVXd(X,V,k,ks,ts,dt,v=NULL,qc=rep(-1,3))
   d1 <- rowSums((Xf%*%V)*Xf)

   ## all diagonal elements...

   oo <- .C("idiagXLLtXt",diag=numeric (100),as.double(L),as.double(unlist(X)),as.integer(k-1L),as.integer(ks-1L),
      m=as.integer(unlist(lapply(X,nrow))),p=as.integer(unlist(lapply(X,ncol))),n=as.integer(100),nx=as.integer(length(X)),
      as.integer(ts-1),as.integer(dt),as.integer(length(dt)),v=as.double(0),qc=as.integer(rep(-1,3)),
      pl=as.integer(9),cl=as.integer(5),as.integer(0:99),as.integer(0:99),as.integer(100),as.integer(2),PACKAGE="mgcv")
   d2 <- oo$diag
   range(d1-d2)
 
   ## not all diagonal elements...

   oo <- .C("idiagXLLtXt",diag=numeric (50),as.double(L),as.double(unlist(X)),as.integer(k-1L),as.integer(ks-1L),
      m=as.integer(unlist(lapply(X,nrow))),p=as.integer(unlist(lapply(X,ncol))),n=as.integer(100),nx=as.integer(length(X)),
      as.integer(ts-1),as.integer(dt),as.integer(length(dt)),v=as.double(0),qc=as.integer(rep(-1,3)),
      pl=as.integer(9),cl=as.integer(5),as.integer(0:49),as.integer(0:49),as.integer(50),as.integer(1),PACKAGE="mgcv")
   d2 <- oo$diag
   range(d1[1:50]-d2)

   ## scattered elements...
   i <- c(2,4,8,25); j <- c(1,7,4,56)
   d1 <- mgcv:::ijXVXd(i,j,X,V,k,ks,ts,dt,v=NULL,qc=rep(0,3))
   
   oo <- .C("idiagXLLtXt",diag=numeric (4),as.double(L),as.double(unlist(X)),as.integer(k-1L),as.integer(ks-1L),
      m=as.integer(unlist(lapply(X,nrow))),p=as.integer(unlist(lapply(X,ncol))),n=as.integer(100),nx=as.integer(length(X)),
      as.integer(ts-1),as.integer(dt),as.integer(length(dt)),v=as.double(0),qc=as.integer(rep(-1,3)),
      pl=as.integer(9),cl=as.integer(5),as.integer(i-1),as.integer(j-1),as.integer(4),as.integer(1),PACKAGE="mgcv")
   d2 <- oo$diag
   d2;d1;range(d2-d1)

 
*/
  double *xl,*dc,*p0,*p1,*p2,*p3,*dwork;
  ptrdiff_t bsj,bs,bsf,i,j,kk,nu,*pwork,space[3];
  int one=1,*ku,*i1p,*i2p,*cs,ncs=-1,*iu,*ip,*iir,*iic,*ii,ck=0,*iwork;
  /* first find the unique rows that need to be extracted - need to find unique
     entries iu in ri,ci and the index of ri and ci in iu...
  */
  iu = (int *)CALLOC((size_t) 2 * *nrc,sizeof(int));
  for (ip=iu,i1p=ri,i2p = ri + *nrc ;i1p<i2p;i1p++,ip++) *ip = *i1p;
  for (i1p=ci,i2p = ci + *nrc ;i1p<i2p;i1p++,ip++) *ip = *i1p;
  iir = ii = (int *)CALLOC((size_t) 2 * *nrc,sizeof(int));
  nu = *nrc * 2; kunique(iu, iir, &nu); iic = iir + *nrc;
  /* there are nu unique values in iu, ri[i] = iu[iir[i]], ci[i] = iu[iic[i]] */
  if (nu<*n) { /* not all rows of XL are required - reduce k accordingly */
    for (i1p=ks,i2p=ks+2 * *nx;i1p<i2p;i1p++) if (ck < *i1p) ck = *i1p; /* find number cols of k */
    ku = (int *)CALLOC((size_t) ck   * nu,sizeof(int));
    for (i=0;i<nu;i++) {
      i1p = ku+i;i2p = k + iu[i];
      for (j=0;j<ck;j++,i1p+=nu,i2p+= *n) *i1p = *i2p; 
    }  
  } else { /* can use k, ri and ci directly */
    iir = ri;iic = ci;
    ku=k;
  }  
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  if (*nthreads > *cl) *nthreads = *cl;
  Xbdspace(space,m,p,n,nx,dt,nt);
  iwork = (int *) CALLOC((size_t) space[0] * *nthreads,sizeof(int));
  pwork = (ptrdiff_t *) CALLOC((size_t)space[1] * *nthreads,sizeof(ptrdiff_t));
  dwork = (double *)CALLOC((size_t)space[2] * *nthreads,sizeof(double));
  
  xl = (double *) CALLOC((size_t) *nthreads * nu,sizeof(double)); /* storage for cols of XL */
  dc = (double *) CALLOC((size_t) *nthreads * *nrc,sizeof(double)); /* storage for components of required elements */
  cs = (int *)CALLOC((size_t) *nt,sizeof(int));
  ncs = *nt; for (i=0;i<ncs;i++) cs[i] = i; /* set explicitly to avoid setting in Xbd in parallel section */
  if (*nthreads>1) {
    bs = *cl / *nthreads;
    while (bs * *nthreads < *cl) bs++;
    while (bs * *nthreads - bs >= *cl) (*nthreads)--;
    bsf = *cl - (bs * *nthreads - bs); 
  } else {
    bsf = bs = *cl;
  }
  #ifdef OPENMP_ON
#pragma omp parallel for private(j,bsj,i,kk,i1p,i2p,p1,p2,p3) num_threads(*nthreads)
  #endif  
  for (j=0;j < *nthreads;j++) {
    if (j == *nthreads - 1) bsj = bsf; else bsj = bs;
    for (i=0;i<bsj;i++) { /* work through this block's columns */
      kk = j * bs + i; /* column being worked on */   
      Xbd(xl + j * nu,L + kk * *pl,X,ku,ks,m,p,&nu,nx,ts,dt,nt,v,qc,&one,cs,&ncs,
	  iwork+j*space[0],pwork+j*space[1],dwork+j*space[2]); /* XL[:,kk] */
      /* Now form dc[i] += XL[k1i[i],kk]*XL[k2i[i],kk] */
      p1=xl + j * nu;p2 = dc + j * *nrc; p3 = p2 + *nrc;
      for (i1p=iir,i2p=iic;p2<p3;p2++,i1p++,i2p++) *p2 += p1[*i1p] * p1[*i2p];
    } 
  } /* parallel loop end */
  /* sum the contributions from the different threads into diag... */
  for (p0=diag,p1=p0+ *nrc,p2=dc;p0<p1;p0++,p2++) *p0 = *p2;
  for (i=1;i< *nthreads;i++) for (p0=diag;p0<p1;p0++,p2++) *p0 += *p2;
  FREE(xl);FREE(dc);FREE(cs);FREE(iu);FREE(ii);if (nu<*n) FREE(ku);
  FREE(iwork);FREE(pwork);FREE(dwork);
} /* diagXLLtXt */

void idiagXLUtXt(double *diag,double *L,double *U,double *X,int *k,int *ks,int *m,int *p, int *n, 
		 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pl,int *cl,
		 int *ri,int *ci,int *nrc,int *nthreads) {
  ptrdiff_t ni,nrci; ni = (ptrdiff_t) *n;nrci = (ptrdiff_t) *nrc;
  diagXLUtXt(diag,L,U,X,k,ks,m,p,&ni,nx,ts,dt,nt,v,qc,pl,cl,ri,ci,&nrci,nthreads);
}
  
void diagXLUtXt(double *diag,double *L,double *U,double *X,int *k,int *ks,int *m,int *p, ptrdiff_t *n, 
	      int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pl,int *cl,
	      int *ri,int *ci,ptrdiff_t *nrc,int *nthreads) {
/* Extracts elements indexed by ri and ci of X(LU'+UL')X'. So on exit diag[i] = (X(LU'+UL')X')[ri[i],ci[i]];

   L and U are pl by cl matrices. ri and ci are nrc vectors. X is packed in standard disrete form, other arguments 
   relating to this.
   
   Parallelization is by splitting the columns of L into nthreads subsets.

   Basic algorithm is to compute XL and XU then the row sums of XL.XU. 
 
   In practice this is done by computing one column of XLand XU  at a time, and only the required rows 
   are computed.

   To test from R it is easiest to create a calling wrapper to deal with the ptrdiff_t *n and defined 
   in mgcv.h and init.c.

   Called using diagXVXd in R
   library(mgcv)
   X <- list(matrix(runif(4*3),4,3),matrix(runif(8*4),8,4),matrix(runif(5*2),5,2))
   k <- cbind(sample(1:4,100,replace=TRUE),sample(1:8,100,replace=TRUE),sample(1:5,100,replace=TRUE))
   ks <- cbind(1:3,2:4); ts <- 1:3;dt <- rep(1,3)
   attr(X,"lpip") <- list(1:3,4:7,8:9)
   Xf <- cbind(X[[1]][k[,1],],X[[2]][k[,2],],X[[3]][k[,3],])
   L <- matrix(runif(45)-.5,9,5); U <- matrix(runif(45)-.5,9,5)
   V <- L %*% t(U)
   d0 <-diagXVXd(X,V,k,ks,ts,dt,v=NULL,qc=rep(-1,3))
   d1 <- rowSums((Xf%*%V)*Xf)

   ## all diagonal elements...

   oo <- .C("idiagXLUtXt",diag=numeric (100),as.double(L),as.double(U),as.double(unlist(X)),as.integer(k-1L),as.integer(ks-1L),
      m=as.integer(unlist(lapply(X,nrow))),p=as.integer(unlist(lapply(X,ncol))),n=as.integer(100),nx=as.integer(length(X)),
      as.integer(ts-1),as.integer(dt),as.integer(length(dt)),v=as.double(0),qc=as.integer(rep(-1,3)),
      pl=as.integer(9),cl=as.integer(5),as.integer(0:99),as.integer(0:99),as.integer(100),as.integer(1),PACKAGE="mgcv")
   d2 <- oo$diag
   range(d1-d2)
 
   ## not all diagonal elements...

   oo <- .C("idiagXLUtXt",diag=numeric (50),as.double(L),as.double(U),as.double(unlist(X)),as.integer(k-1L),as.integer(ks-1L),
      m=as.integer(unlist(lapply(X,nrow))),p=as.integer(unlist(lapply(X,ncol))),n=as.integer(100),nx=as.integer(length(X)),
      as.integer(ts-1),as.integer(dt),as.integer(length(dt)),v=as.double(0),qc=as.integer(rep(-1,3)),
      pl=as.integer(9),cl=as.integer(5),as.integer(0:49),as.integer(0:49),as.integer(50),as.integer(1),PACKAGE="mgcv")
   d2 <- oo$diag
   range(d1[1:50]-d2)

   ## scattered elements...
   i <- c(2,4,8,25); j <- c(1,7,4,56)
   d1 <- mgcv:::ijXVXd(i,j,X,V,k,ks,ts,dt,v=NULL,qc=rep(0,3))
   
   oo <- .C("idiagXLUtXt",diag=numeric (4),as.double(L),as.double(U),as.double(unlist(X)),as.integer(k-1L),as.integer(ks-1L),
      m=as.integer(unlist(lapply(X,nrow))),p=as.integer(unlist(lapply(X,ncol))),n=as.integer(100),nx=as.integer(length(X)),
      as.integer(ts-1),as.integer(dt),as.integer(length(dt)),v=as.double(0),qc=as.integer(rep(-1,3)),
      pl=as.integer(9),cl=as.integer(5),as.integer(i-1),as.integer(j-1),as.integer(4),as.integer(1),PACKAGE="mgcv")
   d2 <- oo$diag
   d2;d1;range(d2-d1)

 
*/
  double *xl,*xu,*dc,*p0,*p1,*p2,*p3,*dwork;
  ptrdiff_t bsj,bs,bsf,i,j,kk,nu,*pwork,space[3];
  int one=1,*ku,*i1p,*i2p,*cs,ncs=-1,*iu,*ip,*iir,*iic,*ii,ck=0,*iwork;
  /* first find the unique rows that need to be extracted - need to find unique
     entries iu in ri,ci and the index of ri and ci in iu...
  */
  iu = (int *)CALLOC((size_t) 2 * *nrc,sizeof(int));
  for (ip=iu,i1p=ri,i2p = ri + *nrc ;i1p<i2p;i1p++,ip++) *ip = *i1p;
  for (i1p=ci,i2p = ci + *nrc ;i1p<i2p;i1p++,ip++) *ip = *i1p;
  iir = ii = (int *)CALLOC((size_t) 2 * *nrc,sizeof(int));
  nu = *nrc * 2; kunique(iu, iir, &nu); iic = iir + *nrc;
  /* there are nu unique values in iu, ri[i] = iu[iir[i]], ci[i] = iu[iic[i]] */
  if (nu<*n) { /* not all rows of XL are required - reduce k accordingly */
    for (i1p=ks,i2p=ks+2 * *nx;i1p<i2p;i1p++) if (ck < *i1p) ck = *i1p; /* find number cols of k */
    ku = (int *)CALLOC((size_t) ck * nu,sizeof(int));
    for (i=0;i<nu;i++) {
      i1p = ku+i;i2p = k + iu[i];
      for (j=0;j<ck;j++,i1p+=nu,i2p+= *n) *i1p = *i2p; 
    }  
  } else { /* can use k, ri and ci directly */
    iir = ri;iic = ci;
    ku=k;
  }  
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  if (*nthreads > *cl) *nthreads = *cl;
  Xbdspace(space,m,p,n,nx,dt,nt);
  iwork = (int *) CALLOC((size_t) space[0] * *nthreads,sizeof(int));
  pwork = (ptrdiff_t *) CALLOC((size_t)space[1] * *nthreads,sizeof(ptrdiff_t));
  dwork = (double *)CALLOC((size_t)space[2] * *nthreads,sizeof(double));
  
  xl = (double *) CALLOC((size_t) *nthreads * nu,sizeof(double)); /* storage for cols of XL */
  xu = (double *) CALLOC((size_t) *nthreads * nu,sizeof(double)); /* storage for cols of XU */
  dc = (double *) CALLOC((size_t) *nthreads * *nrc,sizeof(double)); /* storage for components of required elements */
  cs = (int *)CALLOC((size_t) *nt,sizeof(int));
  ncs = *nt; for (i=0;i<ncs;i++) cs[i] = i; /* avoid resetting in Xbd in parallel section */
  if (*nthreads>1) {
    bs = *cl / *nthreads;
    while (bs * *nthreads < *cl) bs++;
    while (bs * *nthreads - bs >= *cl) (*nthreads)--;
    bsf = *cl - (bs * *nthreads - bs); 
  } else {
    bsf = bs = *cl;
  }
  #ifdef OPENMP_ON
#pragma omp parallel for private(j,bsj,i,kk,i1p,i2p,p0,p1,p2,p3) num_threads(*nthreads)
  #endif  
  for (j=0;j < *nthreads;j++) {
    if (j == *nthreads - 1) bsj = bsf; else bsj = bs;
    for (i=0;i<bsj;i++) { /* work through this block's columns */
      kk = j * bs + i; /* column being worked on */   
      Xbd(xl + j * nu,L + kk * *pl,X,ku,ks,m,p,&nu,nx,ts,dt,nt,v,qc,&one,cs,&ncs,
	  iwork+j*space[0],pwork+j*space[1],dwork+j*space[2]); /* XL[:,kk] */
      Xbd(xu + j * nu,U + kk * *pl,X,ku,ks,m,p,&nu,nx,ts,dt,nt,v,qc,&one,cs,&ncs,
	  iwork+j*space[0],pwork+j*space[1],dwork+j*space[2]); /* XU[:,kk] */
      /* Now form dc[i] += XL[k1i[i],kk]*XL[k2i[i],kk] */
      p0 = xu + j * nu;p1=xl + j * nu;p2 = dc + j * *nrc; p3 = p2 + *nrc;
      for (i1p=iir,i2p=iic;p2<p3;p2++,i1p++,i2p++) *p2 += p1[*i1p] * p0[*i2p] + p0[*i1p] * p1[*i2p] ;
    } 
  } /* parallel loop end */
  /* sum the contributions from the different threads into diag... */
  for (p0=diag,p1=p0+ *nrc,p2=dc;p0<p1;p0++,p2++) *p0 = *p2;
  for (i=1;i< *nthreads;i++) for (p0=diag,p1=p0+ *nrc;p0<p1;p0++,p2++) *p0 += *p2;
  FREE(xu);FREE(xl);FREE(dc);FREE(cs);FREE(iu);FREE(ii);if (nu<*n) FREE(ku);
  FREE(iwork);FREE(pwork);FREE(dwork);
} /* diagXLUtXt */


SEXP CXWyd(SEXP XWyr, SEXP yr, SEXP Xr, SEXP wr, SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP cyr, SEXP tsr, SEXP dtr,
	    SEXP vr,SEXP qcr, SEXP ar_stopr, SEXP ar_rowr, SEXP ar_weightsr,
	    SEXP csr) {
/* .Call wrapper for XWyd allowing R long vector storage for k. Note that 
  this does not allow more than maxint data - that would require re-writting R code 
  to avoid storing k in a matrix (which is only allowed maxint rows).

  n is number of rows of k, nx is length of m or p, nt is the length of ts or dt, 
  ncs and nrs are length of cs/rs.
*/
  double *XWy,*X,*v,*w,*ar_weights,*y;
  int *k,*ks,*m,*p,*ts,*dt,*qc,*cs,nx,nt,ncs,*ar_stop,*ar_row,*cy;
  ptrdiff_t n;
  n = (ptrdiff_t) nrows(kr); XWy = REAL(XWyr);
  X = REAL(Xr);w = REAL(wr);y = REAL(yr);
  k = INTEGER(kr); ks = INTEGER(ksr);
  m = INTEGER(mr); nx = length(mr);
  p = INTEGER(pr);cy=INTEGER(cyr);
  ar_stop = INTEGER(ar_stopr); ar_row = INTEGER(ar_rowr);
  ar_weights = REAL(ar_weightsr);
  ts = INTEGER(tsr); dt = INTEGER(dtr); nt = length(tsr);
  v = REAL(vr);qc = INTEGER(qcr);
  cs = INTEGER(csr);  ncs = length(csr);
  XWyd(XWy,y,X,w,k,ks,m,p,&n,cy,&nx,ts,dt,&nt,v,qc,ar_stop,ar_row,ar_weights,cs,&ncs);
  return(R_NilValue);
} /*  CXWyd */



void XWyd(double *XWy,double *y,double *X,double *w,int *k,int *ks, int *m,int *p, ptrdiff_t *n, int *cy,
         int *nx, int *ts, int *dt, int *nt,double *v,int *qc,
	  int *ar_stop,int *ar_row,double *ar_weights,int *cs, int *ncs) {
/* NOT thread safe. 
   If ncs > 0 then cs contains the subset of terms (blocks of model matrix columns) to include.   */
  double *Wy,*p0,*p1,*p2,*Xy0,*work,*work1,*work2;
  ptrdiff_t i,j,*off,*voff;
  int *tps,maxm=0,maxp=0,one=1,zero=0,*pt,add,q,kk,n_XWy;
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
    if (qc[i]==0) voff[i+1] = voff[i]; else if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else {  /* start of ith Q matrix */
      kk = (int) round(v[voff[i]]); /* number of contrasts in this KP contrast */
      voff[i+1] = voff[i] + kk + 2; 
    }  
    if (maxp < pt[i]) maxp=pt[i];
  }
  if (*ncs<=0) { /* return everything */
    for (j=0;j<*nt;j++) cs[j] = j;
    *ncs = *nt;
  }  
  for (kk=j=0;j<*ncs;j++) { /* get the offsets for the returned terms in the output */
    i = cs[j];tps[i] = kk;
    if (qc[i]==0) kk += pt[i]; /* where cth terms starts in param vector */ 
    else if (qc[i]>0) kk += pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
    else { /* Kronecker product of sum-to-zero contrasts */
      q = (int) round(v[voff[i]]); /* number of contrasts */
      kk += pt[i] - (int) round(v[voff[i]+q+1]); /* subtracting number of constraints */
    }  
  } /* kk is number of rows of XWy, at this point */

  n_XWy = kk;
  
  Xy0 =  (double *) CALLOC((size_t)maxp,sizeof(double));
  work =  (double *) CALLOC((size_t)*n,sizeof(double));
  work1 = (double *) CALLOC((size_t)maxm,sizeof(double));
  work2 = (double *) CALLOC((size_t)maxp * 2,sizeof(double));
  
  Wy = (double *) CALLOC((size_t)*n,sizeof(double)); /* Wy */
  for (j=0;j<*cy;j++) { /* loop over columns of y */
    for (p0=Wy,p1=Wy + *n,p2=w;p0<p1;p0++,y++,p2++) *p0 = *y * *p2; /* apply W to y */
    if (*ar_stop>=0) { /* AR components present (weights are sqrt, therefore) */
      rwMatrix(ar_stop,ar_row,ar_weights,Wy,(int *) n,&one,&zero,work);
      rwMatrix(ar_stop,ar_row,ar_weights,Wy,(int *) n,&one,&one,work); /* transpose of transform applied */
      for (p0=w,p1=w + *n,p2=Wy;p0<p1;p0++,p2++) *p2 *= *p0; /* square weights again */
    }
    /* now loop through terms applying the components of X'...*/
    for (kk=0;kk<*ncs;kk++) {
      i = cs[kk]; /* term to deal with now */ 
      add=0; 
      if (dt[i]>1) { /* it's a tensor */
        //tensorXty(Xy0,work,work1,Wy,X+off[ts[i]],m+ts[i],p+ts[i],dt+i,k+ts[i] * (ptrdiff_t) *n,n);
        for (q=0;q<ks[ts[i] + *nx]-ks[ts[i]];q++) {  /* loop through index columns */
          tensorXty(Xy0,work,work1,Wy,X+off[ts[i]],m+ts[i],p+ts[i],dt+i,k,n,&add,ks+ts[i],&q);
          add=1;
        }
        if (qc[i]!=0) { /* there is a constraint to apply Z'Xy0: form Q'Xy0 and discard first row... */
          /* Q' = I - vv' */
	  Ztb(XWy+tps[i],Xy0,v+voff[i],qc+i,&one,pt+i,work2); // BUG: tps[i] wrong for general constraints - general issue with handling dimensions under general constraint
        } else { /* straight copy */
          for (p0=Xy0,p1=p0+pt[i],p2=XWy+tps[i];p0<p1;p0++,p2++) *p2 = *p0;
        }
      } else { /* it's a singleton */
        //singleXty(XWy+tps[i],work1,Wy,X+off[ts[i]], m+ts[i],p+ts[i], k+ts[i] * (ptrdiff_t) *n,n);
        for (q=ks[ts[i]];q<ks[ts[i] + *nx];q++) { /* loop through index columns */  
          singleXty(XWy+tps[i],work1,Wy,X+off[ts[i]], m+ts[i],p+ts[i], k + q *  *n,n,&add);
          add=1;
        }
      }
    } /* term loop */
    XWy += n_XWy; /* moving on to next column */
  }  
  FREE(Wy); FREE(Xy0); FREE(work); FREE(work1); FREE(work2);
  FREE(pt); FREE(off); FREE(voff); FREE(tps);
} /* XWy */

void XWXd(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, ptrdiff_t *n, int *nx, int *ts, 
          int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,int *ar_row,double *ar_weights) {
/* Based on Wood, Li, Shaddick and Augustin (2017) JASA algorithms. i.e. vector oriented.
   Not currently used in mgcv. Not updated for Kronecker contrast constraints.  

   This version uses open MP by column, within sub-block.

   An alternative in which one thread did each sub-block scaled poorly, largely 
   because the work is very uneven between sub-blocks. 

   Forms Xt'WXt when Xt is divided into blocks of columns, each stored in compact form
   using arguments X and k. 
   * 'X' contains 'nx' blocks, the ith is an m[i] by p[i] matrix containing the unique rows of 
     the ith marginal model matrix. 
   * There are 'nt' model terms. Each term is made up of one or more maginal model matrices.
   * The jth term starts at block ts[j] of X, and has dt[j] marginal matrices. The terms model
     matrix is the row tensor product of its full (n row) marginals.
   * The index vectors converting the unique row matrices to full marginal matrices are in 
     'k', an n-row matrix of integers. Conceptually if Xj and kj represent the jth unique 
      row matrix and index vector then the ith row of the corresponding full marginal matrix 
      is Xj[kj[i],], but things are more complicated when each full term matrix is actually 
      the sum of several matrices (summation convention).
   * To handle the summation convention, each marginal matrix can have several index vectors. 
     'ks' is an nx by 2 matrix giving the columns of k corresponding to the ith marginal 
     model matrix. Specifically columns ks[i,1]:(ks[i,2]-1) of k are the index vectors for the ith 
     marginal. All marginals corresponding to one term must have the same number of index columns.
     The full model matrix for the jth term is constucted by summing over q the full 
     model matrices corresponding to the qth index vectors for each of its marginals.    
   * For example the exression for the full model matrix of the jth term is...
  
     X^full_j = sum_q prod_i X_{ts[j]+i}[k[,ks[i]+q],]  

     - q runs from 0 to ks[i,2] - ks[i,1] - 1; i runs from 0 to dt[j] - 1.
         
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
        Xi = XiB + kk *  *n; Xj = XjB + kk *  *n;
        tempn = tempnB + kk * *n;
        for (i=start[kk];i<start[kk+1];i++) { /* loop over columns of Xb */ 
          /* extract Xb[:,i]... */
          if (ks[ts[b]]==ks[ts[b] + *nx]-1) { /* single index column */
            if (dt[b]>1) { /* tensor */
              for (p0=Xi,p1=p0+*n;p0<p1;p0++) *p0 = 1.0; /* set Xi to 1 */
	      tensorXj(Xi,X+off[ts[b]], m + ts[b] , p + ts[b], dt+b, 
	      		    k, n,&i,ks + ts[b],&zero);
            } else { /* singleton */
              singleXj(Xi,X+off[ts[b]], m + ts[b], k + *n * ks[ts[b]] /* ts[b]*/, n,&i);
            }
          } else { /* need to sum over multiple cols of k */
            for (q = 0;q<ks[ts[b] + *nx]-ks[ts[b]];q++) { /* k column loop */
              if (dt[b]>1) { /* tensor */
                for (p0=Xj,p1=p0+*n;p0<p1;p0++) *p0 = 1.0; /* set Xj to 1 */
                tensorXj(Xj,X+off[ts[b]], m + ts[b] , p + ts[b], dt+b, 
			    k, n,&i,ks+ts[b],&q);
              } else { /* singleton */
                singleXj(Xj,X+off[ts[b]], m + ts[b], k + *n * (q + ks[ts[b]]) /* ts[b]*/, n,&i);
              }
              if (q == 0) for (p0=Xj,p1=p0+*n,p2=Xi;p0<p1;p0++,p2++) *p2 = *p0;
              else for (p0=Xj,p1=p0+*n,p2=Xi;p0<p1;p0++,p2++) *p2 += *p0;
            } /* k column loop */
          }
          /* apply W to Xi - for AR models this is sqrt(W) */
          for (p0=w,p1=w + *n,p2=Xi;p0<p1;p0++,p2++) *p2 *= *p0; 
          if (*ar_stop>=0) { /* AR components present (weights are sqrt, therefore) */
            rwMatrix(ar_stop,ar_row,ar_weights,Xi,(int *) n,&one,&zero,tempn);
            rwMatrix(ar_stop,ar_row,ar_weights,Xi,(int *) n,&one,&one,tempn); /* transpose of transform applied */
            for (p0=w,p1=w + *n,p2=Xi;p0<p1;p0++,p2++) *p2 *= *p0; /* sqrt weights again */
          }
          /* now form Xa'WXb[:,i]... */ 
          add=0; 
          for (q=0;q<ks[ts[a] + *nx]-ks[ts[a]];q++) {
            if (dt[a]>1) { /* tensor */
              tensorXty(xwx + i * pt[a],tempn,temp,Xi,X+off[ts[a]],m+ts[a],p+ts[a],
	      			    dt+a,k, n,&add,ks+ts[a],&q);
            } else { /* singleton */
              singleXty(xwx + i * pt[a],temp,Xi,X+off[ts[a]],m+ts[a],p+ts[a],k + *n * (q + ks[ts[a]]),n,&add);
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
} /* XWXd */




ptrdiff_t XWXijspace(int i,int j,int r,int c,int *k, int *ks, int *m, int *p,int nx,int n,int *ts, int *dt,int nt, int tri) {
/*  computes working memory requirement of XWXijs for given block - called by XWXspace below
*/
  int si,sj,//ri,rj,
    jm,im,ddtj,
    ii,rfac,ddti,tensi,tensj,acc_w,alpha;
  ptrdiff_t nwork=0,mim,mjm; /* avoid integer overflow in large pointer calculations */ 
  si = ks[ts[i]+nx]-ks[ts[i]]; /* number of terms in summation convention for i */
  im = ts[i]+dt[i]-1; /* the index of the final marginal for term i */
  mim = (ptrdiff_t) m[im];
  /* Allocate work space for dXi(n), dXj(n) and initialze pdXj*/
  nwork += 2*n;
  if (dt[i]==1&&dt[j]==1&&m[ts[i]]==n&&m[ts[j]]==n) { /* both sub matrices are dense  */
    // no allocation
  } else if (!tri && i==j && si==1) {/* simplest setup - just accumulate diagonal */ 
    /* Allocate space for wb(m[im]), wbs and wbl*/
    nwork += mim;		  
  } else { /* general case */
    sj = ks[ts[j]+nx]-ks[ts[j]]; /* number of terms in summation convention for j */
    jm = ts[j]+dt[j]-1; /* the index of the final marginal for term j */
    ddti = dt[i]-1;ddtj = dt[j]-1; /* number of marginals excluding final */
    if (ddti) tensi = 1; else tensi = 0; /* is term i a tensor? */
    if (ddtj) tensj = 1; else tensj = 0;
    mjm = (ptrdiff_t) m[jm];
   
    if (n>mjm*mim) acc_w = 1; else acc_w = 0; /* accumulate \bar W or \bar W X_j / \bar W'X_i)? */ 
    if (acc_w) {
      if (p[im]*mim*mjm + p[im]*p[jm]*mjm > mim*mjm*p[jm] + p[im]*p[jm]*mim) rfac=0; else rfac=1;
      /* Allocate storage for W (mim*mjm) */
      nwork += mim*mjm;
    } else {
      /* now establish whether to form left product, D, or right product C */
      if (tensi) ii = 2; else ii = 1;if (tensj) ii++;
      if (tri) alpha = ii*3 + 3; else alpha = ii + 1; /* ~ ops per iteration of accumulation loop */ 
      if (alpha*si*sj*n*p[im]+mjm*p[im]*p[jm]<alpha*si*sj*n*p[jm]+mim*p[im]*p[jm]) rfac = 0; else rfac=1; //rfac = 1; 
      if (mim == n) rfac = 0; else if (mjm == n) rfac = 1; /* make absolutely sure we do not form n by p*m product */
     
    } /* end of accumulation storage allocation */
    
    if (rfac) {
      /* Allocate storge for C mim by p[jm] */
      nwork +=  mim * p[jm];
    } else {
      /* Allocate storage for D mjm by p[im] */
      nwork += mjm * p[im];
    }	
    if (!acc_w &&((rfac && p[jm]>15)||(!rfac && p[im]>15))) {
	if (tri) nwork += 3*n; else nwork += n;
    }
  }        
  return(nwork);
} /* XWXijspace */  

ptrdiff_t XWXspace(int N,int *sb,int *b,int *B,int *R,int *C,int *k, int *ks, int *m, int *p,int *pt,int *pd,int nx,ptrdiff_t n,
		   int *ts, int *dt,int nt, int tri) {
/* Tedious routine to evaluate workspace requirment of XWXijs. Basically does a dummy run through the 
   blocks computing the memory requirement for each and recording the maximum used.
   Avoids over allocating. 
*/
  int j,kk,kb,i,rb,cb,rt,ct,r,c;
  ptrdiff_t nn,nmax=0;
  for (j=0;j<sb[N];j++) { /* the block loop */
    kk = b[j];kb=B[kk];
    rb = R[kb];cb=C[kb]; /* set up allows blocks to be computed in any order by re-arranging B */
    i = kk - sb[kb]; /* sub-block index */
    rt = pt[rb]/pd[rb]; ct = pt[cb]/pd[cb]; /* total rows and cols of sub-blocks */
    /* compute sub-row and column */
    if (sb[kb+1]-sb[kb]<rt*ct) { /* symmetric upper half only needed */ 
      r=0; while (i >= rt - r) { i -= rt - r;r++;}
      c = i + r;
    } else {
      r = i / ct;
      c = i % ct;
    }
    nn = XWXijspace(rb,cb,r,c,k,ks,m,p,nx,n,ts, dt,nt,tri);
    if (nmax<nn) nmax=nn;
  }
  return(nmax);
} /*XWXspace*/

void XVXijs(double *XWX,int i,int j,int r,int c, double *X,int *k, int *ks, int *m, int *p,int nx,ptrdiff_t n,int *ts, int *dt,
	    int nt, double *e,int *a,ptrdiff_t *ma,ptrdiff_t *off,double *work,int *worki,int nxwx) {
/* Forms product X_i'VX_j or the r,cth sub-block of this where X_i and X_j are stored in compact form. 
   V is the matrix of elements V[i,j] = e[i]*e[j] if j in nei(i) and 0 otherwise.

   Neighbourhood structures are stored in a and ma. Specifically a[ma[i-1]+1:ma[i]] are the neighbours of 
   i where ma[i-1] = -1 by convention. 
  
   As with X'WX, X'VX can be partitioned into blocks related to model terms. row block i column block j
   relates to term i and term j. However, if i and or j are tensor product terms, then they can be
   partitioned into sub-blocks. A tensor term has (total number of cols)/(number of cols of final marginal) 
   sub blocks associated with it, so a singleton has one sub-block. This routine computes sub-block r,c of 
   term block i,j (sub-block indexing is within its term block). The reason for working at this finer grain is that
   it makes load balancing easier when parallelizing.

   * off[i] is offset to start of ith matrix (out of nx)  
   If Xk is a tensor product term, let dXk denote row tensor product of all but its final marginal. 
   
  Basically as XWXijs, but with the algorithms for scattered V[i,j] substituted in polace of those for (tri-)
  diagonal W, and no sparse methods used.

*/
  int si,sj,//ri,rj,
    jm,im,kk,ddtj,*Ki,*Kj,pim,pjm,
    ii,jj,akk,rfac,t,s,ddti,tensi,tensj,acc_w,alpha,*Kik,q;
  ptrdiff_t mim,mjm,*map,kp; /* avoid integer overflow in large pointer calculations */ 
  double x,*dXi,*dXj,*Xi,*Xj,done=1.0,dzero=0.0,*Cq,*Dq,*ep,*ep1,
    *C=NULL,*D=NULL,*W=NULL,*p2,*p3;
  char trans = 'T',ntrans = 'N';
  si = ks[ts[i]+nx]-ks[ts[i]]; /* number of terms in summation convention for i */
  im = ts[i]+dt[i]-1; /* the index of the final marginal for term i */
  mim = (ptrdiff_t) m[im];
  /* Allocate work space for dXi(n), dXj(n) and initialze pdXj*/
  dXi = work;work += n;/*pdXj = */dXj = work;work += n;
  if (dt[i]==1&&dt[j]==1&&m[ts[i]]==n&&m[ts[j]]==n) { /* both sub matrices are dense - i.e. what is stored has n rows  */
    jm = ts[j];
    mjm = (ptrdiff_t) m[jm];
    pim = p[im];pjm = p[jm];
    sj = ks[ts[j]+nx]-ks[ts[j]]; /* number of terms in summation convention for j */
    for (ii=0;ii<pim;ii++) for (jj=0;jj<pjm;jj++)  XWX[ii + (ptrdiff_t)nxwx * jj] = 0.0;   
    for (s=0;s<si;s++) for (t=0;t<sj;t++) {
      Ki = k +  (ptrdiff_t)(ks[im]+s) * n; /* index for i */
      Kj = k +  (ptrdiff_t)(ks[jm]+t) * n; /* index for j */
      for (ii=0;ii<pim;ii++) { /* row of XWX */
        Xi = X + off[im] + mim * ii;
        if (i==j) jj=ii; else jj=0;
	for (jj=ii;jj<pjm;jj++) { /* cols of XWX symmetric case */
  	  Xj = X + off[jm] + mjm * jj;  
	  x=0.0;Kik=Ki;
	  for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++) for (;kp <= *map;kp++) {
	      akk = a[kp];  
	      x += *ep * e[akk] *  Xi[*Kik] * Xj[Kj[akk]];
	  }    
	  XWX[ii + (ptrdiff_t) nxwx * jj] += x;
	  if (i==j && jj!=ii) XWX[jj + (ptrdiff_t) nxwx * ii] += x;
        }  
      }
    }  
  } else { /* general case */
    sj = ks[ts[j]+nx]-ks[ts[j]]; /* number of terms in summation convention for j */
    jm = ts[j]+dt[j]-1; /* the index of the final marginal for term j */
    ddti = dt[i]-1;ddtj = dt[j]-1; /* number of marginals excluding final */
    if (ddti) tensi = 1; else tensi = 0; /* is term i a tensor? */
    if (ddtj) tensj = 1; else tensj = 0;
    mjm = (ptrdiff_t) m[jm];
   
    if (n>mjm*mim) acc_w = 1; else acc_w = 0; /* accumulate \bar W or \bar W X_j / \bar W'X_i)? */ 
    if (acc_w) {
      if (p[im]*mim*mjm + p[im]*p[jm]*mjm > mim*mjm*p[jm] + p[im]*p[jm]*mim) rfac=0; else rfac=1; /* which order for final product */
      /* Allocate storage for W (mim*mjm) */
      W = work; work += mim*mjm;
    } else {
      /* now establish whether to form left product, D, or right product C */
      if (tensi) ii = 2; else ii = 1;if (tensj) ii++;
      kk = ma[n-1]/n; /* average nei size */ 
      alpha = ii + kk; /* ~ ops per iteration of accumulation loop */ 
      if (alpha*si*sj*n*p[im]+mjm*p[im]*p[jm]<alpha*si*sj*n*p[jm]+mim*p[im]*p[jm]) rfac = 0; else rfac=1; //rfac = 1; 
      if (mim == n) rfac = 0; else if (mjm == n) rfac = 1; /* make absolutely sure we do not form n by p*m product */
     
    } /* end of accumulation storage allocation */
    
    if (rfac) {
      /* Allocate storage for C mim by p[jm] */
      C = work; work += mim * p[jm];
    } else {
      /* Allocate storage for D mjm by p[im] */
      D = work; work += mjm * p[im]; 
    }	
   
    if (acc_w) for (kk=0;kk<mim*mjm;kk++) W[kk] = 0.0; /* clear W */
    else if (rfac) for (kk=0;kk<mim*p[jm];kk++) C[kk] = 0.0; /* clear C */
    else for (kk=0;kk<mjm*p[im];kk++) D[kk] = 0.0; /* clear D */

    for (s=0;s<si;s++) for (t=0;t<sj;t++) { /* summation convention loop */
      if (tensi) { /* extract col r of dXi according to sth set of index vectors */
	for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddti, 
		 k, &n, &r, ks + ts[i],&s);
      }
      if (tensj) { /* extract col c of dXj according to tth set of index vectors */
	for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
        tensorXj(dXj, X + off[ts[j]], m + ts[j], p + ts[j],&ddtj, 
		 k, &n, &c, ks + ts[j],&t);
      }
      Ki = k +  (ptrdiff_t)(ks[im]+s) * n; /* index for final margin of i */
      Kj = k +  (ptrdiff_t)(ks[jm]+t) * n; /* index for final margin of j */
      if (acc_w) { /* weight accumulation */    
        if (tensi&&tensj) { /* i and j are tensors */
	   Kik=Ki;p2=dXi;p3=dXj;
	   for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++,p2++) for (;kp <= *map;kp++) {
	       akk = a[kp];
	       W[*Kik + mim * Kj[akk]] += *ep * e[akk] * *p2 * p3[akk];
	   }    
	 } else if (tensi) { /* only i is tensor */
	   Kik=Ki;p2=dXi;
	   for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++,p2++) for (;kp <= *map;kp++) {
	       akk = a[kp];
	       W[*Kik + mim * Kj[akk]] += *ep * e[akk] * *p2;
	   }   
	 } else if (tensj) { /* only j is tensor */
	   Kik=Ki;p3=dXj;
	   for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++) for (;kp <= *map;kp++) {
	       akk = a[kp];
	       W[*Kik + mim * Kj[akk]] += *ep * e[akk] * p3[akk];
	   }
	 } else { /* i and j are singletons */  
	   Kik=Ki;
	   for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++) for (;kp <= *map; kp++) {
	       akk = a[kp];
	       W[*Kik + mim * Kj[akk]] += *ep * e[akk];
	   } 
	}	    
      } else if (rfac) { /* use direct accumulation */
          for (q=0;q<p[jm];q++) { 
	    Cq = C + q * mim;Xj = X + off[jm] + q * mjm; /* qth cols of C and X */ 
	    if (tensi&&tensj) { /* X_i and X_j are tensors */
	       Kik=Ki;p2=dXi;p3=dXj;
	       for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++,p2++) for (;kp <= *map;kp++) {
	         akk = a[kp];
	         Cq[*Kik] += *ep * e[akk] * *p2 * p3[akk] * Xj[Kj[akk]];
	       }
	    } else if (tensi) { /* only X_i is tensor */
	       Kik=Ki;p2=dXi;
	       for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++,p2++) for (;kp <= *map; kp++) {
	         akk = a[kp];
	         Cq[*Kik] += *ep * e[akk] * *p2 * Xj[Kj[akk]];
	       }
	    } else if (tensj) { /* only X_j is tensor */
	       Kik=Ki;p3=dXj;
	       for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++) for (;kp <= *map;kp++) {
	         akk = a[kp];
	         Cq[*Kik] += *ep * e[akk]  * p3[akk] * Xj[Kj[akk]];
	       }
	    } else { /* both singletons */
	       Kik=Ki;
	       for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++) for (;kp <= *map;kp++) {
	         akk = a[kp];
	         Cq[*Kik] += *ep * e[akk] * Xj[Kj[akk]];
	       }
	    }
	  } /* q loop */
      } else { /* direct accumulation of left factor D (mjm by p[im]) = \bar W' X_i */
          for (q=0;q<p[im];q++) {
	      Dq = D + q * mjm;Xi = X + off[im] + q * mim; /* qth cols of C and Xi */ 
	      if (tensi&&tensj) { /* X_i and X_j are tensors */
		Kik=Ki;p2=dXi;p3=dXj;
	        for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++,p2++) for (;kp <= *map;kp++) {
	          akk = a[kp];
	          Dq[Kj[akk]] += *ep * e[akk] * *p2 * p3[akk] * Xi[*Kik];
	        }
	      } else if (tensi) { /* only X_i is tensor */
		Kik=Ki;p2=dXi;
	        for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++,p2++) for (;kp <= *map;kp++) {
	          akk = a[kp];
	          Dq[Kj[akk]] += *ep * e[akk] * *p2 * Xi[*Kik];
	        }
	      } else if (tensj) { /* only X_j is tensor */
		Kik=Ki;p3=dXj;
	        for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++) for (;kp <= *map;kp++) {
	          akk = a[kp];
	          Dq[Kj[akk]] += *ep * e[akk] * p3[akk] * Xi[*Kik];
	        }
	      } else { /* both singletons */
		Kik=Ki;
	        for (kp=0,ep=e,ep1=e+n,map=ma;ep<ep1;ep++,map++,Kik++) for (;kp <= *map;kp++) {
	          akk = a[kp];
	          Dq[Kj[akk]] += *ep * e[akk] * Xi[*Kik];
	        }
	      }
	  } /* q loop */
      } /* direct  accumulation of D */
      
  } /* end of summation convention loop */
	if (acc_w) { /* form X_im' \bar W X_j from \bar W */
	  if (rfac) { /* form C = \bar W X_j first */
	  /* dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *A,
                   int *lda, double *B, int *ldb, double *beta,double *C,int *ldc) 
             transa/b = 'T' or 'N' for A/B transposed or not. C = alpha op(A) op(B) + beta C,
             where op() is transpose or not. C is m by n. k is cols of op(A). ldx is rows of X
             in calling routine (to allow use of sub-matrices) */   
            F77_CALL(dgemm)(&ntrans,&ntrans,m+im,p+jm,m+jm,&done,W,m+im,X+off[jm],m+jm,&dzero,C,m+im FCONE FCONE); 
	  } else { /* form D = \bar W' X_i first */
            F77_CALL(dgemm)(&trans,&ntrans,m+jm,p+im,m+im,&done,W,m+im,X+off[im],m+im,&dzero,D,m+jm FCONE FCONE);
	  }  
	} /* if acc_w */
	if (rfac) { /* Xi'C direct to the r,c p[im] by p[jm] block of XWX */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+im,&done,X+off[im],m+im,C,m+im,&dzero,XWX+r*p[im]+c*p[jm]*(ptrdiff_t)nxwx,&nxwx FCONE FCONE);
	} else { /* D'Xj (same block of XWX) */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+jm,&done,D,m+jm,X+off[jm],m+jm,&dzero,XWX+r*p[im]+c*p[jm]*(ptrdiff_t)nxwx,&nxwx FCONE FCONE);
	}  

  } /* general case */
} /* XVXijs */  

void XWXijs(double *XWX,int i,int j,int r,int c, double *X,int *k, int *ks, int *m, int *p,int nx,ptrdiff_t n,int *ts, int *dt,
	    int nt, double *w,double *ws, int tri,ptrdiff_t *off,double *work,int *worki,int nxwx,unsigned long long *ht,SM **sm,SM * SMstack) {
/* Forms product X_i'WX_j or the r,cth sub-block of this where X_i and X_j are stored in compact form. 
   Specifically, X'WX can be partitioned into blocks related to model terms. row block i column block j
   relates to term i and term j. However, if i and or j are tensor product terms, then they can be
   partitioned into sub-blocks. A tensor term has (total number of cols)/(number of cols of final marginal) 
   sub blocks associated with it, so a singleton has one sub-block. This routine computes sub-block r,c of 
   term block i,j (sub-block indexing is within its term block). The reason for working at this finer grain is that
   it makes load balancing easier when parallelizing.

   * W = diag(w) if tri!=0 and is tridiagonal otherwise, with super in ws and sub in ws + n-1;     
   * off[i] is offset to start of ith matrix (out of nx)  
   If Xk is a tensor product term, let dXk denote row tensor product of all but its final marginal. 
   
   Workspace: Let mi and mj index the final marginal of term i and j, then the dim of work needed is
              2*n + 3*m[mi] + I(n>m[mi]*m[mj])*m[mi]*m[mj] + max(m[mi]*p[mj],m[mj]*p[mi]) 
              + 3n 
             SMstack should be an n-vectors if tri==0 and a 3n-vector otherwise. sm is an n vector.

   This version uses sparse accumulation if \bar W too large, based on calling indReduce.
   ht is a 256 vector initialized by SMinihash. sm is length n. SMstack is length n if tri==0, 3n otherwise.

*/
  int si,sj,//ri,rj,
    jm,im,kk,ddt,ddtj,koff,*K,*Ki,*Kj,pim,pjm,
    ii,jj,rfac,t,s,ddti,tensi,tensj,acc_w,alpha,*Kik,*Kjk,*Kik1,*Kjk1,*Kjl1,q;
  ptrdiff_t mim,mjm; /* avoid integer overflow in large pointer calculations */ 
  double x,*wl,*dXi,*dXj,*pdXj,*Xt,*Xi,*Xj,done=1.0,dzero=0.0,*Cq,*Dq,
    *C=NULL,*D=NULL,*W=NULL,*wb,*p0,*p1,*p2,*p3,*pw,*pw1,*pl,*ps,*wi,*wsi,*wli,*wo,*psi,*pwi,*pli;
  char trans = 'T',ntrans = 'N';
  si = ks[ts[i]+nx]-ks[ts[i]]; /* number of terms in summation convention for i */
  if (tri) wl = ws + n - 1; else wl=ws; /* sub-diagonal else only to keep compiler happy */
  im = ts[i]+dt[i]-1; /* the index of the final marginal for term i */
  mim = (ptrdiff_t) m[im];
  /* Allocate work space for dXi(n), dXj(n) and initialze pdXj*/
  dXi = work;work += n;pdXj = dXj = work;work += n;
  if (dt[i]==1&&dt[j]==1&&m[ts[i]]==n&&m[ts[j]]==n) { /* both sub matrices are dense - i.e. what is stored has n rows  */
    jm = ts[j];
    mjm = (ptrdiff_t) m[jm];
    pim = p[im];pjm = p[jm];
    sj = ks[ts[j]+nx]-ks[ts[j]]; /* number of terms in summation convention for j */
    for (ii=0;ii<pim;ii++) for (jj=0;jj<pjm;jj++)  XWX[ii + (ptrdiff_t)nxwx * jj] = 0.0;   
    for (s=0;s<si;s++) for (t=0;t<sj;t++) {
      Ki = k +  (ptrdiff_t)(ks[im]+s) * n; /* index for i */
      Kj = k +  (ptrdiff_t)(ks[jm]+t) * n; /* index for j */
      for (ii=0;ii<pim;ii++) {
        Xi = X + off[im] + mim * ii;
        if (i==j) for (jj=ii;jj<pjm;jj++) { /* symmetric case */
  	  Xj = X + off[jm] + mjm * jj;  
	  if (tri) {
	    Kik = Ki;Kjk = Kj;Kjk1 = Kj+1;pw=w;ps=ws;pl=wl;
	    x = Xi[*Kik]*(Xj[*Kjk] * *pw + Xj[*Kjk1] * *ps);
	    ps++;pw++;Kjl1=Kj;Kjk++;Kik++;Kjk1++;p0 = w + n-1;
            for (;pw < p0;ps++,pw++,pl++,Kik++,Kjk++,Kjk1++,Kjl1++) x += Xi[*Kik]*(*pl * Xj[*Kjl1] + *pw * Xj[*Kjk] + *ps * Xj[*Kjk1]);
            x += Xi[*Kik]*(*pl * Xj[*Kjl1] + *pw * Xj[*Kjk]);			   
	  } else for (x=0.0,Kik=Ki,Kjk=Kj,pw=w,p0=w +n;pw<p0;Kik++,Kjk++,pw++) x += *pw * Xi[*Kik] * Xj[*Kjk];
	  XWX[jj + (ptrdiff_t) nxwx * ii] = XWX[ii + (ptrdiff_t) nxwx * jj] += x;
        } else for (jj=0;jj<pjm;jj++) {
          Xj = X + off[jm] + mjm * jj;
	  if (tri) {
	    Kik = Ki;Kjk = Kj;Kjk1 = Kj+1;pw=w;ps=ws;pl=wl;
	    x = Xi[*Kik]*(Xj[*Kjk] * *pw + Xj[*Kjk1] * *ps);
	    ps++;pw++;Kjl1=Kj;Kjk++;Kik++;Kjk1++;p0 = w + n-1;
            for (;pw < p0;ps++,pw++,pl++,Kik++,Kjk++,Kjl1++,Kjk1++) x += Xi[*Kik]*(*pl * Xj[*Kjl1] + *pw * Xj[*Kjk] + *ps * Xj[*Kjk1]);
            x += Xi[*Kik]*(*pl * Xj[*Kjl1] + *pw * Xj[*Kjk]);			   
	  } else for (x=0.0,Kik=Ki,Kjk=Kj,pw=w,p0=w +n;pw<p0;Kik++,Kjk++,pw++) x += *pw * Xi[*Kik] * Xj[*Kjk];
	  XWX[ii + (ptrdiff_t)nxwx * jj] += x;
        }	  
      }
    }  
  } else if (!tri && i==j && si==1) {/* simplest setup - just accumulate diagonal */
    /* note that if you turn this branch off in debugging then i==j case is forced to general
       code, which does NOT handle fact that only upper triangular blocks computed!! */
    
    /* Allocate space for wb(m[im]), wbs and wbl*/
    wb = work; work += mim;// wbs = work; work += m[i]; wbl = work; work += m[i];
    
    if (dt[i]>1) { /* tensor */
      ddt = dt[i]-1; /* number of marginals, excluding final */
      koff = 0; /* only one index vector per marginal, so no offset */
    }  
   
    if (dt[i]>1) { /* extract col r of dXi */
	for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddt, 
		 k, &n, &r, ks + ts[i],&koff);
    }
     
    if (dt[i]>1) { /* extract col c of dXi */
      if (r!=c) {
        dXj=pdXj;for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
	tensorXj(dXj, X + off[ts[i]], m + ts[i], p + ts[i],&ddt, 
		     k, &n, &c, ks + ts[i],&koff);
      } else dXj = dXi;
    }
    /* Clear work space to zero... */
    for (ii=0;ii<mim;ii++) wb[ii]=0.0;
   
    K = k +  (ptrdiff_t)ks[im] * n; /* index for final margin */
    /* Accumulate the weights ... */
    if (dt[i]>1) {
      for (kk=0;kk<n;kk++) wb[K[kk]] += dXi[kk]*dXj[kk]*w[kk];
    } else { /* singleton */
      for (kk=0;kk<n;kk++) wb[K[kk]] += w[kk];
    }
    /* Now form the Xi'WXi... */
    Xt = X + off[im]; /* final marginal model matrix */
    pim = p[im]; /* cols of Xt */
    if (r!=c) { for (jj=0;jj<pim;jj++) for (ii=0;ii<pim;ii++) { /* block need not be symmetric */
        Xi = Xt + mim * ii; /* iith col of Xt */
        Xj = Xt + mim * jj; /* jjth col of Xt */
        for (x=0.0,kk=0;kk<mim;kk++) x += wb[kk]*Xj[kk]*Xi[kk];
        XWX[c*pim+jj+(r*pim +ii)* (ptrdiff_t) nxwx] = XWX[r*pim+ii+(c*pim +jj)*(ptrdiff_t)nxwx] = x; 
      }
    } else for (ii=0;ii<pim;ii++) for (jj=ii;jj<pim;jj++) { /* diagonal and symmetric */ 
        Xi = Xt + mim * ii; /* iith col of Xt */
        Xj = Xt + mim * jj; /* jjth col of Xt */
        for (x=0.0,kk=0;kk<mim;kk++) x += wb[kk]*Xj[kk]*Xi[kk];
        XWX[r*pim+ ii + (c*pim+jj)*(ptrdiff_t)nxwx] = XWX[r*pim + jj + (c*pim + ii)*(ptrdiff_t)nxwx] =
        XWX[(r*pim+ ii)*(ptrdiff_t)nxwx + c*pim+jj] = XWX[(r*pim + jj)*(ptrdiff_t)nxwx + c*pim + ii] = x;
    }					  
  } else { /* general case */
    sj = ks[ts[j]+nx]-ks[ts[j]]; /* number of terms in summation convention for j */
    jm = ts[j]+dt[j]-1; /* the index of the final marginal for term j */
    ddti = dt[i]-1;ddtj = dt[j]-1; /* number of marginals excluding final */
    if (ddti) tensi = 1; else tensi = 0; /* is term i a tensor? */
    if (ddtj) tensj = 1; else tensj = 0;
    mjm = (ptrdiff_t) m[jm];
   
    if (n>mjm*mim) acc_w = 1; else acc_w = 0; /* accumulate \bar W or \bar W X_j / \bar W'X_i)? */ 
    if (acc_w) {
      if (p[im]*mim*mjm + p[im]*p[jm]*mjm > mim*mjm*p[jm] + p[im]*p[jm]*mim) rfac=0; else rfac=1;
      /* Allocate storage for W (mim*mjm) */
      W = work; work += mim*mjm;
    } else {
      /* now establish whether to form left product, D, or right product C */
      if (tensi) ii = 2; else ii = 1;if (tensj) ii++;
      if (tri) alpha = ii*3 + 3; else alpha = ii + 1; /* ~ ops per iteration of accumulation loop */ 
      if (alpha*si*sj*n*p[im]+mjm*p[im]*p[jm]<alpha*si*sj*n*p[jm]+mim*p[im]*p[jm]) rfac = 0; else rfac=1; //rfac = 1; 
      if (mim == n) rfac = 0; else if (mjm == n) rfac = 1; /* make absolutely sure we do not form n by p*m product */
     
    } /* end of accumulation storage allocation */
    
    if (rfac) {
      /* Allocate storge for C mim by p[jm] */
      C = work; work += mim * p[jm];
    } else {
      /* Allocate storage for D mjm by p[im] */
      D = work; work += mjm * p[im]; 
    }	
   
    if (acc_w) for (kk=0;kk<mim*mjm;kk++) W[kk] = 0.0; /* clear W */
    else if (rfac) for (kk=0;kk<mim*p[jm];kk++) C[kk] = 0.0; /* clear C */
    else for (kk=0;kk<mjm*p[im];kk++) D[kk] = 0.0; /* clear D */

    for (s=0;s<si;s++) for (t=0;t<sj;t++) { /* summation convention loop */
      if (tensi) { /* extract col r of dXi according to sth set of index vectors */
	for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddti, 
		 k, &n, &r, ks + ts[i],&s);
      }
      if (tensj) { /* extract col c of dXj according to tth set of index vectors */
	for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
        tensorXj(dXj, X + off[ts[j]], m + ts[j], p + ts[j],&ddtj, 
		 k, &n, &c, ks + ts[j],&t);
      }
      Ki = k +  (ptrdiff_t)(ks[im]+s) * n; /* index for final margin of i */
      Kj = k +  (ptrdiff_t)(ks[jm]+t) * n; /* index for final margin of j */
      if (acc_w) { /* weight accumulation */    
        if (tensi&&tensj) { /* i and j are tensors */
	  if (tri) {
	     for (Kik=Ki,Kik1=Ki+1,Kjk=Kj,Kjk1=Kj+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p0=dXi,p1=dXj,p2=dXi+1,p3=dXj+1;
		     pw<pw1;ps++,pw++,pl++,Kik++,Kik1++,Kjk++,Kjk1++,p0++,p1++,p2++,p3++) {
		  W[*Kik + mim * *Kjk1] += *ps * *p0 * *p3;
		  W[*Kik1 + mim * *Kjk] += *pl * *p2 * *p1;
		  W[*Kik + mim * *Kjk] += *pw * *p0 * *p1; 
	     }
	     W[*Kik + mim * *Kjk] += *pw * *p0 * *p1;
	   } else
	   for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,p2++,p3++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0 * *p2 * *p3;
	 } else if (tensi) { /* only i is tensor */
           if (tri) {
	      for (Kik=Ki,Kik1=Ki+1,Kjk=Kj,Kjk1=Kj+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p0=dXi,p2=dXi+1;
		     pw<pw1;ps++,pw++,pl++,Kik++,Kik1++,Kjk++,Kjk1++,p0++,p2++) {
		  W[*Kik + mim * *Kjk1] += *ps * *p0 ;
		  W[*Kik1 + mim * *Kjk] += *pl * *p2 ;
		  W[*Kik + mim * *Kjk] += *pw * *p0;
	      }
	      W[*Kik + mim * *Kjk] += *pw * *p0;
	    } else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,p2++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0 * *p2;   
	  } else if (tensj) { /* only j is tensor */
	    if (tri) {
		for (Kik=Ki,Kik1=Ki+1,Kjk=Kj,Kjk1=Kj+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p1=dXj,p3=dXj+1;
		     pw<pw1;ps++,pw++,pl++,Kik++,Kik1++,Kjk++,Kjk1++,p1++,p3++) {
		  W[*Kik + mim * *Kjk1] += *ps *  *p3;
		  W[*Kik1 + mim * *Kjk] += *pl  * *p1;
		  W[*Kik + mim * *Kjk] += *pw  * *p1;
	     }
	     W[*Kik + mim * *Kjk] += *pw  * *p1;
	   } else  for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,p3++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0  * *p3;    
	 } else { /* i and j are singletons */  
	   if (tri) {
		for (Kik=Ki,Kik1=Ki+1,Kjk=Kj,Kjk1=Kj+1,ps=ws,pw=w,pw1=w+n-1,pl=wl;
		     pw<pw1;ps++,pw++,pl++,Kik++,Kik1++,Kjk++,Kjk1++) {
		  W[*Kik + mim * *Kjk1] += *ps;
		  W[*Kik1 + mim * *Kjk] += *pl;
		  W[*Kik + mim * *Kjk] += *pw;
	    }
	    W[*Kik + mim * *Kjk] += *pw;
	  } else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0;
	}	    
      } else { /* \bar W too large to accumulate as dense - use sparse or direct accumulation */
        if ((rfac && p[jm]>15)||(!rfac && p[im]>15)) { 
	  ii = n; /* will contain compressed index length after indReduce call */
	  if (tri) {
	    wi = work;wsi = work+ii;wli = work+2*ii; /* do not increment work - inside s,t, loop!! */
	    if (tensi&&tensj) {
	      ps=ws;pw=w;pl=wl;pwi = wi;psi=wsi;pli=wli;wo=w+ii-1; 
	      for (p0=dXi,p1=dXj,p2=dXi+1,p3=dXj+1;pw<wo;p0++,p1++,ps++,psi++,pw++,pwi++,pli++,pl++,p2++,p3++) {
	        *pwi = *pw * *p0 * *p1;*psi = *ps * *p0 * *p3;*pli = *pl * *p2 * *p1;
	      }  
              *pwi = *pw * *p0 * *p1;
	    } else if (tensi) {
              ps=ws;pw=w;pl=wl;pwi = wi;psi=wsi;pli=wli;wo=w+ii-1; 
	      for (p0=dXi,p2=dXi+1;pw<wo;p0++,p2++,ps++,psi++,pw++,pwi++,pli++,pl++) {
	        *pwi = *pw * *p0;*psi = *ps * *p0;*pli = *pl * *p2;
	      }  
              *pwi = *pw * *p0;
	    } else if (tensj) {
              ps=ws;pw=w;pl=wl;pwi = wi;psi=wsi;pli=wli;wo=w+ii-1; 
	      for (p1=dXj,p3=dXj+1;pw<wo;p3++,p1++,ps++,psi++,pw++,pwi++,pli++,pl++) {
	        *pwi = *pw * *p1;*psi = *ps * *p3;*pli = *pl * *p1;
	      }  
              *pwi = *pw * *p1;
	    } else {
	      /* note: copied because indReduce over-writes input w */
	      ps=ws;pw=w;pl=wl;pwi = wi;psi=wsi;pli=wli;wo=w+ii-1; 
	      for (p1=dXj,p3=dXj+1;pw<wo;p3++,p1++,ps++,psi++,pw++,pwi++,pli++,pl++) {
	        *pwi = *pw;*psi = *ps;*pli = *pl;
	      }  
              *pwi = *pw;
	    }  
	  } else { /* not tri */
	    wi = work; /* do not increment work - inside s,t, loop!! */
	    if (tensi&&tensj) for (p0=wi,wo=wi + ii,p1=w,p2=dXi,p3=dXj;p0<wo;p0++,p1++,p2++,p3++) *p0 = *p1 * *p2 * *p3;
	    else if (tensi) for (p0=wi,wo=wi + ii,p1=w,p2=dXi;p0<wo;p0++,p1++,p2++) *p0 = *p1 * *p2;
	    else if (tensj) for (p0=wi,wo=wi + ii,p1=w,p3=dXj;p0<wo;p0++,p1++,p3++) *p0 = *p1 * *p3;
	    else for (p0=wi,wo=wi + ii,p1=w;p0<wo;p0++,p1++) *p0 = *p1;
	  }
	  /* form C or D using sparse matrix accumulation of \bar W */
	  if (rfac) indReduce(Ki,Kj,wi,tri,&ii,ht,sm,SMstack,C,X+off[jm],mim,p[jm],mjm,0,worki,1);
	  else indReduce(Ki,Kj,wi,tri,&ii,ht,sm,SMstack,D,X+off[im],mjm,p[im],mim,1,worki,1);
        } else if (rfac) { /* not worth using sparse methods (overhead too high) use direct accumulation */
          for (q=0;q<p[jm];q++) { 
	    Cq = C + q * mim;Xj = X + off[jm] + q * mjm; /* qth cols of C and X */ 
	    if (tensi&&tensj) { /* X_i and X_j are tensors */
		if (tri) {
		  for (Kik=Ki,Kjk=Kj,Kjk1=Kj+1,Kik1=Ki+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p0=dXi,p1=dXj,p2=dXi+1,p3=dXj+1;pw<pw1;
		       Kik++,Kik1++,Kjk++,Kjk1++,ps++,pw++,pl++,p0++,p1++,p2++,p3++) {
		    Cq[*Kik] += *p0 * (*ps * *p3 * Xj[*Kjk1] + *pw * *p1 * Xj[*Kjk]);
		    Cq[*Kik1] += *pl * *p2 * *p1 * Xj[*Kjk];
		  }
		  Cq[*Kik] += *pw * *p0 * *p1 * Xj[*Kjk];
		} else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p2++,p3++) Cq[*Kik] += *p0 * *p2 * *p3 * Xj[*Kjk];
	    } else if (tensi) { /* only X_i is tensor */
		if (tri) {
		  for (Kik=Ki,Kjk=Kj,Kjk1=Kj+1,Kik1=Ki+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p0=dXi,p2=dXi+1;pw<pw1;
		       Kik++,Kik1++,Kjk++,Kjk1++,ps++,pw++,pl++,p0++,p2++) {
		    Cq[*Kik] += *p0 * (*ps * Xj[*Kjk1] + *pw  * Xj[*Kjk]);
		    Cq[*Kik1] += *pl * *p2 * Xj[*Kjk];
		  }
		  Cq[*Kik] += *p0 * *pw * Xj[*Kjk];
		} else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,Kik++,Kjk++,p2++) Cq[*Kik] += *p0 * *p2  * Xj[*Kjk];
	    } else if (tensj) { /* only X_j is tensor */
		if (tri) {
		  for (Kik=Ki,Kjk=Kj,Kjk1=Kj+1,Kik1=Ki+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p1=dXj,p3=dXj+1;pw<pw1;
		       Kik++,Kik1++,Kjk++,Kjk1++,ps++,pw++,pl++,p1++,p3++) {
		    Cq[*Kik] += *ps * *p3 * Xj[*Kjk1] + *pw * *p1 * Xj[*Kjk];
		    Cq[*Kik1] += *pl * *p1 * Xj[*Kjk];
		  }
		  Cq[*Kik] += *pw * *p1 * Xj[*Kjk];
		} else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p3++) Cq[*Kik] += *p0 * *p3 * Xj[*Kjk];
	    } else { /* both singletons */
		if (tri) {
		  for (Kik=Ki,Kjk=Kj,Kjk1=Kj+1,Kik1=Ki+1,ps=ws,pw=w,pw1=w+n-1,pl=wl;pw<pw1;
		       Kik++,Kik1++,Kjk++,Kjk1++,pw++,pl++,ps++) {
		    Cq[*Kik] += *ps * Xj[*Kjk1] + *pw * Xj[*Kjk];
		    Cq[*Kik1] += *pl * Xj[*Kjk];
		  }
		  Cq[*Kik] +=  *pw * Xj[*Kjk]; 
		} else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) Cq[*Kik] += *p0 * Xj[*Kjk];
	    }
	  } /* q loop */
	} else { /* direct accumulation of left factor D (mjm by p[im]) = \bar W' X_i */
          for (q=0;q<p[im];q++) {
	      Dq = D + q * mjm;Xi = X + off[im] + q * mim; /* qth cols of C and Xi */ 
	      if (tensi&&tensj) { /* X_i and X_j are tensors */
		if (tri) {
		  for (Kik=Ki,Kjk=Kj,Kjk1=Kj+1,Kik1=Ki+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p0=dXi,p1=dXj,p2=dXi+1,p3=dXj+1;pw<pw1;
		       Kik++,Kik1++,Kjk++,Kjk1++,ps++,pw++,pl++,p0++,p1++,p2++,p3++) {
		    Dq[*Kjk] += *p1 * (*p2 * *pl * Xi[*Kik1] + *p0 * *pw * Xi[*Kik]);
		    Dq[*Kjk1] += *ps * *p0 * *p3 * Xi[*Kik];
		  }
		  Dq[*Kjk] += *p1 * *p0 * *pw * Xi[*Kik];
		} else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p2++,p3++) Dq[*Kjk] += *p0 * *p2 * *p3 * Xi[*Kik];
	      } else if (tensi) { /* only X_i is tensor */
		if (tri) {
		  for (Kik=Ki,Kjk=Kj,Kjk1=Kj+1,Kik1=Ki+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p0=dXi,p2=dXi+1;pw<pw1;
		       Kik++,Kik1++,Kjk++,Kjk1++,ps++,pw++,pl++,p0++,p2++) {
		    Dq[*Kjk] += *pl * *p2  * Xi[*Kik1] + *pw * *p0 * Xi[*Kik];
		    Dq[*Kjk1] += *ps * *p0  * Xi[*Kik];
		  }
		  Dq[*Kjk] += *pw * *p0 * Xi[*Kik];
		} else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,Kik++,Kjk++,p2++) Dq[*Kjk] += *p0 * *p2  * Xi[*Kik];
	      } else if (tensj) { /* only X_j is tensor */
		if (tri) {
		  for (Kik=Ki,Kjk=Kj,Kjk1=Kj+1,Kik1=Ki+1,ps=ws,pw=w,pw1=w+n-1,pl=wl,p1=dXj,p3=dXj+1;pw<pw1;
		       Kik++,Kik1++,Kjk++,Kjk1++,ps++,pw++,pl++,p1++,p3++) {
		    Dq[*Kjk] += *p1  * ( *pl * Xi[*Kik1] + *pw * Xi[*Kik]);
		    Dq[*Kjk1] += *ps * *p3 * Xi[*Kik];
		  }
		  Dq[*Kjk] += *p1  * *pw * Xi[*Kik];
		} else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p3++) Dq[*Kjk] += *p0 * *p3 * Xi[*Kik];
	      } else { /* both singletons */
		if (tri) {
		  for (Kik=Ki,Kjk=Kj,Kjk1=Kj+1,Kik1=Ki+1,ps=ws,pw=w,pw1=w+n-1,pl=wl;pw<pw1;
		       Kik++,Kik1++,Kjk++,Kjk1++,ps++,pw++,pl++) {
		    Dq[*Kjk] += *pl  * Xi[*Kik1] + *pw * Xi[*Kik];
		    Dq[*Kjk1] += *ps  * Xi[*Kik];
		  }
		  Dq[*Kjk] +=  *pw * Xi[*Kik];
		} else for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) Dq[*Kjk] += *p0 * Xi[*Kik];
	      }
	  } /* q loop */
	} /* direct  accumulation of D */
      }
    } /* end of summation convention loop */
	if (acc_w) { /* form X_im' \bar W X_j from \bar W */
	  if (rfac) { /* form C = \bar W X_j first */
	  /* dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *A,
                   int *lda, double *B, int *ldb, double *beta,double *C,int *ldc) 
             transa/b = 'T' or 'N' for A/B transposed or not. C = alpha op(A) op(B) + beta C,
             where op() is transpose or not. C is m by n. k is cols of op(A). ldx is rows of X
             in calling routine (to allow use of sub-matrices) */   
            F77_CALL(dgemm)(&ntrans,&ntrans,m+im,p+jm,m+jm,&done,W,m+im,X+off[jm],m+jm,&dzero,C,m+im FCONE FCONE); 
	  } else { /* form D = \bar W' X_i first */
            F77_CALL(dgemm)(&trans,&ntrans,m+jm,p+im,m+im,&done,W,m+im,X+off[im],m+im,&dzero,D,m+jm FCONE FCONE);
	  }  
	}
	if (rfac) { /* Xi'C direct to the r,c p[im] by p[jm] block of XWX */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+im,&done,X+off[im],m+im,C,m+im,&dzero,XWX+r*p[im]+c*p[jm]*(ptrdiff_t)nxwx,&nxwx FCONE FCONE);
	} else { /* D'Xj (same block of XWX) */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+jm,&done,D,m+jm,X+off[jm],m+jm,&dzero,XWX+r*p[im]+c*p[jm]*(ptrdiff_t)nxwx,&nxwx FCONE FCONE);
	}  

  } /* general case */
} /* XWXijs */  


SEXP CXWXd0(SEXP XWXr, SEXP Xr, SEXP wr, SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP tsr, SEXP dtr,
	    SEXP vr,SEXP qcr, SEXP nthreadsr, SEXP ar_stopr, SEXP ar_weightsr) {
/* .Call wrapper for XWXd0 allowing R long vector storage for k. Note that 
  this does not allow more than maxint data - that would require re-writting R code 
  to avoid storing k in a matrix (which is only allowed maxint rows).

  n is number of rows of k, nx is length of m or p, nt is the length of ts or dt.
*/
  double *XWX,*X,*v,*w,*ar_weights;
  int *k,*ks,*m,*p,*ts,*dt,*qc,*nthreads,nx,nt,*ar_stop;
  ptrdiff_t n;
  n = (ptrdiff_t) nrows(kr); XWX = REAL(XWXr);
  X = REAL(Xr);w = REAL(wr);
  k = INTEGER(kr); ks = INTEGER(ksr);
  m = INTEGER(mr); nx = length(mr);
  p = INTEGER(pr);
  ar_stop = INTEGER(ar_stopr); 
  ar_weights = REAL(ar_weightsr);
  ts = INTEGER(tsr); dt = INTEGER(dtr); nt = length(tsr);
  v = REAL(vr);qc = INTEGER(qcr);
  nthreads = INTEGER(nthreadsr);
  XWXd0(XWX,X,w,k,ks,m,p,&n,&nx,ts,dt,&nt,v,qc,nthreads,ar_stop,ar_weights);
  return(R_NilValue);
} /*  CXWXd0 */

void XWXd0(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, ptrdiff_t  *n, int *nx, int *ts, 
          int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,double *ar_weights) {
/* This version is the original without allowing the selection of sub-blocks

   essentially a driver routine for XWXij implementing block oriented cross products
   
   This version has the looping over sub-blocks, associated with tensor product terms, located in this routine
   to get better load balancing.
   
   Requires XWX to be over-sized on entry - namely n.params + n.terms by n.params + n.terms instead of
   n.params by n.params.


   Forms Xt'WXt when Xt is divided into blocks of columns, each stored in compact form
   using arguments X and k.  
   * 'X' contains 'nx' blocks, the ith is an m[i] by p[i] matrix containing the unique rows of 
     the ith marginal model matrix. 
   * There are 'nt' model terms. Each term is made up of one or more maginal model matrices.
   * The jth term starts at block ts[j] of X, and has dt[j] marginal matrices. The terms model
     matrix is the row tensor product of its full (n row) marginals.
   * The index vectors converting the unique row matrices to full marginal matrices are in 
     'k', an n-row matrix of integers. Conceptually if Xj and kj represent the jth unique 
      row matrix and index vector then the ith row of the corresponding full marginal matrix 
      is Xj[kj[i],], but things are more complicated when each full term matrix is actually 
      the sum of several matrices (summation convention).
   * To handle the summation convention, each marginal matrix can have several index vectors. 
     'ks' is an nx by 2 matrix giving the columns of k corresponding to the ith marginal 
     model matrix. Specifically columns ks[i,1]:(ks[i,2]-1) of k are the index vectors for the ith 
     marginal. All marginals corresponding to one term must have the same number of index columns.
     The full model matrix for the jth term is constucted by summing over q the full 
     model matrices corresponding to the qth index vectors for each of its marginals.    
   * For example the exression for the full model matrix of the jth term is...
  
     X^full_j = sum_q prod_i X_{ts[j]+i}[k[,ks[i]+q],]  

     - q runs from 0 to ks[i,2] - ks[i,1] - 1; i runs from 0 to dt[j] - 1.
 
     * ar_weights[0,2,4,6,...,2*n-2] is the ld of the square root of the tri-diagonal AR weight matrix
       ar_weights[1,3,5,...,2*n-1] is the sub-diagonal.
         
   Tensor product terms may have constraint matrices Z, which post multiply the tensor product 
   (typically imposing approximate sum-to-zero constraints). Actually Z is Q with the first column 
   dropped where Q =  I - vv'. qc[i]==0 for singleton terms.  

*/   
  int *pt, *pd,si,maxp=0,tri,r,c,rb,cb,rt,ct,pa,*tps,*tpsu,ptot,*b,*B,*C,*R,*sb,N,
    kk,kb,tid=0,nxwx=0,*worki,one=1;
  ptrdiff_t *off,*voff,mmp,q,qi=0,i,j;
  double *work,*ws=NULL,*Cost,*cost,*x0,*x1,*p0,*p1,x;
  unsigned long long ht[256];
  SM **sm,*SMstack;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  SMinihash(ht);
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  pd = (int *) CALLOC((size_t)*nt,sizeof(int)); /* storage for last marginal size */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvectors starts */
  tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts */
  tpsu = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the unconstrained term starts */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) {
      if (j==dt[i]-1) pd[i] = p[q]; /* the relevant dimension for deciding product ordering */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
    }
    if (qc[i]==0) voff[i+1] = voff[i]; else if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else {
      si = (int) round(v[voff[i]]); /* number of contrasts in this KP contrast */
      voff[i+1] = voff[i] + si + 2; 
    } /* start of ith v vector */
    if (maxp<pt[i]) maxp=pt[i];
    if (qc[i]==0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else if (qc[i]>0) tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
    else { /* Kronecker product of sum to zero contrasts */ 
      si = (int) round(v[voff[i]]); /* number of contrasts */
      tps[i+1] = tps[i] + pt[i] - (int) round(v[voff[i]+si+1]); /* subtracting number of constraints */
    }  
    tpsu[i+1] = tpsu[i] + pt[i]; /* where ith term starts in unconstrained param vector */ 
  }
  qi =  6 * *n; /* integer work space */

  worki = (int *)CALLOC((size_t)qi * *nthreads,sizeof(int));
  mmp = maxp;mmp = mmp*mmp;
  ptot = tps[*nt]; /* total number of parameters */
  nxwx = tpsu[*nt];
  if (*ar_stop>=0) { /* model has AR component*/
    /* ar_weights[0,2,4,6,...,2*n-2] is the ld of the square root of the tri-diagonal AR weight matrix
       ar_weights[1,3,5,...,2*n-1] is the sub-diagonal. */
    ws = w; /* input weights (pointer) */
    w = (double *)CALLOC((size_t) *n,sizeof(double)); /* need to avoid over-writing weights - so re-allocate */
    for (p0 = w,p1 = w + *n;p0<p1;p0++,ws++) *p0 = sqrt(*ws); /* sqrt weights */
    ws = (double *)CALLOC((size_t) 2 * *n -2,sizeof(double)); /* super and sub diagonals */
    for (i=0;i<*n-1;i++) ws[i+*n-1] = ws[i] = ar_weights[2*i+1] * ar_weights[2*i+2] * w[i+1] * w[i];
    for (i=0;i<*n-1;i++) w[i] *= (ar_weights[2*i+1]*ar_weights[2*i+1]+ar_weights[2*i]*ar_weights[2*i])*w[i]; 
    i = *n-1;w[i] *= w[i]*ar_weights[2*i]*ar_weights[2*i];
    /* so now w contains the leading diagonal of the tri-diagonal weight matrix,
       ws is the super diagonal and ws + n - 1 is the sub-diagonal. */
    tri = 1;
  } else tri = 0;
  sm = (SM **)CALLOC((size_t) *n * *nthreads,sizeof(SM *));
  SMstack = (SM *)CALLOC((size_t) 3 * *n * *nthreads,sizeof(SM));
  N = ((*nt + 1) * *nt)/2;
  C = (int *) CALLOC((size_t)N,sizeof(int)); R = (int *) CALLOC((size_t)N,sizeof(int));
  sb = (int *) CALLOC((size_t)N+1,sizeof(int)); /* at which sub-block does block start */
  Cost = (double *)CALLOC((size_t) N,sizeof(double));
  sb[0] = 0; /* start of first sub-block of block 0 */
  for (kk=r=0;r < *nt;r++) for (c=r;c< *nt;c++,kk++) { /* loop over main blocks */
    
    R[kk]=r;C[kk]=c;
    si = ks[ts[r] + *nx] - ks[ts[r]]; /* number of terms in summation convention for i */
    if (r==c && si==1 && !tri) {
      i = pt[r]/pd[r];
      i = i*(i+1)/2;
    } else {
      i = (pt[r]/pd[r])*(pt[c]/pd[c]);
    }  
    sb[kk+1] = sb[kk] + i; /* next block starts at this sub-block */
    /* compute rough cost per sub block figures */
    if (m[r]*m[c] < *n) { /* weight accumulation */
       Cost[kk] = m[r]*m[c]*(double) pd[c];
       x = m[r]*m[c]*(double) pd[r]; if (x<Cost[kk]) Cost[kk] = x;
    } else { /* direct accumulation */
      Cost[kk] = *n * pd[c];;
      x = *n * pd[r]; if (x < Cost[kk]) Cost[kk] = x;
    }
  }
  b = (int *) CALLOC((size_t)sb[N],sizeof(int));
  B = (int *) CALLOC((size_t)sb[N],sizeof(int));
  cost = (double *)CALLOC((size_t)sb[N],sizeof(double));
  for (kb=0,i=0;i<sb[N];i++) {
    b[i]=i;while (i>=sb[kb+1]) kb++; /* kb is main block */
    rb = R[kb];cb=C[kb];
    cost[i] = Cost[kb];
    B[i] = kb; 
  }  
  revsort(cost,b,sb[N]); /* R reverse sort on cost, to re-order b - see R.h*/
  q = XWXspace(N,sb,b,B,R,C,k,ks,m,p,pt,pd,*nx,*n,ts,dt,*nt,tri); /* compute the maximum workspace required per thread */
  work = (double *)CALLOC((size_t)q * *nthreads,sizeof(double)); /* allocate it */
  /* In what follows rb and cb are the whole term row column indices. r and c are the sub blocks within 
     the cross-product between two terms. The sub blocks arise when we have tensor product terms. The cleaner 
     design in which the sub-blocks are dealt with in XWXij does not load balance so well, hence this design.*/ 
  kb=0;
  #ifdef OPENMP_ON
  #pragma omp parallel for private(j,kb,kk,r,c,rb,cb,rt,ct,tid,i) num_threads(*nthreads) schedule(dynamic)
  #endif
  for (j=0;j<sb[N];j++) { /* the block loop */
    kk = b[j];kb=B[kk];
    rb = R[kb];cb=C[kb]; /* set up allows blocks to be computed in any order by re-arranging B */
    i = kk - sb[kb]; /* sub-block index */
    rt = pt[rb]/pd[rb]; ct = pt[cb]/pd[cb]; /* total rows and cols of sub-blocks */
    /* compute sub-row and column */
    if (sb[kb+1]-sb[kb]<rt*ct) { /* symmetric upper half only needed */ 
      r=0; while (i >= rt - r) { i -= rt - r;r++;}
      c = i + r;
    } else {
      r = i / ct;
      c = i % ct;
    }
   
    #ifdef OPENMP_ON
    tid = omp_get_thread_num(); /* needed for providing thread specific work space to XWXij */
    #endif
    XWXijs(XWX+tpsu[rb] +  (ptrdiff_t) nxwx * tpsu[cb],rb,cb,r,c,X,k,ks,m,p,*nx,*n,ts, dt,*nt,w,ws,
	   tri,off,work + tid * q,worki + tid * (ptrdiff_t) qi,nxwx,ht,
	   sm + tid *  *n,SMstack + 3 * tid * *n ); /* compute r,c block */
    /* NOTE: above will write directly to oversized XWX, then have constraints applied post-hoc. */ 
  } /* block loop */

  FREE(work);
  work = (double *)CALLOC((size_t) 2*nxwx,sizeof(double)); /* working for Ztb */
  
  /* now XWX contains the unconstrained X'WX, but the constraints have to be applied to blocks involving tensor products */
  for (r=0;r < *nt;r++) for (c=r;c< *nt;c++) {
    /* if Xr is tensor, may need to apply constraint */
    if (dt[r]>1&&qc[r]!=0) { /* first term is a tensor with a constraint */
      /* col by col form (I-vv')xwx, dropping first row... */
      /* col by col form Z' xwx where Z is constraint matrix */
      for (j=0;j<pt[c];j++) { /* loop over columns */
        x0 = XWX + tpsu[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx;   /* jth col of constrained block */
	Ztb(x1,x0,v+voff[r],qc+r,&one,pt+r,work);
      }
      pa = pt[r]-1;
    } else {
      pa = pt[r];
      if (tpsu[r]!=tps[r]) { /* still need to shift rows upwards */
        for (j=0;j<pt[c];j++) { /* loop over columns */
	  x0 = XWX + tpsu[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx; /* jth col of raw block */
	  x1 = XWX + tps[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx;   /* jth col of shifted block */
	  for (p1 = x0 + pt[r];x0<p1;x0++,x1++) *x1 = *x0;
	}  
      }	
    }  
    if (dt[c]>1&&qc[c]!=0) { /* Xc term is a tensor with a constraint */
      /* row by row form xwx(I-vv') dropping first col... */
      /* row by row form xwx Z, where Z is constraint matrix */
      for (j=0;j<pa;j++) { /* work down rows */
        x0 = XWX + tps[r] + j + tpsu[c] * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + j + tps[c] * (ptrdiff_t) nxwx;   /* jth col of constrained block */
	Ztb(x1,x0,v+voff[c],qc+c,&nxwx,pt+c,work);
      }
    } else if (tpsu[c]!=tps[c]) { /* still need to shift cols leftwards */
      for (j=0;j<pa;j++) { /* work down rows */
	x0 = XWX + tps[r] + j + tpsu[c] * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + j + tps[c] * (ptrdiff_t) nxwx;   /* jth col of shifted block */
	for (p1 = x0 + nxwx * pt[c];x0<p1;x0 += nxwx,x1 += nxwx) *x1 = *x0;
      }
    }  
  }

  if (ptot<nxwx) row_squash(XWX,ptot,nxwx,ptot); /* drop the now redundant trailing rows */
  up2lo(XWX,ptot); /* copy upper triangle to lower */
  
  FREE(pt);FREE(pd);FREE(off);FREE(voff);FREE(tps);FREE(tpsu);FREE(work); if (tri) {FREE(ws);free(w);}//FREE(xwx);FREE(xwx0);
  FREE(B);FREE(R);FREE(C);FREE(sb);FREE(Cost);FREE(cost);FREE(b);FREE(sm);FREE(SMstack);FREE(worki);
} /* XWXd0 */ 



void XVXd0(double *XWX,double *X,double *e,int *k,int *ks, int *m,int *p, ptrdiff_t  *n, int *nx, int *ts, 
	   int *dt, int *nt,double *v,int *qc,int *nthreads,int *a,ptrdiff_t *ma) {
/* Version of XWXd0 in which (diagonal or tri-diagonal) W is replaced by V such that 
   V[i,j] = e[i]*e[j] if j in nei(i) and 0 otherwise. 


   Forms Xt'VXt when Xt is divided into blocks of columns, each stored in compact form
   using arguments X and k.  
   * 'X' contains 'nx' blocks, the ith is an m[i] by p[i] matrix containing the unique rows of 
     the ith marginal model matrix. 
   * There are 'nt' model terms. Each term is made up of one or more maginal model matrices.
   * The jth term starts at block ts[j] of X, and has dt[j] marginal matrices. The terms model
     matrix is the row tensor product of its full (n row) marginals.
   * The index vectors converting the unique row matrices to full marginal matrices are in 
     'k', an n-row matrix of integers. Conceptually if Xj and kj represent the jth unique 
      row matrix and index vector then the ith row of the corresponding full marginal matrix 
      is Xj[kj[i],], but things are more complicated when each full term matrix is actually 
      the sum of several matrices (summation convention).
   * To handle the summation convention, each marginal matrix can have several index vectors. 
     'ks' is an nx by 2 matrix giving the columns of k corresponding to the ith marginal 
     model matrix. Specifically columns ks[i,1]:(ks[i,2]-1) of k are the index vectors for the ith 
     marginal. All marginals corresponding to one term must have the same number of index columns.
     The full model matrix for the jth term is constucted by summing over q the full 
     model matrices corresponding to the qth index vectors for each of its marginals.    
   * For example the exression for the full model matrix of the jth term is...
  
     X^full_j = sum_q prod_i X_{ts[j]+i}[k[,ks[i]+q],]  

     - q runs from 0 to ks[i,2] - ks[i,1] - 1; i runs from 0 to dt[j] - 1.
 
   * a[ma[i-1]+1:ma[i]] are the neighbours of point i (ma[-1]+1 = 0 by convention).
         
   Tensor product terms may have constraint matrices Z, which post multiply the tensor product 
   (typically imposing approximate sum-to-zero constraints). Actually Z is Q with the first column 
   dropped where Q =  I - vv'. qc[i]==0 for singleton terms.  

*/   
  int *pt, *pd,si,maxp=0,r,c,rb,cb,rt,ct,pa,*tps,*tpsu,ptot,*b,*B,*C,*R,*sb,N,
    kk,kb,tid=0,nxwx=0,*worki,one=1;
  ptrdiff_t *off,*voff,mmp,q,qi=0,i,j;
  double *work,*Cost,*cost,*x0,*x1,*p1,x;
  //unsigned long long ht[256];
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  //SMinihash(ht);
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  pd = (int *) CALLOC((size_t)*nt,sizeof(int)); /* storage for last marginal size */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvectors starts */
  tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts */
  tpsu = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the unconstrained term starts */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) {
      if (j==dt[i]-1) pd[i] = p[q]; /* the relevant dimension for deciding product ordering */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
    }
    if (qc[i]==0) voff[i+1] = voff[i]; else if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else {
      si = (int) round(v[voff[i]]); /* number of contrasts in this KP contrast */
      voff[i+1] = voff[i] + si + 2; 
    } /* start of ith v vector */
    if (maxp<pt[i]) maxp=pt[i];
    if (qc[i]==0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else if (qc[i]>0) tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
    else { /* Kronecker product of sum to zero contrasts */ 
      si = (int) round(v[voff[i]]); /* number of contrasts */
      tps[i+1] = tps[i] + pt[i] - (int) round(v[voff[i]+si+1]); /* subtracting number of constraints */
    }  
    tpsu[i+1] = tpsu[i] + pt[i]; /* where ith term starts in unconstrained param vector */ 
  }
  qi =  6 * *n; /* integer work space */

  worki = (int *)CALLOC((size_t)qi * *nthreads,sizeof(int));
  mmp = maxp;mmp = mmp*mmp;
  ptot = tps[*nt]; /* total number of parameters */
  nxwx = tpsu[*nt];

  N = ((*nt + 1) * *nt)/2;
  C = (int *) CALLOC((size_t)N,sizeof(int)); R = (int *) CALLOC((size_t)N,sizeof(int));
  sb = (int *) CALLOC((size_t)N+1,sizeof(int)); /* at which sub-block does block start */
  Cost = (double *)CALLOC((size_t) N,sizeof(double));
  sb[0] = 0; /* start of first sub-block of block 0 */
  for (kk=r=0;r < *nt;r++) for (c=r;c< *nt;c++,kk++) { /* loop over main blocks */
    
    R[kk]=r;C[kk]=c;
    si = ks[ts[r] + *nx] - ks[ts[r]]; /* number of terms in summation convention for i */

    i = (pt[r]/pd[r])*(pt[c]/pd[c]);
   
    sb[kk+1] = sb[kk] + i; /* next block starts at this sub-block */
    /* compute rough cost per sub block figures */
    if (m[r]*m[c] < *n) { /* weight accumulation */
       Cost[kk] = m[r]*m[c]*(double) pd[c];
       x = m[r]*m[c]*(double) pd[r]; if (x<Cost[kk]) Cost[kk] = x;
    } else { /* direct accumulation */
      Cost[kk] = *n * pd[c];;
      x = *n * pd[r]; if (x < Cost[kk]) Cost[kk] = x;
    }
  }
  b = (int *) CALLOC((size_t)sb[N],sizeof(int));
  B = (int *) CALLOC((size_t)sb[N],sizeof(int));
  cost = (double *)CALLOC((size_t)sb[N],sizeof(double));
  for (kb=0,i=0;i<sb[N];i++) {
    b[i]=i;while (i>=sb[kb+1]) kb++; /* kb is main block */
    rb = R[kb];cb=C[kb];
    cost[i] = Cost[kb];
    B[i] = kb; 
  }  
  revsort(cost,b,sb[N]); /* R reverse sort on cost, to re-order b - see R.h*/
  i=0;
  q = XWXspace(N,sb,b,B,R,C,k,ks,m,p,pt,pd,*nx,*n,ts,dt,*nt,i); /* compute the maximum workspace required per thread */
  work = (double *)CALLOC((size_t)q * *nthreads,sizeof(double)); /* allocate it */
  /* In what follows rb and cb are the whole term row column indices. r and c are the sub blocks within 
     the cross-product between two terms. The sub blocks arise when we have tensor product terms. The cleaner 
     design in which the sub-blocks are dealt with in XWXij does not load balance so well, hence this design.*/ 
  kb=0;
  #ifdef OPENMP_ON
  #pragma omp parallel for private(j,kb,kk,r,c,rb,cb,rt,ct,tid,i) num_threads(*nthreads) schedule(dynamic)
  #endif
  for (j=0;j<sb[N];j++) { /* the block loop */
    kk = b[j];kb=B[kk];
    rb = R[kb];cb=C[kb]; /* set up allows blocks to be computed in any order by re-arranging B */
    i = kk - sb[kb]; /* sub-block index */
    rt = pt[rb]/pd[rb]; ct = pt[cb]/pd[cb]; /* total rows and cols of sub-blocks */
    /* compute sub-row and column */
    if (sb[kb+1]-sb[kb]<rt*ct) { /* symmetric upper half only needed */ 
      r=0; while (i >= rt - r) { i -= rt - r;r++;}
      c = i + r;
    } else {
      r = i / ct;
      c = i % ct;
    }
   
    #ifdef OPENMP_ON
    tid = omp_get_thread_num(); /* needed for providing thread specific work space to XWXij */
    #endif
    XVXijs(XWX+tpsu[rb] +  (ptrdiff_t) nxwx * tpsu[cb],rb,cb,r,c,X,k,ks,m,p,*nx,*n,ts, dt,*nt,e,a,ma,
	   off,work + tid * q,worki + tid * (ptrdiff_t) qi,nxwx); /* compute r,c block */
    /* NOTE: above will write directly to oversized XWX, then have constraints applied post-hoc. */ 
  } /* block loop */

  FREE(work);
  work = (double *)CALLOC((size_t) 2*nxwx,sizeof(double)); /* working for Ztb */
  
  /* now XWX contains the unconstrained X'WX, but the constraints have to be applied to blocks involving tensor products */
  for (r=0;r < *nt;r++) for (c=r;c< *nt;c++) {
    /* if Xr is tensor, may need to apply constraint */
    if (dt[r]>1&&qc[r]!=0) { /* first term is a tensor with a constraint */
      /* col by col form (I-vv')xwx, dropping first row... */
      /* col by col form Z' xwx where Z is constraint matrix */
      for (j=0;j<pt[c];j++) { /* loop over columns */
        x0 = XWX + tpsu[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx;   /* jth col of constrained block */
	Ztb(x1,x0,v+voff[r],qc+r,&one,pt+r,work);
      }
      pa = pt[r]-1;
    } else {
      pa = pt[r];
      if (tpsu[r]!=tps[r]) { /* still need to shift rows upwards */
        for (j=0;j<pt[c];j++) { /* loop over columns */
	  x0 = XWX + tpsu[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx; /* jth col of raw block */
	  x1 = XWX + tps[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx;   /* jth col of shifted block */
	  for (p1 = x0 + pt[r];x0<p1;x0++,x1++) *x1 = *x0;
	}  
      }	
    }  
    if (dt[c]>1&&qc[c]!=0) { /* Xc term is a tensor with a constraint */
      /* row by row form xwx(I-vv') dropping first col... */
      /* row by row form xwx Z, where Z is constraint matrix */
      for (j=0;j<pa;j++) { /* work down rows */
        x0 = XWX + tps[r] + j + tpsu[c] * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + j + tps[c] * (ptrdiff_t) nxwx;   /* jth col of constrained block */
	Ztb(x1,x0,v+voff[c],qc+c,&nxwx,pt+c,work);
      }
    } else if (tpsu[c]!=tps[c]) { /* still need to shift cols leftwards */
      for (j=0;j<pa;j++) { /* work down rows */
	x0 = XWX + tps[r] + j + tpsu[c] * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + j + tps[c] * (ptrdiff_t) nxwx;   /* jth col of shifted block */
	for (p1 = x0 + nxwx * pt[c];x0<p1;x0 += nxwx,x1 += nxwx) *x1 = *x0;
      }
    }  
  }

  if (ptot<nxwx) row_squash(XWX,ptot,nxwx,ptot); /* drop the now redundant trailing rows */
  up2lo(XWX,ptot); /* copy upper triangle to lower */
  
  FREE(pt);FREE(pd);FREE(off);FREE(voff);FREE(tps);FREE(tpsu);FREE(work);
  FREE(B);FREE(R);FREE(C);FREE(sb);FREE(Cost);FREE(cost);FREE(b);FREE(worki);
} /* XVXd0 */ 

SEXP CXVXd0(SEXP XWXr, SEXP Xr, SEXP er, SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP tsr, SEXP dtr,
	    SEXP vr,SEXP qcr, SEXP nthreadsr, SEXP ar, SEXP mar) {
/* .Call wrapper for XVXd0 allowing R long vector storage for k and a. Note that 
  this does not allow more than maxint data - that would require re-writting R code 
  to avoid storing k in a matrix (which is only allowed maxint rows).

  n is number of rows of k, nx is length of m or p, nt is the length of ts or dt.
*/
  double *XWX,*X,*v,*mad=NULL,*e;
  int *k,*ks,*m,*p,*ts,*dt,*qc,*nthreads,nx,nt,*a,*mai=NULL,nma,i;
  ptrdiff_t n,*ma;
  n = (ptrdiff_t) nrows(kr); XWX = REAL(XWXr);
  X = REAL(Xr);
  k = INTEGER(kr); ks = INTEGER(ksr);
  m = INTEGER(mr); nx = length(mr);
  p = INTEGER(pr); e = REAL(er);
  a= INTEGER(ar);nma = length(mar);
  ma = (ptrdiff_t *)CALLOC((size_t)nma,sizeof(ptrdiff_t)); 
  if (isReal(mar)) { /* case where a is a long vector so mar[i] with be passed as double */
    mad=REAL(mar);
    for (i=0;i<nma;i++) ma[i] = (ptrdiff_t) mad[i];
  } else { /* case where mar[i] all fit in integers */
    mai=INTEGER(mar);
    for (i=0;i<nma;i++) ma[i] = (ptrdiff_t) mai[i];
  }  

  ts = INTEGER(tsr); dt = INTEGER(dtr); nt = length(tsr);
  v = REAL(vr);qc = INTEGER(qcr);
  nthreads = INTEGER(nthreadsr);
  XVXd0(XWX,X,e,k,ks,m,p,&n,&nx,ts,dt,&nt,v,qc,nthreads,a,ma);
  FREE(ma);
  return(R_NilValue);
} /*  CXVXd0 */



SEXP CXWXd1(SEXP XWXr, SEXP Xr, SEXP wr, SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP tsr, SEXP dtr,
	    SEXP vr,SEXP qcr, SEXP nthreadsr, SEXP ar_stopr, SEXP ar_weightsr,
	    SEXP rsr,SEXP csr) {
/* .Call wrapper for XWXd1 allowing R long vector storage for k. Note that 
  this does not allow more than maxint data - that would require re-writting R code 
  to avoid storing k in a matrix (which is only allowed maxint rows).

  n is number of rows of k, nx is length of m or p, nt is the length of ts or dt, 
  ncs and nrs are length of cs/rs.
*/
  double *XWX,*X,*v,*w,*ar_weights;
  int *k,*ks,*m,*p,*ts,*dt,*qc,*nthreads,*cs,*rs,nx,nt,ncs,nrs,*ar_stop;
  ptrdiff_t n;
  n = (ptrdiff_t)nrows(kr); XWX = REAL(XWXr);
  X = REAL(Xr);w = REAL(wr);
  k = INTEGER(kr); ks = INTEGER(ksr);
  m = INTEGER(mr); nx = length(mr);
  p = INTEGER(pr);
  ar_stop = INTEGER(ar_stopr);
  ar_weights = REAL(ar_weightsr);
  ts = INTEGER(tsr); dt = INTEGER(dtr); nt = length(tsr);
  v = REAL(vr);qc = INTEGER(qcr);
  nthreads = INTEGER(nthreadsr);
  cs = INTEGER(csr); rs = INTEGER(rsr);
  nrs = length(rsr); ncs = length(csr);
  if (*cs<0) ncs = 0;
  if (*rs<0) nrs = 0;
  XWXd1(XWX,X,w,k,ks,m,p,&n,&nx,ts,dt,&nt,v,qc,nthreads,ar_stop,ar_weights,
	rs,cs,&nrs,&ncs);
  return(R_NilValue);
} /*  CXWXd1 */



void XWXd1(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, ptrdiff_t *n, int *nx, int *ts, 
	   int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,double *ar_weights,
	   int *rs, int *cs, int *nrs, int *ncs) {
/* essentially a driver routine for XWXij implementing block oriented cross products

   This version allows selection of only some subsets of term blocks:
   * rs is the length nrs array of required row indices.
   * cs is the length ncs array of required col indices.
   * if nrs and ncs are <= zero then the full set or nt r/c is used (arrays should then be dim nt). 
     If only one is > zero, then its array denotes the rows and cols required.   
   - the selected blocks are returned as a dense compact matrix (rather than being inserted into the full X'WX 
     in the appropriate places, for example.) 
   
   This version has the looping over sub-blocks, associated with tensor product terms, located in this routine
   to get better load balancing.
   
   Requires XWX to be over-sized on entry - namely n.params + n.terms by n.params + n.terms instead of
   n.params by n.params.

   Forms Xt'WXt when Xt is divided into blocks of columns, each stored in compact form
   using arguments X and k. 
   * 'X' contains 'nx' blocks, the ith is an m[i] by p[i] matrix containing the unique rows of 
     the ith marginal model matrix. 
   * There are 'nt' model terms. Each term is made up of one or more maginal model matrices.
   * The jth term starts at block ts[j] of X, and has dt[j] marginal matrices. The terms model
     matrix is the row tensor product of its full (n row) marginals.
   * The index vectors converting the unique row matrices to full marginal matrices are in 
     'k', an n-row matrix of integers. Conceptually if Xj and kj represent the jth unique 
      row matrix and index vector then the ith row of the corresponding full marginal matrix 
      is Xj[kj[i],], but things are more complicated when each full term matrix is actually 
      the sum of several matrices (summation convention).
   * To handle the summation convention, each marginal matrix can have several index vectors. 
     'ks' is an nx by 2 matrix giving the columns of k corresponding to the ith marginal 
     model matrix. Specifically columns ks[i,1]:(ks[i,2]-1) of k are the index vectors for the ith 
     marginal. All marginals corresponding to one term must have the same number of index columns.
     The full model matrix for the jth term is constucted by summing over q the full 
     model matrices corresponding to the qth index vectors for each of its marginals.    
   * For example the exression for the full model matrix of the jth term is...
  
     X^full_j = sum_q prod_i X_{ts[j]+i}[k[,ks[i]+q],]  

     - q runs from 0 to ks[i,2] - ks[i,1] - 1; i runs from 0 to dt[j] - 1.
         
   Tensor product terms may have constraint matrices Z, which post multiply the tensor product 
   (typically imposing approximate sum-to-zero constraints). Actually Z is Q with the first column 
   dropped where Q =  I - vv'. qc[i]==0 for singleton terms.  


*/   
  int *pt, *pd,ri,ci,si,maxp=0,nxwx=0,tri,r,c,rb,cb,rt,ct,pa,*tpsr,*tpsur,*tpsc,*tpsuc,ptot,
    *b,*B,*C,*R,*sb,N,kk,kb,tid=0,*worki,symmetric=1,one=1;
  ptrdiff_t *off,*voff,mmp,q,i,j,qi=0;
  double *work,*ws=NULL,*Cost,*cost,*x0,*x1,*p0,*p1,x;
  unsigned long long ht[256];
  SM **sm,*SMstack;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  SMinihash(ht);
  /* check row/col subset arrays... */
  if (*nrs <= 0) {
    if (*ncs > 0) {
      *nrs = *ncs;rs=cs;
    } else {
      for (i=0;i<*nt;i++) rs[i] = cs[i] = i;
    }  
  } else {
    if (*ncs>0) symmetric = 0; else {
      *ncs = *nrs; cs = rs;
    }  
  }  
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  pd = (int *) CALLOC((size_t)*nt,sizeof(int)); /* storage for last marginal size */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvectors starts */
  tpsr = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts, row set */
  tpsur = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the unconstrained term starts, row set */
  tpsc = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts, col set */
  tpsuc = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the unconstrained term starts, col set */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) { /* loop through the marginals for this model term */
      if (j==dt[i]-1) pd[i] = p[q]; /* the relevant dimension for deciding product ordering */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
    } 
    if (qc[i]==0) voff[i+1] = voff[i]; else if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else {
      ri = (int) round(v[voff[i]]); /* number of contrasts in this KP contrast */
      voff[i+1] = voff[i] + ri + 2; 
    } /* start of ith v vector */
    if (maxp<pt[i]) maxp=pt[i];
  }
  for (kk=i=j=0;j<*nrs;j++) {
    r = rs[j];
    tpsr[r] = kk;tpsur[r] = i;
    if (qc[r]==0) kk += pt[r]; /* where rth terms starts in param vector */ 
    else if (qc[r]>0) kk += pt[r] - 1; /* there is a tensor constraint to apply - reducing param count*/
    else { /* Kronecker product of sum-to-zero contrasts */
      ri = (int) round(v[voff[r]]); /* number of contrasts */
      kk += pt[r] - (int) round(v[voff[r]+ri+1]); /* subtracting number of constraints */
    }	
    i += pt[r]; /* where rth term starts in unconstrained param vector */ 
  }
  ptot = kk; /* rows of computed XWX post constraint */
  nxwx = i; /* rows of computed XWX */
  for (kk=i=j=0;j<*ncs;j++) {
    c = cs[j];
    tpsc[c] = kk;tpsuc[c] = i;
    if (qc[c]==0) kk += pt[c]; /* where cth terms starts in param vector */ 
    else if (qc[c]>0) kk += pt[c] - 1; /* there is a tensor constraint to apply - reducing param count*/
    else { /* Kronecker product of sum-to-zero contrasts */
      ri = (int) round(v[voff[c]]); /* number of contrasts */
      kk += pt[c] - (int) round(v[voff[c]+ri+1]); /* subtracting number of constraints */
    }	
    i +=  pt[c]; /* where cth term starts in unconstrained param vector */ 
  }
    
  qi =  6 * *n; /* integer work space */
  worki = (int *)CALLOC((size_t)qi * *nthreads,sizeof(int));
  mmp = maxp;mmp = mmp*mmp;
 
  if (*ar_stop>=0) { /* model has AR component*/
    /* ar_weights[0,2,4,6,...,2*n-2] is the ld of the square root of the tri-diagonal AR weight matrix
       ar_weights[1,3,5,...,2*n-1] is the sub-diagonal. */
    ws = w; /* input weights (pointer) */
    w = (double *)CALLOC((size_t) *n,sizeof(double)); /* need to avoid over-writing weights - so re-allocate */
    for (p0 = w,p1 = w + *n;p0<p1;p0++,ws++) *p0 = sqrt(*ws); /* sqrt weights */
    ws = (double *)CALLOC((size_t) 2 * *n -2,sizeof(double)); /* super and sub diagonals */
    for (i=0;i<*n-1;i++) ws[i+*n-1] = ws[i] = ar_weights[2*i+1] * ar_weights[2*i+2] * w[i+1] * w[i];
    for (i=0;i<*n-1;i++) w[i] *= (ar_weights[2*i+1]*ar_weights[2*i+1]+ar_weights[2*i]*ar_weights[2*i])*w[i]; 
    i = *n-1;w[i] *= w[i]*ar_weights[2*i]*ar_weights[2*i];
    /* so now w contains the leading diagonal of the tri-diagonal weight matrix,
       ws is the super diagonal and ws + n - 1 is the sub-diagonal. */
    tri = 1;
  } else tri = 0;
  /* memory allocation for sparse accumulation */
  sm = (SM **)CALLOC((size_t) *n * *nthreads,sizeof(SM *));
  SMstack = (SM *)CALLOC((size_t) 3 * *n * *nthreads,sizeof(SM));
  if (symmetric) N = ((*nrs + 1) * *nrs)/2; else N = *nrs * *ncs;
  C = (int *) CALLOC((size_t)N,sizeof(int)); R = (int *) CALLOC((size_t)N,sizeof(int));
  sb = (int *) CALLOC((size_t)N+1,sizeof(int)); /* at which sub-block does block start */
  Cost = (double *)CALLOC((size_t) N,sizeof(double));
  sb[0] = 0; /* start of first sub-block of block 0 */
  for (kk=ri=0;ri < *nrs;ri++) {
    if (symmetric) ci=ri; else ci = 0;
    for (;ci< *ncs;ci++,kk++) { /* loop over main blocks */
      r = rs[ri];c=cs[ci]; /* selected main block row/col */
      R[kk]=r;C[kk]=c;
      si = ks[ts[r] + *nx] - ks[ts[r]]; /* number of terms in summation convention for i */
      if (symmetric && r==c && si==1 && !tri) {
        i = pt[r]/pd[r];
        i = i*(i+1)/2;
      } else {
        i = (pt[r]/pd[r])*(pt[c]/pd[c]);
      }  
      sb[kk+1] = sb[kk] + i; /* next block starts at this sub-block */
      /* compute rough cost per sub block figures */
      if (m[r]*m[c] < *n) { /* weight accumulation */
         Cost[kk] = m[r]*m[c]*(double) pd[c];
         x = m[r]*m[c]*(double) pd[r]; if (x<Cost[kk]) Cost[kk] = x;
      } else { /* direct accumulation */
        Cost[kk] = *n * pd[c];;
        x = *n * pd[r]; if (x < Cost[kk]) Cost[kk] = x;
      }
    }  
  }
  b = (int *) CALLOC((size_t)sb[N],sizeof(int)); /* index of all subblocks for sorting to improve load balance (expensive first) */
  B = (int *) CALLOC((size_t)sb[N],sizeof(int)); /* index the main block to which each subblock belongs */
  cost = (double *)CALLOC((size_t)sb[N],sizeof(double)); /* approx sub-block cost */
  for (kb=0,i=0;i<sb[N];i++) { /* loop over all sub-blocks */
    b[i]=i;while (i>=sb[kb+1]) kb++; /* kb is main block */
    rb = R[kb];cb=C[kb]; /* which main block row and column are we at? */
    cost[i] = Cost[kb];
    B[i] = kb; /* record the main block we are in */
  }  
  revsort(cost,b,sb[N]); /* R reverse sort on cost, to re-order b - see R.h*/
  q = XWXspace(N,sb,b,B,R,C,k,ks,m,p,pt,pd,*nx,*n,ts,dt,*nt,tri); /* compute the maximum workspace required per thread */
  work = (double *)CALLOC((size_t)q * *nthreads,sizeof(double)); /* allocate it */
  /* In what follows rb and cb are the whole term row column indices. r and c are the sub blocks within 
     the cross-product between two terms. The sub blocks arise when we have tensor product terms. The cleaner 
     design in which the sub-blocks are dealt with in XWXij does not load balance so well, hence this design.*/ 
  kb=0;
  #ifdef OPENMP_ON
  #pragma omp parallel for private(j,kb,kk,r,c,rb,cb,rt,ct,tid,i) num_threads(*nthreads) schedule(dynamic)
  #endif
  for (j=0;j<sb[N];j++) { /* the all sub-blocks loop */
    kk = b[j];kb=B[kk]; /* kb is main block */
    rb = R[kb];cb=C[kb]; /* set up allows blocks to be computed in any order by re-arranging B */
    i = kk - sb[kb]; /* sub-block index - i.e. where is this subblock within its block */
    rt = pt[rb]/pd[rb]; ct = pt[cb]/pd[cb]; /* total rows and cols of sub-blocks */
    /* compute sub-row and column */
    if (symmetric && sb[kb+1]-sb[kb]<rt*ct) { /* symmetric upper half only needed */ 
      r=0; while (i >= rt - r) { i -= rt - r;r++;}
      c = i + r;
    } else {
      r = i / ct;
      c = i % ct;
    }
   
    #ifdef OPENMP_ON
    tid = omp_get_thread_num(); /* needed for providing thread specific work space to XWXij */
    #endif
    XWXijs(XWX+tpsur[rb] +  (ptrdiff_t) nxwx * tpsuc[cb],rb,cb,r,c,X,k,ks,m,p,*nx,*n,ts, dt,*nt,w,ws,
	   tri,off,work + tid * q,worki + tid * (ptrdiff_t) qi,nxwx,ht,
	   sm + tid * *n,SMstack + 3 * tid *  *n ); /* compute r,c block */
    /* NOTE: above will write directly to oversized XWX, then have constraints applied post-hoc. */ 
  } /* block loop */

  FREE(work);
  work = (double *)CALLOC((size_t) 2*nxwx,sizeof(double)); /* working for Ztb */

  /* now XWX contains the unconstrained X'WX, but the constraints have to be applied to blocks involving tensor products */
  for (ri=0;ri < *nrs;ri++) { /* loop over required block rows */
    if (symmetric) ci=ri; else ci = 0;
    for (;ci< *ncs;ci++) { /* and over required block cols */
      /* if Xr is tensor, may need to apply constraint */
      r = rs[ri];c = cs[ci];
      if (dt[r]>1&&qc[r]!=0) { /* first term is a tensor with a constraint */
      /* col by col form (I-vv')xwx, dropping first row... */
      /* col by col form Z' xwx where Z is constraint matrix */
        for (j=0;j<pt[c];j++) { /* loop over columns */
          x0 = XWX + tpsur[r] + (tpsuc[c]+j) * (ptrdiff_t) nxwx; /* jth col of raw block */
	  x1 = XWX + tpsr[r] + (tpsuc[c]+j) * (ptrdiff_t) nxwx;   /* jth col of constrained block */
	  Ztb(x1,x0,v+voff[r],qc+r,&one,pt+r,work);
        }
        pa = pt[r]-1; /* number of block rows after constraint */ 
      } else { /* block may still need to be shifted because of prior constraints */
        pa = pt[r];
        if (tpsur[r]!=tpsr[r]) { /* still need to shift rows upwards */
          for (j=0;j<pt[c];j++) { /* loop over columns */
	    x0 = XWX + tpsur[r] + (tpsuc[c]+j) * (ptrdiff_t) nxwx; /* jth col of raw block */
	    x1 = XWX + tpsr[r] + (tpsuc[c]+j) * (ptrdiff_t) nxwx;   /* jth col of shifted block */
	    for (p1 = x0 + pt[r];x0<p1;x0++,x1++) *x1 = *x0;
	  }  
        }	
      }  
      if (dt[c]>1&&qc[c]!=0) { /* Xc term is a tensor with a constraint */
        /* row by row form xwx(I-vv') dropping first col... */
	/* row by row form xwx Z, where Z is constraint matrix */
        for (j=0;j<pa;j++) { /* work down rows */
          x0 = XWX + tpsr[r] + j + tpsuc[c] * (ptrdiff_t) nxwx; /* jth col of raw block */
	  x1 = XWX + tpsr[r] + j + tpsc[c] * (ptrdiff_t) nxwx;   /* jth col of constrained block */
	  Ztb(x1,x0,v+voff[c],qc+c,&nxwx,pt+c,work);
        }
      } else if (tpsuc[c]!=tpsc[c]) { /* still need to shift cols leftwards */
        for (j=0;j<pa;j++) { /* work down rows */
	  x0 = XWX + tpsr[r] + j + tpsuc[c] * (ptrdiff_t) nxwx; /* jth col of raw block */
	  x1 = XWX + tpsr[r] + j + tpsc[c] * (ptrdiff_t) nxwx;   /* jth col of shifted block */
	  for (p1 = x0 + nxwx * pt[c];x0<p1;x0 += nxwx,x1 += nxwx) *x1 = *x0;
        }
      }
    } /* ci loop */
  } /* constraint loop */

  if (ptot<nxwx) row_squash(XWX,ptot,nxwx,ptot); /* drop the now redundant trailing rows */
  if (symmetric) up2lo(XWX,ptot); /* copy upper triangle to lower */
  
  FREE(pt);FREE(pd);FREE(off);FREE(voff);FREE(tpsr);FREE(tpsur);FREE(tpsc);FREE(tpsuc);
  FREE(work); if (tri) {FREE(ws);FREE(w);}
  FREE(B);FREE(R);FREE(C);FREE(sb);FREE(Cost);FREE(cost);FREE(b);FREE(sm);FREE(SMstack);FREE(worki);
} /* XWXd1 */ 


void fill_lt(double *A,int k) {
/* fill lower triangle of symmetric A from upper triangle */
  int i;
  double *Al,*Au,*Au0,*Al0,*Au1;
  //for (i=0;i<k;i++) for (j=0;j<i;j++) A[i + (j-1)*k] = A[j+(i-1)*k];
  for (Al0=Au0=A,i=0;i<k;i++,Au0+=k,Al0++) for (Au=Au0,Au1=Au+i,Al=Al0;Au<Au1;Au++,Al+=k) *Al = *Au;
} /* fill_lt */  

SEXP CNCV(SEXP NCVr, SEXP NCV1r, SEXP NCV2r, SEXP Gr, SEXP rsdr, SEXP betar, SEXP beta1r, SEXP wr,SEXP spr,
	    SEXP Xr,SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP tsr, SEXP dtr,
	    SEXP vr,SEXP qcr, SEXP nthreadsr, SEXP nei, SEXP Sr) {
/* wrapper function for ncvd() suitable for .Call-ing from R. */
  double *NCV,*NCV1,*NCV2,*w,*X,*v,*beta,*G,*rsd,*sp,**S,*db,*mar,*mdr;
  int *a,*d,*k,*ks,*m,nx,*p,*ts,*dt,nt,*qc,*nthreads,nprot=0,nS,ns,nsp,*sr,*soff,i,j,q,pg,nn,ck;
  SEXP K,M,I,MI,Sl,Si;
  ptrdiff_t n,*ma,*md,*p1,*p2;
  n = (ptrdiff_t)nrows(kr);
  ck = (int) ncols(kr);
  pg = (int) nrows(Gr);
  G = REAL(Gr); X = REAL(Xr);w = REAL(wr);
  beta = REAL(betar);rsd= REAL(rsdr);db = REAL(beta1r);
  k = INTEGER(kr); ks = INTEGER(ksr);
  m = INTEGER(mr); nx = length(mr);
  p = INTEGER(pr);sp = REAL(spr);nsp=length(spr);
  ts = INTEGER(tsr); dt = INTEGER(dtr); nt = length(tsr);
  v = REAL(vr);qc = INTEGER(qcr);
  nthreads = INTEGER(nthreadsr);
  NCV = REAL(NCVr);NCV1 = REAL(NCV1r);NCV2 = REAL(NCV2r);
  /* unpack neighbourhood list, nei */
  K = getListEl(nei,"a"); // was k
  K = PROTECT(coerceVector(K,INTSXP));nprot++;
  a = INTEGER(K);
  M = getListEl(nei,"ma"); // was m
  M = PROTECT(coerceVector(M,REALSXP));nprot++;
  mar = REAL(M);nn = length(M);
  ma = (ptrdiff_t *)CALLOC(nn,sizeof(ptrdiff_t));
  for (p1=ma,p2=p1+nn;p1<p2;p1++,mar++) *p1 = (ptrdiff_t) *mar;
  I = getListEl(nei,"d"); // was i
  I = PROTECT(coerceVector(I,INTSXP));nprot++;
  d = INTEGER(I);
  MI = getListEl(nei,"md"); // was mi
  MI = PROTECT(coerceVector(MI,REALSXP));nprot++;
  mdr = REAL(MI);
  md = (ptrdiff_t *)CALLOC(nn,sizeof(ptrdiff_t));
  for (p1=md,p2=p1+nn;p1<p2;p1++,mdr++) *p1 = (ptrdiff_t) *mdr;
  /* now get the smoothing penalty matrix structure, Sr */
  nS = length(Sr);
  S = (double **)CALLOC(nsp,sizeof(double *));
  soff = (int *)CALLOC(nsp,sizeof(int));
  sr = (int*) CALLOC(nsp,sizeof(int));
  for (ns=i=0;i<nS;i++) {
    Sl = VECTOR_ELT(Sr, i);
    soff[ns] = asInteger(getListEl(Sl,"start")) - 1;
    sr[ns] = asInteger(getListEl(Sl,"stop")) - soff[ns];
    Si = getListEl(Sl,"S0");
    q = length(Si);
    for (j=0;j<q;j++,ns++) {
      K = VECTOR_ELT(Si, j);
      K = PROTECT(coerceVector(K,REALSXP));nprot++;
      S[ns] = REAL(K);
      if (j>0) { soff[ns] = soff[ns-1]; sr[ns] = sr[ns-1];}
    }
    // Following was used when X and G were in incompatible parameterizations
    // with S compatible with G but not X...
    //    if (q>1) for (j=0;j<q;j++,ns++) {
    //  K = VECTOR_ELT(Si, j);sr[ns] = ncols(K); /* uses fact that multi-S terms always have unpenalized last */
    //  K = PROTECT(coerceVector(K,REALSXP));nprot++;
    //  S[ns] = REAL(K);
    //  if (j>0) { soff[ns] = soff[ns-1]; /*sr[ns] = sr[ns-1];*/} 
    //} else {
    //  sr[ns] = -sr[ns];
    //  /* diagonal of Si - one or zero */
    //  S[ns] = REAL(PROTECT(coerceVector(getListEl(Sl,"ind"),REALSXP)));nprot++;
    //  ns++;
    //} /* signals identity penalty */ 
  }
  // NOTE: should throw error if ns != nsp
  ncvd(NCV,NCV1,NCV2,beta,db,G,rsd,w,&pg,&nn,a,ma,d,md,X,k,&ck,ks,m,p,&n,
       S,ns,sr,soff,sp,&nx,ts,dt,&nt,v,qc,nthreads);
  FREE(S);FREE(soff);FREE(sr);FREE(ma);FREE(md);
  UNPROTECT(nprot);return(R_NilValue);
} /* CNCV */


void tad(double *M,int p) {
/* Forms M + t(M) in situ, where M is p by p */
  int i;
  double *Md,*Ms,*M1;
  for (i=0;i<p;i++) {
    Md = M+i; Ms = M + i*p; M1 = Ms + i;
    for (; Ms<=M1; Md+= p, Ms++) { *Md += *Ms;*Ms = *Md; }
  }  
}  /* tad */

static int ucomp(const void *a, const void *b) {
  /* comparison function for qsort as used in upair. a and b are pointers to integer pointers. Each points in particular to an
     integer triplet to sort, on the basis of the first pair (the third element is an index tracker).
  */   
  int k;
  k = *(*(int **)a) - *(*(int **)b); /* difference between first elements pointed to */
  if (k==0) k = *(*(int **)a + 1) - *(*(int **)b + 1); /* difference between second elements pointed to, if needed */
  return(k);
} /* ucomp */  

void upair(int *k1,int *k2, int *ind,ptrdiff_t *n) {
/* finds the unique integer pairs in k1, k2, returning the number of these in *n (which is array length on entry).
   The pairs are returned in the early elements of k1 and k2 and ind[i] gives the elements of these unique pair arrays
   giving the ith input pair. i.e. k1[i], k2[i] on input is located in k1[ind[i]],k2[ind[i]] on output.
   - could be optimized once tested.
*/
  int **a,*ap,*a0,i,j;
  a = (int **)CALLOC(*n,sizeof(int *));
  a0 = ap = (int *)CALLOC(*n * 3,sizeof(int));
  for (i=0;i<*n;i++,ap += 3) { /* a[i] points to an integer triplet (k1[i],k2[i],i) */
    a[i] = ap;ap[0] = k1[i];
    ap[1] = k2[i];ap[2] = i;
  }  
  qsort(a,*n,sizeof(int *),ucomp); /* sort the triplet pointers to order k1,k2 */
  /* read sorted back into k1, k2 and ii */
  for (i=0;i<*n;i++,ap += 3) {
    ap = a[i];
    k1[i] = *ap;k2[i] = ap[1];ind[i] = ap[2];
  }
  for (i=0;i<*n;i++) a0[i] = ind[i]; /* need a copy for re-arranging into ind */
  /* now compress unique into k1 and k2, while recording index into ind */
  j = 0;
  ind[*a0] = j;
  for (i=1;i<*n;i++) {
    if (k1[i]!=k1[i-1]||k2[i]!=k2[i-1]) { /* k1,k2 change from previous */
      j++; k1[j] = k1[i]; k2[j] = k2[i];
    }
    ind[a0[i]] = j; /* record that original index maps to unique row j */
  }
  *n = j+1; /* unique pair count */
  FREE(a0);FREE(a);
} /* upair */  

void thread_lim(int i0,int i1,int *thlim,int nth) {
/* divide [i0,i1] into nth subintervals [thlim[i],thlim[i+1]-1], i=0..nth-1. */ 
  int di,ri,i;
  ri = (i1-i0+1); /* number of i values */
  di = ri/nth;
  while (di*nth<ri) di++; /* ensure nth*di more than covers range */
  for (i=0;i<=nth;i++) {
    thlim[i] = di*i + i0;
    if (thlim[i]>i1) thlim[i]=i1+1;
  } 
} /* thlim */  

void ncvd(double *NCV,double *NCV1,double *NCV2,double *beta,double *db, double *G,double *rsd,
	  double *w,int *pg,int *nn,int *a,ptrdiff_t *ma,int *d,ptrdiff_t *md,
	  double *X,int *k, int *ck, int *ks,int *m,int *p, ptrdiff_t *n,double **S,int ns,int *sr,
	  int *soff,double *sp,
	  int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *nthreads) {
/* computes the NCV criterion and its first two derivatives for the working linear regression of 
   a discretized smooth regression model. On entry G = (X'WX + S_t)^{-1} and rsd are the current 
   full model residuals.  nn is the number of neighbourhood blocks.
   a is the array of data rows to drop, ma is the  nn-vector of block ends (a[ma[i-1]+1:ma[i]]
   is index rows to drop for the ith neighbourhood). 
   d is the array of data rows to predict, md is the nn-vector of block ends (d[md[i-1]+1:md[i]]
   is index of predictions to be made for ith neighbourhood).
   NOTE: both a and d are assumed to have row indices sorted into ascending order within 
   neighbourhoods: see mgcv.r:onei for a suitable fast routine to do this up front. This is
   needed to facilitate matching of dropped and predicted points. 
   S[i] is ith of ns penalty matrices, dimension sr[i] by sr[i] first element penalized soff[i]
   if sr[i] < 0 then S[i] is an identity matrix -sr[i] by -sr[i], but is not actually supplied.

   The first nn entries of 'rsd' are modified on output, containing the first weighted cross 
   validated residual for each prediction neighbourhood. If prediction neighbourhoods are 
   all of size one, these are the only weighted cross validated residuals for each prediction 
   neighbourhood, of course. 

   Created below, k1 and k2 are the indexing nk-vectors required to extract the required 
   influence matrix blocks. mk gives the block-ends.

   
   NOTE: never more than maxint data, but various arrays can be larger (e.g. a and d), hence need for
         various ptrdiff_t indexing arrays and indices.

   Set up for OMP but so far tests suggest no benefit at all. schedule(dynamic) tried too, as was
   a thread per neighbourhood. Appears that it's not worth the setup overhead, basically.

   nthreads passed to diagXVXt can make a big difference, but too
   many threads harms performance.
   
   On exit rsd contains first weighted cross-validated residual for each neighbourhood.  
*/
  double *A,*A2,**A1,*Aaa,*IAaa,*Ap,*Aaap,*IA,*IAp,*ba,*bp,*rp,one=1.0,zero=0.0,*dp0,*dp1,*dp2,
    *dp3,xx,*L,**beta1,**mu1,**P,*wij,*wp,**ba1,**IAdAIA,*B,*ba1p,**L1,*IAdAIAp,***beta2,*mu2,
    *V2,*U,*ba2,*NCVth,*NCV1th,*NCV2th,*ba2p,*Arp,*Alp,*IAdAIArp,*work,***L2,***U2,*Vl,*xbdwork,*rsd1;
  int ncs=0,dum,i1,i,j,M,q,r,info=1,*k1p,*k2p,*a1p,*a2p,tid=0,maxM,maxM2,*xbiwork,
    *p1,ione=1,l,*cs,srl,*ind,*a1,*a2,nchunk,*ichunk,ch,ch_start,*piv,ncvth=1,*thlim,*l1c,**l2c,*iwork;
  ptrdiff_t nb,nu,nu_size=0,*mk,ii,jj,kk,kk0,target_nk,max_chunk,max_nb,*xbpwork,space[3];
  char uplo='U',ntrans='N',side='L';
  cs =  (int *)CALLOC(*nt,sizeof(int));
  thlim = (int *)CALLOC(ncvth+1,sizeof(int)); /* i range limits for threads */
  if (ncvth>1) {
    NCVth = (double *)CALLOC(ncvth,sizeof(double));
    NCV1th = (double *)CALLOC(ncvth*ns,sizeof(double));
    NCV2th = (double *)CALLOC(ncvth*ns*ns,sizeof(double));
  } else {
    NCVth = NCV;NCV1th = NCV1; NCV2th = NCV2;
  }

  Xbdspace(space,m,p,n,nx,dt,nt);
  xbiwork = (int *) CALLOC((size_t) space[0],sizeof(int));
  xbpwork = (ptrdiff_t *) CALLOC((size_t)space[1],sizeof(ptrdiff_t));
  xbdwork = (double *)CALLOC((size_t)space[2],sizeof(double));
  /* first create P[l], V[l], and derivatives of beta w.r.t. log smoothing parameters */ 
  P = (double **)CALLOC(ns,sizeof(double *));
  L1 = (double **)CALLOC(ns,sizeof(double *)); //!!
  beta1 = (double **)CALLOC(ns,sizeof(double *)); beta2 = (double ***)CALLOC(ns,sizeof(double **));
  L2 = (double ***)CALLOC(ns,sizeof(double **));U2 = (double ***)CALLOC(ns,sizeof(double **)); // !!
  l1c = (int *)CALLOC(ns,sizeof(int *)); /* cols in V[l] square roots */ //!!
  l2c = (int **)CALLOC(ns,sizeof(int **)); // !!
  A1 = (double **)CALLOC(ns,sizeof(double *));
  mu1 = (double **)CALLOC(ns,sizeof(double *));
  rsd1 = (double *)CALLOC((size_t) *nn,sizeof(double *));
  IAdAIA = (double **)CALLOC(ns,sizeof(double *));ba1 = (double **)CALLOC(ns,sizeof(double *));
  L = (double *)CALLOC((size_t) *pg * *pg,sizeof(double)); // !!
  Vl = (double *)CALLOC((size_t) *pg * *pg,sizeof(double)); // !!
  V2 = (double *)CALLOC((size_t) *pg * *pg,sizeof(double));
  U = (double *)CALLOC((size_t) *pg * *pg,sizeof(double)); // !!
  mu2 = (double *)CALLOC(*n,sizeof(double));
  for (dp1=L,dp2=L+ *pg * *pg,dp3=G;dp1<dp2;dp1++,dp3++) *dp1 = *dp3; /* copy G to L */
  work = (double *)CALLOC((size_t) *pg * 2,sizeof(double)); iwork = (int *)CALLOC((size_t) *pg ,sizeof(int)); // !!

  mtrf(L,L,pg,pg,&ione,&zero,work,iwork); /* LL' = G */
  
  for (l=0;l<ns;l++) { /* loop over smoothing parameters computing P[l], V[l], beta1[l], beta2[l][r] and mu1[l] */
    srl = sr[l]; if (srl<0) srl = -srl; 
    P[l] = (double *)CALLOC(srl * *pg,sizeof(double));
    beta2[l] = (double **)CALLOC(ns,sizeof(double *));
    L2[l] = (double **)CALLOC(ns,sizeof(double *));
    U2[l] = (double **)CALLOC(ns,sizeof(double *)); 
    l2c[l] = (int *)CALLOC(ns,sizeof(int *)); 
    mu1[l] = (double *)CALLOC(*n,sizeof(double));
    if (sr[l]>0) {
      F77_CALL(dgemm)(&ntrans,&ntrans,pg,sr+l,sr+l,&one,G+soff[l] * *pg,pg,S[l],sr+l,&zero,P[l],pg FCONE FCONE); /* P[l] = [G]^l S_l */
    } else { /* S[l] = I so P[l] = [G]^l */
      for (dp1=G+soff[l] * *pg,dp2=P[l],j=0;j<srl;j++) {
	xx = S[l][j];
	for (dp3=dp2 + *pg;dp2<dp3;dp2++,dp1++) *dp2 = *dp1 * xx;
      }	
    }  
    beta1[l] = db + l * *pg;
    xx = -sp[l];
    F77_CALL(dgemv)(&ntrans,pg,&srl,&xx,P[l],pg,beta+soff[l],&ione,&zero,beta1[l],&ione FCONE); /* beta1[l] = dbeta/d rho_l */
    ncs=0;xx = sp[l];
    Xbd(mu1[l],beta1[l],X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,&ione,cs,&ncs,xbiwork,xbpwork,xbdwork); /* dmu/d rho_l note limited thread safety */
  
    for (r=0;r<=l;r++) {

      if (sr[r]>sr[l]) i=sr[r]; else i=sr[l];
      for (dp1=V2,dp2=V2 + i*i;dp1<dp2;dp1++) *dp1 = 0.0;
      xx = sp[r]*sp[l];
      F77_CALL(dgemm)(&ntrans,&ntrans,sr+l,sr+r,sr+l,&xx,S[l],sr+l,P[r]+soff[l],pg,&zero,V2,&i FCONE FCONE); /* nz block of V2 is Sl [G]^r_l Sr */
      xx = 1e-14; // tol
      if (l==r) q = 1; else q = 0;
      j = -1; mtrf(V2,U,&i,&j,&q,&xx,work,iwork); /* Let V2_o be input V2. On exit V2 = L where LU'= V2_o or LL' = V2_o if l==r  */
      L2[l][r] = L2[r][l] = (double *)CALLOC(*pg * j,sizeof(double)); // !!
   
      F77_CALL(dgemm)(&ntrans,&ntrans,pg,&j,sr+l,&one,G+soff[l] * *pg,pg,V2,&i,&zero,L2[l][r],pg FCONE FCONE);
      dp0=L2[l][r];dp1=dp0+j * *pg;dp3=V2;	
      if (!q) { // i!=j
	U2[l][r] = U2[r][l] = (double *)CALLOC(*pg * j,sizeof(double));
	F77_CALL(dgemm)(&ntrans,&ntrans,pg,&j,sr+r,&one,G+soff[r] * *pg,pg,U,&i,&zero,U2[l][r],pg FCONE FCONE);
      }
      l2c[r][l]=l2c[l][r] = j;
      
      beta2[r][l] = beta2[l][r] =  (double *)CALLOC(*pg,sizeof(double));
      /* get d2beta/drho_l drho_r */
      xx = -sp[l];
      F77_CALL(dgemv)(&ntrans,pg,&srl,&xx,P[l],pg,beta1[r]+soff[l],&ione,&zero,beta2[l][r],&ione FCONE);
      xx = -sp[r];
      q = sr[r]; if (q<0) q = -q;
      F77_CALL(dgemv)(&ntrans,pg,&q,&xx,P[r],pg,beta1[l]+soff[r],&ione,&one,beta2[l][r],&ione FCONE);
      if (l==r) for (bp=beta1[r],dp1=beta2[l][r],dp2=dp1 + *pg;dp1<dp2;bp++,dp1++) *dp1 += *bp;
    }

    // alternative computation based on root of S[l]...
    xx = sp[l];
    for (dp1 = S[l],dp2= S[l] + sr[l]*sr[l],dp3=Vl;dp1<dp2;dp1++,dp3++) *dp3 = *dp1 *xx;
    j = -1; mtrf(Vl,Vl,sr+l,&j,&ione,&zero,work,iwork); l1c[l] = j; /* Vl Vl' = S[l]*sp[l] */
    L1[l] = (double *)CALLOC(*pg * j,sizeof(double)); /* L1[l]L1[l]' = sp[l] P[l] [G]_l */
    F77_CALL(dgemm)(&ntrans,&ntrans,pg,&j,sr+l,&one,G+soff[l] * *pg,pg,Vl,sr+l,&zero,L1[l],pg FCONE FCONE);
    
  } /* end of P[l], V[l] etc. loop */
 
  *NCV=0.0;
  /* create indices for extracting A_aa blocks of influence matrix */
  mk = (ptrdiff_t *)CALLOC(*nn,sizeof(ptrdiff_t)); /* ends of a1,a2 blocks defining each A_aa */
  for (ii=kk = -1,i=0;i<*nn;i++) { /* create mk */
    q = ma[i]-kk; kk = ma[i];
    ii = mk[i] = ii + q*(q+1)/2; /* mk[i] is end of block i */
  }

  /* get total storage requirements for all CV residuals + largest A_aa */
  for (M=nb=j=kk=i=0;i<*nn;i++) {
    dum = mk[i]-kk+1; kk = mk[i]+1; if (j<dum) j = dum; /* largest A_aa */
    dum = (int)(sqrt(8*dum+1)-1)/2;
    if (j<dum) j = dum; /* largest dimension of A_aa */
    nb += dum; /* total storage for adjusted residuals - but must be ma[*nn]+1 */
  }
  Aaa = (double *)CALLOC(ncvth*j*j,sizeof(double));
  IAaa = (double *)CALLOC(ncvth*j*j,sizeof(double));
  B = (double *)CALLOC(ncvth*j*j,sizeof(double));
  rp = (double *)CALLOC(j*ncvth,sizeof(double)); /* storage for residuals left out */
  wij = (double *)CALLOC(j*ncvth,sizeof(double)); /* storage for corresponding weights */
  piv = (int *)CALLOC(j*ncvth,sizeof(int)); /* pivot sequence for LU */
  ba2 = (double *)CALLOC(j*ncvth,sizeof(double)); /* storage for derivs of ba wrt rho_l */
  maxM=j;maxM2=j*j;
  /* Note that CV residuals are stored for whole drop-fold, rather than just the predict
     fold in order to avoid repeated recomputation */

  /* at this stage nk and nb are known, so it is time to decide on how to split up the
     computation into blocks to save memory */

  target_nk = *n * 2; /* target longest vector length */
  max_chunk = 0; /* largest chunk of split up nk-vectors */
  for (j=1,ii=-1,i=0;i<*nn;i++) { /* loop to split nk-vectors indexed by mk */ 
    if (mk[i]-ii>target_nk) { /* start a new chunk */
      if (mk[i-1]!=ii&&i>0) { /* aim for under target */
	if (mk[i-1]-ii>max_chunk) max_chunk = mk[i-1]-ii; /* record maximum chunk size */
	ii = mk[i-1]; j++;
      } else { /* neighbourhood so large we have to be over target */
	if (mk[i]-ii>max_chunk) max_chunk = mk[i]-ii; /* record maximum chunk size */
	ii=mk[i]; if (i < *nn-1) j++; /* only a new chunk if not at end */
      }	
    }  
  } /* j is number of chunks */
  nchunk = j;
  if (nchunk==1) max_chunk = mk[*nn-1]+1;
  ichunk = (int *)CALLOC(nchunk,sizeof(int));
  max_nb=0;
  for (j=0,jj=ii=-1,i=0;i<*nn;i++) { /* determine elements of mk indexing chunk ends */
    if (mk[i]-ii>target_nk) { /* start a new chunk */
      if (mk[i-1]!=ii&&i>0) { /* aim for under target */
	ii = mk[i-1];
	if (ma[i-1]-jj>max_nb) max_nb = ma[i-1]-jj;
	jj = ma[i-1];
	ichunk[j] = i-1; j++;
      } else { /* chunk so large we have to be over target */
	ii = mk[i];ichunk[j]=i; j++;
	if (ma[i]-jj>max_nb) max_nb = ma[i]-jj;
	jj = ma[i];
      }	
    }  
  }
  ichunk[j] = *nn-1;
  if (nchunk==1) max_nb = ma[*nn-1]+1;
  /* now mk[ichunk[j]] is end of jth chunk and max_nb is largest storage needed for
     CV residuals by a chunk. max_chunk is the largest vector needed to index
     the required elements A_ij for a chunk. ichunk[j] indexes the last neighbourhood in
     chunk j.
  */
  
  /* allocate required memory not dependent on knowing nu */
  ba = (double *)CALLOC(max_nb,sizeof(double)); /* storage adjusted residuals for each neighbourhood */
  
  for (l=0;l<ns;l++) {
    ba1[l] = (double *)CALLOC(max_nb,sizeof(double)); /* storage for derivs of ba wrt rho_l */
  } /* l loop - allocation */
 
  /* obtain the unique i,j elements of A required to produce the A_aa submatrix for each
     neighbourhood along with an index vector for extracting each neighbourhood's lower
     triangle from the unique elements */

  IA = (double *)CALLOC(max_chunk,sizeof(double)); /* storage for (I-A_aa)^-1 - same format as A without unique removal */
  for (l=0;l<ns;l++) IAdAIA[l]  = (double *)CALLOC(max_chunk,sizeof(double)); /* storage for (I-A_aa)^{-1} dA_aa/rho_l (I-A_aa)^&{-1} */
  a1 = (int *)CALLOC(max_chunk,sizeof(int));a2 = (int *)CALLOC(max_chunk,sizeof(int));
  ind = (int *)CALLOC(max_chunk,sizeof(int));

  for (ch=0;ch<nchunk;ch++) { /* loop over chunks */
    if (ch==0) {
      ch_start = kk = 0;
      nu = mk[ichunk[ch]]+1;
    } else {
      ch_start=ichunk[ch-1]+1; /* neighbourhood starting this chunk */ 
      kk = ma[ch_start-1]+1; /* start index of first drop set */
      nu = mk[ichunk[ch]] - mk[ichunk[ch-1]] ; /* number of A_ij elements needed */
    }
   
    for (k1p=a1,k2p=a2,i=ch_start;i<=ichunk[ch];i++) { /* loop over neighbourhoods */
      M = ma[i]-kk+1; /* size of neighbourhood i */
      a1p = a + kk; 
      for (j=0;j<M;j++,a1p++) for (a2p = a+kk,q=0;q<=j;q++,k1p++,k2p++,a2p++) {
	*k1p = *a1p; /* a[kk+j] */
	*k2p = *a2p; /* a[kk+q] */
      }
      kk = ma[i] + 1; 
    }
    
    upair(a1,a2,ind,&nu);
    /* now first nu elements of a1 and a2 are unique pairs, while a1[ind[i]], a2[ind[i]] is
       original a1[i],a2[i]. This enables everything to be based only on unique elements
       of A etc, rather than repeating element extractions many times for each neighbourhood
       so mk now contains block ends for ind elements defining each A_aa
    */
   
    if (ch&&nu>nu_size) { /* need to extend storage */
      nu_size = nu;
      for (l=0;l<ns;l++) A1[l] = (double *)REALLOC(A1[l],nu*sizeof(double)); /* storage for dA_aa/d rho_l packed as A */ 
      A = (double *)REALLOC(A,nu*sizeof(double));A2 = (double *)REALLOC(A2,nu*sizeof(double));
    } else { /* set up initial storage */
      nu_size = nu;
      for (l=0;l<ns;l++) A1[l] = (double *)CALLOC(nu,sizeof(double)); /* storage for dA_aa/d rho_l packed as A */ 
      A = (double *)CALLOC(nu,sizeof(double));A2 = (double *)CALLOC(nu,sizeof(double));
    }

    /* extract the unique elements making up the required influence matrix blocks */
    j=1;
    diagXLLtXt(A,L,X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,pg,pg,
	       a1,a2,&nu,nthreads);
    /* Get NCV and (I-A_aa)^-1... */
    if (ch==0) kk0=0; else kk0 = mk[ichunk[ch-1]]+1;
    thread_lim(ch_start,ichunk[ch],thlim,ncvth);
#ifdef OPENMP_ON
#pragma omp parallel for private(i,j,q,tid,ii,jj,kk,bp,IAp,Aaap,i1,dum,wp,dp1,dp2,xx,p1,info,M) num_threads(ncvth) 
#endif
    for (tid=0;tid<ncvth;tid++)
    for (i=thlim[tid];i<thlim[tid+1];i++) { /* loop over neighbourhoods */
      // dpotri, DPOTRF and DPSTRF (piv)
      if (i==0) bp=ba; else {
        if (ch_start==0) bp = ba + ma[i-1]+1; else bp = ba + ma[i-1]- ma[ch_start-1];
      }	
      if (i==0) {
	ii=jj=kk=0;
	IAp = IA;
      } else {
        kk = mk[i-1]+1; /* start index of neighbourhood in A */
        ii = ma[i-1]+1; /* start index of drop neighbourhood in a */
        jj = md[i-1]+1; /* start index of predict neighbourhood in d */
	IAp = IA + kk - kk0; 
      }	
      M = (int) ma[i]-ii+1;//(int)(sqrt(8*(mk[i]-kk+1)+1)-1)/2; /* current Aaa is M by M */
   
      /* Fill upper triangle of I-Aaa */
      for (p1 = ind+kk-kk0,Aaap = Aaa+maxM2*tid,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,Aaap++,p1++) *Aaap = -A[*p1]; /* fill upper triangle of current Aaa*/
      fill_lt(Aaa+maxM2*tid,M); /* fill in lower triangle */
      for (wp=wij+maxM*tid,j=0;j<M;j++) wp[j] = w[a[ii+j]]; /* w_aa */
      for (Aaap=Aaa+maxM2*tid,wp=wij+maxM*tid,dp1=wp+M;wp<dp1;wp++) for (dp2=Aaap+M;Aaap<dp2;Aaap++) *Aaap *= *wp; /*A_aa W_a */
      for (Aaap=Aaa+maxM2*tid,j=0;j<M;j++,Aaap += M+1) *Aaap += 1; /* Aaa = I-A_aa W_a */
      if (M==1) { /* unit neighbourhood - trivial computation */
        *IAp = *(wij+maxM*tid) / *(Aaa+maxM2*tid);  /*W_a/(1-A_aa)*/
        *bp = *IAp * rsd[a[ii]]; /* weighted adjusted residual */
        IAp++;
      } else { /* bigger neighbourhood - need matrix inversion */
        /* now invert I-Aaa */
        F77_CALL(dgetrf)(&M,&M,Aaa+maxM2*tid,&M,piv+tid*maxM,&info); /* LAPACK LU*/
        F77_CALL(dgetri)(&M,Aaa+maxM2*tid,&M,piv+tid*maxM,rp+tid*maxM,&M,&info); /* LAPACK LU to inv */
	for (Aaap=Aaa+maxM2*tid,j=0;j<M;j++) for (wp=wij+tid*maxM,dp1=wp+M;wp<dp1;Aaap++,wp++) *Aaap *= *wp; /*W_a A_aa */
        /* compute adjusted residuals ba = (I-Aaa)^-1 rsd_aa */
        for (dp1=rp+tid*maxM,j=0;j<M;j++,dp1++) *dp1 = rsd[a[ii+j]]; /* rsd_aa */
        F77_CALL(dsymv)(&uplo,&M,&one,Aaa+tid*maxM2,&M,rp+tid*maxM,&ione,&zero,bp,&ione FCONE); /* weighted adjusted residuals */
        /* repack W_a A_aa into IA as A is packed, for use in derivative computation */
        for (Aaap = Aaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,IAp++,Aaap++) *IAp = *Aaap; 
       
      }
      /* Now compute NCV itself  */
      q = md[i]-jj+1; /* size of predict neighbourhood */
      for (dum=j=0;j<q;j++) { /* loop over predict indices */
        i1 = d[jj+j]; /* need to match predict index to a dropped index */
        while (i1 > a[ii+dum] && dum<M) dum++; // NOTE: need to throw an error if dum==M
	xx = bp[dum];
        if (xx!=0.0) NCVth[tid] += xx*xx/w[i1];
	if (j==0) rsd1[i] = xx; /* store first (and only in cases where useful) weighted cross validated residual */
      }
   
    } /* neighbourhood, i loop */
  
    /* Compute the derivatives... */

    for (l=0;l<ns;l++) { /* loop over smoothing parameters */
      srl = sr[l]; if (srl<0) srl = -srl; 
  
      j=1;diagXLLtXt(A1[l],L1[l],X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,pg,l1c+l,
		     a1,a2,&nu,nthreads);
      for (dp1=A1[l],dp2=dp1+nu;dp1<dp2;dp1++) *dp1 = - *dp1;
      if (ch==0) kk0=0; else kk0= mk[ichunk[ch-1]]+1;
      
#ifdef OPENMP_ON
#pragma omp parallel for private(i,j,q,tid,ii,jj,kk,bp,ba1p,IAp,Ap,Aaap,IAdAIAp,i1,dum,dp1,dp2,xx,p1,M,side) num_threads(ncvth)
#endif
       for (tid=0;tid<ncvth;tid++)
       for (i=thlim[tid];i<thlim[tid+1];i++) { /* loop over neighbourhoods */
	if (i==0) {
	  bp=ba;ba1p=ba1[l];
	} else {
          if (ch_start==0) {
	    bp = ba + ma[i-1]+1;
	    ba1p = ba1[l] + ma[i-1]+1;
	  } else {
	    bp = ba + ma[i-1]- ma[ch_start-1];
	    ba1p = ba1[l] + ma[i-1] - ma[ch_start-1];
	  }  
        }	
        if (i==0) {
	  ii=jj=kk=0;
	  IAp = IA;IAdAIAp = IAdAIA[l];
        } else {
          kk = mk[i-1]+1; /* start index of neighbourhood in A */
          ii = ma[i-1]+1; /* start index of drop neighbourhood in a */
          jj = md[i-1]+1; /* start index of predict neighbourhood in d */
	  IAp = IA + kk - kk0;
	  IAdAIAp = IAdAIA[l] + kk - kk0;
        }	

	M = (int) ma[i]-ii+1; //(int)(sqrt(8*(mk[i]-kk+1)+1)-1)/2; /* current Aaa is M by M */
  
        for (p1 = ind+kk-kk0,Ap=A1[l],Aaap = Aaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,Aaap++,p1++) *Aaap = Ap[*p1]; /* current A1_aa to upper tri Aaa */
        fill_lt(Aaa+tid*maxM2,M); /* fill in lower triangle */
       
	for (Aaap = IAaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,IAp++,Aaap++) *Aaap = *IAp; /* current W_a(I-A_aa W_a)^{-1} to upper tri IAaa */
        side = 'L';
        F77_CALL(dsymm)(&side,&uplo,&M,&M,&one,IAaa+tid*maxM2,&M,Aaa+tid*maxM2,&M,&zero,B+tid*maxM2,&M FCONE FCONE); /* B = W_a(I-A_aa W_a)^{-1} dA */
        side = 'R';
        F77_CALL(dsymm)(&side,&uplo,&M,&M,&one,IAaa+tid*maxM2,&M,B+tid*maxM2,&M,&zero,Aaa+tid*maxM2,&M FCONE FCONE); /* Aaa =  W_a(I-A_aa W_a)^{-1} dA W_a(I-A_aa W_a)^{-1} */
        for (dp2=rp+tid*maxM,dp1=mu1[l],j=0;j<M;j++,dp2++) *dp2 = -dp1[a[ii+j]]; /* -d mu_a / d rho_l */
        F77_CALL(dsymv)(&uplo,&M,&one,IAaa+tid*maxM2,&M,rp+tid*maxM,&ione,&zero,ba1p,&ione FCONE); /* ba1[l] = -W_a(I-A_aa W_a)^{-1}dmu_a/drho_l */ 
        F77_CALL(dgemv)(&ntrans,&M,&M,&one,B+tid*maxM2,&M,bp,&ione,&one,ba1p,&ione FCONE); /* ba1p = db^a/drho_l */
        /* store W_a (I-A_aa W_a)^{-1} dA W_a (I-A_aa W_a)^{-1} upper tri for future use */
        /* repack into AI as A is packed, for use in derivative computation */
        for (Aaap = Aaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,IAdAIAp++,Aaap++) {
          *IAdAIAp = *Aaap; /* store - no need to clear only upper tri written actually accessed */
        }
        /* compute the NCV derivative w.r.t. rho_l... */
        q = md[i]-jj+1; /* size of predict neighbourhood */
        for (dum=j=0;j<q;j++) { /* loop over predict indices */
          i1 = d[jj+j]; /* need to match predict index to a dropped index */
          while (i1 > a[ii+dum] && dum<M) dum++; // NOTE: need to throw an error if dum==M
	  xx = ba1p[dum] * bp[dum];
          if (xx!=0.0) NCV1th[tid*ns+l] += 2*xx/w[i1]; /* deriv of NCV wrt rho_l */
        }
       
      } /* neighbourhood i loop */
    
      /* The second derivative loop... */
      for (r=0;r<=l;r++) {   
        ncs=0;
        Xbd(mu2,beta2[l][r],X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,&ione,cs,&ncs,xbiwork,xbpwork,xbdwork); /* d2mu/drho_l drho_r note limited thread safety */
	
      	if (l==r) {
	  diagXLLtXt(A2,L2[l][r],X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,
			  pg,l2c[l]+r,a1,a2,&nu,nthreads);
	  for (dp1=A1[l],dp2=dp1+nu,dp3=A2;dp1<dp2;dp1++,dp3++) *dp3 += *dp3 + *dp1;
        } else diagXLUtXt(A2,L2[l][r],U2[l][r],X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,
			  pg,l2c[l]+r,a1,a2,&nu,nthreads);

	/* Compute the second derivative of NCV w.r.t. rho_l and rho_r */
       
        if (ch==0) kk0=0; else kk0= mk[ichunk[ch-1]]+1; 
#ifdef OPENMP_ON
#pragma omp parallel for private(i,j,q,tid,ii,jj,kk,bp,ba1p,ba2p,IAp,Alp,Arp,Aaap,IAdAIAp,IAdAIArp,i1,dum,dp1,dp2,xx,p1,M) num_threads(ncvth) 
#endif
	for (tid=0;tid<ncvth;tid++)
        for (i=thlim[tid];i<thlim[tid+1];i++) { /* loop over neighbourhoods */
	  if (i==0) {
	    bp=ba;ba1p=ba1[l];ba2p=ba1[r];
	  } else {
            if (ch_start==0) {
	      bp = ba + ma[i-1]+1;
	      ba1p = ba1[l] + ma[i-1]+1;
	      ba2p = ba1[r] + ma[i-1]+1;
	    } else {
	      bp = ba + ma[i-1]- ma[ch_start-1];
	      ba1p = ba1[l] + ma[i-1] - ma[ch_start-1];
	      ba2p = ba1[r] + ma[i-1] - ma[ch_start-1];
	    }  
          }	
          if (i==0) {
	    ii=jj=kk=0;
	    IAp = IA;IAdAIAp = IAdAIA[l];IAdAIArp = IAdAIA[r];
          } else {
            kk = mk[i-1]+1; /* start index of neighbourhood in A */
            ii = ma[i-1]+1; /* start index of drop neighbourhood in a */
            jj = md[i-1]+1; /* start index of predict neighbourhood in d */
	    IAp = IA + kk - kk0;
	    IAdAIAp = IAdAIA[l] + kk - kk0;
	    IAdAIArp = IAdAIA[r] + kk - kk0;
          }	

	  M = (int) ma[i]-ii+1;// (sqrt(8*(mk[i]-kk+1)+1)-1)/2; /* current Aaa is M by M */
	  /*  need second derivative of b_a wrt rho_l and rho_r... */
          for (Aaap = IAaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,IAp++,Aaap++) *Aaap =  *IAp; /* current W_a(I-A_aa W_a)^{-1} to upper tri IAaa */
          for (p1 = ind+kk-kk0,Aaap = Aaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,Aaap++,p1++) *Aaap = A2[*p1]; /* Aaa =  d2A_aa/drho_l drho_r */

	  F77_CALL(dsymv)(&uplo,&M,&one,Aaa+tid*maxM2,&M,bp,&ione,&zero,B+tid*maxM2,&ione FCONE); /* B =  d2 A_aa/drho_l drho_r b^a */
	  F77_CALL(dsymv)(&uplo,&M,&one,IAaa+tid*maxM2,&M,B+tid*maxM2,&ione,&zero,ba2+tid*maxM,&ione FCONE); /*  W_a(I-A_aa W_a)^{-1} d2 A_aa/drho_l drho_r b^a */
	  for (dp2=rp+tid*maxM,dp1=mu2,j=0;j<M;j++,dp2++) *dp2 = -dp1[a[ii+j]]; /* rp = -d2mu/d rho_l d rho_r */
	  F77_CALL(dsymv)(&uplo,&M,&one,IAaa+tid*maxM2,&M,rp+tid*maxM,&ione,&one,ba2+tid*maxM,&ione FCONE); /* ba2 has 2nd deriv based terms - now add first deriv based */

	  for (p1 = ind+kk-kk0,Alp=A1[l],Aaap = Aaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,Aaap++,p1++) *Aaap = Alp[*p1]; /* Aaa = dA_aa/drho_l */
          F77_CALL(dsymv)(&uplo,&M,&one,Aaa+tid*maxM2,&M,bp,&ione,&zero,B+tid*maxM2,&ione FCONE); /* B = dA_aa/drho_l b^a */

	  for (Aaap = Aaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,IAdAIArp++,Aaap++) *Aaap =  *IAdAIArp; 
          F77_CALL(dsymv)(&uplo,&M,&one,Aaa+tid*maxM2,&M,B+tid*maxM2,&ione,&one,ba2+tid*maxM,&ione FCONE);
	  /* ...added W_a(I-A_aa W_a)^{-1}dA_aa/drho_r W_a(I-A_aa W_a)^{-1} dA_aa/drho_l b^a */

	  for (dp2=rp+tid*maxM,dp1=mu1[l],j=0;j<M;j++,dp2++) *dp2 = -dp1[a[ii+j]]; /* rp = -dmu/d rho_l */
	  F77_CALL(dsymv)(&uplo,&M,&one,Aaa+tid*maxM2,&M,rp+tid*maxM,&ione,&one,ba2+tid*maxM,&ione FCONE);  /* -W_a(I-A_aa W_a)^{-1}dA_aa/drho_r W_a(I-A_aa W_a)^{-1} d mu_a / drho_l */

          /* now r <-> l for last two terms added to ba2 */
	  for (p1 = ind+kk-kk0,Arp=A1[r],Aaap = Aaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,Aaap++,p1++) *Aaap = Arp[*p1]; /* Aaa = dA_aa/drho_r */
	  F77_CALL(dsymv)(&uplo,&M,&one,Aaa+tid*maxM2,&M,bp,&ione,&zero,B+tid*maxM2,&ione FCONE); /* B = dA_aa/drho_r b^a */

	  for (Aaap = Aaa+tid*maxM2,j=0;j<M;j++,Aaap += M-j) for (q=0;q<=j;q++,IAdAIAp++,Aaap++) *Aaap =  *IAdAIAp; 
          F77_CALL(dsymv)(&uplo,&M,&one,Aaa+tid*maxM2,&M,B+tid*maxM2,&ione,&one,ba2+tid*maxM,&ione FCONE);
	  /* ...added W_a(I-A_aa W_a)^{-1}dA_aa/drho_l W_a(I-A_aa W_a)^{-1} dA_aa/drho_r b^a */

	  for (dp2=rp+tid*maxM,dp1=mu1[r],j=0;j<M;j++,dp2++) *dp2 = -dp1[a[ii+j]]; /* rp = -dmu/d rho_r */
	  F77_CALL(dsymv)(&uplo,&M,&one,Aaa+tid*maxM2,&M,rp+tid*maxM,&ione,&one,ba2+tid*maxM,&ione FCONE);  /* -W_a(I-A_aa W_a)^{-1}dA_aa/drho_l W_a(I-A_aa W_a)^{-1} d mu_a / drho_r */

	  q = md[i]-jj+1; /* size of predict neighbourhood */
          for (dum=j=0;j<q;j++) { /* loop over predict indices */
            i1 = d[jj+j]; /* need to match predict index to a dropped index */
            while (i1 > a[ii+dum]&&dum<M) dum++; // NOTE: need to throw an error if dum==M
	    xx = ba1p[dum] * ba2p[dum] + bp[dum]*ba2[dum+tid*maxM];
            if (xx!=0.0) NCV2th[tid*ns*ns+l+r*ns] += 2 * xx/w[i1]; /* deriv of NCV wrt rho_l */
          }
	
        } /* neighbourhood i loop */
      } /* r sp loop */   
    } /* l sp loop */
  } /* chunk loop */
  if (ncvth>1) {
    for (i=0;i<ncvth;i++) {
      *NCV += NCVth[i];
      for (j=0;j<ns;j++) NCV1[j] += NCV1th[i*ns+j];
      for (j=0;j<ns;j++) for (q=0;q<=j;q++) NCV2[j+q*ns] += NCV2th[i*ns*ns+j+q*ns];
    }  
  }
  for (j=0;j<ns;j++) for (q=0;q<=j;q++) NCV2[q+j*ns]=NCV2[j+q*ns];
  for (j=0;j<*nn;j++) rsd[j] = rsd1[j]; /* return first weighted cv residual for each nei */
  /* free memory */
  for (l=0;l<ns;l++) {
    for (r=0;r<=l;r++) {
      FREE(beta2[l][r]);FREE(L2[l][r]);if (r<l) FREE(U2[l][r]);
    }  
    FREE(beta2[l]);FREE(L2[l]);FREE(U2[l]);
    FREE(P[l]);FREE(L1[l]);FREE(A1[l]);FREE(ba1[l]);
    FREE(mu1[l]);FREE(IAdAIA[l]);
    FREE(l2c[l]);
  }
  FREE(l1c);FREE(l2c);FREE(L);FREE(V2);FREE(U);FREE(Vl);FREE(P);FREE(beta1);FREE(L1);FREE(A1);
  FREE(ba1);FREE(mu1);FREE(rsd1);FREE(IAdAIA);FREE(piv);FREE(work);FREE(iwork);
  FREE(B);FREE(cs);FREE(mk);FREE(A);FREE(Aaa);FREE(IAaa);FREE(IA);FREE(ba);FREE(rp);
  FREE(beta2);FREE(L2);FREE(U2);FREE(mu2);FREE(wij);FREE(A2);FREE(ba2);FREE(ind);FREE(a1);FREE(a2);FREE(thlim);
  if (ncvth>1) {
    FREE(NCV1th);FREE(NCV2th);
    FREE(NCVth);
  }
  FREE(xbiwork);FREE(xbpwork);FREE(xbdwork);
} /* ncvd */





/* PENDING discrete NCV implementation notes:
   - need to adjust nei for any dropped data in bam!! 
*/   

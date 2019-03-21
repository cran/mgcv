/* Convenient C wrappers for calling LAPACK and LINPACK + other matrix routines using same 
   packing format....

   See R-x/include/R_ext/Lapack.h...

   parallel support (using openMP) offered by...
   * mgcv_pqr - parallel QR based on row blocking ok for few cores and n>>p but
                somewhat involved. Advantage is that just one thread used per 
                core, so threading overhead minimal. Disadvantage is that an extra 
                single thread QR is used to combine everything.
   * mgcv_piqr - pivoted QR that simply parallelizes the 'householder-to-unfinished-cols'
                 step. Storage exactly as standard LAPACK pivoted QR.
   * mgcv_Rpiqr - wrapper for above for use via .call 
   * mgcv_pmmult - parallel matrix multiplication.
   * mgcv_Rpbsi - parallel inversion of upper triangular matrix.
   * Rlanczos - parallel on leading order cost step (but note that standard BLAS seems to 
                use Strassen for square matrices.)   

*/

/* dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *A,
         int *lda, double *B, int *ldb, double *beta,double *C,int *ldc) 
   transa/b = 'T' or 'N' for A/B transposed or not. C = alpha op(A) op(B) + beta C,
   where op() is transpose or not. C is m by n. k is cols of op(A). ldx is rows of X
   in calling routine (to allow use of sub-matrices) */

/* dsyrk(char *uplo, char *trans,int *n, int *k,double *a, double *A, int *lda,
         double *b, double *C,int *ldc)
   uplo = 'U' or 'L' for upper or lower tri of C used. trans = 'N' for C = aAA' + bC  
   and A n by k, or 'T' for C =  aA'A + bC and A k by n; C is n by n. lda and ldc
   are actual number of rows in A and C respectively (allows use on submatrices). */

/* dswap(int *n,double *x,int *dx,double *y,int *dy)
   Swaps n elements of vectors x and y. Spacing between elements is dx and dy. */

/* dgemv(char *trans,int *m, int *n,double a,double *A,int *lda,double *x,int *dx,
         double *b, double *y,int *dy)
   trans='T' to transpose A, 'N' not to. A is m by n. Forms y = a*A'x + b*y, or
   y = a*Ax + b*y. lda is number of actual rows in A (to allow for sub-matrices)
   dx and dy are increments of x and y indices.  */


/* dtrmm(char *side,char *uplo,char *transa, char *diag,int *m, int *n, double *alpha,
         double *A,int *lda,double *B,int *ldb)
   B = alpha*op(A)*B or B = alpha*B*op(A) for side = 'L' or 'R'. op is identity of transpose
   depending on whether transa = 'N' or 'T'. B is m by n, upper/lower triangular if uplo = 'U'
   or 'L'. B is unit diagonal if diag = 'U' rather than 'N'. lda and ldb are physical rows 
   of A and B. */

#include "mgcv.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
//#include <R_ext/Linpack.h> /* only needed for pivoted chol - see note in mgcv_chol */
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifdef OPENMP_ON
#include <omp.h>
#endif

/*#include <dmalloc.h>*/

void mgcv_omp(int *a) {
#ifdef OPENMP_ON
  *a=1;
#else
  *a=0;
#endif  
}

void rpmat(double *A,int n) {

  int i,j;
  for (i=0;i<n;i++) {
    Rprintf("\n");
    for (j=0;j<n;j++) Rprintf("%7.2g  ",A[i + n * j]);
  }
  Rprintf("\n");
} /* rpmat */ 

void dump_mat(double *M,int *r,int*c,const char *path) {
  /* dump r by c matrix M to path - intended for debugging use only */
  FILE *mf;
  mf = fopen(path,"wb");
  if (mf == NULL) { 
    Rprintf("\nFailed to open file\n");
    return;
  }
  fwrite(r,sizeof(int),1,mf); fwrite(c,sizeof(int),1,mf);
  fwrite(M,sizeof(double),*r * *c,mf);
  fclose(mf);
}

void read_mat(double *M,int *r,int*c,char *path) {
/* routine to facilitate reading dumped matrices back into R - debugging use only
   
   e.g. (actually path doesn't work here)
   oo <- .C("read_mat",as.double(0),r=as.integer(0),c=as.integer(0),
             as.character("/home/sw283/tmp/badmat.dat"),PACKAGE="mgcv")
   oo <- .C("read_mat",M=as.double(rep(0,oo$c*oo$r)),r=as.integer(oo$r),c=as.integer(oo$c),
             as.character("/home/sw283/tmp/badmat.dat"),PACKAGE="mgcv")
   M <- matrix(oo$M,oo$r,oo$c)
*/
 size_t j;
 FILE *mf;
 mf = fopen("/home/sw283/tmp/badmat.dat","rb"); 
 if (mf == NULL) { 
    Rprintf("\nFailed to open file\n");
    return;
 }
 if (*r < 1) { /* dimension query */
   j=fread(r,sizeof(int),1,mf); j=fread(c,sizeof(int),1,mf);
 } else {
   j=fread(r,sizeof(int),1,mf); j=fread(c,sizeof(int),1,mf);
   j=fread(M,sizeof(double),*r * *c,mf);
   if (j!= *r * *c)  Rprintf("\nfile dim problem\n");
 }
 fclose(mf);
}

void up2lo(double * A, int n) {
/* copies upper triangle of n by n matrix A to lower triangle */
  double *pu,*pl,*plf;
  int i;
  for (i=0;i<n;i++)
    for (plf = A + i * (ptrdiff_t) n, pu=pl=plf + i,plf=plf+n,pl++,pu+=n;pl<plf;pl++,pu+=n) *pl = *pu;

} /* up2lo */  

void row_squash(double *X,int rnew,int rold,int col) {
/* rnew<rold. Copies rold by col matrix in X to rnew by col matrix in X on output.
   Does not clear beyond new matix to zero.
*/   
  int c;
  double *Xnew,*Xold,*X1;
  Xnew = X;
  for (c=0;c<col;c++)
    for (Xold = X + c * (ptrdiff_t) rold,X1=Xold+rnew;Xold<X1;Xold++,Xnew++) *Xnew = *Xold;
} /* drop_rows */  

void tile_ut(int n,int *m,int *K, int *C,int *R,int *B) {
/* Split n by n upper triangular matrix into blocks distributed on m threads.     
   Routine may reduce m if request pointless otherwise. Tiling done on assumption 
   that work on leading diagonal blocks is half the cost of work on off-diagonals.
   There are N = m * (m+1)/2 blocks. The ith block is block row R[i] and block col
   C[i]. Block i is rows K[R[i]]:K[R[i]+1]-1 and cols K[C[i]]:K[C[i]+1].
   The jth thread deals with blocks B[j]:B[j+1]-1. K and B are dimension m + 1;
   C and R are dimension N. n and m are input params. m and the rest are output. 

   The basic idea is that if there are nt threads then we divide the computation into m blocks of rows 
   and m blocks of columns (of approximately equal dimension). Leading diagonal blocks have about half the 
   cost of off diagonal blocks, by symmetry, so when dividing up the blocks between threads this must be taken 
   into account. Generally each core gets m/2 full blocks. If m is odd, then each core gets one leading
   diagonal block and (m-1)/2 off-diagonal blocks. If m is even then the first m/2 cores each get 2 leading 
   diagonal blocks plus (m-2)/2 off diagonal blocks. The remaining m/2 cores each get m/2 off diagonal blocks.
*/
  double x,dn=0.0;
  int i,N,kk,l,nb,r,c;
  (*m) ++;
  while (dn<1 && *m>1) { (*m)--;dn = n/(double) *m; } /* m is number of threads to use */
  N = (*m * (*m + 1))/2; /* total number of blocks to process */
  /* get K such that i,jth block has rows K[i] to K[i+1]-1 and cols K[j] to K[j+1]-1 */ 
  K[0]=0;
  for (x=0.0,i=1;i < *m;i++) { x+=dn; K[i] = (int) floor(x);}
  K[*m] = n;
  if (*m % 2) { /* odd number of threads */
    kk=l=nb=B[0]=C[0]=R[0]=0;i=1;
    for (r=0;r < *m;r++) for (c=r+1;c < *m;c++,i++) {
      if (kk==(*m - 1)/2) { /* each new block gets one leading diagonal block */
	  kk=0;l++;
	  C[i]=R[i]=l;
	  nb++;B[nb]=i;i++;
      }
      C[i] = c;R[i]=r;kk++;
    }	
  } else { /* even number of threads */
    nb=l=kk=B[0]=i=0;
    for (r=0;r < *m;r++) for (c=r+1;c < *m;c++) {
	if (kk == *m/2) {kk=0;nb++;B[nb]=i;}
	if (kk==0 && l < *m) { /* allocate 2 leading diagonal blocks */
          R[i]=C[i]=l;l++;i++;R[i]=C[i]=l;l++;i++;kk++;
	  if (kk == *m/2) { kk = 0;nb++;B[nb]=i;}
	}
        R[i]=r;C[i]=c;kk++;i++;
    }	
  }
  B[*m] = N;
} /* tile_ut */  

/* dtrmm(char *side,char *uplo,char *transa, char *diag,int *m, int *n, double *alpha,
         double *A,int *lda,double *B,int *ldb)
   B = alpha*op(A)*B or B = alpha*B*op(A) for side = 'L' or 'R'. op is identity of transpose
   depending on whether transa = 'N' or 'T'. B is m by n, upper/lower triangular if uplo = 'U'
   or 'L'. B is unit diagonal if diag = 'U' rather than 'N'. lda and ldb are physical rows 
   of A and B. */

/* dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *A,
         int *lda, double *B, int *ldb, double *beta,double *C,int *ldc) 
   transa/b = 'T' or 'N' for A/B transposed or not. C = alpha op(A) op(B) + beta C,
   where op() is transpose or not. C is m by n. k is cols of op(A). ldx is rows of X
   in calling routine (to allow use of sub-matrices) */

void pdtrmm(int *n,int *q,double *alpha, double *A,int *lda,double *D,int *ldd,int *nt,int *iwork,double *work) {
/* iwork is dim 3 *nt * (nt + 1)/2 + 2 * nt + 2. work is dim q*(n+nt)*(nt+1)/2. D is n by q. A is n by n upper triangular.
   D = alpha * A * D.  

*/
  int i,j,m,N,*K,*C,*R,*B,*off,r,c,nr,nc,ldt;
  double *p0,*p1,*p2,*p3,*p4,*p5,zero=0.0;
  char side ='L',nope='N',up='U';
  m = *nt;N = (m * (m + 1))/2; /* total number of blocks to process (upper bound)*/
  K = iwork;iwork += m+1;C = iwork;iwork += N;
  R = iwork; iwork += N;B = iwork;iwork += m+1;off = iwork;
  tile_ut(*n,&m,K,C,R,B); /* set up tiling and allocation of tiles to threads */
  N = (m * (m + 1))/2; /* total number of blocks to process (actual) */
  off[0] = 0;
  for (i=1;i<N;i++) off[i] = off[i-1] + K[R[i-1]+1] - K[R[i-1]]; /* start row of each temporary block (in work) */
 
  ldt = off[N-1] + K[R[N-1]+1] - K[R[N-1]]; /* rows of work */

  #ifdef OPENMP_ON
  #pragma omp parallel for private(i,j,r,c,nr,p0,p1,p2,p3,p4,p5,nc) num_threads(m)
  #endif
  for (i=0;i<m;i++) {
    for (j=B[i];j<B[i+1];j++) {
      r = R[j]; c = C[j]; nr = K[r+1] - K[r];
      if (r==c) { /* diagonal block */
	/* copy D[K[c]:K[c+1]-1,0:k-1] to work[off[j]:off[j]+nr-1,0:k-1] */
	for (p0=D+K[c],p1=work+off[j],p2 = D + (ptrdiff_t) *ldd * *q;p0<p2;p0 += *ldd,p1 += ldt)
	  for (p3=p0,p4=p3+nr,p5=p1;p3<p4;p3++,p5++) *p5 = *p3;
	
        F77_CALL(dtrmm)(&side,&up,&nope,&nope,&nr,q,alpha,A + K[r] + K[c] * (ptrdiff_t) *lda,lda,work+off[j],&ldt);
	
      } else { /* a general block */
  	nc = K[c+1] - K[c];
	F77_CALL(dgemm)(&nope,&nope,&nr,q,&nc,alpha,A + K[r] + K[c] * (ptrdiff_t) *lda,lda,D + K[c],ldd,&zero,work+off[j],&ldt);
      }    
    }
  }
  /* now sew it all back together */
  for (p0 = D,p2 = D + (ptrdiff_t) *ldd * *q;p0<p2;p0 += *ldd) for (p1=p0,p3=p0 + *n;p1<p3;p1++) *p1 = 0.0;
  for (i=0;i<N;i++) {
    r = R[i];nr = K[r+1]-K[r]; /* add work[off[i]:off[i]+nr-1,0:q-1] to  D[K[r]:K[r+1]-1,0:q-1] */
    for (p0=work + off[i],p1 = D+K[r];p1<p2;p0 += ldt,p1 += *ldd)
      for (p3=p0,p4=p1,p5=p4 + nr;p4<p5;p4++,p3++) *p4 += *p3;
  }  
  
} /* pdtrmm */  


void pdsyrk(int *n,int *k,double *a,double *A,int *lda,double *b,double *D, int *ldd,int *work,int *nt) {
/* Forms D = aA'A + bD in upper triangle of D using nt threads. A is n by n.
   A is k by n. lda and ldd are physical dimensions. work is length nt * (nt + 3) + 2. 

   The basic idea is that if there are nt threads then we divide the computation into nt blocks of rows 
   and nt blocks of columns (of approximately equal dimension). Leading diagonal blocks have about half the 
   cost of off diagonal blocks, by symmetry, so when dividing up the blocks between threads this must be taken 
   into account. Generally each core gets nt/2 full blocks. If nt is odd, then each core gets one leading
   diagonal block and (nt-1)/2 off-diagonal blocks. If nt is even then the first nt/2 cores each get 2 leading 
   diagonal blocks plus (nt-2)/2 off diagonal blocks. The remaining nt/2 cores each get nt/2 off diagonal blocks.

   Leading diagonal blocks are computed using dsyrk and off diagonals using dgemm. For a tuned BLAS the speed up
   from multi-threading is modest. For the reference BLAS scaling is better, but of course the tunded BLAS is 
   still better.  

*/ 
  int i,j,N,m,*K,*C,*R,*B,l,r,c,nb;
  char uplo = 'U',trans = 'T',ntrans = 'N';
  m = *nt;N = (m * (m + 1))/2; /* total number of blocks to process (upper bound)*/
  K = work;work += m+1;C = work;work += N;
  R = work; work += N;B = work;
  tile_ut(*n,&m,K,C,R,B); /* set up tiling and allocation of tiles to threads */
 
  #ifdef OPENMP_ON
  #pragma omp parallel for private(i,j,r,c,nb,l) num_threads(m)
  #endif
  for (i=0;i<m;i++) {
    for (j=B[i];j<B[i+1];j++) {
      r=R[j];c=C[j];
      if (r==c) { /* symmetric block, use  dsyrk */
	nb = K[r+1]-K[r]; /* block dimension */
        F77_CALL(dsyrk)(&uplo,&trans,&nb,k,a,A + (ptrdiff_t) *lda * K[r],lda,b,D + (ptrdiff_t) *ldd * K[c] + K[r],ldd);
      } else { /* general block, use dgemm */
        nb = K[r+1]-K[r]; l = K[c+1]-K[c]; /* block dimension is nb by l */
	F77_CALL(dgemm)(&trans,&ntrans,&nb,&l,k,a,A + (ptrdiff_t)*lda * K[r],lda,A + *lda * K[c],lda,b,
			D + (ptrdiff_t) *ldd * K[c] + K[r],ldd);
	
      }	
    }  
  }  
  
} /* pdsyrk */  


void mgcv_tensor_mm(double *X,double *T,int *d,int *m,int *n) {
/* Code for efficient production of row tensor product model matrices.
   X contains rows of matrices to be producted. Contains m matrices,
   d[i] is number of columns in ith matrix (which are stored in ascending 
   order). Each column has n rows. T is the target matrix, with n rows
   and \prod_i d[i] columns.   
*/
  ptrdiff_t start, i,j,k, tp=1, xp=0,pd;
  double *Xj,*Xj1,*Xk, *Tk,*p,*p1,*p2;
  /*Rprintf("m = %d  n = %d  d = ",*m,*n);
    for (i=0;i<*m;i++) Rprintf(" %d,",d[i]);*/
  /* compute total columns in X, xp, and total columns in T, tp ... */
  for (i=0;i<*m;i++) { xp += d[i];tp *= d[i];} 
  Xk = X + (xp-d[*m-1]) * *n; /* start of last matrix in X */
  Tk = T + (tp-d[*m-1]) * *n; /* start of last (filled) block in T */
  /* initialize by putting final matrix in X into end block of T... */
  p = Xk;p1 = Xk + *n * (ptrdiff_t) d[*m-1];p2 = Tk;
  for (;p<p1;p++,p2++) *p2 = *p;
  pd = d[*m-1]; /* cols in filled (end) block of T */
  for (i = *m - 2;i>=0;i--) { /* work down through matrices stored in X */
    Xk -= *n * (ptrdiff_t) d[i]; /* start of ith X matrix */
    Xj = Xk;
    start = tp - pd * d[i]; /* start column of target block in T */
    p = T + start * *n; /* start location in T */
    for (j=0;j<d[i];j++) { /* work through columns of ith X matrix */
      p1 = Tk; /* start of T block to be multiplied by jth col of ith X matrix*/
      Xj1 = Xj + *n; /* end of col j of ith X matrix */
      /* following multiplies end T block by jth col of ith X matrix, storing 
         result in appropriate T block... */
      for (k=0;k<pd;k++) for (p2 = Xj;p2<Xj1;p2++,p1++,p++) *p = *p1 * *p2;
      Xj += *n; /* start of col j+1 of ith X matrix */
    }
    pd *= d[i]; /* number of cols in filled in T block */
    Tk = T + (tp - pd) * *n; /* start of filled in T block */
  }  
} /* mgcv_tensor_mm */

void mgcv_tmm(SEXP x,SEXP t,SEXP D,SEXP M, SEXP N) {
  /* wrapper for calling mgcv_tensor_mm using .Call */
  double *X,*T;
  int *d,*m,*n;
  X = REAL(x);T=REAL(t);
  d = INTEGER(D);
  m = INTEGER(M);
  n = INTEGER(N);
  mgcv_tensor_mm(X,T,d,m,n);
}

void mgcv_mmult0(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n)
/* This code doesn't rely on the BLAS...
 
   Not thread safe as C modified and restored internally...

   Forms r by c product of B and C, transposing each according to bt and ct.
   n is the common dimension of the two matrices, which are stored in R 
   default column order form. Algorithm is inner loop optimized in each case 
   (i.e. inner loops only update pointers by 1, rather than jumping).
   
  
   r<-1000;c<-1000;n<-1000
   A<-matrix(0,r,c);B<-matrix(rnorm(r*n),n,r);C<-matrix(rnorm(n*c),c,n)
   system.time(
   A<-matrix(.C("mgcv_mmult",as.double(A),as.double(B),as.double(C),as.integer(1),as.integer(1),
          as.integer(r),as.integer(c),as.integer(n))[[1]],r,c))
   range(as.numeric(t(B)%*%t(C)-A))
   A<-matrix(0,r,c);B<-matrix(rnorm(r*n),n,r);C<-matrix(rnorm(n*c),n,c)
   system.time(
   A<-matrix(.C("mgcv_mmult",as.double(A),as.double(B),as.double(C),as.integer(1),as.integer(0),
          as.integer(r),as.integer(c),as.integer(n))[[1]],r,c))
   range(as.numeric(t(B)%*%C-A))
   A<-matrix(0,r,c);B<-matrix(rnorm(r*n),r,n);C<-matrix(rnorm(n*c),c,n)
   system.time(
   A<-matrix(.C("mgcv_mmult",as.double(A),as.double(B),as.double(C),as.integer(0),as.integer(1),
          as.integer(r),as.integer(c),as.integer(n))[[1]],r,c))
   range(as.numeric(B%*%t(C)-A))
   A<-matrix(0,r,c);B<-matrix(rnorm(r*n),r,n);C<-matrix(rnorm(n*c),n,c)
   system.time(
   A<-matrix(.C("mgcv_mmult",as.double(A),as.double(B),as.double(C),as.integer(0),as.integer(0),
          as.integer(r),as.integer(c),as.integer(n))[[1]],r,c))
   range(as.numeric(B%*%C-A))  
*/

{ double xx,*bp,*cp,*cp1,*cp2,*cp3,*ap,*ap1;
  int cr,i,j;
  if (*bt)
  { if (*ct) /* A=B'C' */
    { /* this one is really awkward: have to use first row of C' as working storage
         for current A row; move most slowly through B */
      for (i=0;i<*r;i++) /* update A *row-wise*, one full sweep of C per i */    
      { cp1 = C + *c;
        /* back up row 1 of C' (in A), and initialize current row of A (stored in row 1 of C') */
        xx = *B;
        for (ap=A,cp=C;cp<cp1;cp++,ap+= *r) { *ap = *cp; *cp *= xx;}  
        B++;cp2=cp1;
        for (j=1;j< *n;j++,B++) 
	    for (xx= *B,cp=C;cp<cp1;cp++,cp2++) *cp += xx * *cp2;      
        /* row i of A is now in row 1 of C', need to swap back */
        for (ap=A,cp=C;cp<cp1;cp++,ap+= *r) { xx = *ap; *ap = *cp; *cp = xx;}
        A++;
      } 
    } else {/* A=B'C - easiest case: move most slowly through A*/
      cr= *n;cp2 = C + *c * (ptrdiff_t) cr;
      for (ap=A,cp1=C;cp1< cp2;cp1+=cr) for (bp=B,i=0;i< *r;i++,ap++)  
      { for (xx=0.0,cp=cp1,cp3=cp1+ *n;cp< cp3;cp++,bp++) xx += *cp * *bp; /* B[k+br*i]*C[k+cr*j];*/
        *ap=xx;
      }
    }
  } else
  { if (*ct) /* A=BC' - update A column-wise, move most slowly through C (but in big steps)*/
    { cp = C;
      for (j=0;j< *c;j++) /* go through columns of A, one full B sweep per j */
      { ap1 = A + *r;C=cp;xx = *C;bp=B;
        for (ap=A;ap<ap1;ap++,bp++) *ap = xx * *bp ;
        C += *c;
        for (i=1;i< *n;i++,C+= *c) 
	for (xx=*C,ap=A;ap<ap1;ap++,bp++) *ap += xx * *bp;
        A=ap1;cp++;
      }
    } else /* A=BC - update A column-wise, moving most slowly through C */
    { for (j=0;j< *c;j++) /* go through columns of A, one full B sweep per j */
      { ap1 = A + *r;xx = *C;bp=B;
        for (ap=A;ap<ap1;ap++,bp++) *ap = xx * *bp ;
        C++;
        for (i=1;i< *n;i++,C++) 
	for (xx=*C,ap=A;ap<ap1;ap++,bp++) *ap += xx * *bp;
        A=ap1;
      }
    }
  }
} /* end mgcv_mmult0*/

void mgcv_mmult(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n) {
  /* Forms A = B C transposing B and C according to bt and ct.
     A is r by c. Common dimension is n;
     BLAS version A is c (result), B is a, C is b, bt is transa ct is transb 
     r is m, c is n, n is k.
     Does nothing if r,c or n <= zero. 
  */
  char transa='N',transb='N';
  int lda,ldb,ldc;
  double alpha=1.0,beta=0.0;
  if (*r<=0||*c<=0||*n<=0) return;
  if (B==C) {
    if (*bt&&(!*ct)&&(*r==*c)) { getXtX(A,B,n,r);return;} 
    else if (*ct&&(!*bt)&&(*r==*c)) { getXXt(A,B,c,n);return;}
  }
  if (*bt) { /* so B is n by r */
    transa = 'T';  
    lda = *n;
  } else  lda = *r; /* B is r by n */

  if (*ct) { /* C is c by n */ 
    transb = 'T';
    ldb = *c;
  } else ldb = *n; /* C is n by c */
  
  ldc = *r;
  F77_CALL(dgemm)(&transa,&transb,r,c,n, &alpha,
		B, &lda,C, &ldb,&beta, A, &ldc);
} /* end mgcv_mmult */

SEXP mgcv_madi(SEXP a, SEXP b,SEXP ind,SEXP diag) {
/* Performs 
     a[ind,ind] <- a[ind,ind] + b
   or, if diag != 0,
     diag(a)[ind] <- diag(a)[ind] + b
*/
  int dia,*ii,ij,n,i,j,k;
  double *A,*B;
  SEXP kr;
  dia = asInteger(diag);n = nrows(a); k = length(ind);
  ind = PROTECT(coerceVector(ind,INTSXP));
  b = PROTECT(coerceVector(b,REALSXP));
  a = PROTECT(coerceVector(a,REALSXP));
  ii = INTEGER(ind);
  A = REAL(a);
  B = REAL(b);
   if (dia==0) {
    for (j=0;j<k;j++) { /* B is a matrix */
      ij = (ii[j]-1)*n-1;
      for (i=0;i<k;i++,B++) A[ii[i] + ij] += *B;
    }  
  } else if (dia > 0) { /* B is a vector */
    for (i=0;i<k;i++,B++) {
      ij = (ii[i]-1)*(n+1); A[ij] += *B;
    }
  } else { /* B is a scalar */
    for (i=0;i<k;i++) {
      ij = (ii[i]-1)*(n+1);
        A[ij] += *B;
     }
  }
  
  PROTECT(kr=allocVector(REALSXP,1));
  REAL(kr)[0] = 1.0; 
  UNPROTECT(4);
  return(kr); /* dummy so that mtrace doesn't cause a crash on return! */
} /* mgcv_madi */

SEXP mgcv_pmmult2(SEXP b, SEXP c,SEXP bt,SEXP ct, SEXP nthreads) {
/* parallel matrix multiplication using .Call interface */
  double *A,*B,*C;
  int r,col,n,nt,Ct,Bt;
  #ifdef OPENMP_ON
  int m;
  #endif

  SEXP a; 
 
  nt = asInteger(nthreads);
  Bt = asInteger(bt);
  Ct = asInteger(ct);
 
  if (Bt) {
    r = ncols(b);
    n = nrows(b);
  } else {
    r = nrows(b);
    n = ncols(b);
  }
 
  if (Ct) col = nrows(c);
  else col = ncols(c);
  /* Rprintf("nt = %d, Bt = %d, Ct = %d, r = %d, c = %d, n = %d\n",nt,Bt,Ct,r,col,n);*/
  B = REAL(b);C=REAL(c);
  a = PROTECT(allocMatrix(REALSXP,r,col));
  A = REAL(a);
  #ifdef OPENMP_ON
  m = omp_get_num_procs(); /* detected number of processors */
  if (nt > m || nt < 1) nt = m; /* no point in more threads than m */
  /*Rprintf("\n open mp %d cores, %d used\n",m,nt);*/
  #else
  /*Rprintf("\n no openmp\n");*/
  nt = 1; /* no openMP support - turn off threading */
  #endif
  mgcv_pmmult(A,B,C,&Bt,&Ct,&r,&col,&n,&nt);
  UNPROTECT(1);
  return(a);
} /* mgcv_pmmult2 */

int mgcv_bchol0(double *A,int *piv,int *n,int *nt,int *nb) {
/* Lucas (2004) "LAPACK-Style Codes for Level 2 and 3 Pivoted Cholesky Factorizations" 
   block pivoted Choleski algorithm 5.1. Note some misprints in paper, noted below. 
   nb is block size, nt is number of threads, A is symmetric
   +ve semi definite matrix and piv is pivot sequence. 

   This version is BLAS free.
*/  
  int i,j,k,l,q,r=-1,*pk,*pq,jb,n1,m,N,*a,b;
  double tol=0.0,*dots,*pd,*p1,*Aj,*Aq0,*Aj0,*Aj1,*Ajn,*Ail,xmax,x,*Aq,*Ajj,*Aend;
  dots = (double *)CALLOC((size_t) *n,sizeof(double));
  for (pk = piv,i=0;i < *n;pk++,i++) *pk = i; /* initialize pivot record */
  jb = *nb; /* block size, allowing final to be smaller */
  n1 = *n + 1;
  Ajn = A;
  m = *nt;if (m<1) m=1;if (m>*n) m = *n; /* threads to use */
  a = (int *)CALLOC((size_t) (*nt+1),sizeof(int)); /* thread block cut points */
  a[m] = *n;
  for (k=0;k<*n;k+= *nb) {
    if (*n - k  < jb) jb = *n - k ; /* end block */ 
    for (pd = dots + k,p1 = dots + *n;pd<p1;pd++) *pd = 0;
    for (j=k;j<k+jb;j++,Ajn += *n) {
      pd = dots + j;Aj = Ajn + j; Aj1 = Aj - 1;
      xmax = -1.0;q=j;p1 = dots + *n;
      if (j>k) for (;pd<p1;pd++,Aj1 += *n) *pd += *Aj1 * *Aj1; /* dot product update */
      for (l=j,pd = dots + j;pd<p1;pd++,Aj += n1,l++) {   
        x = *Aj - *pd; 
        if (x>xmax) { xmax = x;q=l;} /* find the pivot */
      } 
      if (j==0) tol = *n * xmax * DOUBLE_EPS;
      Aq = A + *n * q + q;
      // Rprintf("\n n = %d k = %d j = %d  q = %d,  A[q,q] = %g  ",*n,k,j,q,*Aq);
      if (*Aq - dots[q]<tol) {r = j;break;} /* note Lucas (2004) has 'dots[q]' missing */
      /* swap dots... */
      pd = dots + j;p1 = dots + q;
      x = *pd;*pd = *p1;*p1 = x;
      /* swap pivots... */
      pk = piv + j;pq = piv +q;
      i = *pk;*pk = *pq;*pq = i;
      /* swap rows ... */
      Aend = A + *n * *n;
      Aj = Ajn + j; Aq = Ajn + q;
      for (;Aj<Aend;Aj += *n,Aq += *n) {
        x = *Aj;*Aj = *Aq;*Aq = x;
      }
      /* swap cols ... */
      Aj = Ajn; Aq = A + *n * q;Aend = Aj + *n;
      for (;Aj < Aend;Aj++,Aq++)  {
        x = *Aj;*Aj = *Aq;*Aq = x;
      }
      /* now update */
      Ajj = Ajn + j;  
      // Rprintf(" %g  %g",*Aj,*pd);
      *Ajj = sqrt(*Ajj - *pd); /* sqrt(A[j,j]-dots[j]) */      
      Aend = A + *n * *n;
      if (j > k&&j < *n) { /* Lucas (2004) has '1' in place of 'k' */
        Aj = Ajn + *n;
        Aq = Aj + k; /* Lucas (2004) has '1' in place of 'k' */
        Aj += j;        
        Aj1 = Ajn + k; /* Lucas (2004) has '1' in place of 'k' */
        for (;Aj<Aend;Aj += *n,Aq += *n) {
	  for (pd = Aj1,p1=Aq;pd < Ajj;pd++,p1++) *Aj -= *pd * *p1;
          //*Aj -= x;
        }  
      }
      if (j < *n) {
        Aj = Ajj; x = *Aj;Aj += *n;
        for (;Aj<Aend;Aj += *n) *Aj /= x;
      }    
    } /* j loop */
    if (r > 0) break;
    /* now the main work - updating the trailing factor... */
  
    if (k + jb < *n) {
      /* create the m work blocks for this... */
      N = *n - j; /* block to be processed is N by N */
      if (m > N) { m = N;a[m] = *n; } /* number of threads to use must be <= r */
      *a = j; /* start of first block */
      x = (double) N;x = x*x / m;
      /* compute approximate optimal split... */
      for (i=1;i < m;i++) a[i] = round(N - sqrt(x*(m-i)))+j;
      for (i=1;i <= m;i++) { /* don't allow zero width blocks */
          if (a[i]<=a[i-1]) a[i] = a[i-1]+1;
      }     
      #ifdef OPENMP_ON
#pragma omp parallel private(b,i,l,Aj,Aend,Aq,Aj1,Ail,Aj0,Aq0) num_threads(m)
      #endif 
      { /* start parallel section */
        #ifdef OPENMP_ON
        #pragma omp for
        #endif
        for (b=0;b<m;b++)
	for (i=a[b];i<a[b+1];i++) { 
          Aj0 = A + i * *n;
          Aq0 = Aj0 + k;
          Aend = Aj0 + j;
          Ail = Aj1 = Aj0 + i;
          Aj0 += k; 
	  for (l=i;l<*n;l++,Aq0 += *n,Aj1++,Ail += *n) {
	    // Aj = A + i * *n;Aend = Aj + j;Aj1 = Aj + l;
            //Aj+=k;Aq = A + l * *n + k; /* Lucas (2004) has '1' in place of 'k' */
            Aj = Aj0;
            Aq = Aq0;
            for (;Aj < Aend;Aj++,Aq++) *Aj1 -= *Aq * *Aj;
            // A[i + *n * l ] = *Aj1;
            *Ail = *Aj1;
          }
	}
      } /* end parallel section */
    } /* if (k + jb < *n) */
  } /* k loop */
  if (r<0) r = *n;
  FREE(dots);
  
  for (Ajn=A,j=0;j<*n;j++,Ajn += *n) {
    Aj = Ajn;Aend = Aj + *n;
    if (j<r) Aj += j+1; else Aj += r;
    for (;Aj<Aend;Aj++) *Aj = 0.0;
  }
  FREE(a);
  return(r);
} /* mgcv_bchol0 */

int mgcv_bchol(double *A,int *piv,int *n,int *nt,int *nb) {
/* Lucas (2004) "LAPACK-Style Codes for Level 2 and 3 Pivoted Cholesky Factorizations" 
   block pivoted Choleski algorithm 5.1. Note some misprints in paper, noted below. 
   nb is block size, nt is number of threads, A is symmetric
   +ve semi definite matrix and piv is pivot sequence. 

   This version calls BLAS routines at least as efficiently as the LAPACK version (it's slightly faster). 
   A parallel version of dsyrk is then called for the key trailing block update.
*/  
  int i,j,k,l,q,r=-1,*pk,*pq,jb,n1,m,N,*a,one=1,*work;
  double tol=0.0,*dots,*pd,*p1,*Aj,*Aj1,*Ajn,xmax,x,*Aq,*Ajj,*Aend,alpha=-1.0,beta=1.0;
  char uplo = 'L',trans='N';
  dots = (double *)CALLOC((size_t) *n,sizeof(double));
  work = (int *)CALLOC((size_t) *nt * (*nt+3) + 2,sizeof(int));
  for (pk = piv,i=0;i < *n;pk++,i++) *pk = i; /* initialize pivot record */
  jb = *nb; /* block size, allowing final to be smaller */
  n1 = *n + 1;
  Ajn = A;
  m = *nt;if (m<1) m=1;if (m>*n) m = *n; /* threads to use */
  a = (int *)CALLOC((size_t) (*nt+1),sizeof(int)); /* thread block cut points */
  a[m] = *n;
  for (k=0;k<*n;k+= *nb) {
    if (*n - k  < jb) jb = *n - k ; /* end block */ 
    for (pd = dots + k,p1 = dots + *n;pd<p1;pd++) *pd = 0;
    for (j=k;j<k+jb;j++,Ajn += *n) {
      pd = dots + j;Aj = Ajn + j; Aj1 = Aj - 1;
      xmax = -1.0;q=j;p1 = dots + *n;
      if (j>k) for (;pd<p1;pd++,Aj1 += *n) *pd += *Aj1 * *Aj1; /* dot product update (upper triangle only) */
      for (l=j,pd = dots + j;pd<p1;pd++,Aj += n1,l++) {   
        x = *Aj - *pd; 
        if (x>xmax) { xmax = x;q=l;} /* find the pivot q >= j (leading diag only used)*/
      } 
      if (j==0) tol = *n * xmax * DOUBLE_EPS;
      Aq = A + *n * q + q; 
      if (*Aq - dots[q]<tol) {r = j;break;} /* note Lucas (2004) has 'dots[q]' missing */
      /* swap dots... */
      pd = dots + j;p1 = dots + q;
      x = *pd;*pd = *p1;*p1 = x;
      /* swap pivots... */
      pk = piv + j;pq = piv +q;
      i = *pk;*pk = *pq;*pq = i;

      /* dswap(int *n,double *x,int *dx,double *y,int *dy)
         Swaps n elements of vectors x and y. Spacing between elements is dx and dy. 
      */

      /* row-col exchange j and q - only update above the diagonal, but rows j and q 
         are exchanged only from col j, while cols j and q are fully exchanged. 
         This is because we are working through cols of factor and cols before j are
         finished. 
      */
  
      Aj = Ajn + j; 
      x =  *Aj; *Aj = *Aq;
      *Aq = x; /* A[j,j] <-> A[q,q] */
      /* A[j,j+1:q-1] <-> A[j+1:q-1,q] */
      N = q-j-1;
      if (N>0) {
        Aj += *n; /* A[j,j+1] */
        Aq = A + q * *n + j + 1; /* A[j+1,q] */
        F77_CALL(dswap)(&N,Aj,n,Aq,&one);
      }	
      /* A[q,q+1:n-1] <-> A[j,q+1:n-1] */
      N = *n-q-1;
      if (N>0) {
        Aq = A + (q+1) * *n + q; /* A[q,q+1] */
        Aj = A + (q+1) * *n + j; /* A[j,q+1] */
        F77_CALL(dswap)(&N,Aj,n,Aq,n);
      }	
      /* A[0:(j-1),j] <-> A[0:(j-1),q] --- note we have to start at row 0 (not k) or
         already complete rows will not end up correctly pivoted */
      N = j;
      if (N>0) {
        Aq = A + q * *n;
        Aj = Ajn;
        F77_CALL(dswap)(&N,Aj,&one,Aq,&one);
      }
      /* now update (only accesses upper triangle) */

     /* dgemv(char *trans,int *m, int *n,double *a,double *A,int *lda,double *x,int *dx,
         double *b, double *y,int *dy)
         trans='T' to transpose A, 'N' not to. A is m by n. Forms y = a*A'x + b*y, or
         y = a*Ax + b*y. lda is number of actual rows in A (to allow for sub-matrices)
         dx and dy are increments of x and y indices.  */
     
      Ajj = Ajn + j;  
      *Ajj = sqrt(*Ajj - *pd); /* sqrt(A[j,j]-dots[j]) */      
      Aend = A + *n * *n;
      if (j > k&&j < *n) { /* Lucas (2004) has '1' in place of 'k' */
        /* A[j,j+1:n-1] += -A[k:j-1,j]*A[k:j-1,j+1:n-1] */ 
	trans='T';N = *n - j - 1;i=j-k;
        F77_CALL(dgemv)(&trans,&i,&N,&alpha,A+(j+1) * *n + k,n,A + j * *n + k,&one,&beta,A + (j+1) * *n+j,n);
      }
      if (j < *n) {
        Aj = Ajj; x = *Aj;Aj += *n;
        for (;Aj<Aend;Aj += *n) *Aj /= x; /* upper triangle access only */ 
      }    
    } /* j loop */
    if (r > 0) break;
    /* now the main work - updating the trailing factor... 
       A[j:n-1,j:n-1] += - A[k:j-1,j:n-1]'A[k:j-1,j:n-1]
    */
    /* dsyrk(char *uplo, char *trans,int *n, int *k,double *a, double *A, int *lda,
         double *b, double *C,int *ldc)
       uplo = 'U' or 'L' for upper or lower tri of C used. trans = 'N' for C = aAA' + bC  
       and A n by k, or 'T' for C =  aA'A + bC and A k by n; C is n by n. lda and ldc
       are actual number of rows in A and C respectively (allows use on submatrices).
    */
    if (k + jb < *n) {
      N = *n - j ; /* block to be processed is N by N */
      i = j - k; /* number of rows in */
      trans = 'T';uplo='U';//alpha = -1.0;beta=1.0;
      if (0) F77_CALL(dsyrk)(&uplo,&trans,&N, &i, &alpha,A+j * *n + k,n,&beta,A + j * *n + j,n);
      else pdsyrk(&N,&i,&alpha,A+j * *n + k,n,&beta,A + j * *n + j,n,work,nt);
      //rpmat(A,*n);
    }   
  } /* k loop */
  if (r<0) r = *n;
  FREE(dots);
  
  for (Ajn=A,j=0;j<*n;j++,Ajn += *n) {
    Aj = Ajn;Aend = Aj + *n;
    if (j<r) Aj += j+1; else Aj += r;
    for (;Aj<Aend;Aj++) *Aj = 0.0;
  }
  FREE(a);FREE(work);
  return(r);
} /* mgcv_bchol */

int mgcv_pchol(double *A,int *piv,int *n,int *nt) {
/* Obtain pivoted Choleski factor of n by n matrix A using 
   algorithm 4.2.4 of Golub and van Loan 3 (1996) and using 
   openMP to parallelize on nt cores. 
   Note that 4.2.4 is only correct if the pivot instructions are 
   replaced by those below, which assume that only the lower triangle 
   of A is to be accessed. GvL4 (2013) Alg 4.2.2 is very similar, and works 
   as is, with a much simpler correction of the pivoting, but in that 
   case both triangles of A are accessed and the cost doubles. Modifying
   to use only the lower triangle again gives the following. 

   This version is almost identical to LAPACK with default BLAS in speed 
   (4% slower at n=4000), as a single thread algorithm.

*/
  int i,j,k,r,q,n1,*pk,*pq,kn,qn,*a,N,m,b;
  double *Aj,*Ak,*Aq,*Aend,x,Ajk,Akk,thresh=0.0;
  if (*nt < 1) *nt = 1;
  if (*nt > *n) *nt = *n;
  m = *nt;  
  a = (int *)CALLOC((size_t) (*nt+1),sizeof(int));
  a[0] = 0;a[m] = *n; /* ... initialize column block splitting array */ 
  r = 0;n1 = *n + 1;
  for (pk = piv,i=0;i < *n;pk++,i++) *pk = i; /* initialize pivot record */
  for (pk=piv,k=0;k< *n;k++,pk++) {
    kn = k * *n;
    /* find largest element of diag(A), from k onwards */
    Ak = A + kn + k;x = *Ak;q=k;Ak+=n1;
    for (i=k+1;i < *n;i++,Ak+=n1) if (*Ak>x) {x = *Ak;q=i;}
    qn = q * *n;
    if (k==0) thresh = *n * x * DOUBLE_EPS;
    if (x>thresh) { /* A[q,q] =x > 0 */
      r++;
      /* piv[k] <-> piv[q] */
      pq = piv + q;i = *pq; *pq = *pk;*pk = i;
      
      /* A[k,k] <-> A[q,q] */
      Ak = A + kn + k;Aq = A + qn + q;
      x = *Ak;*Ak = *Aq;*Aq = x;
      /* A[k+1:q-1,k] <-> A[q,k+1:q-1] */
      Ak++; Aend = Aq; Aq = A + q + kn + *n;
      for (;Aq<Aend;Ak++,Aq += *n) {
        x = *Ak;*Ak = *Aq;*Aq = x;
      }
      /* A[q,1:k-1] <-> A[k,1:k-1] */
      Ak = A + k;Aend=Ak + kn;Aq = A + q;
      for (;Ak < Aend;Ak += *n,Aq += *n) {x = *Aq;*Aq = *Ak;*Ak = x;}
      /* A[q+1:n,k] <-> A[q+1:n,q] */
      Ak = A + kn; Aq = A + qn+q+1;Aend = Ak + *n;Ak+=q+1;
      for (;Ak<Aend;Ak++,Aq++) {x = *Aq;*Aq = *Ak;*Ak = x;}
      /* kth column scaling */
      Ak = A + kn + k;*Ak = sqrt(*Ak);Akk = *Ak;Ak++;
      for (;Ak<Aend;Ak++) *Ak /= Akk;
      /* finally the column updates...*/
      /* create the nt work blocks for this... */
      N = *n - k - 1; /* block to be processed is r by r */
      if (m > N) { m = N;a[m] = *n; } /* number of threads to use must be <= r */
      (*a)++;
      x = (double) N;x = x*x / m;
      /* compute approximate optimal split... */
      for (i=1;i < m;i++) a[i] = round(N - sqrt(x*(m-i)))+k+1;
      for (i=1;i <= m;i++) { /* don't allow zero width blocks */
        if (a[i]<=a[i-1]) a[i] = a[i-1]+1;
      }    
      /* check load balance... */
      // for (i=0;i<m;i++) Rprintf("%d  ",(a[i+1]-a[i])*(*n-(a[i+1]-1+a[i])/2));Rprintf("\n");
      #ifdef OPENMP_ON
      #pragma omp parallel private(b,j,Aj,Aend,Ak,Ajk) num_threads(m)
      #endif 
      { /* begin parallel section */
        #ifdef OPENMP_ON
        #pragma omp for
        #endif
        for (b=0;b<m;b++)
	  for (j=a[b];j<a[b+1];j++) {  
	    Aj = A + *n * j + j;
            Ak = A + kn; Aend = Ak + *n;Ak += j;Ajk = *Ak;
            for (;Ak < Aend;Ak++,Aj++) *Aj -= *Ak * Ajk;
          }
      } /* end parallel section */
    } else break; /* if x not greater than zero */
  } /* k loop */ 
  /* Wipe redundant columns... */ 
  Aend = A + *n * *n; /* 1 beyond end of A */  
  for (Ak = A + r * *n;Ak<Aend;Ak++) *Ak = 0.0;
  
  /* transpose and wipe lower triangle ... */
  /* first compute block split */
  a[0] = 0;a[*nt] = *n;
  x = (double) *n;x = x*x / *nt;
  /* compute approximate optimal split... */
  for (i=1;i < *nt;i++) a[i] = round(*n - sqrt(x*(*nt-i)));
  for (i=1;i <= *nt;i++) { /* don't allow zero width blocks */
    if (a[i]<=a[i-1]) a[i] = a[i-1]+1;
  }
  #ifdef OPENMP_ON
  #pragma omp parallel private(b,i,Ak,Aend,Aq) num_threads(*nt)
  #endif 
  { /* start parallel */
    #ifdef OPENMP_ON
    #pragma omp for
    #endif
    for (b=0;b < *nt;b++)
      for (i=a[b];i<a[b+1];i++) {
        Ak = A + i * *n; Aend = Ak + *n;Ak += i;
        Aq = Ak + *n; /* row i pointer A[i,i+1] */
        Ak ++; /* col i pointer A[i+1,i] */
        for (;Ak<Aend;Ak++,Aq += *n) { /* copy and wipe */
          *Aq = *Ak;*Ak=0.0;
        }
      }
  } /* end parallel */
  /* wipe upper triangle code 
    for (i = 0 ; i < r;i++) {
    Ak = A + *n * i;Aend = Ak + i;
    for (;Ak<Aend;Ak++) *Ak = 0.0; 
    }*/
  
  FREE(a);
  return(r); 
} /* mgcv_pchol */


SEXP mgcv_Rpchol(SEXP Amat,SEXP PIV,SEXP NT,SEXP NB) {
/* routine to Choleski decompose n by n  matrix Amat with pivoting 
   using routine mgcv_bchol to do the work. Uses NT parallel
   threads. NB is block size.

   Designed for use with .call rather than .C
   Returns detected rank, and overwrites Amat with its pivoted
   Choleski factor. PIV contains pivot vector.

*/
  int n,nt,*piv,r,*rrp,nb;
  double *A;
  SEXP rr;
  nb = asInteger(NB);
  nt = asInteger(NT);
  n = nrows(Amat);
  A = REAL(Amat);
  piv = INTEGER(PIV);
  // r = mgcv_pchol(A,piv,&n,&nt); 
  r = mgcv_bchol(A,piv,&n,&nt,&nb);
  /* should return rank (r+1) */
  rr = PROTECT(allocVector(INTSXP, 1));
  rrp = INTEGER(rr);
  *rrp = r;
  UNPROTECT(1);
  return(rr);
} /* mgcv_Rpchol */


int bpqr(double *A,int n,int p,double *tau,int *piv,int nb,int nt) {
/* block pivoted QR based on 
    Quintana-Orti, G.; Sun, X. & Bischof, C. H. (1998)
    "A BLAS-3 version of the QR factorization with column pivoting" 
    SIAM Journal on Scientific Computing, 19: 1486-1494
   but note quite alot of correction of algorithm as stated in paper (there 
   are a number of indexing errors there, and down-date cancellation 
   strategy is only described in words). 
*/ 
  int jb,pb,i,j,k=0,m,*p0,nb0,q,one=1,ok_norm=1,*mb,*kb,rt,nth;
  double *cn,*icn,x,*a0,*a1,*F,*Ak,*Aq,*work,tol,xx,done=1.0,dmone=-1.0,dzero=0.0; 
  char trans='T',nottrans='N';
#ifdef OMP_REPORT
  Rprintf("bpqr...");
#endif
  tol = pow(DOUBLE_EPS,.8);
  mb = (int *)CALLOC((size_t) nt,sizeof(int));
  kb = (int *)CALLOC((size_t) nt,sizeof(int));  
  for (p0=piv,i=0;i<p;i++,p0++) *p0 = i; /* initialize pivot index */
  work = (double *)CALLOC((size_t) nb,sizeof(double));
  /* create column norm vectors */
  cn =(double *)CALLOC((size_t) p,sizeof(double)); /* version for downdating */
  icn =(double *)CALLOC((size_t) p,sizeof(double)); /* inital for judging cancellation error */
  for (a0=A,i=0;i<p;i++) { /* cn[i]=icn[i] = || A[:,i] ||^2 */
    for (x=0.0,a1 = a0 + n;a0<a1;a0++) x += *a0 * *a0;
    cn[i] = icn[i] = x;   
  }
  if (p<nb) nb = p;
  nb0 = nb; /* target block size nb */  
  jb=0; /* start column of current block */
  pb = p; /* columns left to process */
  F = (double *)CALLOC((size_t) p * nb0,sizeof(double));
  while (jb < p) {
    nb = p-jb;if (nb>nb0) nb = nb0;/* attempted block size */
    for (a0=F,a1=F+nb*pb;a0<a1;a0++) *a0 = 0.0; /* F[1:pb,1:nb] = 0 - i.e. clear F */
    for (j=0;j<nb;j++) { /* loop through cols of this block */
      k = jb + j; /* index of current column of A to process */
      a0 = cn + k;x = *a0;q=k;a0++;
      for (i=k+1;i<p;i++,a0++) if (*a0>x) { x = *a0;q=i; } /* find pivot col q */
      if (q!=k) { /* then pivot */
        i = piv[q];piv[q]=piv[k];piv[k] = i;
        x = cn[q];cn[q]=cn[k];cn[k] = x;
        x = icn[q];icn[q]=icn[k];icn[k] = x;
        Aq = A + q * (ptrdiff_t) n;Ak = A + k * (ptrdiff_t) n;a1 = Aq + n;
        for (;Aq<a1;Aq++,Ak++) { x = *Aq; *Aq = *Ak; *Ak = x;} /* A[:,k] <-> A[:,q] */
        Aq = F + q - jb;Ak = F + j;a1 = F + nb * (ptrdiff_t) pb;
        for (;Aq<a1;Aq+=pb,Ak+=pb) { x = *Aq; *Aq = *Ak; *Ak = x;} /* F[q-jb,:] <-> F[j,:] */        
      }
      /* update the pivot column: A[k:n-1,k] -= A[k:n-1,jb:k-1]F[j,0:j-1]' using
         BLAS call to DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
                    y := alpha*A*x + beta*y (or A' if TRANS='T')*/
      m = n-k;Ak = A + (ptrdiff_t)n * k + k;
      if (j) {
        q = m ; /* total number of rows to split between threads */
        rt = q/nt;if (rt*nt < q) rt++; /* rows per thread */
        nth = nt; while (nth>1&&(nth-1)*rt>q) nth--; /* reduce number of threads if some empty */
        kb[0] = k; /* starting row */
        for (i=0;i<nth-1;i++) { mb[i] = rt;kb[i+1]=kb[i]+rt;}
        mb[nth-1]=q-(nth-1)*rt;
        #ifdef OPENMP_ON
        #pragma omp parallel private(i) num_threads(nth)
        #endif 
        { /* start of parallel section */
          #ifdef OPENMP_ON
          #pragma omp for
          #endif
          for (i=0;i<nth;i++) {
	    /* Only 10th argument changed on exit... */ 
            F77_CALL(dgemv)(&nottrans, mb+i, &j,&dmone,A+jb*(ptrdiff_t)n+kb[i],&n,F+j,&pb,&done,
	                    A + (ptrdiff_t)n*k + kb[i], &one);
            //F77_CALL(dgemv)(&nottrans, &m, &j,&dmone,A+jb*n+k,&n,F+j,&pb,&done,Ak, &one);
          }
        }
      }
      /* now compute the Householder transform for A[,k], by calling 
         dlarfg (N, ALPHA, X, INCX, TAU), N is housholder dim */
      xx = *Ak; /* A[k,k] on entry and transformed on exit */
      Ak++; /* on exit all but first element of v will be in Ak, v(1) = 1 */
      F77_CALL(dlarfg)(&m,&xx,Ak,&one,tau+k);
      Ak--;*Ak = 1.0; /* for now, have Ak point to whole of v */
      /* F[j+1:pb-1,j] = tau[k] * A[k:n-1,k+1:p-1]'v */ 
      
      if (k<p-1) { /* up to O(np) step - most expensive after block */
        // i=p-k-1;
        q = p - k - 1 ; /* total number of rows to split between threads */
        rt = q/nt;if (rt*nt < q) rt++; /* rows per thread */
        nth = nt; while (nth>1&&(nth-1)*rt>q) nth--; /* reduce number of threads if some empty */
        kb[0] = j+1;
        for (i=0;i<nth-1;i++) { mb[i] = rt;kb[i+1]=kb[i]+rt;}
        mb[nth-1]=q-(nth-1)*rt;
	//  #ifdef OPENMP_ON
        //#pragma omp parallel private(i) num_threads(nth)
        //#endif 
        { /* start of parallel section */
          #ifdef OPENMP_ON
          #pragma omp parallel for private(i) num_threads(nth)
          #endif
          for (i=0;i<nth;i++) {
            //#pragma flush(trans,m,mb,tau,k,A,kb,jb,n,Ak,one,dzero,F,pb)
	  F77_CALL(dgemv)(&trans, &m, mb+i,tau+k,A+(kb[i]+jb) * (ptrdiff_t) n+k,&n,
                          Ak,&one,&dzero,F+kb[i]+(ptrdiff_t)j*pb, &one);
          }
          //F77_CALL(dgemv)(&trans, &m, &i,tau+k,A+(k+1)*n+k,&n,Ak,&one,&dzero,F+j+1+j*pb, &one);
        } /* end of parallel section */
      } 
      /* F[0:pb-1,j] -= tau[k] F[0:pb-1,0:j-1] A[k:n-1,jb:k-1]'v */
      if (j>0) {
        q = j ; /* total number of rows to split between threads */
        rt = q/nt;if (rt*nt < q) rt++; /* rows per thread */
        nth = nt; while (nth>1&&(nth-1)*rt>q) nth--; /* reduce number of threads if some empty */
        kb[0] = jb; /* starting row */
        for (i=0;i<nth-1;i++) { mb[i] = rt;kb[i+1]=kb[i]+rt;}
        mb[nth-1]=q-(nth-1)*rt;
        #ifdef OPENMP_ON
        #pragma omp parallel private(i) num_threads(nth)
        #endif 
        { /* start of parallel section */
          #ifdef OPENMP_ON
          #pragma omp for
          #endif
          for (i=0;i<nth;i++) {
  	    F77_CALL(dgemv)(&trans, &m, mb+i,&dmone,A+kb[i]* (ptrdiff_t)n+k,&n,Ak,&one,&dzero,work+kb[i]-jb, &one);
          // F77_CALL(dgemv)(&trans, &m, &j,&dmone,A+jb*n+k,&n,Ak,&one,&dzero,work, &one);
          }
        }
        q = pb ; /* total number of rows to split between threads */
        rt = q/nt;if (rt*nt < q) rt++; /* rows per thread */
        nth = nt; while (nth>1&&(nth-1)*rt>q) nth--; /* reduce number of threads if some empty */
        kb[0] = 0; /* starting row */
        for (i=0;i<nth-1;i++) { mb[i] = rt;kb[i+1]=kb[i]+rt;}
        mb[nth-1]=q-(nth-1)*rt;
        for (i=0;i<nth;i++) { 
	  F77_CALL(dgemv)(&nottrans, mb+i, &j,tau+k,F+kb[i],&pb,work,&one,&done,F+(ptrdiff_t)j*pb+kb[i], &one);
	  // F77_CALL(dgemv)(&nottrans, &pb, &j,tau+k,F,&pb,work,&one,&done,F+j*pb, &one);
        }
      }
      /* update pivot row A[k,k+1:p-1] -= A[k,jb:k]F(j+1:pb-1,1:j)' */
      if (k<p-1) {
        q = pb-j-1 ; /* total number of cols to split between threads */
        rt = q/nt;if (rt*nt < q) rt++; /* colss per thread */
        nth = nt; while (nth>1&&(nth-1)*rt>q) nth--; /* reduce number of threads if some empty */
        kb[0] = j+1; /* starting col in F, and jb to this to get start in A */
        for (i=0;i<nth-1;i++) { mb[i] = rt;kb[i+1]=kb[i]+rt;}
        mb[nth-1]=q-(nth-1)*rt;
        q=j+1; 
        #ifdef OPENMP_ON
        #pragma omp parallel private(i) num_threads(nth)
        #endif 
        { /* start of parallel section */
          #ifdef OPENMP_ON
          #pragma omp for
          #endif
          for (i=0;i<nth;i++) {
	    F77_CALL(dgemv)(&nottrans, mb+i, &q,&dmone,F+kb[i],&pb,A+(ptrdiff_t)jb*n+k,
	                    &n,&done,A+(kb[i]+jb)*(ptrdiff_t)n+k, &n); 
          }
        }
        //m=pb-j-1;i=j+1;
        //F77_CALL(dgemv)(&nottrans, &m, &i,&dmone,F+j+1,&pb,A + jb * n + k,&n,&done,Ak+n, &n);
      }
      *Ak = xx; /* restore A[k,k] */
      /* Now down date the column norms */
      ok_norm = 1;
      if (k < p-1) {
        Ak += n;
        a0 = cn + k + 1;a1=cn + p;Aq = icn + k + 1;
        for (;a0<a1;a0++,Aq++,Ak+=n) { 
          *a0 -= *Ak * *Ak;
          if (*a0<*Aq*tol) ok_norm = 0; /* cancellation error problem in downdate */ 
        }  
      } /* down dates complete */
      if (!ok_norm) { 
        j++;
        nb = j;
        break; /* a column norm has suffered serious cancellation error. Need to break out and recompute */
      }
    } /* for j */
    j--; /* make compatible with current k */

    /* now the block update - about half the work is here*/    
    if (k<p-1) {
      /* A[k+1:n,k+1:p] -= A[k+1:n,jb:k]F[j+1:pb,0:nb-1]'
         dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)*/ 
      m = n - k - 1 ;
      rt = m/nt;if (rt*nt < m) rt++; /* rows per thread */
      nth = nt; while (nth>1&&(nth-1)*rt>m) nth--; /* reduce number of threads if some empty */
      kb[0] = k+1;
      for (i=0;i<nth-1;i++) { mb[i] = rt;kb[i+1]=kb[i]+rt;}
      mb[nth-1]=m-(nth-1)*rt;
      rt = p - k - 1;  
      //Rprintf("nth = %d  nt = %d\n",nth,nt);   
      #ifdef OPENMP_ON
      #pragma omp parallel private(i,Ak,Aq) num_threads(nth)
      #endif 
      { /* start of parallel section */
        #ifdef OPENMP_ON
        #pragma omp for
        #endif
        for (i=0;i<nth;i++) {
          Ak = A + (k+1)*(ptrdiff_t) n + kb[i];Aq = A + jb * (ptrdiff_t) n + kb[i];
          /* Argument 12 changed on exit... */
          F77_CALL(dgemm)(&nottrans,&trans,mb+i,&rt,&nb,&dmone,Aq,&n,F+j+1,&pb,&done,Ak,&n);
          // Ak = A + (k+1)*n + k + 1;Aq = A + jb * n + k + 1;
          // F77_CALL(dgemm)(&nottrans,&trans,&m,&rt,&nb,&dmone,Aq,&n,F+j+1,&pb,&done,Ak,&n);
        }
      } /* end of parallel section */
    }

    if (!ok_norm) { /* recompute any column norms that had cancelled badly */ 
      for (i = k+1;i<p; i++) { 
        if (cn[i]<icn[i]*tol) { /* do the recompute */
          x=0.0; Ak = A + i * (ptrdiff_t)n; 
          Aq = Ak + n; // Aq += n; <-- this was old code: Aq not initialised!
          Ak += k+1;
          for (;Ak<Aq;Ak++) x += *Ak * *Ak;
          cn[i] = icn[i] = x;
        } 
      } 
    }
    pb -= nb;
    jb += nb;
  } /* end while (jb<p) */
  FREE(F); FREE(mb); FREE(kb);
  FREE(cn);
  FREE(icn);
  FREE(work);
#ifdef OMP_REPORT
  Rprintf("done\n");
#endif
  return(p); /* NOTE: not really rank!! */
} /* bpqr */

int mgcv_piqr(double *x,int n, int p, double *beta, int *piv, int nt) {
/* Do the work for parallel QR: 
   Key is to inline the expensive step that is parallelized, noting
   that householder ops are not in BLAS, but are LAPACK auxilliary 
   routines (so no loss in re-writing)
   
   Parallelization here is based on splitting the application of householder
   rotations to the 'remaining columns' between cores using openMP. 
*/

  int i,k,r,nh,j,one=1,cpt,nth,cpf,ii;
  double *c,*p0,*p1,*p2,xx,tau,*work,zz,*v,*z,*z1,br;
  /* const char side = 'L';*/

  #ifndef OPENMP_ON
  nt = 1;
  #endif
  #ifdef OMP_REPORT
  Rprintf("mgcv_piqr...");
  #endif
  c =(double *)CALLOC((size_t) p,sizeof(double)); 
  work =(double *)CALLOC((size_t) (p*nt),sizeof(double));
  k=0;tau=0.0;
  /* find column norms O(np) */
  for (p0=x,i=0;i<p;i++) {
    piv[i] = i;
    for (xx=0.0,p1=p0 + n;p0<p1;p0++) xx += *p0 * *p0;
    c[i] = xx;
    if (xx>tau) {tau=xx;k=i;}
  }
  r = -1;
  nh = n; /* householder length */
  
  while (tau > 0) {
    r++;
    i=piv[r]; piv[r] = piv[k];piv[k] = i;
    /* swap r with k O(n) */
    xx = c[r];c[r] = c[k];c[k] = xx;
    for (p0 = x + n * r, p1 = x + (ptrdiff_t)n * k,p2 = p0 + n;p0<p2;p0++,p1++) {
          xx = *p0; *p0 = *p1; *p1 = xx;
    }
    /* now generate the householder reflector for column r O(n)*/
    p0 = x + r * n + r; /* first element of column */
    p1 = p0 + 1; /* remaining elements of column */
    xx = *p0; /* contains first element of column to be worked on */
    F77_CALL(dlarfg)(&nh,&xx,p1,&one,beta+r);
    /* now xx contains first element of rotated column - i.e. leading 
       diagonal element of R[r,r]. Elements of column r of x below the 
       diagonal now contain elements of the housholder vector apart from 
       the first, which is 1. So if v' = [1,x[(r+1):n,r]'] then 
       H = I - beta[r]*v*v' */

    /* next apply the rotation to the remaining columns of x O(np) */
       
    *p0 = 1.0; /* put 1 into leading diagonal element of x so that v = x[r:n,r] */
    j=p-r-1;
    
    /* now distribute the j columns between nt threads */
    if (j) {
          cpt = j / nt; /* cols per thread */
          if (cpt * nt < j) cpt++; 
          nth = j/cpt; /* actual number of threads */
          if (nth * cpt < j) nth++;
          cpf = j - cpt * (nth-1); /* columns on final block */
    } else nth=cpf=cpt=0;

    br = beta[r];
    //#ifdef OPENMP_ON
    //#pragma omp parallel private(i,j,p1,v,z,z1,zz,ii) num_threads(nt)
    //#endif
    { j = cpt;
      if (j) {        
        #ifdef OPENMP_ON
        #pragma omp parallel for private(i,j,p1,v,z,z1,zz,ii) num_threads(nt)
        #endif
        for (i=0;i<nth;i++) {
	    if (i == nth-1) j = cpf; else j = cpt;
            p1 = p0 + n + n * cpt *i; 
            z1=p1+nh;
            for (ii =0;ii<j;ii++,p1+=n,z1+=n) {
            /* apply reflection I - beta v v' */ 
            for (zz=0.0,v=p0,z=p1;z<z1;z++,v++) zz += *z * *v * br;
            for (z=p1,v=p0;z<z1;z++,v++) *z -= zz * *v;
	   }
	   /*F77_CALL(dlarfx)(&side, &nh, &j, p0, beta+r, p0 + n + n * cpt * i , &n , work + p * i);*/
        }
      } /* if (j) */
    } /* end parallel */ 
    nh--;
    *p0 = xx; /* now restore leading diagonal element of R to x[r,r] */
    /* update c, get new k... */
    k = r + 1;
    for (tau=0.0,p0+=n,i=r+1;i<p;i++,p0+=n) { 
          c[i] -= *p0 * *p0;
          if (c[i]>tau) { 
            tau = c[i];
            k=i;
          }
     }
     if (r==n-1) tau = 0.0;
    
    } /* end while (tau > 0) */
  
  FREE(c); FREE(work); 
  #ifdef OMP_REPORT
  Rprintf("done\n");
  #endif
  return(r+1);
} /* mgcv_piqr */

SEXP mgcv_Rpiqr(SEXP X, SEXP BETA,SEXP PIV,SEXP NT, SEXP NB) {
/* routine to QR decompose N by P matrix X with pivoting.
   Work is done by bpqr.

   Designed for use with .call rather than .C
   Return object is as 'qr' in R.

*/
  int n,p,nt,*piv,r,*rrp,nb;
  double *x,*beta;
  SEXP rr;
  nt = asInteger(NT);nb = asInteger(NB);
  n = nrows(X);
  p = ncols(X);
  x = REAL(X);beta = REAL(BETA);
  piv = INTEGER(PIV);
  r = bpqr(x,n,p,beta,piv,nb,nt); /* block version */
  /* should return rank (r+1) */
  rr = PROTECT(allocVector(INTSXP, 1));
  rrp = INTEGER(rr);
  *rrp = r;
 
  UNPROTECT(1);
  return(rr);

} /* mgcv_piqr */

void mgcv_pmmult(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n,int *nt) {
  /* 
     Forms r by c product, A, of B and C, transposing each according to bt and ct.
     n is the common dimension of the two matrices, which are stored in R 
     default column order form.
   
     This version uses openMP parallelization. nt is number of threads to use. 
     The strategy is rather simple, and this routine is really only useful when 
     B and C have numbers of rows and columns somewhat higher than the number of 
     threads. Assumes number of threads already set on entry and nt reset to 1
     if no openMP support.

     BLAS version A is c (result), B is a, C is b, bt is transa ct is transb 
     r is m, c is n, n is k.
     Does nothing if r,c or n <= zero. 
  */
  char transa='N',transb='N';
  int lda,ldb,ldc,cpt,cpf,c1,i,nth;
  double alpha=1.0,beta=0.0;
  if (*r<=0||*c<=0||*n<=0) return;
  #ifdef OMP_REPORT
  Rprintf("mgcv_pmmult...");
  #endif
  if (B==C) { /* this is serial, unfortunately. note case must be caught as B can be re-ordered!  */
    if (*bt&&(!*ct)&&(*r==*c)) { getXtX(A,B,n,r);return;} 
    else if (*ct&&(!*bt)&&(*r==*c)) { getXXt(A,B,c,n);return;}
  }
  #ifndef OPENMP_ON
  *nt = 1;
  #endif
  if (*nt == 1) {
    mgcv_mmult(A,B,C,bt,ct,r,c,n); /* use single thread version */
    return;
  }
  if (*bt) { /* so B is n by r */
    transa = 'T';  
    lda = *n;
  } else  lda = *r; /* B is r by n */

  if (*ct) { /* C is c by n */ 
    transb = 'T';
    ldb = *c;
  } else ldb = *n; /* C is n by c */

  ldc = *r;

  if (*ct) { /* have to split on B, which involves re-ordering */
    if (*bt) { /* B'C': can split on columns of n by r matrix B, but (r by c) A then needs re-ordering */
      cpt = *r / *nt; /* cols per thread */
      if (cpt * *nt < *r) cpt++; 
      nth = *r/cpt;
      if (nth * cpt < *r) nth++;
      cpf = *r - cpt * (nth-1); /* columns on final block */
      #ifdef OPENMP_ON
      #pragma omp parallel private(i,c1) num_threads(nth)
      #endif
      { /* open parallel section */
        //c1 = cpt;
        #ifdef OPENMP_ON
        #pragma omp for
        #endif
        for (i=0;i<nth;i++) {
          if (i == nth-1) c1 = cpf; else c1 = cpt;
          /* note integer after each matrix is its leading dimension */
          if (c1>0) F77_CALL(dgemm)(&transa,&transb,&c1,c,n, &alpha,
				    B + i * (ptrdiff_t) cpt * *n, n ,C, c,&beta, A + i * (ptrdiff_t) cpt * *c, &c1);
        }
      } /* parallel section ends */
      /* now re-order the r by c matrix A, which currently contains the sequential
         blocks corresponding to each cpt rows of A */
      row_block_reorder(A,r,c,&cpt,bt); /* bt used here for 'reverse' as it contains a 1 */
    } else { /* BC':worst case - have to re-order r by n mat B and then reverse re-ordering of B and A at end */
      cpt = *r / *nt; /* cols per thread */
      if (cpt * *nt < *r) cpt++; 
      nth = *r/cpt;
      if (nth * cpt < *r) nth++;
      cpf = *r - cpt * (nth-1); /* columns on final block */
      /* re-order cpt-row blocks of B into sequential cpt by n matrices (in B) */ 
      row_block_reorder(B,r,n,&cpt,bt); /* bt contains a zero - forward mode here */
      #ifdef OPENMP_ON
      #pragma omp parallel private(i,c1) num_threads(nth)
      #endif
      { /* open parallel section */
        //c1 = cpt;
        #ifdef OPENMP_ON
        #pragma omp for
        #endif
        for (i=0;i<nth;i++) {
          if (i == nth-1) c1 = cpf;else c1=cpt;
          if (c1>0) F77_CALL(dgemm)(&transa,&transb,&c1,c,n, &alpha,
		B + i * (ptrdiff_t) cpt * *n, &c1,C,c,&beta, A + i * (ptrdiff_t) cpt * *c, &c1);
        }
      } /* parallel ends */
      /* now reverse the re-ordering */
      row_block_reorder(B,r,n,&cpt,ct);
      row_block_reorder(A,r,c,&cpt,ct);
    }
  } else { /* can split on columns of n by c matrix C, which avoids re-ordering */
    cpt = *c / *nt; /* cols per thread */
    if (cpt * *nt < *c) cpt++;
    nth = *c/cpt;
    if (nth * cpt < *c) nth++;
    cpf = *c - cpt * (nth-1); /* columns on final block */
    #ifdef OPENMP_ON
    #pragma omp parallel private(i,c1) num_threads(*nt)
    #endif
    { /* open parallel section */
      //c1 = cpt;
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (i=0;i< nth;i++) {
        if (i == nth-1) c1 = cpf;else c1=cpt; /* how many columns in this block */
        if (c1>0) F77_CALL(dgemm)(&transa,&transb,r,&c1,n, &alpha,
		B, &lda,C + i * (ptrdiff_t) *n * cpt, &ldb,&beta, A + i * (ptrdiff_t) *r * cpt, &ldc);
      }
    } /* end parallel */
  }
  #ifdef OMP_REPORT
  Rprintf("done\n");
  #endif
} /* end mgcv_pmmult */

void pcrossprod(double *B, double *A,int *R, int *C,int *nt,int *nb) {
/* B=A'A if t==0. A is R by C. nb^2 is the target number of elements in a block. 
   nt is the number of threads to use. B is C by C.
   30/4 memorial edition 
*/
  int M,N,nf,nrf,kmax,kk,i,r,c,k,bn,an,cn; 
  ptrdiff_t as,bs,cs;
  char uplo = 'U',trans='T',ntrans='N'; 
  double alpha=1.0,beta=1.0;
  M = ceil(((double) *C)/ *nb);
  N = ceil(((double) *R)/ *nb);
  if (M==1) { /* perform single threaded crossprod */
    beta = 0.0;
    F77_CALL(dsyrk)(&uplo,&trans,C,R,&alpha,A,R,&beta,B,C);
  } else {
    nf = *C - (M-1) * *nb;  /* cols in last col block of A */
    nrf = *R - (N-1) * *nb;  /* rows in last row block of A */
    kmax = (M+1)*M/2; /* number of blocks in upper triangle */ 
    #ifdef OPENMP_ON
    #pragma omp parallel for private(kk,i,r,c,bn,bs,k,as,an,beta,cs,cn) num_threads(*nt)
    #endif
    for (kk=0;kk<kmax;kk++) { /* work through blocks on upper triangle */
      i=kk;r=0;while (i >= M-r) { i -= M - r; r++;}; c = r + i; /* convert kk to row/col */
      if (r==M-1) bn = nf; else bn = *nb; /* (row) B block size */
      bs = r * (ptrdiff_t) *nb; /* (row) B block start */
      if (c==r) { /* diagonal block */
        for (k=0;k<N;k++) { /* work through row blocks of A */
          as = k * (ptrdiff_t) *nb; /* A row block start */ 
          if (k==N-1) an = nrf; else an = *nb; /* A row block size */ 
          /* uplo 'U' or 'L' for upper/lower tri. trans = 'T' for A'A.*/
          if (k) beta=1.0; else beta = 0.0; /* multiplier for B block */
          F77_CALL(dsyrk)(&uplo,&trans,&bn, /* cols of A block */
			&an, /* rows of A block */ 
                        &alpha,A + as + (ptrdiff_t) *R * bs, /* start of A block */
                        R, /* R is physical number of rows in A */
                        &beta,B + bs + (ptrdiff_t) *C * bs,C);
        }
      } else {
	cs = c * (ptrdiff_t) *nb; /* b col start */
        if (c==M-1) cn = nf; else cn = *nb; /* b col block size */
        for (k=0;k<N;k++) { /* work through row blocks of A */
          as = k * (ptrdiff_t) *nb; /* A row block start */ 
          if (k==N-1) an = nrf; else an = *nb; /* A row block size */ 
          /* uplo 'U' or 'L' for upper/lower tri. trans = 'T' for A'A.*/
          if (k) beta=1.0; else beta = 0.0; /* multiplier for B block */
          F77_CALL(dgemm)(&trans,&ntrans,
                          &bn,&cn,&an, 
                          &alpha,A + *R * bs + as, 
                          R,A + *R * cs + as,R,&beta, B + bs + *C * cs, C);
        }
      }
    }
  }
  /* fill in lower block */
  for (r=0; r < *C;r++) for (c=0;c < r;c++) B[r + *C * c] = B[c + *C * r];
} /* pcrossprod */


SEXP mgcv_Rpcross(SEXP A, SEXP NT,SEXP NB) {
/* .Call wrapper for pcrossproduct - forms A'A using 
   NT threads and blocksize NB.
   Return object a matrix.
*/
  int c,r,nt,nb;
  double *A0,*B0;
  SEXP B;
  nt = asInteger(NT); /* number of threads */
  nb = asInteger(NB); /* block size */
  r = nrows(A);c = ncols(A);
  A0 = REAL(A);  
  B = PROTECT(allocMatrix(REALSXP,c,c));
  B0 = REAL(B);
  
  pcrossprod(B0,A0,&r,&c,&nt,&nb);
   
  UNPROTECT(1);
  return(B);
} /* mgcv_Rpcross */


void getXtX0(double *XtX,double *X,int *r,int *c)
/* form X'X (nearly) as efficiently as possible - BLAS free
   thread safety: modifies XtX. Nothing else.
*/
{ double *p0,*p1,*p2,*p3,*p4,x;
  int i,j;
  for (p0=X,i=0;i<*c;i++,p0 += *r) 
  for (p1=X,j=0;j<=i;j++,p1 += *r) {
    for (x=0.0,p2=p0,p3=p1,p4=p0 + *r;p2<p4;p2++,p3++) x += *p2 * *p3;    
    XtX[i + j * *c] = XtX[j + i * *c] = x;
  }
}

void getXtX(double *XtX,double *X,int *r,int *c)
/* form X'X (nearly) as efficiently as possible - uses BLAS
   Thread safety: Only XtX changed on exit.
*/
{ double alpha=1.0,beta=0.0;
  int i,j;
  char uplo = 'L',trans='T';
  F77_CALL(dsyrk)(&uplo,&trans,c, r, &alpha,X,r,&beta,XtX,c);
  /* fill in upper triangle from lower */
  for (i=0;i<*c;i++) 
  for (j=0;j<i;j++)  XtX[j + i * *c] = XtX[i + j * *c];
 
}

void getXXt(double *XXt,double *X,int *r,int *c)
/* form XX' (nearly) as efficiently as possible - uses BLAS*/
{ double alpha=1.0,beta=0.0;
  ptrdiff_t ii,jj;
  char uplo = 'L',trans='N';
  F77_CALL(dsyrk)(&uplo,&trans,r, c, &alpha,X,r,&beta,XXt,r);
  /* fill in upper triangle from lower */
  for (ii=0;ii<*r;ii++) 
  for (jj=0;jj<ii;jj++)  XXt[jj + ii * *r] = XXt[ii + jj * *r];
 
}


void getXtWX0(double *XtWX, double *X,double *w,int *r,int *c,double *work)
/* Original BLAS free version...
   forms X'WX as efficiently as possible, where W = diag(w)
   and X is an r by c matrix stored column wise. 
   work should be an r-vector (longer is no problem).
   Thread safety: work and XtWX modified; rest not.
*/ 
{ int i,j;
  double *p,*p1,*p2,*pX0,*pX1,xx;
  pX0=X;
  for (i=0;i< *c;i++) { 
    p2 = work + *r;
    for (p=w,p1=work;p1<p2;p++,p1++,pX0++) *p1 = *pX0 * *p; 
    for (pX1=X,j=0;j<=i;j++) {
      for (xx=0.0,p=work;p<p2;p++,pX1++) xx += *p * *pX1;
      XtWX[i * *c + j] = XtWX[j * *c + i] = xx;
    }
  }
}

void getXtWX(double *XtWX, double *X,double *w,int *r,int *c,double *work)
/* BLAS version...
   forms X'WX as efficiently as possible, where W = diag(w)
   and X is an r by c matrix stored column wise. 
   work should be an r-vector (longer is no problem).
   Thread safety: work and XtWX modified; rest unmodified.
*/ 
{ int i,j,one=1;
  char trans='T';
  double *p,*p1,*p2,*pX0,xx=0.0,*w2,alpha=1.0,beta=0.0;
  pX0=X;
  w2 = XtWX; /* use first column as work space until end */
  for (i=0;i< *c;i++) { 
    p2 = work + *r;
    /* form work=WX[,i]... */
    for (p=w,p1=work;p1<p2;p++,p1++,pX0++) *p1 = *pX0 * *p;  
    /* Now form X[,1:i]'work ... */
    j = i+1; /* number of columns of X to use */
    F77_CALL(dgemv)(&trans, r, &j,&alpha,X, r,work,&one,&beta,w2, &one);
    if (i==0) xx = w2[0]; /* save the 0,0 element of XtWX (since its in use as workspace) */
    else for (j=0;j<=i;j++) XtWX[i * *c + j] = w2[j];
  }
  /* now fill in the other triangle */
  if (*r * *c > 0) *XtWX = xx;
  for (i=0;i< *c;i++) for (j=0;j<i;j++) XtWX[j * *c + i] =  XtWX[i * *c + j];
}


void getXtMX(double *XtMX,double *X,double *M,int *r,int *c,double *work)
/* BLAS free version (currently no BLAS version as this is only used on 
   small matrices).
   forms X'MX as efficiently as possible, where M is a symmetric matrix
   and X is an r by c matrix. X and M are stored column wise. 
   work should be an r-vector (longer is no problem).
*/

{ int i,j;
  double *p,*p1,*p2,*pX0,*pX1,xx,*pM;
  pX0=X;
  for (i=0;i< *c;i++) { 
    /* first form MX[:,i] */
    p2 = work + *r;pM=M;
    for (p1=work;p1<p2;pM++,p1++) *p1 = *pX0 * *pM;
    pX0++;
    for (j=1;j< *r;j++,pX0++) 
    for (p1=work;p1<p2;pM++,p1++) *p1 += *pX0 * *pM;
    /* now form ith row and column of X'MX */
    for (pX1=X,j=0;j<=i;j++) {
      for (xx=0.0,p=work;p<p2;p++,pX1++) xx += *p * *pX1;
      XtMX[i * *c + j] = XtMX[j * *c + i] = xx;
    }
  }
} /* getXtMX */

void vcorr(double *dR,double *Vr,double *Vb,int *p,int *M) {
/* dR contains M upper p by p matrices. The dR matrices are upper
   triangular, but the zeros are included. Vr is an M by M matrix.
   Vb is a p by p matrix. The kth dR contains dR/d\rho_k the derivative 
   of the Choleski factor of the covariance matrix of some parameters
   \beta w.r.t. \rho_k. Vr is the covariance matrix of \rho. If
   b = \sum_k dR'/d\rho_k z (\rho_k - \hat \rho_k) where z ~ N(0,I_p)
   then Vb is its cov matrix 
   If *M < 0 then R is not transposed, i.e. 
   b = \sum_k dR/d\rho_k z (\rho_k - \hat \rho_k)
*/ 
  double *Vi,*ViV,*p0,*p1,*p2,zero=0.0,one=1.0,x;
  int i,j,k;
  char trans = 'N';
  k = *p * *M; if (k<0) k = -k;
  Vi = (double *)CALLOC((size_t) k,sizeof(double));
  ViV = (double *)CALLOC((size_t) k,sizeof(double));
  if (*M>0) { /* R is transposed */
    for (i=0;i<*p;i++) {   
      for (p0=Vi,k=0;k<*M;k++) {
        /* Vi is i by M */
        p1 = dR + k * *p * *p + i * *p; /* start of col i of kth dR */
        p2 = p1 + i + 1; /* first zero in col i of kth dR */ 
        for (;p1<p2;p1++,p0++) *p0 = *p1;
      }    
      /* create ViV = Vi Vr (i by M) */
      k = i + 1;
      F77_CALL(dgemm)(&trans,&trans,&k,M,M, &one,
  		Vi, &k, Vr , M, &zero, ViV, &k);
      for (j=i;j<*p;j++) {
        for (x=0.0,p0=ViV,k=0;k<*M;k++) {
          /* Vj is i by M */
          p1 = dR + k * *p * *p + j * *p; /* start of col j of kth dR */
          p2 = p1 + i + 1; /* one past last relevant row in col j of kth dR */ 
          for (;p1<p2;p1++,p0++) x += *p0 * *p1;
        }
        Vb[i + *p * j] = Vb[j + *p * i] = x;
      }
    }
  } else { /* R not transposed */
    *M = - *M;
    for (i=0;i<*p;i++) {   
      for (p0=Vi,k=0;k<*M;k++) {
        /* Vi is p-i by M */
        p1 = dR + k * *p * *p + i ; /* start of row i of kth dR */
        p2 = p1 + *p * *p; /* 1 past row end */
        p1 += i * *p; /* skip the starting zeroes */ 
        for (;p1<p2;p1+= *p,p0++) *p0 = *p1;
      }    
      /* create ViV = Vi Vr (p-i by M) */
      k = *p - i;
      F77_CALL(dgemm)(&trans,&trans,&k,M,M, &one,
  		Vi, &k, Vr , M, &zero, ViV, &k);
      for (j=i;j<*p;j++) {
        for (x=0.0,p0=ViV,k=0;k<*M;k++) {
          p0 += (j-i);
          /* Vj is p-i by M */
          p1 = dR + k * *p * *p + j; /* start of row j of kth dR */
          p2 = p1 + *p * *p; /* one past last relevant row in col j of kth dR */ 
          p1 += j * *p; /* skip the starting zeroes */
          for (;p1<p2;p1 += *p,p0++) x += *p0 * *p1;
        }
        Vb[i + *p * j] = Vb[j + *p * i] = x;
      }
    }
  }
  FREE(Vi);FREE(ViV);
} /* vcorr */

void dchol(double *dA, double *R, double *dR,int *p) {
/* R is the Choleski factor of A, s.t. R'R = A (no pivoting). dA is dA/dx.
   This function computes dR/dx and returns it in dR. Lower triangles not accessed.
   FLOP cost double choleski, but no square roots.
*/
  int i,j,k;
  double x,*pdRi,*pdRj,*pRi,*pRj,*pend;
  for (i=0;i < *p;i++) for (j=i;j < *p;j++) {
    k = i * *p;pdRi = dR + k; pRi = R + k;
    k = j * *p; pdRj = dR + k; pRj = R + k;pend = pRj + i; 
    for (x=0.0;pRj<pend;pRi++,pRj++,pdRi++,pdRj++) x += *pRi * *pdRj + *pRj * *pdRi;
    k = i + j * *p;
    if (j>i) dR[k] = (dA[k] - x - R[k]*dR[i + i * *p])/R[i + i * *p]; 
    else  dR[k] = (dA[k] - x)*.5/R[i + i * *p]; 
  }
} /* dchol */



void mgcv_chol_down(SEXP r,SEXP ru,SEXP N,SEXP K, SEXP UT) {
/* wrapper for calling chol_down using .Call */
  double *R,*Rup;
  int *n,*k,*ut;
  R = REAL(r);Rup=REAL(ru);
  n = INTEGER(N);
  k = INTEGER(K);
  ut = INTEGER(UT);
  chol_down(R,Rup,n,k,ut);
}

inline double hypot(double x, double y) {
/* stable computation of sqrt(x^2 + y^2) */  
  double t;
  x = fabs(x);y=fabs(y);
  if (y>x) { t = x;x = y; y = t;}
  if (x==0) return(y); else t = y/x;
  return(x*sqrt(1+t*t));
} /* hypot */   

void chol_down(double *R,double *Rup, int *n,int *k,int *ut) {
/* R is an n by n choleski factor of an n by n matrix A. We want the downdated
   factor for A[-k,-k] returned in n-1 by n-1 matrix Rup. 
   If ut!=0 then R'R = A, with R and Rup upper triangular. 
   Otherwise RR'=A, with R and Rup lower triangular. The latter update 
   is more Cache friendly with less effort, since the update is then from 
   the right and operates columnwise. However, code below is column oriented in both 
   cases, storing the givens rotations as they are computed to facilitate
   column orientation when R is upper triangular.
   Calls from R should ideally be made from a wrapper called from .Call, 
  since otherwise copying can be the dominant cost. 

*/
  int i,j,n1;
  double x,*Ri1,*Ri,*Rj,c,s,*Re,*ca,*sa,*sp,*cp;
  n1 = *n-1;
  if (*ut) { /* upper trianglar col oriented computation */
      ca = R + 2;sa = ca + *n; /* Givens storage */
      for (i=0;i<n1;i++) { /* work across columns */
	Ri = Rup + i * n1; /* destination column */
	if (i < *k) {
	  Rj = R + i * *n;
	  Re = Rj + i;
	} else {
	  Rj = R + (i+1) * *n; /* source column */
	  Re = Rj + *k;
        }
	for (;Rj<=Re;Ri++,Rj++) *Ri=*Rj; /* copy column starts */
	/* now the updates */
	if (i >= *k) { /* cols from beyond dropped col k need updating */
          /* first the stored rotations */
          Re = Rup + i * n1 + i;
	  Ri1=Ri;Ri--;
	  for (cp=ca,sp=sa;Ri<Re;cp++,sp++,Ri++,Ri1++,Rj++) {
            *Ri1 = - *sp * *Ri  + *cp * *Rj;
            *Ri = *cp * *Ri + *sp * *Rj;
	  }  
	  /* Now the final rotation */
	  x = sqrt(*Ri * *Ri + *Rj * *Rj);
          c = *Ri / x; s = *Rj / x;*Ri = x;
	  if (i<n1-1) { *cp = c;*sp=s;} /* store rotation */
	}  
      }
      /* now clear Givens storage */
      Re = R + *n;
      for (;ca<Re;ca++,sa++) *ca = *sa = 0.0;
  } else { /* lower triangular */
    /* First copy the input rows before row k */
    for (i=0;i< *k;i++) {
      for (Ri=Rup+i*n1,Rj=R+i * *n,Re=Ri+*k;Ri<Re;Ri++,Rj++) *Ri = *Rj;
    }
    /* now copy cols 0:k from row k+1 */
    for (i=0;i<=*k;i++) {
      for (Ri=Rup+i*n1,Re=Ri+n1,Ri+= *k,Rj=R+i * *n+ *k+1;Ri<Re;Ri++,Rj++) *Ri = *Rj;
    }
    for (i = *k;i<n1;i++) { /* work down diagonal */
      Ri = Rup + i  + i * n1; /* rotate into here */
      Rj = R + (i+1) + (i+1) * *n;  /* zeroing this */
      /* x = sqrt(*Ri * *Ri + *Rj * *Rj);*/
      x = hypot(*Ri,*Rj);
      c = *Ri / x; s = *Rj / x;*Ri = x;
      Rj++;Ri++;Ri1=Ri+n1;Re = Rup + (i+1) * n1;
      for (;Ri<Re;Rj++,Ri++,Ri1++) {
	*Ri1 = -s * *Ri  + c * *Rj;
        *Ri = c * *Ri + s * *Rj;
      }     
    }  
  }  
 
} /* chol_down */  

void chol_up(double *R,double *u, int *n,int *up,double *eps) {
/* Rank 1 update of a cholesky factor. Works as follows:
 
   [up=1] R'R + uu' = [u,R'][u,R']' =  [u,R']Q'Q[u,R']', and then uses Givens rotations to
   construct Q such that Q[u,R']' = [0,R1']'. Hence R1'R1 = R'R + uu'. The construction 
   operates from first column to last.

   [up=0] uses an almost identical sequence, but employs hyperbolic rotations
   in place of Givens. See Golub and van Loan (2013, 4e 6.5.4)

   Givens rotations are of form [c,-s] where c = cos(theta), s = sin(theta). 
                                [s,c]

   Currently does not check that result of down date is positive definite. 
 
   Assumes R upper triangular, and that it is OK to use first two columns 
   below diagonal as temporary strorage for Givens rotations (the storage is 
   needed to ensure algorithm is column oriented). 
*/ 
  double c0,s0,*c,*s,z,*x,z0,*c1;
  int j,j1,n1;
  n1 = *n - 1;
  if (*up) for (j1=-1,j=0;j<*n;j++,u++,j1++) { /* loop over columns of R */
    z = *u; /* initial element of u */
    x = R + *n * j; /* current column */
    c = R + 2;s = R + *n + 2; /* Storage for first n-2 Givens rotations */ 
    for (c1=c+j1;c<c1;c++,s++,x++) { /* apply previous Givens */
      z0 = z;
      z = *c * z - *s * *x;
      *x = *s * z0 + *c * *x;
    }
    if (j) {
      /* apply last computed Givens */
      z0 = z;
      z = c0 * z - s0 * *x;
      *x = s0 * z0 + c0 * *x;
      x++;
      if (j<n1) {*c = c0; *s = s0;} /* store if needed for further columns */
    }  
    
    /* now construct the next rotation u[j] <-> R[j,j] */
    z0 = hypot(z,*x); /* sqrt(z^2+R[j,j]^2) */
    c0 = *x/z0; s0 = z/z0;  /* need to zero z */
    /* now apply this rotation and this column is finished (so no need to update z) */
    *x = s0 * z + c0 * *x;
  } else  for (j1=-1,j=0;j<*n;j++,u++,j1++) { /* loop over columns of R for down-dating */
    z = *u; /* initial element of u */
    x = R + *n * j; /* current column */
    c = R + 2;s = R + *n + 2; /* Storage for first n-2 hyperbolic rotations */ 
    for (c1=c+j1;c<c1;c++,s++,x++) { /* apply previous hyperbolic */
      z0 = z;
      z =  *c * z - *s * *x;
      *x = -*s * z0 + *c * *x;
    }
    if (j) {
      /* apply last computed hyperbolic */
      z0 = z;
      z = c0 * z - s0 * *x;
      *x = -s0 * z0 + c0 * *x;
      x++;
      if (j<n1) {*c = c0; *s = s0;} /* store if needed for further columns */
    }  
    
    /* now construct the next hyperbolic rotation u[j] <-> R[j,j] */
    z0 = z / *x; /* sqrt(z^2+R[j,j]^2) */
    if (fabs(z0)>=1) { /* downdate not +ve def */
      //Rprintf("j = %d  d = %g ",j,z0);
      if (*n>1) R[1] = -2.0;return; /* signals error */
    }
    if (z0 > 1 - *eps) z0 = 1 - *eps;
    c0 = 1/sqrt(1-z0*z0);s0 = c0 * z0;
    /* now apply this rotation and this column is finished (so no need to update z) */
    *x = -s0 * z + c0 * *x;
  }

 /* now zero c and s storage */
  c = R + 2;s = R + *n + 2;
  for (x = c + *n - 2;c<x;c++,s++) *c = *s = 0.0;
} /* chol_up */  

void mgcv_chol_up(SEXP r,SEXP U,SEXP N,SEXP UP,SEXP EPS) {
/* wrapper for calling chol_down using .Call */
  double *R,*u,*eps;
  int *n,*up;
  R = REAL(r);u=REAL(U);
  eps = REAL(EPS);
  n = INTEGER(N);
  up = INTEGER(UP);
  chol_up(R,u,n,up,eps);
} /* mgcv_chol_up */

void mgcv_chol(double *a,int *pivot,int *n,int *rank)
/* a stored in column order, this routine finds the pivoted choleski decomposition of matrix a 
   library(mgcv)
   X<-matrix(rnorm(16),4,4);D<-X%*%t(X)
   rank<-0
   er<-.C("mgcv_chol",as.double(D),as.integer(rep(0,4)),as.integer(4),as.integer(rank))
   rD<-matrix(er[[1]],4,4);piv<-er[[2]]
   chol(D,pivot=TRUE);rD;piv
   n<-length(piv);ind<-1:n;
   ind[piv]<-1:n
   rD<-rD[,ind]
   L<-mroot(D)
   D;t(rD)%*%rD;L%*%t(L)
*/
{ double *work,*p1,*p2,*p,*p3,tol=-1.0;
  int info=1;
  char uplo='U';
  work=(double *)CALLOC((size_t) *n * 2,sizeof(double));
  //F77_CALL(dchdc)(a,n,n,work,pivot,&info,rank); /* LINPACK */
  F77_CALL(dpstrf)(&uplo,n,a,n,pivot,rank,&tol,work,&info); /* LAPACK */
  /* zero stuff below the leading diagonal */
  for (p2=a + *n,p1=a+1,p3=a+ *n * *n;p2<p3;p1+= *n+1,p2+= *n) for (p=p1;p<p2;p++) *p=0.0;
  FREE(work);
}

void mroot(double *A,int *rank,int *n)
/* finds the minimum rank or supplied rank square root of n by n matrix  A by pivoted choleski 
   decomposition returned in A is B such that B'B = A (this is different from R routine mroot)

   R testing code....
   library(mgcv)
   X<-matrix(rnorm(12),4,3);D<-X%*%t(X)
   rank<-0
   er<-.C("mroot",as.double(D),as.integer(rank),as.integer(4))
   rD<-matrix(er[[1]],er[[2]],4)
   D;t(rD)%*%rD
   
*/
{ int *pivot,erank,i,j;
  double *pi,*pj,*p0,*p1,*B;
  pivot=(int *)CALLOC((size_t)*n,sizeof(int));
  mgcv_chol(A,pivot,n,&erank);
  if (*rank<=0) *rank = erank;
  /* unscramble the returned choleski factor */ 
  B=CALLOC((size_t)(*n * *n),sizeof(double));
  /* copy A to B and zero A */
  for (p0=A,p1=B,i=0;i< *n;i++,p0+= *n,p1+= *n)
  for (pi=p0,pj=p1;pi<=p0+i;pi++,pj++) { *pj = *pi; *pi=0.0;}
  /* now copy B into A undoing pivoting */
  for (p0=B,i=0;i< *n;p0+= *n,i++) 
  for (j=pivot[i],pj=p0,pi=A + (j-1)* *n;pj<=p0+i;pi++,pj++) *pi = *pj;
  /* remove trailing rows from choleski factor */
  for (pi=A,p0=A,i=0;i< *n;i++,p0+= *n)
  for (pj=p0;pj<p0+ *rank;pj++,pi++) *pi = *pj; 
  FREE(pivot);FREE(B);
}


void mgcv_svd(double *x,double *u,double *d,int *r,int *c)
/* call LA_PACK svd routine to form x=UDV'. Assumes that V' not wanted, but that full square U is required.
*/
{ const char jobu='A',jobvt='N';
  int lda,ldu,ldvt=1,lwork;
  int info;
  double *vt=NULL,work1,*work;
  ldu=lda= *r;
  lwork=-1;
  /* workspace query */
  F77_CALL(dgesvd)(&jobu,&jobvt, r, c, x, &lda, d, u, &ldu, vt,&ldvt,
  		   &work1, &lwork, &info);
  lwork=(int)floor(work1);
  if (work1-lwork>0.5) lwork++;
  work=(double *)CALLOC((size_t)lwork,sizeof(double));
  /* actual call */
  F77_CALL(dgesvd)(&jobu,&jobvt, r, c, x, &lda, d, u, &ldu, vt,&ldvt,
  		   work, &lwork, &info);
  FREE(work);
}

void mgcv_svd_full(double *x,double *vt,double *d,int *r,int *c)
/* call LA_PACK svd routine to form x=UDV'. U returned in x. V' returned in vt.
   assumed r >= c. U is r by c. D is length c. V is c by c.

# Here is R test code.....
library(mgcv)
n<-4;q<-3
X<-matrix(rnorm(n*q),n,q)
um<-.C("mgcv_svd_full",as.double(X),double(q*q),double(q),as.integer(n),as.integer(q),
        PACKAGE="mgcv")
er<-La.svd(X)
matrix(um[[1]],n,q);er$u
um[[3]];er$d
matrix(um[[2]],q,q);er$v
*/
{ const char jobu='O',jobvt='A';
  int lda,ldu,ldvt,lwork;
  int info;
  double work1,*work,*u=NULL;
  ldu=lda= *r;ldvt = *c;
  lwork=-1;
  /* workspace query */
  F77_CALL(dgesvd)(&jobu,&jobvt, r, c, x, &lda, d, u, &ldu, vt,&ldvt,
  		   &work1, &lwork, &info);
  lwork=(int)floor(work1);
  if (work1-lwork>0.5) lwork++;
  work=(double *)CALLOC((size_t)lwork,sizeof(double));
  /* actual call */
  F77_CALL(dgesvd)(&jobu,&jobvt, r, c, x, &lda, d, u, &ldu, vt,&ldvt,
  		   work, &lwork, &info);
  FREE(work);
}


void mgcv_td_qy(double *S,double *tau,int *m,int *n, double *B,int *left,int *transpose)
/* Multiplies m by n matrix B by orthogonal matrix returned from mgcv_tri_diag and stored in 
   S, tau. B is overwritten with result. 
    
   Note that this is a bit inefficient if really only a few rotations matter!

   Calls LAPACK routine dormtr 
*/
{ char trans='N',side='R',uplo='U';
  int nq,lwork=-1,info;
  double *work,work1;
  if (*left) { side = 'L';nq = *m;} else nq = *n;
  if (*transpose) trans = 'T';
  /* workspace query ... */
  F77_CALL(dormtr)(&side,&uplo,&trans,m,n,S,&nq,tau,B,m,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)CALLOC((size_t)lwork,sizeof(double));
  /* actual call ... */
  F77_CALL(dormtr)(&side,&uplo,&trans,m,n,S,&nq,tau,B,m,work,&lwork,&info);
  FREE(work);
}

void tri_chol(double *ld,double *sd,int *n,int *info) {
/* compute LDL' decomposition of n by n symm tridiagonal matrix
   with leading diagonal ld and sub/sup diagonals sd. 
   Returns D in ld and sub-diagonal of L in sd (leading diagonal 
   of L is all ones). info is 0 for success, -k if kth argument illegal
   and k<n if not +ve definite, k==n if positive semi-definite and could 
   just complete.
*/  
  F77_CALL(dpttrf)(n,ld,sd,info); 
}

void band_chol(double *B,int *n,int *k,int *info) {
/* k by n matrix B contains k non-zero diagonal of a +ve definite 2k-1 banded matrix
   The leading diagonal is stored in the first row and is of length n.
   Subsequent diagonal are zero padded at the end, to also be of length n. 
   Note the unusual `by row' storage.
   On exit the diagonals of the factor are returned in the same way.
*/  
  char uplo='L';
  int kd;
  kd = *k - 1;  
  F77_CALL(dpbtrf)(&uplo,n,&kd,B,k,info);
  // DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
}

void mgcv_tri_diag(double *S,int *n,double *tau)
/* Reduces symmetric n by n matrix S to tridiagonal form T 
   by similarity transformation. Q'SQ = T, where Q is an
   orthogonal matrix. Only the upper triangle of S is actually
   used.

   On exit the diagonal and superdiagonal of T are written in 
   the corresponding position in S. The elements above the first 
   superdiagonal, along with tau, store the householder reflections 
   making up Q. 

   Note that this is not optimally efficient if actually only a 
   few householder rotations are needed because S does not have 
   full rank.

   The routine calls dsytrd from LAPACK.
     
*/
{ int lwork=-1,info;
  double *work,work1,*e,*d;
  char uplo='U';
  d = (double *)CALLOC((size_t)*n,sizeof(double));
  e = (double *)CALLOC((size_t)(*n-1),sizeof(double));
  /* work space query... */
  F77_CALL(dsytrd)(&uplo,n,S,n,d,e,tau,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)CALLOC((size_t)lwork,sizeof(double));
  /* Actual call... */
  F77_CALL(dsytrd)(&uplo,n,S,n,d,e,tau,work,&lwork,&info);
  FREE(work);FREE(d);FREE(e);
}


void mgcv_pbsi0(double *R, int *n,int *nt) {
/* Form inverse of upper triangular R (n by n) in situ.
   Block oriented, using BLAS/LAPACK calls. 
*/
  int j,info,jb,nb=50;
  char left = 'L',right = 'R',up='U',ntrans='N',diag='N';
  double d1 = 1.0,m1= -1.0;
  for (j=0;j< *n;j+=nb) {
    jb = *n - j; if (jb>nb) jb=nb; /* block size */
    if (j) {
      /* A[0:j-1,j:j+jb-1] <- A[0:j-1,0:j-1] A[0:j-1,j:j+jb-1] O(j^2 jb)*/
      F77_CALL(dtrmm)(&left,&up,&ntrans,&diag,&j,&jb,&d1,R,n,R + j * *n,n);
      /* A[0:j-1,j:j+jb-1] <- A[0:j-1,j:j+jb-1] A[j:j+jb-1,j:j+jb-1]^{-1} O(j jb^2) */
      F77_CALL(dtrsm)(&right,&up,&ntrans,&diag,&j,&jb,&m1,R + *n * j + j,n,R + j * *n,n);
    }
    /* invert A[j:j+jb-1,j:j+jb-1] O(jb^3)*/
    F77_CALL(dtrti2)(&up,&diag,&jb,R + *n * j + j,n,&info); 
  }  

} /* mgcv_pbsi0 */

void mgcv_pbsi1(double *R, int *n,int *nt) {
/* Form inverse of upper triangular R (n by n) in situ.
   Block oriented, using BLAS/LAPACK calls. 
*/
  int j,info,jb,nb=50,*iwork;
  char left = 'L',right = 'R',up='U',ntrans='N',diag='N';
  double d1 = 1.0,m1= -1.0,*work;
  iwork = (int *)CALLOC((size_t)3 * (*nt *(*nt+1))/2 + 2 * *nt + 2,sizeof(int));
  work = (double *)CALLOC((size_t) nb * ((ptrdiff_t) *n + *nt + 1) * (*nt + 1)/2,sizeof(double));
  for (j=0;j< *n;j+=nb) {
    jb = *n - j; if (jb>nb) jb=nb; /* block size */
    if (j) {
      /* A[0:j-1,j:j+jb-1] <- A[0:j-1,0:j-1] A[0:j-1,j:j+jb-1] O(j^2 jb)*/
      if (0) F77_CALL(dtrmm)(&left,&up,&ntrans,&diag,&j,&jb,&d1,R,n,R + j * *n,n);
                                         else pdtrmm(&j,&jb,&d1,R,n,R + j * *n,n,nt,iwork,work);
      /* A[0:j-1,j:j+jb-1] <- A[0:j-1,j:j+jb-1] A[j:j+jb-1,j:j+jb-1]^{-1} O(j jb^2) */
      F77_CALL(dtrsm)(&right,&up,&ntrans,&diag,&j,&jb,&m1,R + *n * j + j,n,R + j * *n,n);
    }
    /* invert A[j:j+jb-1,j:j+jb-1] O(jb^3)*/
    F77_CALL(dtrti2)(&up,&diag,&jb,R + *n * j + j,n,&info); 
  }  
  FREE(work);FREE(iwork);
} /* mgcv_pbsi1 */

void mgcv_pbsi(double *R,int *r,int *nt) {
/* parallel back substitution inversion of upper triangular matrix using nt threads
   i.e. Solve of R Ri = I for Ri. Idea is to work through columns of I
   exploiting fact that full solve is never needed. Results are stored
   in lower triangle of R and an r-dimensional array, and copied into R 
   before exit (working back trough columns over-writing columns of R
   as we go is not an option without waiting for threads to return in 
   right order - which can lead to blocking by slowest). 

   This version avoids BLAS calls which have function call overheads. 
   In single thread mode it is 2-3 times faster than a BLAS call to dtrsm.
   It uses a load balancing approach, splitting columns up to threads in advance,
   so that thread initialisation overheads are minimised. This is important, 
   e.g. the alternative of sending each column to a separate thread removes 
   almost all advantage of parallel computing, because of the thread initiation 
   overheads. 

   Only disadvantage of this approach is that all cores are assumed equal. 

*/
  int i,j,k,r1,*a,b;
  double x,*d,*dk,*z,*zz,*z1,*rr,*Rjj,*r2;
  #ifdef OMP_REPORT
  Rprintf("mgcv_pbsi...");
  #endif
  d = (double *)CALLOC((size_t) *r,sizeof(double));
  if (*nt < 1) *nt = 1;
  if (*nt > *r) *nt = *r; /* no point having more threads than columns */
  /* now obtain block start columns, a. a[i] is start column of block i.  */
  a = (int *)CALLOC((size_t) (*nt+1),sizeof(int));
  a[0] = 0;a[*nt] = *r;
  x = (double) *r;x = x*x*x / *nt;
  /* compute approximate optimal split... */
  for (i=1;i < *nt;i++) a[i] = round(pow(x*i,1/3.0));
  for (i=*nt-1;i>0;i--) { /* don't allow zero width blocks */
    if (a[i]>=a[i+1]) a[i] = a[i+1]-1;
  }
  r1 = *r + 1;
  #ifdef OPENMP_ON
  #pragma omp parallel private(b,i,j,k,zz,z,z1,rr,Rjj,dk) num_threads(*nt)
  #endif
  { /* open parallel section */
    #ifdef OPENMP_ON
    #pragma omp for
    #endif
    for (b=0;b< *nt;b++) { /* b is thread/block index */
      for (i=a[b];i<a[b+1];i++) { /* i is column index */	
        k = *r - i -1; /* column of R in which to store result */
        Rjj = rr = R + *r * (ptrdiff_t) i + i; /* Pointer to R[i,i] */
        dk = d + k;
        *dk = 1 / *rr; /* Ri[i,i] */
        z = R + *r * (ptrdiff_t)k + k + 1; /* pointer to solution storage, below l.d. of R */
        rr -= i; /* pointer moved back to start of R[,i] */
        for (zz=z,z1 = z+i;zz<z1;zz++,rr++) *zz = *rr * *dk; /* sum_0_{i-1} Rii * zz[i] */
        for (j=i-1;j>=0;j--) {
          Rjj -= r1;
          dk = z + j;
          *dk /= - *Rjj;
          for (zz=z,z1=z+j,rr=Rjj-j;zz<z1;zz++,rr++) *zz += *rr * *dk; 
        }
      } /* end of for i (columns for block b) */
    } /* end of block, b, loop */
  } /* end parallel section */

  /* now copy the solution back into original R */  
  /* different block structure is required here */
  x = (double) *r;x = x*x / *nt;
  /* compute approximate optimal split... */
  for (i=1;i < *nt;i++) a[i] = round(sqrt(x*i));
  for (i=*nt-1;i>0;i--) { /* don't allow zero width blocks */
    if (a[i]>=a[i+1]) a[i] = a[i+1]-1;
  }

  #ifdef OPENMP_ON
#pragma omp parallel private(b,i,k,zz,rr,r2) num_threads(*nt)
  #endif 
  { /* open parallel section */
    #ifdef OPENMP_ON
    #pragma omp for
    #endif 
    for (b=0;b<*nt;b++) {
      for (i=a[b];i<a[b+1];i++) {
        k = *r - i - 1;
        zz = R + i * (ptrdiff_t) *r + i; /* leading diagonal element */
        *zz = d[k];
        zz -= i; /* shift back to start of column */
        /* copy rest of column and wipe below diagonal storage */
        for (rr = R + k * (ptrdiff_t)*r + k + 1, r2 = rr + i;rr<r2;rr++,zz++) {*zz = *rr;*rr = 0.0;}
      }
    } 
  } /* end parallel section */
  FREE(d);FREE(a);
  #ifdef OMP_REPORT
  Rprintf("done\n");
  #endif
} /* mgcv_pbsi */



void mgcv_Rpbsi(SEXP A, SEXP NT) {
/* parallel back sub inversion of upper triangular matrix A 
   designed for use with .call, but using mgcv_pbsi to do the work. 
*/
  int nt,r;
  double *R;
  nt = asInteger(NT);
  r = nrows(A);
  R = REAL(A);
  mgcv_pbsi1(R,&r,&nt);
} /* mgcv_Rpbsi */

void mgcv_PPt1(double *A,double *R,int *r,int *nt) {
/* Computes A=RR', where R is r by r upper triangular. BLAS version.
   Modified and stripped down from LAPACK dlauum.
   Currently uses nt threads to split dominant term by row. 

   n <- 4000;a <- rep(1:n,n);b <- rep(1:n,each=n);which(a<=b) -> ii;a[ii]+(b[ii]-1)*n->ii ## upper
   library(mgcv);R <- matrix(0,n,n);R[ii] <- runif(n*(n+1)/2)
   system.time(A <- mgcv:::pRRt(R,2))
   system.time(A2 <- tcrossprod(R));range(A-A2);plot(A,A2,pch=".")

*/  
  int i,j,k,ik,nb=50,ib,m,*s;
  char right='R',up='U',trans='T',no='N';
  double one=1.0,di,x;
  // first pointlessly copy R to A
  s = (int *)CALLOC((size_t) *nt + 1,sizeof(int));
  for (i=0;i<*r;i++) for (j=i;j<*r;j++) A[i + j * *r] = R[i + j * *r]; 
  for (i=0;i< *r;i+=nb) {
    ib= *r-i; if (ib>nb) ib = nb;
    /* A[0:i-1,i:i+ib-1] = A[0:i-1,i:i+ib-1]A'[i:i+ib-1,i:i+ib-1] O(ib^2 i) */
    F77_CALL(dtrmm)(&right,&up,&trans,&no,&i,&ib,&one,A + *r * i + i,r,A + *r * i, r);
    /* A[i:i+ib-1,i:i+ib-1] A[i:i+ib-1,i:i+ib-1]' (self overwrite)*/
    F77_CALL(dlauu2)(&up,&ib,A + *r * i + i,r,&j);
    if (i+ib < *r) {
      j = *r - i - ib;
      m = *nt;while (5*m>i && m>1) m--; /* don't use a pointless number of threads */
      di = i / (double) m; /* rows per thread */
      s[0]=0; for (x=0.0,k=1;k<m;k++) { x+= di;s[k] = (int)floor(x);} s[m] = i; /* row starts */ 
      /* A[0:i-1,i:i+ib-1] +=  A[0:i-1,i+ib:n-1] A'[i:i+ib-1,i+ib:n-1] O(i j ib) */
      #ifdef OPENMP_ON
      #pragma omp parallel for private(k,ik) num_threads(m)
      #endif
      for (k=0;k<m;k++) { /* row block loop */
	ik = s[k+1] - s[k]; /* number of rows to process */
        F77_CALL(dgemm)(&no,&trans,&ik,&ib,&j,&one,A + s[k] + (i+ib)* *r,r,A + i + (i+ib) * *r,
	                r,&one,A +  i * *r + s[k],r);
      }
      //F77_CALL(dgemm)(&no,&trans,&i,&ib,&j,&one,A+ (i+ib)* *r,r,A + i + (i+ib) * *r,
      //	              r,&one,A +  i * *r,r); // single thread version
      /* A[i:i+ib-1,i:i+ib-1] += A[i:i+ib-1,i:i+j-1] A[i:i+ib-1,i:i+j-1] O(j ib^2)*/
      F77_CALL(dsyrk)(&up,&no,&ib,&j,&one,A+i + (i+ib) * *r,r,&one,A+ i * *r + i, r);
    }
  }
  FREE(s);
  /* now copy upper triangle to lower */
  for (i=0;i<*r;i++) for (j=i+1;j<*r;j++) A[j + i * *r] =  A[i + j * *r];
} /* mgcv_PPt1 */



void mgcv_PPt(double *A,double *R,int *r,int *nt) {
/* Computes A=RR', where R is r by r upper triangular.
   Computation uses *nt cores. */
  double x,*ru,*rl,*r1,*ri,*rj,*Aji,*Aij;
  int *a,i,j,b;
  #ifdef OMP_REPORT
  Rprintf("mgcv_PPt...");
  #endif
  if (*nt < 1) *nt = 1;
  if (*nt > *r) *nt = *r; /* no point having more threads than columns */
  a = (int *)CALLOC((size_t) (*nt+1),sizeof(int));
  a[0] = 0;a[*nt] = *r;
  /* It is worth transposing R into lower triangle */
  x = (double) *r;x = x*x / *nt;
  /* compute approximate optimal split... */
  for (i=1;i < *nt;i++) a[i] = round(*r - sqrt(x*(*nt-i)));
  for (i=1;i <= *nt;i++) { /* don't allow zero width blocks */
    if (a[i]<=a[i-1]) a[i] = a[i-1]+1;
  }

  #ifdef OPENMP_ON
  #pragma omp parallel private(b,i,ru,rl,r1) num_threads(*nt)
  #endif
  { /* open parallel section */
    #ifdef OPENMP_ON
    #pragma omp for
    #endif 
    for (b=0;b<*nt;b++) {
      for (i=a[b];i<a[b+1];i++) { /* work through rows */
        ru = rl = R + i * (ptrdiff_t)*r + i; /* diagonal element R[i,i] */
        ru += *r; /* move to R[i,i+1] */
        /* now copy R[i,i:r] to R[i:r,i]... */
        for (r1=rl + *r - i,rl++;rl<r1;rl++,ru += *r) *rl = *ru;
      }
    } /* end of for b block loop */
  } /* end parallel section */

  /* Now obtain approximate load balancing block starts for product... */
  x = (double) *r;x = x*x*x / *nt;
  /* compute approximate optimal split... */
  for (i=1;i < *nt;i++) a[i] = round(*r - pow(x*(*nt-i),1/3.0));
  for (i=1;i <= *nt;i++) { /* don't allow zero width blocks */
    if (a[i]<=a[i-1]) a[i] = a[i-1]+1;
  }
  #ifdef OPENMP_ON
  #pragma omp parallel private(x,b,i,j,ru,rl,r1,rj,ri,Aij,Aji) num_threads(*nt)
  #endif
  { /* start parallel block */
    #ifdef OPENMP_ON
    #pragma omp for
    #endif
    for (b=0;b< *nt;b++) {
      for (i=a[b];i<a[b+1];i++) { /* loop over this block's rows */
        Aij = Aji = A + i * (ptrdiff_t)*r + i;
        rl = ru =  R + i * (ptrdiff_t)*r + i;
        r1 = R + (i+1) * (ptrdiff_t)*r; /* upper limit for ri */
        for (j=i;j < *r;j++,Aij += *r,Aji++,ru++,rl+= *r + 1) {
          //rj = rl;// R + j * *r + j;
          //ri = ru; // ri=R + i * *r + j;
          //for (x=0.0,k=j;k<*r;k++,rj++,ri++) x += *rj * *ri;
          for (x=0.0,rj=rl,ri=ru;ri < r1;rj++,ri++) x += *rj * *ri;
          *Aij = *Aji = x;  
        }
      }
    } /* block loop end */
  } /* end parallel block */
  /* wipe lower triangle of input R */
  x = (double) *r;x = x*x / *nt;
  /* compute approximate optimal split... */
  for (i=1;i < *nt;i++) a[i] = round(*r - sqrt(x*(*nt-i)));
  for (i=1;i <= *nt;i++) { /* don't allow zero width blocks */
    if (a[i]<=a[i-1]) a[i] = a[i-1]+1;
  }
  #ifdef OPENMP_ON
  #pragma omp parallel private(b,i,rl,r1) num_threads(*nt)
  #endif
  { /* start parallel block */
    #ifdef OPENMP_ON
    #pragma omp for
    #endif 
    for (b=0;b<*nt;b++) {
      for (i=a[b];i<a[b+1];i++) { /* work through rows */
        rl = R + i * (ptrdiff_t)*r + i; /* diagonal element R[i,i] */
        /* now wipe sub diagonal elements...*/
        for (r1=rl + *r - i,rl++;rl<r1;rl++) *rl = 0.0;
      }
    } /* end of for b block loop */
  } /* end parallel block */ 
  FREE(a); 
  #ifdef OMP_REPORT
  Rprintf("done\n");
  #endif
} /* mgcv_PPt */

void mgcv_RPPt(SEXP a,SEXP r, SEXP NT) {
/* Form a = rr' where r is upper triangular. 
   designed for use with .Call, but using mgcv_PPt to do the work. 
*/
  int nt,n;
  double *R,*A;
  nt = asInteger(NT);
  n = nrows(a);
  A = REAL(a);
  R = REAL(r);
  mgcv_PPt1(A,R,&n,&nt);
} /* mgcv_Rpbsi */


void mgcv_backsolve0(double *R,int *r,int *c,double *B,double *C, int *bc) 
/* BLAS free version
   Finds C = R^{-1} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard back 
   substitution loop.
*/  
{ int i,j,k;
  double x,*pR,*pC;
  for (j=0;j<*bc;j++) { /* work across columns of B & C */
    for (i = *c-1;i>=0;i--) { /* work up each column of B & C */
      x = 0.0;
      /* for (k=i+1;k<*c;k++) x += R[i + *r * k] * C[k + j * *c]; ...following replaces...*/
      pR = R + i + (i+1) * (ptrdiff_t)*r;pC = C + j * (ptrdiff_t)*c + i + 1;
      for (k=i+1;k<*c;k++,pR+= *r,pC++) x += *pR * *pC;      
      C[i + (ptrdiff_t)j * *c] = (B[i + (ptrdiff_t)j * *c] - x)/R[i + (ptrdiff_t)*r * i];
    }
  }
}
void mgcv_forwardsolve0(double *R,int *r,int *c,double *B,double *C, int *bc) 
/* BLAS free version
   Finds C = R^{-T} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard forward 
   substitution loop.
*/  
{ int i,j,k;
  double x;
  for (j=0;j<*bc;j++) { /* work across columns of B & C */
    for (i = 0;i< *c;i++) { /* work down each column of B & C */
      x=0.0;
      for (k=0;k<i;k++) x+= C[k + (ptrdiff_t)j * *c] * R[k + (ptrdiff_t)i * *r]; 
      C[i + (ptrdiff_t)j * *c] = (B[i + (ptrdiff_t)j * *c] - x) / R[i + (ptrdiff_t)i * *r]; 
    }
  }
}

void mgcv_backsolve(double *R,int *r,int *c,double *B,double *C, int *bc,int *right) 
/* BLAS version
   If *right==0:
   Finds C = R^{-1} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard back 
   substitution loop.
   Otherwise:
   Finds C=BR^{-1} where B is bc by c
*/  
{ double *pR,*pC,alpha=1.0;
  int n,m;
  char side='L',uplo='U',transa='N',diag='N';
  if (*right) { side = 'R';m = *bc,n= *c;} else {m = *c;n= *bc;}
  for (pC=C,pR=pC+ *bc * (ptrdiff_t)*c;pC<pR;pC++,B++) *pC = *B; /* copy B to C */
  F77_CALL(dtrsm)(&side,&uplo,&transa, &diag,&m, &n, &alpha,R, r,C,&m);
} /* mgcv_backsolve */


 void mgcv_forwardsolve(double *R,int *r,int *c,double *B,double *C, int *bc,int *right) 
/* BLAS version
   If *right == 0:
   Finds C = R^{-T} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard forward 
   substitution loop.
   Otherwise:
   Finds C = B R^{-T} where B is bc by c.
*/  
{ double *pR,*pC,alpha=1.0;
  int m,n;
  char side='L',uplo='U',transa='T',diag='N';
  if (*right) { side='R';m = *bc;n= *c; } else {m = *c;n= *bc;}
  for (pC=C,pR=pC+ *bc * (ptrdiff_t) *c;pC<pR;pC++,B++) *pC = *B; /* copy B to C */
  F77_CALL(dtrsm)(&side,&uplo,&transa, &diag,&m, &n, &alpha,R, r,C,&m);
} /* mgcv_forwardsolve */

void mgcv_pforwardsolve(double *R,int *r,int *c,double *B,double *C, int *bc,int *nt) 
/* parallel forward solve, using nt threads.
   Finds C = R^{-T} B where R is the c by c matrix stored in the **upper** triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard forward 
   substitution loop. NOTE: R stored in UPPER not lower triangle.  
*/  
{ double *pR,*pC,alpha=1.0;
  int cpt,cpf,nth,i,cp;
  char side='L',uplo='U',transa='T',diag='N';
  cpt = *bc / *nt; /* cols per thread */
  if (cpt * *nt < *bc) cpt++;
  nth = *bc/cpt;
  if (nth * cpt < *bc) nth++;
  cpf = *bc - cpt * (nth-1); /* columns on final block */ 
  for (pC=C,pR=pC+ *bc * (ptrdiff_t)*c;pC<pR;pC++,B++) *pC = *B; /* copy B to C */
  //Rprintf("r = %d, c= %d bc = %d, nth= %d, cpt = %d, cpf = %d",*r,*c,*bc,nth,cpt,cpf);
  #ifdef OPENMP_ON
  #pragma omp parallel for private(i,cp) num_threads(nth)
  #endif
  for (i=0;i<nth;i++) {
    if (i==nth-1) cp = cpf; else cp = cpt;
    F77_CALL(dtrsm)(&side,&uplo,&transa, &diag,c, &cp, &alpha,R, r,C + i * (ptrdiff_t) cpt * *c,c);
  }
} /* mgcv_pforwardsolve */

void mgcv_pbacksolve(double *R,int *r,int *c,double *B,double *C, int *bc,int *nt) 
/* parallel back solve, using nt threads.
   Finds C = R^{-T} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). 
*/  
{ double *pR,*pC,alpha=1.0;
  int cpt,cpf,nth,i,cp;
  char side='L',uplo='U',transa='N',diag='N';
  cpt = *bc / *nt; /* cols per thread */
  if (cpt * *nt < *bc) cpt++;
  nth = *bc/cpt;
  if (nth * cpt < *bc) nth++;
  cpf = *bc - cpt * (nth-1); /* columns on final block */ 
  for (pC=C,pR=pC+ *bc * (ptrdiff_t)*c;pC<pR;pC++,B++) *pC = *B; /* copy B to C */
  //Rprintf("r = %d, c= %d bc = %d, nth= %d, cpt = %d, cpf = %d",*r,*c,*bc,nth,cpt,cpf);
  #ifdef OPENMP_ON
  #pragma omp parallel for private(i,cp) num_threads(nth)
  #endif
  for (i=0;i<nth;i++) {
    if (i==nth-1) cp = cpf; else cp = cpt;
    F77_CALL(dtrsm)(&side,&uplo,&transa, &diag,c, &cp, &alpha,R, r,C + i * (ptrdiff_t) cpt * *c,c);
  }
} /* mgcv_pbacksolve */

SEXP mgcv_Rpforwardsolve(SEXP R, SEXP B,SEXP NT) {
/* .Call wrapper for mgcv_pforwardsolve.
   Return object a matrix.

*/
  int c,r,nt,bc;
  double *R0,*b,*C0;
  SEXP C;
  nt = asInteger(NT);
  r = nrows(R);c = ncols(R);
  R0 = REAL(R);
  bc = ncols(B);
  b = REAL(B);  
  C = PROTECT(allocMatrix(REALSXP,c,bc));
  C0 = REAL(C);

  mgcv_pforwardsolve(R0,&r,&c,b,C0,&bc,&nt);
   
  UNPROTECT(1);
  return(C);
} /* mgcv_Rpforwardsolve */

SEXP mgcv_Rpbacksolve(SEXP R, SEXP B,SEXP NT) {
/* .Call wrapper for mgcv_pbacksolve.
   Return object a matrix.

*/
  int c,r,nt,bc;
  double *R0,*b,*C0;
  SEXP C;
  nt = asInteger(NT);
  r = nrows(R);c = ncols(R);
  R0 = REAL(R);
  bc = ncols(B);
  b = REAL(B);  
  C = PROTECT(allocMatrix(REALSXP,c,bc));
  C0 = REAL(C);

  mgcv_pbacksolve(R0,&r,&c,b,C0,&bc,&nt);
   
  UNPROTECT(1);
  return(C);
} /* mgcv_Rpbacksolve */



void row_block_reorder(double *x,int *r,int *c,int *nb,int *reverse) {
/* x contains an r by c matrix, which is to be split into k = ceil(r/nb)
   row-wise blocks each containing nb rows, except for the last.
  
   In forward mode, this splits the matrix blocks so they are stored one 
   after another in x.

   This routine re-orders these k submatrices, so that they occur one after 
   another in x. This is done largely in-situ, with the aid of 2 indexing
   arrays, and some extra storage required to deal with fact that last block
   may have fewer rows than rest.

   The elements of one column belonging to one block are referred to as a 
   segment.
*/
  ptrdiff_t *a,*s,k,nbf=0,i,j,q,si,ti,ns,ns_main,ns_extra;
  double *x0,*x1,tmp,*extra=NULL;
  k = *r / *nb; /* number of blocks */
  if (*r > *nb * k) {
    nbf = *r - *nb * k;k++; /* nbf number of rows in final block */
  }
  /* first task is to pad the end block segments, so that all segments
     have equal length, otherwise efficient segment swapping is not
     possible. This requires spilling over into extra storage. */
  ns = k * *c; /* total number of segments */
  if (nbf) { /* only do this if final block shorter than rest */
    ns_main = (*r * *c) / *nb; /* full segments fitting in x */
    ns_extra = ns - ns_main; /* segments requiring extra storage */
    extra = (double *) CALLOC((size_t) (*nb * ns_extra),sizeof(double));
    x0 = extra + *nb * ns_extra - 1; /* end of extra */ 
    x1 = x + *r * *c -1 ; /* end of x */

    if (*reverse) { /* blocks back into single matrix */
      /* expand end segments out into extra storge */ 
      for (i=ns-1;i >= ns_main;i--) {
        x0 -= *nb - nbf; /* skip padding in target */
        for (j=0;j<nbf;j++,x0--,x1--) *x0 = *x1; 
      }
      x0 = x + ns_main * *nb -1; /* target moved to main block */
      for (;i >= ns - *c;i--) {
        x0 -= *nb - nbf; /* skip padding in target */
        for (j=0;j<nbf;j++,x0--,x1--) *x0 = *x1; 
      }
    } else { /* forward mode - matrix to separated blocks */
      /* first copy the end segments into extra */
      for (i=ns-1;i>=ns_main;i--) { /* work down through segments */
        if ((i+1)%k) { /* not a short segment */     
          for (j = 0;j < *nb;j++,x0--,x1--) *x0 = *x1;
        } else {
          x0 -= (*nb - nbf); /* skip padding in target */
          for (j = 0;j < nbf;j++,x0--,x1--) *x0 = *x1; /* fill rest from source */
        }
      }
      /* now copy segments into x with padding ... */
      x0 = x + ns_main * *nb - 1; /* end of main block segment storage */
      for (;i>=0;i--) { /* continue down through segments */
        if ((i+1)%k) { /* not a short segment */     
          for (j = 0;j < *nb;j++,x0--,x1--) *x0 = *x1;
        } else {
          x0 -= (*nb - nbf);
          for (j = 0;j < nbf;j++,x0--,x1--) *x0 = *x1;
        }
      }
    } /* end of forward mode padding */
  } else { /* segments already equal length */
    ns_main = ns;ns_extra=0; /* all segments fit into x */
  }
  /* now re-arrange row-block wise... */

  /* a[i] is original segment now in segment i... */ 
  a = (ptrdiff_t *) CALLOC((size_t) (k * *c),sizeof(ptrdiff_t));
  /* s[i] is segment now containing original segment i */
  s = (ptrdiff_t *) CALLOC((size_t) (k * *c),sizeof(ptrdiff_t));
  for (i=0;i<k * *c;i++) s[i] = a[i] = i;
  ti=0; /* target segment */
  for (i=0;i<k;i++) for (j=0;j<*c;j++,ti++) {
      if (*reverse) si = s[(ti%k) * *c + ti/k]; /* source segment in reverse mode */
      else si = s[i + j * k]; /* source segment forward mode */
      if (ti < ns_main) x0 = x + ti * *nb; /* target segment in main storage */
      else x0 = extra + (ti - ns_main) * *nb; /* target segment in extra storage */ 
      if (si < ns_main) x1 = x + si * *nb;
      else x1 = extra + (si - ns_main) * *nb;
      /* swap segment contents... */
      for (q=0;q < *nb;q++,x0++,x1++) { tmp = *x0;*x0 = *x1;*x1 = tmp;}
      /*  update indices... */
      q = a[si]; a[si] = a[ti]; a[ti] = q;
      s[a[si]] = si;s[a[ti]]=ti;
  }  

  /* now it only remains to delete padding from end block */
  if (nbf) {
    if (*reverse) { 
      si = k; /* start segment */
      x0 = x + *r; /* start target */
      x1 = x + k * *nb; /* start source */
      for (;si<ns;si++) {
        if (si == ns_main) x1 = extra; /* source shifted to extra */
        if ((si+1)%k) { /* its a full segment */
          for (i=0; i < *nb;i++,x0++,x1++) *x0 = *x1;
        } else { /* short segment */
          for (i=0;i<nbf;i++,x0++,x1++) *x0 = *x1; /* only copy info, not padding */
          x1 += *nb - nbf; /* move source past padding */
        }
      }  
    } else { /* forward mode */
      x1 = x0 = x + (k-1) * *c * *nb; /* address of first target (and first source) */
      si = *c * (k-1); /* source segment index */
      for (i=0;i < *c;i++,si++,x1+=(*nb-nbf)) {
        if (si == ns_main) x1 = extra;  /* reset source blocks to extra */
        for (j=0;j<nbf;j++,x0++,x1++) *x0 = *x1;   
      }
    }
  } 

  FREE(a); FREE(s);if (nbf) FREE(extra);
} /* row_block_reorder */


int get_qpr_k(int *r,int *c,int *nt) {
/* For a machine with nt available cores, computes k, the best number of threads to 
   use for a parallel QR.
*/
  #ifdef OPENMP_ON
  int k;double kd,fkd,x,ckd;
  kd = sqrt(*r/(double)*c); /* theoretical optimum */
  if (kd <= 1.0) k = 1; else
    if (kd > *nt) k = *nt; else {
      fkd = floor(kd);ckd = ceil(kd);
      if (fkd>1) x =  *r / fkd + fkd * *c; else x = *r;
      if (*r / ckd + ckd * *c < x) 
        k = (int)ckd; else k = (int)fkd;
    }
  return(k);
  #else
  return(1); /* can only use 1 thread if no openMP support */ 
  #endif
}

void getRpqr(double *R,double *x,int *r, int *c,int *rr,int *nt) {
/* x contains qr decomposition of r by c matrix as computed by mgcv_pqr 
   This routine simply extracts the c by c R factor into R. 
   R has rr rows, where rr == c if R is square. 

*/
  int i,j,n,rows;
  double *Rs;
  Rs = x;n = *r;
  rows = *c; if (rows > *rr) rows = *rr;
  for (i=0;i<rows;i++) for (j=0;j<*c;j++) if (i>j) R[i + *rr * j] = 0; else
	R[i + *rr * j] = Rs[i + n * j];
} /* getRpqr */

void getRpqr0(double *R,double *x,int *r, int *c,int *rr,int *nt) {
/* x contains qr decomposition of r by c matrix as computed by mgcv_pqr 
   This routine simply extracts the c by c R factor into R. 
   R has rr rows, where rr == c if R is square. 

   This version matches mgcv_pqrqy0, which is inferior to current code.
*/
  int i,j,k,n;
  double *Rs;
  k = get_qpr_k(r,c,nt); /* number of blocks used */
  if (k==1) { /* actually just a regular serial QR */ 
    Rs = x;n = *r;
  } else {
    n = k * *c; /* rows of R */ 
    Rs = x + *r * *c; /* source R */
  }
  for (i=0;i<*c;i++) for (j=0;j<*c;j++) if (i>j) R[i + *rr * j] = 0; else
	R[i + *rr * j] = Rs[i + n * j];
} /* getRpqr0 */

void mgcv_pqrqy0(double *b,double *a,double *tau,int *r,int *c,int *cb,int *tp,int *nt) {
/* Applies factor Q of a QR factor computed in parallel to b. 
   If b is physically r by cb, but if tp = 0 it contains a c by cb matrix on entry, while if
   tp=1 it contains a c by cb matrix on exit. Unused elments of b on entry assumed 0. 
   a and tau are the result of mgcv_pqr  

   This version matches mgcv_pqr0, which scales less well than current code. 
*/
  int i,j,k,l,left=1,n,nb,nbf,nq,TRUE=1,FALSE=0;
  double *x0,*x1,*Qb;  
  #ifdef OMP_REPORT
  Rprintf("mgcv_pqrqy0...");
  #endif
  k = get_qpr_k(r,c,nt); /* number of blocks in use */
  if (k==1) { /* single block case */ 
    if (*tp == 0 ) {/* re-arrange so b is a full matrix */
      x0 = b + *r * *cb -1; /* end of full b (target) */
      x1 = b + *c * *cb -1; /* end of used block (source) */
      for (j= *cb;j>0;j--) { /* work down columns */
        /*for (i = *r;i > *c;i--,x0--) *x0 = 0.0;*/ /* clear unused */
        x0 -= *r - *c; /* skip unused */
        for (i = *c;i>0;i--,x0--,x1--) { 
          *x0 = *x1; /* copy */
          if (x0!=x1) *x1 = 0.0; /* clear source */
        }
      }
    } /* if (*tp) */
    mgcv_qrqy(b,a,tau,r,cb,c,&left,tp);
    if (*tp) { /* need to strip out the extra rows */
      x1 = x0 = b;
      for (i=0;i < *cb;i++,x1 += *r - *c) 
	for (j=0;j < *c;j++,x0++,x1++) *x0 = *x1; 
    }
    return;
  }
  
  /* multi-block case starts here */

  nb = (int)ceil(*r/(double)k); /* block size - in rows */
  nbf = *r - (k-1)*nb; /* end block size */ 
  Qb = (double *)CALLOC((size_t) (k * *c * *cb),sizeof(double));
  nq = *c * k;
  if (*tp) { /* Q'b */
    /* first the component Q matrices are applied to the blocks of b */
    if (*cb > 1) { /* matrix case - repacking needed */
      row_block_reorder(b,r,cb,&nb,&FALSE);
    }
    #ifdef OPENMP_ON
    #pragma omp parallel private(i,j,l,n,x1) num_threads(k)
    #endif
    { /* open parallel section */
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (i=0;i<k;i++) {
        if (i == k-1) n = nbf; else n = nb; /* number of rows in this block */
        x1 = b+ i * nb * *cb; /* start of current block of b */
        mgcv_qrqy(x1,a + i * nb * *c,tau + i * *c,&n,cb,c,&left,tp);
        /* extract the upper c by cb block into right part of Qb */
        for (j = 0; j < *c;j++) for (l=0;l< *cb;l++) Qb[j + i * *c + nq * l] = x1[j + n * l];    
      }
    } /* end of parallel */      
    /* Now multiply Qb by the combined Q' */
    mgcv_qrqy(Qb,a + *r * *c,tau + k * *c,&nq,cb,c,&left,tp);
    /* and finally copy the top c rows of Qb into b */
    x0 = b; /* target */
    x1 = Qb; /* source */
    for (i=0;i<*cb;i++) { 
      for (j=0;j<*c;j++,x0++,x1++) *x0 = *x1;
      x1 += (k-1) * *c; /* skip rows of Qb */ 
    }
  } else {  /* Qb */
    /* copy c by cb matrix at start of b into first c rows of Qb */
    x0=Qb; /* target */
    x1=b;
    for (i=0;i<*cb;i++) {
      for (j=0;j<*c;j++,x0++,x1++) {
        *x0 = *x1;
        *x1=0.0; /* erase or left in below */
      }
      x0 += (k-1) * *c; /* skip rows of Qb */ 
    }
    /* now apply the combined Q factor */
    mgcv_qrqy(Qb,a + *r * *c,tau + k * *c,&nq,cb,c,&left,tp);
    /* the blocks of Qb now have to be copied into separate blocks in b */
    #ifdef OPENMP_ON
    #pragma omp parallel private(i,j,l,n,x0,x1) num_threads(k)
    #endif
    { /* open parallel section */
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (i=0;i<k;i++) {
        if (i == k-1) n = nbf; else n = nb; /* number of rows in this block */
        x0 = b + nb * *cb * i; /* target block start */
        x1 = Qb + i * *c; /* source row start */
        for (j=0;j < *cb;j++) { /* across the cols */
          for (l=0;l< *c;l++,x0++,x1++) *x0 = *x1;
          x0 += n - *c; /* to start of next column */
          x1 += nq - *c; /* ditto */ 
        }
        x1 = b+ i * nb * *cb; /* start of current block of b */
        mgcv_qrqy(x1,a + i * nb * *c,tau + i * *c,&n,cb,c,&left,tp);
      }
    } /* end of parallel section */
    /* now need to unscramble separate blocks back into b as regular matrix */
    if (*cb>1) row_block_reorder(b,r,cb,&nb,&TRUE);
  } 
  #ifdef OMP_REPORT
  Rprintf("done\n");
  #endif
  FREE(Qb);
}  /* mgcv_pqrqy0 */


void mgcv_pqrqy(double *b,double *a,double *tau,int *r,int *c,int *cb,int *tp,int *nt) {
/* Applies factor Q of a QR factor computed by mgcv_pqr to b. 
   b is physically r by cb, but if tp = 0 it contains a c by cb matrix on entry, while if
   tp=1 it contains a c by cb matrix on exit. Unused elments of b on entry assumed 0. 
   a and tau are the result of mgcv_pqr.

   Note that in multi-threaded mode this uses mgcv_pqr0, which is thread safe, but level 2.
   mgcv_pqr is level 3 but not thread safe.  mgcv_pqrqy itself is not thread safe - i.e. 
   this should not be called from a parallel section (which would be dumb anyway). 

*/
  int i,j,ki,k,left=1,nth;
  double *x0,*x1,*aii,*p0;  
  #ifdef OMP_REPORT
  Rprintf("mgcv_pqrqy...");
  #endif
  //Rprintf("pqrqy %d ",*nt);
  if (*tp == 0 ) {/* re-arrange so b is a full matrix */
    x0 = b + *r * *cb -1; /* end of full b (target) */
    x1 = b + *c * *cb -1; /* end of used block (source) */
    for (j= *cb;j>0;j--) { /* work down columns */
      /*for (i = *r;i > *c;i--,x0--) *x0 = 0.0;*/ /* clear unused */
      x0 -= *r - *c; /* skip unused */
      for (i = *c;i>0;i--,x0--,x1--) { 
        *x0 = *x1; /* copy */
        if (x0!=x1) *x1 = 0.0; /* clear source */
      }
    }
  } /* if (*tp) */
  if (*cb==1 || *nt==1) mgcv_qrqy(b,a,tau,r,cb,c,&left,tp);
  else { /* split operation by columns of b */ 
    /* set leading diagonal elements of a to 1 and store them */
    aii = (double *)CALLOC((size_t)*c,sizeof(double));
    for (k=*r+1,x0=aii,x1=aii + *c,p0=a;x0<x1;x0++,p0 += k) { *x0 = *p0;*p0 = 1.0; } 
    nth = *nt;if (nth > *cb) nth = *cb;
    k = *cb/nth; 
    if (k*nth < *cb) k++; /* otherwise last thread is rate limiting */
    if (k*(nth-1) >= *cb) nth--; /* otherwise last thread has no work */
    #ifdef OPENMP_ON
    #pragma omp parallel private(i,j,ki) num_threads(nth)
    #endif
    { /* open parallel section */
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (i=0;i<nth;i++) {
        j = i * k;
        if (i==nth-1) ki = *cb - j; else ki = k;
        /* note this is level 2 at present, but thread safe */
        mgcv_qrqy0(b + j * *r,a,tau,r,&ki,c,&left,tp);
      }
    } /* end parallel section */
    for (k=*r+1,x0=aii,x1=aii + *c,p0=a;x0<x1;x0++,p0 += k) *p0= *x0;
    FREE(aii); 
  } 
  if (*tp) { /* need to strip out the extra rows */
    x1 = x0 = b;
    for (i=0;i < *cb;i++,x1 += *r - *c) 
    for (j=0;j < *c;j++,x0++,x1++) *x0 = *x1; 
  }  
  /* restore A[i,i]... */
 
  #ifdef OMP_REPORT
  Rprintf("done\n");
  #endif
  return;
}  /* mgcv_pqrqy */

void mgcv_pqr(double *x,int *r, int *c,int *pivot, double *tau, int *nt) {
/* Function called from mgcv to request a multi-threaded pivoted QR decomposition.
   Basically a wrapper to call the actual code, allowing painless changing of the 
   underlying technology.

   Currently Block Pivoted QR scales best from the codes available. 
   Hard coded block size (30) is not ideal. 
*/
  //Rprintf("pqr %d ",*nt);
  if (*nt==1) mgcv_qr(x,r,c,pivot,tau); else { /* call bpqr */
    /* int bpqr(double *A,int n,int p,double *tau,int *piv,int nb,int nt)*/
    bpqr(x,*r,*c,tau,pivot,30,*nt); /* 30 is hard coded block size here */ 
  }
} /* mgcv_pqr */


void mgcv_pqr0(double *x,int *r, int *c,int *pivot, double *tau, int *nt) {
/* parallel qr decomposition using up to nt threads. 
   * x is an r*c+nt*c^2 array. On entry first r*c entries is matrix to decompose.
   * pivot is a c vector and tau is a (nt+1)*c vector.
   On exit x, pivot and tau contain decompostion information in packed form
   for use by getRpqr and mgcv_pqrqy
   Routine first computes k, the optimal number of threads to use from 1 to nt.

   strategy is to split into row blocks, QR each and then combine the decompositions
   - good for n>>p, not so good otherwise. 

   - this is old code, which is uniformly less efficient than replacement.
*/

  int i,j,k,l,*piv,nb,nbf,n,TRUE=1,FALSE=0,nr; 
  double *R,*R1,*xi; 
  #ifdef OMP_REPORT
  Rprintf("mgcv_pqr0...");
  #endif
  k = get_qpr_k(r,c,nt);/* number of threads to use */
 
  if (k==1) mgcv_qr(x,r,c,pivot,tau); else { /* multi-threaded version */
    nb = (int)ceil(*r/(double)k); /* block size */
    nbf = *r - (k-1)*nb; /* end block size */
    /* need to re-arrange row blocks so that they can be split between qr calls */
    row_block_reorder(x,r,c,&nb,&FALSE); 
    piv = (int *)CALLOC((size_t) (k * *c),sizeof(int));
    R = x + *r * *c ; /* pointer to combined unpivoted R matrix */
    nr = *c * k; /* number of rows in R */
    #ifdef OPENMP_ON
    #pragma omp parallel private(i,j,l,n,xi,R1) num_threads(k)
    #endif
    { /* open parallel section */
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (i=0;i<k;i++) { /* QR the blocks */
        if (i==k-1) n = nbf; else n = nb; 
        xi = x + nb * i * *c; /* current block */
        mgcv_qr(xi,&n,c,piv + i * *c,tau + i * *c);
        /* and copy the R factors, unpivoted into a new matrix */
        R1 = (double *)CALLOC((size_t)(*c * *c),sizeof(double));
        for (l=0;l<*c;l++) for (j=l;j<*c;j++) R1[l + *c * j] = xi[l + n * j]; 
        /* What if final nbf is less than c? - can't happen  */
        pivoter(R1,c,c,piv + i * *c,&TRUE,&TRUE); /* unpivoting the columns of R1 */
        for (l=0;l<*c;l++) for (j=0;j<*c;j++) R[i * *c +l + nr *j] = R1[l+ *c * j];
        FREE(R1);
      }
    } /* end of parallel section */
    FREE(piv);
    n = k * *c;
    mgcv_qr(R,&n,c,pivot,tau + k * *c); /* final pivoted QR */
  } 
  #ifdef OMP_REPORT
  Rprintf("done\n");
  #endif
} /* mgcv_pqr0 */

void mgcv_qr(double *x, int *r, int *c,int *pivot,double *tau)
/* call LA_PACK to get pivoted QR decomposition of x
   tau is an array of length min(r,c)
   pivot is array of length c, zeroed on entry, pivoting order on return.
   On exit upper triangle of x is R. Below upper triangle plus tau 
   represent reflectors making up Q.
   pivoting is always performed (not just when matrix is rank deficient), so
   leading diagonal of R is in descending order of magnitude.
   library(mgcv)
   r<-4;c<-3
   X<-matrix(rnorm(r*c),r,c)
   pivot<-rep(1,c);tau<-rep(0,c)
   um<-.C("mgcv_qr",as.double(X),as.integer(r),as.integer(c),as.integer(pivot),as.double(tau))
   qr.R(qr(X));matrix(um[[1]],r,c)[1:c,1:c]
*/
{ int info,lwork=-1,*ip;
  double work1,*work;
  /* workspace query */
  /* Args: M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO */
  F77_CALL(dgeqp3)(r,c,x,r,pivot,tau,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)CALLOC((size_t)lwork,sizeof(double));
   /* actual call */
  F77_CALL(dgeqp3)(r,c,x,r,pivot,tau,work,&lwork,&info); 
  FREE(work);
  /*if (*r<*c) lwork= *r; else lwork= *c;*/ 
  for (ip=pivot;ip < pivot + *c;ip++) (*ip)--;
  /* ... for 'tis C in which we work and not the 'cursed Fortran... */
  
} /* end mgcv_qr */

void mgcv_qr2(double *x, int *r, int *c,int *pivot,double *tau)
/* call LA_PACK to get  QR decomposition of x
   tau is an array of length min(r,c)
   pivot is array of length c, zeroed on entry, pivoting order on return.
   On exit upper triangle of x is R. Below upper triangle plus tau 
   represent reflectors making up Q.
   pivoting is not performed in this case, but the pivoting index is returned anyway. 
   library(mgcv)
   r<-4;c<-3
   X<-matrix(rnorm(r*c),r,c)
   pivot<-rep(1,c);tau<-rep(0,c)
   um<-.C("mgcv_qr",as.double(X),as.integer(r),as.integer(c),as.integer(pivot),as.double(tau))
   qr.R(qr(X));matrix(um[[1]],r,c)[1:c,1:c]
*/
{ int info,*ip,i;
  double *work;
  work=(double *)CALLOC((size_t)*r,sizeof(double));
   /* actual call */
  /* Args: M, N, A, LDA, TAU, WORK, INFO */
  F77_CALL(dgeqr2)(r,c,x,r,tau,work,&info); 
  FREE(work);
  /*if (*r<*c) lwork= *r; else lwork= *c;*/ 
  for (i=0,ip=pivot;ip < pivot + *c;ip++,i++) *ip = i;
  /* ... pivot index equivalent to no pivoting */
  
} /* end mgcv_qr2 */

void mgcv_qrqy0(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp) {
/* mgcv_qrqy is not thread safe, because of the behaviour of dormqr (and similar functions).
   This version uses dlarf for a thread safe routine, but this is then level 2, *and requires
   modification of a before entry*.

   Applies k reflectors of Q of a QR decomposition to r by c matrix b.
   Apply Q from left if left!=0, right otherwise.
   Transpose Q only if tp!=0.
   Information about Q has been returned from mgcv_qr, and is stored in tau and 
   *on and* below the leading diagonal of a. In fact the leading diagonal elements must 
   be set to 1 before entry to this routine. 
    
   SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*/
  char side='L';
  int lda,lwork=-1,incv=1,ri,i0,i1,ii,i;
  double *work,*v;
 
  if (! *left) { 
    side='R';lda = *c;lwork = *r;
  } else { 
    lda= *r;lwork = *c;
  }
  /* calloc is thread safe, CALLOC ambiguous */
  work=(double *)calloc((size_t)lwork,sizeof(double));
  if ((*left && ! *tp)||(! *left && *tp)) { /* kth H applied first */
    i0 = *k - 1;i1=-1;ii=-1;
  } else { /* 1st H applied first */
    i0 = 0;i1 = *k;ii=1;
  }
  for (i=i0;i!=i1;i+=ii) {
    v = a + lda * i + i; /* start of v */
    ri = *r - i; /* number of rows in sub-block to which this applies */
    F77_CALL(dlarf)(&side,&ri,c,v,&incv,tau+i,b+i,r,work);
  }
  free(work);
  } /* mgcv_qrqy0 */

void mgcv_qrqy(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp)
/* applies k reflectors of Q of a QR decomposition to r by c matrix b.
   Apply Q from left if left!=0, right otherwise.
   Transpose Q only if tp!=0.
   Information about Q has been returned from mgcv_qr, and is stored in tau and 
   below the leading diagonal of a.
   library(mgcv)
   r<-4;c<-3
   X<-matrix(rnorm(r*c),r,c)
   qrx<-qr(X)
   pivot<-rep(1,c);tau<-rep(0,c)
   um<-.C("mgcv_qr",a=as.double(X),as.integer(r),as.integer(c),as.integer(pivot),tau=as.double(tau))
   y<-1:4;left<-1;tp<-0;cy<-1
   er<-.C("mgcv_qrqy",as.double(y),as.double(um$a),as.double(um$tau),as.integer(r),as.integer(cy),as.integer(c),
        as.integer(left),as.integer(tp),PACKAGE="mgcv")
   er[[1]];qr.qy(qrx,y)

   dormqr is not thread safe (with feasible memory use, at least).
   A block threadsafe version could be built using dlarft, dlarfb.
   A level 2 thread safe version could be built using dlarf (as in dorm2r)
   
*/

{ char side='L',trans='N';
  int lda,lwork=-1,info;
  double *work,work1;
 
  if (! *left) { side='R';lda = *c;} else lda= *r;
  if ( *tp) trans='T'; 
  /* workspace query */
  
  F77_CALL(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)CALLOC((size_t)lwork,sizeof(double));
  /* actual call */
  F77_CALL(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,work,&lwork,&info); 
  FREE(work);
}


void update_qr(double *Q,double *R,int *n, int *q,double *lam, int *k)
/* Let X=QR where X is n by q and R is q by q upper triangular and Q is n by q.
   A single element extra row, x, is to be appended to X and Q and R updated 
   accordingly. x is zero except for kth element lam.

   Let Q* be the full orthogonal matrix of which Q is the upper left portion, then 
   
   [X] = [Q* 0][R]
   [x]   [0  1][0]
               [x]

   The rhs of the above can be bought into standard QR form by application of givens rotations 
   from the left to the augmented R matrix and the corresponding inverse rotations from the right 
   to the augmented Q* matrix. The rotations from the right applied to the augmented Q* have the 
   effect of rotating columns of Q into the final column of the augmented matrix and vice-versa.

   Since the columns between Q and the final column are not involved, only Q and R need to be updated 
   here, the rest of Q* being irrelevant. This routine does not augment the Q by an extra row, it is 
   assumed that the calling function only requires the update of the input rows.

   All matrices are assumed to be packed in column order. Some very minor further optimizations could 
   be added (e.g. using fact that working and most of x are zero at first iteration), but it's really 
   unlikely to yield a worthwhile saving.

   Some R code for testing the routine:

library(mgcv)
n<-4;q<-3
X<-matrix(rnorm(n*q),n,q)
#X[,q-1]<-X[,q]
qrx<-qr(X,tol=0)
Q<-qr.Q(qrx);R<-qr.R(qrx);lam<-1
um<-.C("update_qr",as.double(Q),as.double(R),as.integer(n),as.integer(q),
        as.double(lam),as.integer(q-1),PACKAGE="mgcv")
R1<-matrix(um[[2]],q,q);Q1<-matrix(um[[1]],n,q)
Xa<-matrix(0,n+1,q)
Xa[1:n,]<-X;Xa[n+1,q]<-lam
qr.R(qr(Xa,tol=0))   

*/

{ double *x,*work,c,s,r,x0,x1,m,*xip,*xjp,*riip,*rijp,*Qp,*wp;
  x=(double *)CALLOC((size_t)*q,sizeof(double)); 
  work=(double *)CALLOC((size_t)*n,sizeof(double)); /* working extra column of Q */
  x[*k] = *lam;
  /* conceptually i runs from k to q in the following loop */
  for (Qp=Q+ *k * *n,riip=R+ *k * *q + *k,xip=x+ *k ;xip< x+ *q;xip++,riip+= *q+1)
  { /* rotate x[i] into R[i,i], using over/underflow proof rotator */
    x0= *xip; /* x[i] */
    x1= *riip; /* R[1 * *q + i] */
    m = fabs(x0);s=fabs(x1); if (s>m) m=s;
    x1/=m;x0/=m;
    r=sqrt(x0*x0+x1*x1);
    c=x1/r;s=x0/r;
    *riip=m*r;/* *xip=0.0; but never neaded*/    
    /* conceptually j runs from i+1 to q in the following loop */
    for (rijp=riip + *q,xjp=xip+1;xjp<x+ *q;xjp++,rijp+= *q) /* apply rotation to x[j], R[i,j] */
    { x1= *rijp;    /* R[j* *q+i] */
      /* *xjp == x[j] */
      *rijp = c*x1 - s* *xjp;
      *xjp = s*x1 + c* *xjp;
    }
    /* apply transpose of rotation from right to Q */
    /* conceptually j runs from 0 to n in following loop */
    for (wp=work;wp<work+ *n;wp++,Qp++)
    { x1 = *Qp;   /*Q[i* *n+j]*/
      /* *wp == work[j] */
      *Qp = c*x1-s* *wp;
      *wp = s*x1+c* *wp;
    }
  }
  FREE(x);FREE(work);
}

void mgcv_symeig(double *A,double *ev,int *n,int *use_dsyevd,int *get_vectors,
                 int *descending)

/* gets eigenvalues and (optionally) vectors of symmetric matrix A. 
   
   Two alternative underlying LAPACK routines can be used, 
   either dsyevd (slower, robust) or dsyevr (faster, seems less robust). 
   Vectors returned in columns of A, values in ev (ascending).
   
   Note: R 2.15.2 upgraded to a patched version of LAPACK 
         3.4.1 which seems to have broken uplo='U' - non-orthogonal 
         eigen-vectors are possible with that option.

   ******************************************************
   *** Eigenvalues are returned  in *ascending* order ***
   *** unless descending is set to be non-zero        ***
   ******************************************************

   Testing R code....
   library(mgcv);m <- 34
   n<-46;A<-matrix(rnorm(n*m),n,m);A<-A%*%t(A);
   er<-eigen(A)
   d<-array(0,ncol(A));n <- ncol(A)
   um<-.C("mgcv_symeig",as.double(A),as.double(d),as.integer(ncol(A)),
           as.integer(0),as.integer(1),as.integer(1),PACKAGE="mgcv")
   U <- matrix(um[[1]],n,n)
   er$vectors;U
   er$values;um[[2]]
*/  

{ char jobz='V',uplo='L',range='A'; 
  double work1,*work,dum1=0,abstol=0.0,*Z,*dum2,x,*p,*p1,*p2,*Acopy;
  int lwork = -1,liwork = -1,iwork1,info,*iwork,dumi=0,n_eval=0,*isupZ,i,j,k,debug=0;
  if (debug && *get_vectors) { /* need a copy to dump in case of trouble */
    Acopy = (double *)CALLOC((size_t)(*n * *n),sizeof(double));
    for (p2=Acopy,p=A,p1=A+ *n * *n;p<p1;p++,p2++) *p2 = *p;
  }

  if (*get_vectors) jobz='V'; else jobz='N';
  if (*use_dsyevd)
  { F77_CALL(dsyevd)(&jobz,&uplo,n,A,n,ev,&work1,&lwork,&iwork1,&liwork,&info);
    lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
    work=(double *)CALLOC((size_t)lwork,sizeof(double));
    liwork = iwork1;iwork= (int *)CALLOC((size_t)liwork,sizeof(int));
    F77_CALL(dsyevd)(&jobz,&uplo,n,A,n,ev,work,&lwork,iwork,&liwork,&info);
    FREE(work);FREE(iwork);
    if (*descending) for (i=0;i<*n/2;i++) {
        /* work in from left and right swapping cols */
	p = A + i * *n; /* start of left col */ 
        p1 = A + *n * (*n - 1 - i); /* start of right col */
        for (p2 = p + *n;p<p2;p++,p1++) { /* do the swap */
          x = *p;*p=*p1;*p1=x;
        }
    }
  } else
  { Z=(double *)CALLOC((size_t)(*n * *n),sizeof(double)); /* eigen-vector storage */
    isupZ=(int *)CALLOC((size_t)(2 * *n),sizeof(int)); /* eigen-vector support */
    F77_CALL(dsyevr)(&jobz,&range,&uplo,
		   n,A,n,&dum1,&dum1,&dumi,&dumi,
		   &abstol,&n_eval,ev, 
    		   Z,n,isupZ, &work1,&lwork,&iwork1,&liwork,&info);
    lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
    work=(double *)CALLOC((size_t)lwork,sizeof(double));
    liwork = iwork1;iwork= (int *)CALLOC((size_t)liwork,sizeof(int));
    F77_CALL(dsyevr)(&jobz,&range,&uplo,
		   n,A,n,&dum1,&dum1,&dumi,&dumi,
		   &abstol,&n_eval,ev, 
    		     Z,n,isupZ, work,&lwork,iwork,&liwork,&info);
    FREE(work);FREE(iwork);
    
    /* if (*descending) for (i=0;i<*n/2;i++) { 
      x = ev[i]; ev[i] = ev[*n-i-1];ev[*n-i-1] = x;
      } - now below*/

    if (*get_vectors) {  /* copy vectors back into A */
     p1 = A; 
     if (*descending) { /* need to reverse order */
        dum2 = Z + *n * (*n-1);
        for (work=dum2;work>=Z;work -= *n)
	  for (p=work;p<work + *n;p++,A++) *A = *p; 
      } else { /* simple copy */
        dum2 = Z + *n * *n;
        for (work=Z;work<dum2;work++,A++) *A = *work;
      }
    }
    A = p1;
    FREE(Z);FREE(isupZ);
  }
  if (*descending) for (i=0;i<*n/2;i++) { /* reverse the eigenvalues */
       x = ev[i]; ev[i] = ev[*n-i-1];ev[*n-i-1] = x;
  }
  if (debug && *get_vectors) { /* are the eigenvectors really orthogonal?? */
    p = (double *)CALLOC((size_t)(*n * *n),sizeof(double));
    getXtX(p,A,n,n); /* cross prod of eigenvec matrix - should be I */
    x=0.0;k=0;
    for (i=0;i<*n;i++) for (j=0;j<i;j++)  if (fabs(p[i + *n * j])>1e-14) { 
	  x += fabs(p[i + *n * j]);k++;  
    }
    Rprintf("**\n");
    j=k;
    if (k) Rprintf("Non orthogonal eigenvectors %d %g\n",k,x/k);
    x=0.0;k=0;
    for (i=0;i<*n;i++) if (fabs(p[i + *n * i]-1)>1e-14) { 
	x += fabs(p[i + *n * i]-1);k++;
    }
    if (k) Rprintf("Eigenvectors not normalized %d %g\n",k,x/k);
    if (k+j>0) dump_mat(Acopy,n,n,"/home/sw283/tmp/badmat.dat");
    FREE(p);FREE(Acopy);
  }

} /* mgcv_symeig */

void mgcv_trisymeig(double *d,double *g,double *v,int *n,int getvec,int descending) 
/* Find eigen-values and vectors of n by n symmetric tridiagonal matrix 
   with leading diagonal d and sub/super diagonals g. 
   eigenvalues returned in d, and eigenvectors in columns of v, if
   getvec!=0. If *descending!=0 then eigenvalues returned in descending order,
   otherwise ascending. eigen-vector order corresponds.  

   Routine is divide and conquer followed by inverse iteration. 

   dstevd could be used instead, with just a name change.
   dstevx may be faster, but needs argument changes.
*/ 
		    
{ char compz;
  double *work,work1,x,*dum1,*dum2;
  int ldz=0,info,lwork=-1,liwork=-1,*iwork,iwork1,i,j;

  if (getvec) { compz='I';ldz = *n;} else { compz='N';ldz=0;}

  /* workspace query first .... */
  F77_CALL(dstedc)(&compz,n,
		   d, g, /* lead and su-diag */
		   v, /* eigenvectors on exit */  
                   &ldz, /* dimension of v */
		   &work1, &lwork,
		   &iwork1, &liwork, &info);

   lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
   work=(double *)CALLOC((size_t)lwork,sizeof(double));
   liwork = iwork1;
   iwork= (int *)CALLOC((size_t)liwork,sizeof(int));

   /* and the actual call... */
   F77_CALL(dstedc)(&compz,n,
		   d, g, /* lead and su-diag */
		   v, /* eigenvectors on exit */  
                   &ldz, /* dimension of v */
		   work, &lwork,
		   iwork, &liwork, &info);

   if (descending) { /* need to reverse eigenvalues/vectors */
     for (i=0;i<*n/2;i++) { /* reverse the eigenvalues */
       x = d[i]; d[i] = d[*n-i-1];d[*n-i-1] = x;
       dum1 = v + *n * i;dum2 = v + *n * (*n-i-1); /* pointers to heads of cols to exchange */
       for (j=0;j<*n;j++,dum1++,dum2++) { /* work down columns */
         x = *dum1;*dum1 = *dum2;*dum2 = x;
       }
     }
   }

   FREE(work);FREE(iwork);
   *n=info; /* zero is success */
} /* mgcv_trisymeig */



void Rlanczos(double *A,double *U,double *D,int *n, int *m, int *lm,double *tol,int *nt) {
/* Faster version of lanczos_spd for calling from R.
   A is n by n symmetric matrix. Let k = m + max(0,lm).
   U is n by k and D is a k-vector.
   m is the number of upper eigenvalues required and lm the number of lower.
   If lm<0 then the m largest magnitude eigenvalues (and their eigenvectors)
   are returned 

   Matrices are stored in R (and LAPACK) format (1 column after another).

   If nt>1 and there is openMP support then the routine computes the O(n^2) inner products in 
   parallel.

   ISSUE: 1. Currently all eigenvectors of Tj are found, although only the next unconverged one
             is really needed. Might be better to be more selective using dstein from LAPACK. 
          2. Basing whole thing on dstevx might be faster
          3. Is random start vector really best? Actually Demmel (1997) suggests using a random vector, 
             to avoid any chance of orthogonality with an eigenvector!
          4. Could use selective orthogonalization, but cost of full orth is only 2nj, while n^2 of method is
             unavoidable, so probably not worth it.  
*/
  int biggest=0,f_check,i,k,kk,ok,l,j,vlength=0,ni,pi,converged,incx=1,ri,ci=0,cir,one=1;
  double **q,*v=NULL,bt,xx,yy,*a,*b,*d,*g,*z,*err,*p0,*p1,*zp,*qp,normTj,eps_stop,max_err,alpha=1.0,beta=0.0;
  unsigned long jran=1,ia=106,ic=1283,im=6075; /* simple RNG constants */
  const char uplo='U',trans='T';
  #ifdef OMP_REPORT
  Rprintf("Rlanczos");
  #endif


  #ifndef OPENMP_ON
  *nt = 1; /* reset number of threads to 1 if openMP not available  */ 
  #endif
  if (*nt > *n) *nt = *n; /* don't use more threads than columns! */

  eps_stop = *tol; 

  if (*lm<0) { biggest=1;*lm=0;} /* get m largest magnitude eigen-values */
  f_check = (*m + *lm)/2; /* how often to get eigen_decomp */
  if (f_check<10) f_check =10;
  kk = (int) floor(*n/10); if (kk<1) kk=1;  
  if (kk<f_check) f_check = kk;

  q=(double **)CALLOC((size_t)(*n+1),sizeof(double *));

  /* "randomly" initialize first q vector */
  q[0]=(double *)CALLOC((size_t)*n,sizeof(double));
  b=q[0];bt=0.0;
  for (i=0;i < *n;i++)   /* somewhat randomized start vector!  */
  { jran=(jran*ia+ic) % im;   /* quick and dirty generator to avoid too regular a starting vector */
    xx=(double) jran / (double) im -0.5; 
    b[i]=xx;xx=-xx;bt+=b[i]*b[i];
  } 
  bt=sqrt(bt); 
  for (i=0;i < *n;i++) b[i]/=bt;

  /* initialise vectors a and b - a is leading diagonal of T, b is sub/super diagonal */
  a=(double *)CALLOC((size_t) *n,sizeof(double));
  b=(double *)CALLOC((size_t) *n,sizeof(double));
  g=(double *)CALLOC((size_t) *n,sizeof(double));
  d=(double *)CALLOC((size_t) *n,sizeof(double)); 
  z=(double *)CALLOC((size_t) *n,sizeof(double));
  err=(double *)CALLOC((size_t) *n,sizeof(double));
  for (i=0;i< *n;i++) err[i]=1e300;
  /* figure out how to split up work if nt>1 */
  if (*nt>1) {
    ci = *n / *nt; /* cols per thread */
    cir = *n - ci * (*nt - 1); /* cols for final thread */  
    if (cir>ci) { /* final thread has more work than normal thread - redo */
      ci++; /* up work per thread by one */
      *nt = (int)ceil(*n/ci); /* drop number of threads */
      cir = *n - ci * (*nt - 1);  /* recompute cols for final thread */
    }
    if (cir == 0) { (*nt)--;cir=ci; } /* no cols left for final thread so drop it */
  }
  //Rprintf("nt = %d, ci = %d, cir = %d\n",*nt,ci,cir);
  /* The main loop. Will break out on convergence. */
  for (j=0;j< *n;j++) {
    /* form z=Aq[j]=A'q[j], the O(n^2) step ...  */
 
    /*blas free version ... for (Ap=A,zp=z,p0=zp+*n;zp<p0;zp++) 
      for (*zp=0.0,qp=q[j],p1=qp+*n;qp<p1;qp++,Ap++) *zp += *Ap * *qp;*/

    /*  BLAS versions y := alpha*A*x + beta*y, */
    if (*nt>1) { /* use parallel computation for the z = A q[j] */
      #ifdef OPENMP_ON
      #pragma omp parallel private(i,ri) num_threads(*nt)
      #endif
      { 
	#ifdef OPENMP_ON
	#pragma omp for
	#endif
        for (i=0;i<*nt;i++) {
          if (i < *nt-1) ri = ci; else ri = cir; /* number of cols of A to process */
          /* note that symmetry, A' == A, is exploited here, (rows a:b of A are same 
             as cols a:b of A, but latter are easier to access as a block) */  
          F77_CALL(dgemv)(&trans,n,&ri,&alpha,A+i * ci * *n,n,q[j],
                        &one,&beta,z+i*ci,&one);
        }
      } /* end parallel */
    } else F77_CALL(dsymv)(&uplo,n,&alpha,
		A,n,
		q[j],&incx,
		&beta,z,&incx);
    /* Now form a[j] = q[j]'z.... */
    for (xx=0.0,qp=q[j],p0=qp+*n,zp=z;qp<p0;qp++,zp++) xx += *qp * *zp;
    a[j] = xx;
 
    /* Update z..... */
    if (!j)
    { /* z <- z - a[j]*q[j] */
      for (zp=z,p0=zp+*n,qp=q[j];zp<p0;qp++,zp++) *zp -= xx * *qp;
    } else
    { /* z <- z - a[j]*q[j] - b[j-1]*q[j-1] */
      yy=b[j-1];      
      for (zp=z,p0=zp + *n,qp=q[j],p1=q[j-1];zp<p0;zp++,qp++,p1++) *zp -= xx * *qp + yy * *p1;    
  
      /* Now stabilize by full re-orthogonalization.... */
      
      for (i=0;i<=j;i++) 
      { /* form xx= z'q[i] */
        /*for (xx=0.0,qp=q[i],p0=qp + *n,zp=z;qp<p0;zp++,qp++) xx += *zp * *qp;*/
        xx = -F77_CALL(ddot)(n,z,&incx,q[i],&incx); /* BLAS version */
        /* z <- z - xx*q[i] */
        /*for (qp=q[i],zp=z;qp<p0;qp++,zp++) *zp -= xx * *qp;*/
        F77_CALL(daxpy)(n,&xx,q[i],&incx,z,&incx); /* BLAS version */
      } 
      
      /* exact repeat... */
      for (i=0;i<=j;i++) 
      { /* form xx= z'q[i] */
        /* for (xx=0.0,qp=q[i],p0=qp + *n,zp=z;qp<p0;zp++,qp++) xx += *zp * *qp; */
        xx = -F77_CALL(ddot)(n,z,&incx,q[i],&incx); /* BLAS version */
        /* z <- z - xx*q[i] */
        /* for (qp=q[i],zp=z;qp<p0;qp++,zp++) *zp -= xx * *qp; */
        F77_CALL(daxpy)(n,&xx,q[i],&incx,z,&incx); /* BLAS version */
      } 
      /* ... stabilized!! */
    } /* z update complete */
    

    /* calculate b[j]=||z||.... */
    for (xx=0.0,zp=z,p0=zp+*n;zp<p0;zp++) xx += *zp * *zp;
    b[j]=sqrt(xx); 
  
     /* get q[j+1]      */
    if (j < *n-1)
    { q[j+1]=(double *)CALLOC((size_t) *n,sizeof(double));
      for (xx=b[j],qp=q[j+1],p0=qp + *n,zp=z;qp<p0;qp++,zp++) *qp = *zp/xx; /* forming q[j+1]=z/b[j] */
    }

    /* Now get the spectral decomposition of T_j.  */

    if (((j>= *m + *lm)&&(j%f_check==0))||(j == *n-1))   /* no  point doing this too early or too often */
    { for (i=0;i<j+1;i++) d[i]=a[i]; /* copy leading diagonal of T_j */
      for (i=0;i<j;i++) g[i]=b[i]; /* copy sub/super diagonal of T_j */   
      /* set up storage for eigen vectors */
      if (vlength) FREE(v); /* free up first */
      vlength=j+1; 
      v = (double *)CALLOC((size_t)(vlength*vlength),sizeof(double));
 
      /* obtain eigen values/vectors of T_j in O(j^2) flops */
    
      kk = j + 1;
      mgcv_trisymeig(d,g,v,&kk,1,1);
      /* ... eigenvectors stored one after another in v, d[i] are eigenvalues */

      /* Evaluate ||Tj|| .... */
      normTj=fabs(d[0]);if (fabs(d[j])>normTj) normTj=fabs(d[j]);

      for (k=0;k<j+1;k++) /* calculate error in each eigenvalue d[i] */
      { err[k]=b[j]*v[k * vlength + j]; /* bound on kth e.v. is b[j]* (jth element of kth eigenvector) */
        err[k]=fabs(err[k]);
      }
      /* and check for termination ..... */
      if (j >= *m + *lm)
      { max_err=normTj*eps_stop;
        if (biggest) { /* getting m largest magnitude eigen values */
	  /* only one convergence test is sane here:
             1. Find the *m largest magnitude elements of d. (*lm is 0)
             2. When all these have converged, we are done.
          */   
          pi=ni=0;converged=1;
          while (pi+ni < *m) if (fabs(d[pi])>= fabs(d[j-ni])) { /* include d[pi] in largest set */
              if (err[pi]>max_err) {converged=0;break;} else pi++;
	    } else { /* include d[j-ni] in largest set */
              if (err[ni]>max_err) {converged=0;break;} else ni++;
            }
   
          if (converged) {
            *m = pi;
            *lm = ni;
            j++;break;
          }
        } else /* number of largest and smallest supplied */
        { ok=1;
          for (i=0;i < *m;i++) if (err[i]>max_err) ok=0;
          for (i=j;i > j - *lm;i--) if (err[i]>max_err) ok=0;
          if (ok) 
          { j++;break;}
        }
      }
    }
  }
  /* At this stage, complete construction of the eigen vectors etc. */
  
  /* Do final polishing of Ritz vectors and load va and V..... */
  /*  for (k=0;k < *m;k++) // create any necessary new Ritz vectors 
  { va->V[k]=d[k];
    for (i=0;i<n;i++) 
    { V->M[i][k]=0.0; for (l=0;l<j;l++) V->M[i][k]+=q[l][i]*v[k][l];}
  }*/

  /* assumption that U is zero on entry! */

  for (k=0;k < *m;k++) /* create any necessary new Ritz vectors */
  { D[k]=d[k];
    for (l=0;l<j;l++)
    for (xx=v[l + k * vlength],p0=U + k * *n,p1 = p0 + *n,qp=q[l];p0<p1;p0++,qp++) *p0 += *qp * xx;
  }

  for (k= *m;k < *lm + *m;k++) /* create any necessary new Ritz vectors */
  { kk=j-(*lm + *m - k); /* index for d and v */
    D[k]=d[kk];
    for (l=0;l<j;l++)
    for (xx=v[l + kk * vlength],p0=U + k * *n,p1 = p0 + *n,qp=q[l];p0<p1;p0++,qp++) *p0 += *qp * xx;
  }
 
  /* clean up..... */
  FREE(a);
  FREE(b);
  FREE(g);
  FREE(d);
  FREE(z);
  FREE(err);
  if (vlength) FREE(v);
  for (i=0;i< *n+1;i++) if (q[i]) FREE(q[i]);FREE(q);  
  *n = j; /* number of iterations taken */
  #ifdef OMP_REPORT
  Rprintf("done\n");
  #endif
} /* end of Rlanczos */


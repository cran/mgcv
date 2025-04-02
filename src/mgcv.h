/* main method routines */

/* See http://developer.r-project.org/blosxom.cgi/R-devel/2019/08/29#n2019-08-29
   for what USE_FC_LEN_T is doing and for why see
   https://developer.r-project.org/Blog/public/2019/05/15/gfortran-issues-with-lapack/index.html 
   
   In a nutshell, the mechanism used to call BLAS/LAPACK from C (by everyone, not just R) is not 
   technically supported by the Fortran standard. Fortran needs to know how long strings are (they 
   are not null terminated) so the lengths are passed as hidden extra function arguments. BLAS/LAPACK
   only ever uses single character strings, so it never needs to access the string lengths and it is 
   then no problem that they are missing (they are at the end of the call arguments), so they are simply 
   not passed in the C call. This was no problem until Gfortran decided to optimize the process of calling 
   a function with the same argument list as the calling function. Basically it passed the call stack of 
   the calling function to the called function assuming that it contained the string lengths - as it 
   didn't this caused stack corruption. 

   The solution is to pass the unit string lengths explicitly using FCONE defined in Rconfig.h if 
   USE_FC_LEN_T is defined. This mechanism is needed since it is compiler specific what type is used 
   to pass the string lengths (what happens then if BLAS/LAPACK and R are compiled using different 
   extra argument types is unclear to me, but no problems of this sort are currently known in any case 
   to get an actual problem the LAPACK/BLAS compiler would have to be using a different number of bytes 
   to the R compiler). 

   In practice when calling BLAS/LAPACK macro FCONE has to be added to the end of the call as
   many times as there are character arguments to the call. mat.c has many examples.

*/

#define USE_FC_LEN_T
#include <Rinternals.h>
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
/* If we are compiling with a version of R before FCONE and the explicit supplying of extra arguments 
   was introduced, then FCONE has to be defined */ 
#ifndef FCONE
#define FCONE
#endif

/* Most compilers with openMP support supply 
   a pre-defined compiler macro _OPENMP. Following 
   facilitates selective turning off (by testing value 
   or defining multiple versions OPENMP_ON1, OPENMP_ON2...)  */

#if defined _OPENMP
#define OPENMP_ON 1 
#endif
/* ... note also that there is no actual *need* to protect #pragmas with 
  #ifdef OPENMP_ON, since C ignores undefined pragmas, but failing 
  to do so may produce alot of compilation warnings if openMP is not supported. 
  In contrast functions from omp.h must be protected, and there is 
  non-avoidable use of these in the mgcv code. */

//#define OMP_REPORT // define to have all routines using omp report on start and end.

/* sed -i 's/old-text/new-text/g' *.c
   is quite useful!!
*/

/* For safe memory handling from R... */
#define CALLOC R_chk_calloc
#define FREE R_chk_free
#define REALLOC R_chk_realloc

/* BUT, this can mess up valgrinding for memory error checking - problems are 
   sometimes missed because standard allocation is being circumvented. Then errors can 
   corrupt R memory management without detection and trigger nothing until R
   messes up internally because of corruption, which then makes it look as if
   R is generating the problem. Hence better to reset for checking. Also sizing
   errors in .C often generate no obvious valgrind error.*/
//#define CALLOC calloc
//#define FREE free
//#define REALLOC realloc

void *R_chk_calloc1(size_t nmemb,size_t size);

void magic(double *y,double *X,double *sp0,double *def_sp,double *S,double *H,double *L,
	   double *lsp0,double *gamma,double *scale, int *control,int *cS,double *rank_tol,
	   double *tol,double *b,double *rV,double *norm_const,int *n_score,int *nt);

void gdi1(double *X,double *E,double *Es,double *rS,double *U1,
	  double *sp,double *z,double *w,double *wf,double *alpha,double *mu,double *eta, double *y,
	 double *p_weights,double *g1,double *g2,double *g3,double *g4,double *V0,
	  double *V1,double *V2,double *V3,double *beta,double *b1,double *w1,double *D1,double *D2,
         double *P0, double *P1,double *P2,double *trA,
         double *trA1,double *trA2,double *rV,double *rank_tol,double *conv_tol, int *rank_est,
	 int *n,int *q, int *M,int *Mp,int *Enrow,int *rSncol,int *deriv,
	  int *REML,int *fisher,int *fixed_penalty,int *nthreads,double *dVkk);     

void gdi2(double *X,double *E,double *Es,double *rS,double *U1,
	  double *sp,double *theta,double *z,double *w,double *wz,double *wf,
          double *Dth,double *Det,double *Det2,double *Dth2,double *Det_th,
          double *Det2_th,double *Det3,double *Det_th2,
          double *Det4, double *Det3_th, double *Det2_th2,
          double *beta,double *b1,double *w1,double *D1,double *D2,double *P,double *P1,double *P2,
          double *ldet, double *ldet1,double *ldet2,double *rV,
          double *rank_tol,int *rank_est,
	  int *n,int *q, int *M,int *n_theta, int *Mp,int *Enrow,int *rSncol,int *deriv,
	  int *fixed_penalty,int *nt,int *type,double *dVkk);

void pls_fit1(double *y,double *X,double *w,double *wy,double *E,double *Es,int *n,int *q,int *rE,double *eta,
	      double *penalty,double *rank_tol,int *nt,int *use_wy);

void get_detS2(double *sp,double *sqrtS, int *rSncol, int *q,int *M, int * deriv, 
               double *det, double *det1, double *det2, double *d_tol,
               double *r_tol,int *fixed_penalty); /* stable determinant of sum evaluation */

void get_stableS(double *S,double *Qf,double *sp,double *sqrtS, int *rSncol, int *q,int *M, int * deriv, 
               double *det, double *det1, double *det2, double *d_tol,
		 double *r_tol,int *fixed_penalty);

/* cox model routines */

void coxpred(double *X,double *t,double *beta,double *off,double *Vb,double *a,double *h,double *q,
             double *tr,int *n,int *p, int *nt,double *s,double *se);
void coxpp(double *eta,double *X,int *r, int *d,double *h,double *q,double *km,
	   int *n,int *p, int *nt);
void coxlpl(double *eta,double *X,int *r, int *d,double *tr, 
            int *n,int *p, int *nt,double *lp,double *g,double *H,
            double *d1beta,double *d1H,double *d2beta,
            double *d2H,int *n_sp,int *deriv);

/* MVN smooth additive */
/*void mvn_ll(double *y,double *X,double *XX,double *beta,int *n,int *lpi,
            int *m,double *ll,double *lb,double *lbb,double *dbeta,
            double *dH,int *deriv,int *nsp,int *nt);*/

SEXP mvnll(SEXP Y,SEXP x,SEXP xx,SEXP BETA,SEXP LPI, SEXP LL, SEXP LB,
	   SEXP LBB, SEXP DBETA, SEXP Dh, SEXP DERIV,SEXP NSP, SEXP NT);

/* discretized covariate methods */
//void XWXd(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, int *n, int *nx, 
//          int *ts, int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,
//          int *ar_row,double *ar_weights);
void XWXd0(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, ptrdiff_t *n, int *nx, 
	   int *ts, int *dt, int *nt,double *v,int *qc,int *nthreads,
          int *ar_stop,double *ar_weights);
SEXP CXWXd0(SEXP XWXr, SEXP Xr, SEXP wr, SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP tsr, SEXP dtr,
	    SEXP vr,SEXP qcr, SEXP nthreadsr, SEXP ar_stopr, SEXP ar_weightsr);
SEXP CXVXd0(SEXP XWXr, SEXP Xr, SEXP er, SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP tsr, SEXP dtr,
	    SEXP vr,SEXP qcr, SEXP nthreadsr, SEXP ar, SEXP mar);
void XWXd1(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, ptrdiff_t *n, int *nx, int *ts, 
	   int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,double *ar_weights,
	   int *rs, int *cs, int *nrs, int *ncs);
SEXP CXWXd1(SEXP XWXr, SEXP Xr, SEXP wr, SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP tsr, SEXP dtr,
	    SEXP vr,SEXP qcr, SEXP nthreadsr, SEXP ar_stopr, SEXP ar_weightsr,
	    SEXP rsr, SEXP csr);
void XWyd(double *XWy,double *y,double *X,double *w,int *k, int *ks, int *m,int *p, ptrdiff_t *n,int *cy, 
	  int *nx, int *ts, int *dt, int *nt,double *v,int *qc,
          int *ar_stop,int *ar_row,double *ar_weights,int *cs,int *ncs);
SEXP CXWyd(SEXP XWyr, SEXP yr, SEXP Xr, SEXP wr, SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP cyr, SEXP tsr,
	   SEXP dtr,SEXP vr,SEXP qcr, SEXP ar_stopr, SEXP ar_rowr, SEXP ar_weightsr,SEXP csr);
void Xbdspace(ptrdiff_t *space,int *m,int *p, ptrdiff_t *n, int *nx, int *dt, int *nt);
void Xbd(double *f,double *beta,double *X,int *k, int *ks, int *m,int *p, ptrdiff_t  *n, 
	 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *bc,int *cs,int *ncs,
	 int *iwork, ptrdiff_t *pwork,double *dwork);
SEXP CXbd(SEXP fr, SEXP betar, SEXP Xr, SEXP kr, SEXP ksr, SEXP mr, SEXP pr,
	  SEXP tsr, SEXP dtr,SEXP vr,SEXP qcr,SEXP bcr,SEXP csr);
void diagXVXt(double *diag,double *V,double *X,int *k1,int *k2,int *ks,int *m,int *p, ptrdiff_t *n, 
	      int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pv,int *cv,int *nthreads,int *cs,int *ncs,int *rs,int *nrs);
void idiagXLLtXt(double *diag,double *L,double *X,int *k,int *ks,int *m,int *p, int *n, 
		 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pl,int *cl,
		 int *ri,int *ci,int *nrc,int *nthreads);
void diagXLLtXt(double *diag,double *L,double *X,int *k,int *ks,int *m,int *p, ptrdiff_t *n, 
	        int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pl,int *cl,
		int *ri,int *ci,ptrdiff_t *nrc,int *nthreads);
void idiagXLUtXt(double *diag,double *L,double *U,double *X,int *k,int *ks,int *m,int *p, int *n, 
		 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pl,int *cl,
		 int *ri,int *ci,int *nrc,int *nthreads);
void diagXLUtXt(double *diag,double *L,double *U,double *X,int *k,int *ks,int *m,int *p, ptrdiff_t *n, 
	      int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pl,int *cl,
		int *ri,int *ci,ptrdiff_t *nrc,int *nthreads);
SEXP CdiagXVXt(SEXP DIAG, SEXP Vp, SEXP x, SEXP K, SEXP KS, SEXP M, SEXP P, SEXP TS, SEXP DT,
	       SEXP vp,SEXP QC, SEXP NTHREADS, SEXP CS, SEXP RS);
SEXP CijXVXt(SEXP DIAG, SEXP Vp, SEXP x, SEXP K,SEXP K1, SEXP KS, SEXP M, SEXP P, SEXP TS, SEXP DT,
	     SEXP vp,SEXP QC, SEXP NTHREADS, SEXP CS, SEXP RS);
SEXP CNCV(SEXP NCVr, SEXP NCV1r, SEXP NCV2r, SEXP Gr, SEXP rsdr, SEXP betar, SEXP beta1r, SEXP wr,SEXP spr,
	  SEXP Xr,SEXP kr, SEXP ksr, SEXP mr, SEXP pr, SEXP tsr, SEXP dtr,
	  SEXP vr,SEXP qcr, SEXP nthreadsr, SEXP nei, SEXP Sr);

/* various service routines */
void davies(double *lb,double *nc,int *n,int *r,double *sigma,double *c,int *lim,
	    double *acc,double *trace,int *ifault);
void tweedious(double *w,double *w1,double *w2, double *w1p,double *w2p,double *w2pp, 
	       double *y,double *eps,int *n,
               double *th,double *rho,double *a, double *b);
void tweedious2(double *w,double *w1,double *w2, double *w1p,double *w2p,double *w2pp, 
	       double *y,double *eps,int *n,
               double *th,double *rho,double *a, double *b);
void psum(double *y, double *x,int *index,int *n);
void rwMatrix(int *stop,int *row,double *w,double *X,int *n,int *p,int *trans,double *work);
void in_out(double *bx, double *by, double *break_code, double *x,double *y,int *in, int *nb, int *n);
void Rlanczos(double *A,double *U,double *D,int *n, int *m, int *lm,double *tol,int *nt);
void RuniqueCombs(double *X,int *ind,int *r, int *c);
void  RPCLS(double *Xd,double *pd,double *yd, double *wd,double *Aind,double *bd,double *Afd,double *Sd,int *off,int *dim,double *theta, int *m,int *nar,int *active);
void RMonoCon(double *Ad,double *bd,double *xd,int *control,double *lower,double *upper,int *n);
void MinimumSeparation(double *x,int *n, int *d,double *t,int *m,double *dist);
void rksos(double *x,int *n,double *eps);
void pivoter(double *x,int *r,int *c,int *pivot, int *col, int *reverse);

/* Routines for linear algebra with direct access to linpack and lapack */
void row_squash(double *X,int rnew,int rold,int col);
void up2lo(double * A, int n);
void band_chol(double *B,int *n,int *k,int *info);
void tri_chol(double *ld,double *sd,int *n,int *info);
void mgcv_omp(int *a);
void mgcv_chol(double *a,int *pivot,int *n,int *rank);
void mgcv_qrqy(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp);
void mgcv_qrqy0(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp);
void mgcv_backsolve(double *R,int *r,int *c,double *B,double *C, int *bc, int *right);
void mgcv_forwardsolve(double *R,int *r,int *c,double *B,double *C, int *bc, int *right);
void mgcv_qr(double *x, int *r, int *c,int *pivot,double *tau);
void mgcv_qr2(double *x, int *r, int *c,int *pivot,double *tau);
void update_qr(double *Q,double *R,int *n, int *q,double *lam, int *k);
extern void mgcv_mmult(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n);
void mgcv_pmmult(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n,int *nt);
SEXP mgcv_pmmult2(SEXP b, SEXP c,SEXP bt,SEXP ct, SEXP nthreads);
void mgcv_mmult0(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n);
void mgcv_svd_full(double *x,double *vt,double *d,int *r,int *c);
void mgcv_symeig(double *A,double *ev,int *n,int *use_dsyevd, int *get_vectors,int *descending);
void mroot(double *A,int *rank,int *n);
void R_cond(double *R,int *r,int *c,double *work,double *Rcondition);
void mgcv_td_qy(double *S,double *tau,int *m,int *n, double *B,int *left,int *transpose);
void mgcv_tri_diag(double *S,int *n,double *tau);
void mgcv_trisymeig(double *d,double *g,double *v,int *n,int *getvec,int *descending);
void mtrf(double *A, double *B,int *n,int *rank,int *psd,double *tol,double *work,int *iwork);
void getXtWX(double *XtWX, double *X,double *w,int *r,int *c,double *work);
void getXtX(double *XtX,double *X,int *r,int *c);
void getXtMX(double *XtMX,double *X,double *M,int *r,int *c,double *work);
void getXXt(double *XXt,double *X,int *r,int *c);
void read_mat(double *M,int *r,int*c, char *path);
void row_block_reorder(double *x,int *r,int *c,int *nb,int *reverse);
void mgcv_pqr(double *x,int *r, int *c,int *pivot, double *tau, int *nt);
void getRpqr(double *R,double *x,int *r, int *c,int *rr,int *nt);
void mgcv_pqrqy(double *b,double *a,double *tau,int *r,int *c,int *cb,int *tp,int *nt);
SEXP wdiag(SEXP a,SEXP IND,SEXP B);
SEXP dpdev(SEXP a);
SEXP mgcv_Rpiqr(SEXP X, SEXP BETA,SEXP PIV,SEXP NT,SEXP NB);
SEXP mgcv_tmm(SEXP x,SEXP t,SEXP D,SEXP M, SEXP N);
SEXP mgcv_Rpbsi(SEXP A, SEXP NT);
SEXP mgcv_RPPt(SEXP a,SEXP r, SEXP NT);
SEXP mgcv_Rpchol(SEXP Amat,SEXP PIV,SEXP NT,SEXP NB);
void dchol(double *dA, double *R, double *dR,int *p);
void chol_down(double *R,double *Rup,int *n,int *k,int *ut);
SEXP mgcv_chol_down(SEXP r,SEXP ru,SEXP N,SEXP K, SEXP UT);
SEXP mgcv_chol_up(SEXP r,SEXP U,SEXP N,SEXP UP,SEXP EPS);
void vcorr(double *dR,double *Vr,double *Vb,int *p,int *M);
SEXP mgcv_Rpforwardsolve(SEXP R, SEXP B,SEXP NT);
SEXP mgcv_Rpbacksolve(SEXP R, SEXP B,SEXP NT);
SEXP mgcv_Rpcross(SEXP A, SEXP NT,SEXP NB);
SEXP mgcv_madi(SEXP a, SEXP b,SEXP ind,SEXP diag);
SEXP mrow_sum(SEXP x,SEXP M, SEXP K);
SEXP ncv(SEXP x, SEXP hi, SEXP W1, SEXP W2, SEXP DB, SEXP DW, SEXP rS, SEXP IND, SEXP MI,SEXP M,
	 SEXP K, SEXP BETA, SEXP SP, SEXP ETA, SEXP DETA,SEXP DLET,SEXP DERIV);
SEXP Rncv(SEXP x, SEXP r, SEXP W1, SEXP W2, SEXP DB, SEXP DW, SEXP rS, SEXP IND, SEXP MI, SEXP M, SEXP K,
	  SEXP BETA, SEXP SP, SEXP ETA,SEXP DETA,SEXP DLET,SEXP DERIV,SEXP EPS,SEXP NT);
SEXP ncvls(SEXP x,SEXP JJ,SEXP h,SEXP hi,SEXP dH,SEXP L1, SEXP L2,SEXP L3,SEXP IND, SEXP MI, SEXP M, SEXP K,SEXP BETA,
	   SEXP ETACV,SEXP DETACV,SEXP DETA,SEXP DB,SEXP DERIV);
SEXP Rncvls(SEXP x,SEXP JJ,SEXP R1,SEXP dH,SEXP L1, SEXP L2,SEXP L3,SEXP IND, SEXP MI, SEXP M, SEXP K,SEXP BETA,
	    SEXP ETACV,SEXP DETACV,SEXP DETA,SEXP DB,SEXP DERIV,SEXP EPS,SEXP NT);
SEXP nei_cov(SEXP v,SEXP d, SEXP d1, SEXP M, SEXP K);
void ncvd(double *NCV,double *NCV1,double *NCV2,double *beta,double *db, double *G,double *rsd, double *w,int *pg,
	  int *nn,int *a,ptrdiff_t *ma,int *d,ptrdiff_t *md,double *X,int *k,int *ck, int *ks,int *m,int *p, ptrdiff_t *n,
	  double **S,
	  int ns,int *sr,int *soff,double *sp,int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *nthreads);
void chol_up(double *R,double *u, int *n,int *up,double *eps);
void minres(double *R, double *u,double *b, double *x, int *p,int *m,double *work);
SEXP QRdrop(SEXP q,SEXP r,SEXP K);
SEXP QRadd(SEXP q,SEXP r,SEXP A);


/* sparse matrix routines */
SEXP isa1p(SEXP L,SEXP S,SEXP NT);
SEXP stmm(SEXP X); /* row Kronecker product */
SEXP sdiagXVXt(SEXP X, SEXP V, SEXP LT, SEXP RT);
SEXP sXbd(SEXP X,SEXP BETA,SEXP LT);
SEXP sXyd(SEXP X,SEXP Y,SEXP LT);
SEXP sXWXd(SEXP X,SEXP W,SEXP LT, SEXP RT,SEXP NT);
SEXP AddBVB(SEXP A,SEXP bt, SEXP vbt);
SEXP spdev(SEXP A);
SEXP getListEl(SEXP list, const char *str);



/* basis constructor/prediction routines*/

void crspl(double *x,int *n,double *xk, int *nk,double *X,double *S, double *F,int *Fsupplied);
void predict_tprs(double *x, int *d,int *n,int *m,int *k,int *M,double *Xu,int *nXu,
                  double *UZ,double *by,int *by_exists,double *X);
void construct_tprs(double *x,int *d,int *n,double *knt,int *nk,int *m,int *k,double *X,double *S,
                    double *UZ,double *Xu,int *nXu,double *C);
void gen_tps_poly_powers(int *pi,int *M,int *m, int *d);
void boundary(int *G, double *d, double *dto, double *x0, double *y0, double *dx, double *dy,
              int *nx, int *ny, double *x, double *y,double *break_code, int *n, int *nb);
void gridder(double *z,double *x,double *y,int *n,double *g, int *G,int *nx, int *ny,double *x0, 
             double *y0,double *dx,double *dy,double NA_code);
void pde_coeffs(int *G,double *x,int *ii,int *jj,int *n,int *nx,int *ny,double *dx,double *dy);

/* sparse smooth related routines */
typedef struct { /* defines structure for kd-tree box */
  double *lo,*hi;    /* box defining co-ordinates */
  int parent,child1,child2, /* indices of parent and 2 offspring */
      p0,p1;         /* indices of first and last point in box */
} box_type; 



typedef struct {
  box_type *box;
  int *ind, /* index of points in coordinate matrix which tree relates to */
      *rind, /* where is ith row of X in ind? */
      n_box, /* number of boxes */
      d, /* dimension */
    n; /* number of points that tree relates to */
  double huge; /* number indicating an open boundary */
} kdtree_type;

void k_newn_work(double *Xm,kdtree_type kd,double *X,double *dist,int *ni,int*m,int *n,int *d,int *k);
void k_nn(double *X,double *dist,double *a,int *ni,int *n,int *d,int *k,int *get_a);
//void Rkdtree(double *X,int *n, int *d,int *idat,double *ddat);
SEXP Rkdtree(SEXP x);
//void Rkdnearest(double *X,int *idat,double *ddat,int *n,double *x, int *m, int *ni, double *dist,int *k);
SEXP Rkdnearest(SEXP kdr,SEXP Xr, SEXP xr,SEXP k);
//void Rkradius(double *r,int *idat,double *ddat,double *X,double *x,int *m,int *off,int *ni,int *op);
SEXP Rkradius(SEXP kdr,SEXP Xr, SEXP xr,SEXP rr,SEXP offr);
double xidist(double *x,double *X,int i,int d, int n);
int closest(kdtree_type *kd, double *X,double *x,int n,int *ex,int nex);
void kd_tree(double *X,int *n, int *d,kdtree_type *kd);
void free_kdtree(kdtree_type kd);

void tri2nei(int *t,int *nt,int *n,int *d,int *off);
void nei_penalty(double *X,int *n,int *d,double *D,int *ni,int *ii,int *off,
		 int *m,int *a_weight,double *kappa);
void sspl_construct(double *lambda,double *x,double *w,double *U,double *V,
             double *diagA,double *lb,int *n,double *tol);
void sspl_mapply(double *y,double *x,double *w,double *U,double *V,int *n,int *nf,double *tol,int *m);

/* just for testing */
void Zb(double *b1,double *b0,double *v,int *qc, int *p,double *w);
void Ztb(double *b1,double *b0,double *v,int *qc,int *di, int *p,double *w);
  

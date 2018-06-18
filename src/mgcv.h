/* main method routines */
#include <Rinternals.h>
#include <Rconfig.h>
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
// For safe memory handling from R...
#define CALLOC R_chk_calloc
#define FREE R_chk_free
// Can reset to check for memory errors...
//#define CALLOC calloc
//#define FREE free
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
void mvn_ll(double *y,double *X,double *XX,double *beta,int *n,int *lpi,
            int *m,double *ll,double *lb,double *lbb,double *dbeta,
            double *dH,int *deriv,int *nsp,int *nt);

/* discretized covariate methods */
void XWXd(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, int *n, int *nx, 
          int *ts, int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,
          int *ar_row,double *ar_weights);
void XWyd(double *XWy,double *y,double *X,double *w,int *k, int *ks, int *m,int *p, int *n, 
	  int *nx, int *ts, int *dt, int *nt,double *v,int *qc,
          int *ar_stop,int *ar_row,double *ar_weights);
void Xbd(double *f,double *beta,double *X,int *k, int *ks, int *m,int *p, int *n, 
	 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *bc);
void diagXVXt(double *diag,double *V,double *X,int *k,int *ks,int *m,int *p, int *n, 
	      int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pv,int *nthreads);

/* various service routines */

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
void  RPCLS(double *Xd,double *pd,double *yd, double *wd,double *Aind,double *bd,double *Afd,double *Sd,int *off,int *dim,double *theta, int *m,int *nar);
void RMonoCon(double *Ad,double *bd,double *xd,int *control,double *lower,double *upper,int *n);
/*void MinimumSeparation(double *gx,double *gy,int *gn,double *dx,double *dy, int *dn,double *dist);*/
void MinimumSeparation(double *x,int *n, int *d,double *t,int *m,double *dist);
void rksos(double *x,int *n,double *eps);
void pivoter(double *x,int *r,int *c,int *pivot, int *col, int *reverse);

/* Routines for linear algebra with direct access to linpack and lapack */ 
void band_chol(double *B,int *n,int *k,int *info);
void tri_chol(double *ld,double *sd,int *n,int *info);
void mgcv_omp(int *a);
void mgcv_chol(double *a,int *pivot,int *n,int *rank);
void mgcv_svd(double *x,double *u, double *d,int *r,int *c);
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
void mgcv_trisymeig(double *d,double *g,double *v,int *n,int getvec,int descending); 
void getXtWX(double *XtWX, double *X,double *w,int *r,int *c,double *work);
void getXtX(double *XtX,double *X,int *r,int *c);
void getXtMX(double *XtMX,double *X,double *M,int *r,int *c,double *work);
void getXXt(double *XXt,double *X,int *r,int *c);
void read_mat(double *M,int *r,int*c, char *path);
void row_block_reorder(double *x,int *r,int *c,int *nb,int *reverse);
void mgcv_pqr(double *x,int *r, int *c,int *pivot, double *tau, int *nt);
void getRpqr(double *R,double *x,int *r, int *c,int *rr,int *nt);
void mgcv_pqrqy(double *b,double *a,double *tau,int *r,int *c,int *cb,int *tp,int *nt);
SEXP mgcv_Rpiqr(SEXP X, SEXP BETA,SEXP PIV,SEXP NT,SEXP NB);
void mgcv_tmm(SEXP x,SEXP t,SEXP D,SEXP M, SEXP N);
void mgcv_Rpbsi(SEXP A, SEXP NT);
void mgcv_RPPt(SEXP a,SEXP r, SEXP NT);
SEXP mgcv_Rpchol(SEXP Amat,SEXP PIV,SEXP NT,SEXP NB);
void dchol(double *dA, double *R, double *dR,int *p);
void vcorr(double *dR,double *Vr,double *Vb,int *p,int *M);
SEXP mgcv_Rpforwardsolve(SEXP R, SEXP B,SEXP NT);
SEXP mgcv_Rpbacksolve(SEXP R, SEXP B,SEXP NT);
SEXP mgcv_Rpcross(SEXP A, SEXP NT,SEXP NB);
SEXP mgcv_madi(SEXP a, SEXP b,SEXP ind,SEXP diag);


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



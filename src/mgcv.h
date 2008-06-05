/* main method routines */
void mgcv(double *yd,double *Xd,double *Cd,double *wd,double *Sd,double *pd, double *sp,
          int *offd,int *dimd,int *md,int *nd,int *qd,int *rd,double *sig2d, double *Vpd, 
          double *edf, double *conv_tol, int *ms_max_half,double *ddiag,int *idiag,double *sdiag, 
          int *direct_mesh,double *min_edf,double *gcvubre,double *target_edf,int *fixed_sp,double *hat);

void update_beta(double *X,double *Sr,double *rS,double *theta,double *w,
		 double *w1, double *z,double *z1,int *Srncol,int *rSncol,
		 int *m, int *n,int *q, int *get_trA,int *deriv,
		 double *rank_tol,double *beta, double *trA, double *beta1,
		 double *trA1,double *rV,int *rank_est);

void magic(double *y,double *X,double *sp0,double *def_sp,double *S,double *H,double *L,
	   double *lsp0,double *gamma,double *scale, int *control,int *cS,double *rank_tol,
	   double *tol,double *b,double *rV,double *norm_const,int *n_score);

void gdi(double *X,double *E,double *rS,
         double *sp,double *z,double *w,double *mu, double *eta, double *y,
         double *p_weights,double *g1,double *g2,double *g3,double *V0,
         double *V1,double *V2,double *beta,double *D1,double *D2,
         double *P0,double *P1,double *P2,double *trA,
         double *trA1,double *trA2,double *rV,double *rank_tol,double *conv_tol, int *rank_est,
	 int *n,int *q, int *M,int *Encol,int *rSncol,int *deriv,int *use_svd);

void pls_fit(double *y,double *X,double *w,double *E,int *n,int *q,int *cE,double *eta,
             double *penalty,double *rank_tol);

/* various service routines */
void RQT(double *A,int *r,int *c);
void RuniqueCombs(double *X,int *ind,int *r, int *c);
void  RPCLS(double *Xd,double *pd,double *yd, double *wd,double *Aind,double *bd,double *Afd,double *Hd,double *Sd,int *off,int *dim,double *theta, int *m,int *nar);
void RMonoCon(double *Ad,double *bd,double *xd,int *control,double *lower,double *upper,int *n);
void mgcv_AtA(double *AA,double *A,int *q,int *n);
void MinimumSeparation(double *gx,double *gy,int *gn,double *dx,double *dy, int *dn,double *dist);
/* Test routines for direct access to linpack and lapack */
void mgcv_chol(double *a,int *pivot,int *n,int *rank);
void mgcv_svd(double *x,double *u, double *d,int *r,int *c);
void mgcv_qrqy(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp);
void mgcv_qr(double *x, int *r, int *c,int *pivot,double *tau);
void update_qr(double *Q,double *R,int *n, int *q,double *lam, int *k);
void mgcv_mmult(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n);
void mgcv_svd_full(double *x,double *vt,double *d,int *r,int *c);
void mgcv_symeig(double *A,double *ev,int *n,int *use_dsyevd);
void mroot(double *A,int *rank,int *n);
void R_cond(double *R,int *r,int *c,double *work,double *Rcondition);


/* basis constructor/prediction routines*/

void construct_cr(double *x,int *nx,double *k,int *nk,double *X,double *S,double *C,int *control);
void predict_tprs(double *x, int *d,int *n,int *m,int *k,int *M,double *Xu,int *nXu,
                  double *UZ,double *by,int *by_exists,double *X);
void construct_tprs(double *x,int *d,int *n,double *knt,int *nk,int *m,int *k,double *X,double *S,
                    double *UZ,double *Xu,int *nXu,double *C);

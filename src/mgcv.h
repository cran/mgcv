void mgcv(double *yd,double *Xd,double *Cd,double *wd,double *Sd,double *pd, double *sp,
          int *offd,int *dimd,int *md,int *nd,int *qd,int *rd,double *sig2d, double *Vpd, 
          double *edf, double *conv_tol, int *ms_max_half,double *ddiag,int *idiag,double *sdiag, 
          int *direct_mesh,double *min_edf,double *gcvubre,double *target_edf,int *fixed_sp,double *hat);
void RGAMsetup(double *Xd,double *Cd,double *Sd,double *UZd,double *Xud,int *xu,double *xpd,
             int *offd,double *xd,int  *md,int  *nd,int *dfd,int *nsdfd,int *dim,int *s_type, int *p_order,
             int *by_exists, double *byd,double *knots, int *n_knots);
void RGAMpredict(int *xu,double *Xud,double *UZd,double *xpd,int *nsdf,int *dim,int *s_type,int *df,
                 int *p_order,int *m,int *n,double *xd,int *np,double *p, double *Vpd,double *etad,
                 double *sed,double *X,int *control,int *by_exists, double *by);
void RQT(double *A,int *r,int *c);
void RuniqueCombs(double *X,int *r, int *c);
void  RPCLS(double *Xd,double *pd,double *yd, double *wd,double *Aind,double *bd,double *Afd,double *Hd,double *Sd,
	    int *off,int *dim,double *theta, int *m,int *nar);
void RMonoCon(double *Ad,double *bd,double *xd,int *control,double *lower,double *upper,int *n);
void mgcv_AtA(double *AA,double *A,int *q,int *n);
/* Test routines for direct access to linpack and lapack */
void mgcv_chol(double *a,int *pivot,int *n,int *rank);
void mgcv_svd(double *x,double *u, double *d,int *r,int *c);
void mgcv_qrqy(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp);
void mgcv_qr(double *x, int *r, int *c,int *pivot,double *tau);
void update_qr(double *Q,double *R,int *n, int *q,double *lam, int *k);
void mgcv_mmult(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n);
void mgcv_svd_full(double *x,double *vt,double *d,int *r,int *c);
void mgcv_symeig(double *A,double *ev,int *n);
void mroot(double *A,int *rank,int *n);


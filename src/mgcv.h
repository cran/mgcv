void mgcv(double *yd,double *Xd,double *Cd,double *wd,double *Sd,
     double *pd, double *sp,int *offd,int *dimd,int *md,
     int *nd,int *qd,int *rd,double *sig2d, double *Vpd, double *edf, double *conv_tol, int *ms_max_half);
void RGAMsetup(double *Xd,double *Cd,double *Sd,double *UZd,double *Xud,int *xu,double *xpd,
             int *offd,double *xd,int  *md,int  *nd,int *dfd,int *nsdfd,int *dim,int *s_type, int *p_order);
void RGAMpredict(int *xu,double *Xud,double *UZd,double *xpd,int *nsdf,int *dim,int *s_type,int *df,
                 int *p_order,int *m,int *n,double *xd,int *np,double *p, double *Vpd,double *etad,
                 double *sed,double *X,int *control);
void RQT(double *A,int *r,int *c);
void RuniqueCombs(double *X,int *r, int *c);
void  RPCLS(double *Xd,double *pd,double *yd, double *wd,double *Aind,double *bd,double *Afd,double *Hd,double *Sd,
	    int *off,int *dim,double *theta, int *m,int *nar);
void RMonoCon(double *Ad,double *bd,double *xd,int *control,double *lower,double *upper,int *n);



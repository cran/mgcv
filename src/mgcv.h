void mgcv(double *yd,double *Xd,double *Cd,double *wd,double *Sd,
     double *pd, double *sp,int *offd,int *dimd,int *md,
     int *nd,int *qd,int *rd,double *sig2d, double *Vpd, double *edf);
void RGAMsetup(double *Xd,double *Cd,double *Sd,double *xpd,
             int *offd,double *xd,int  *md,int  *nd,int *dfd,int *nsdfd);
void RGAMpredict(double *xpd,int *nsdf,int *df,int *m,double *xd,int *np,double *p, double *Vpd,double *etad,double *sed,int *control);
void RQT(double *A,int *r,int *c);


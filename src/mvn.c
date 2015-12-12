/* (c) Simon N Wood. 2014. Released under GPL2. 
  
  likelihood and derivative evaluation for multivariate Gaussian 
  additive models.

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <R_ext/BLAS.h>
#include "mgcv.h"

void mvn_ll(double *y,double *X,double *XX,double *beta,int *n,int *lpi, /* note zero indexing */
            int *m,double *ll,double *lb,double *lbb,double *dbeta,
            double *dH,int *deriv,int *nsp,int *nt) {
/* inputs:
    * 'y' is an m by n matrix, each column of which is a m-dimensional observation of a 
      multivariate normal r.v.
    * 'X' is a sequence of model matrices. The first (0th) model matrix runs from columns 0 to lpi[0]-1,
      the jth from cols lpi[j-1] to lpi[j]-1. lpi indexing starts from 0!!
    * XX is the pre-computed X'X matrix.
    * 'beta' is a parameter vector corresponding to X. The m*(m+1)/2 elements starting at lpi[m] are the 
      parameters of the Choleki factor of the precision matrix.
    * nt is number of threads to use.     
  
    outputs:
    * 'll' is the evaluated log likelihood.
    * 'lb' is the grad vector 
*/
  double *R,*theta,ldetR,*Xl,*bl,oned=1.0,zerod=0.0,*p,*p1,*p2,*p3,xx,zz,yy,*yty,
    *mu,*Rymu,rip,*dtheta,*db,*deriv_theta,*yX,*yRX;
  int i,j,k,l,pl,one=1,bt,ct,nb,*din,ntheta,ncoef,*rri,*rci,ri,rj,ril,rjl,rik,rjk,rij,rjj,q,r;
  const char not_trans='N';
  ntheta = *m * (*m+1)/2;ncoef = lpi[*m-1];
 
  nb = ncoef + ntheta; /* number of coefficients overall */
  /* Create the Choleski factor of the precision matrix */
  R = (double *)CALLOC((size_t)*m * *m,sizeof(double));
  theta = beta + lpi[*m-1]; /* parameters of R */
  ldetR = 0.0; /* log|R| */
  rri = (int *)CALLOC((size_t)ntheta,sizeof(int)); /* theta to R row index */
  rci = (int *)CALLOC((size_t)ntheta,sizeof(int)); /* theta to R col index */
  deriv_theta = (double *)CALLOC((size_t)ntheta,sizeof(double)); /* exp(theta) or 1*/
  for (k=0,i=0;i<*m;i++) { /* fill out R */
    deriv_theta[k] = exp(theta[k]);
    R[i + *m * i] = deriv_theta[k];ldetR += theta[k];
    rri[k]=rci[k]=i;k++; 
    for (j=i+1;j<*m;j++) { 
      R[i + *m * j] = theta[k];
      deriv_theta[k] = 1.0;
      rri[k]=i;rci[k]=j;k++;
    }
  }  
  /* obtain y - mu */
  mu  = (double *)CALLOC((size_t)*n,sizeof(double));
  for (l=0;l<*m;l++) { /* loop through components */
    if (l==0) { Xl = X;pl = lpi[0];bl=beta;} /* Xl is lth model matrix with pl columns, coef vec bl */ 
    else { Xl = X + *n * lpi[l-1];pl = lpi[l]-lpi[l-1];bl = beta + lpi[l-1];}   
    F77_CALL(dgemv)(&not_trans,n,&pl,&oned,Xl,n, bl, &one,&zerod, mu, &one); /* BLAS call for mu = Xl bl */
    /* now subtract mu from relevant component of y */
    for (p=mu,p1= mu + *n,p2=y+l;p<p1;p++,p2 += *m) *p2 -= *p;
  }
  FREE(mu);
  /* so y now contains y-mu */
  
  /* R(y-mu) is required repeatedly... */
  Rymu =  (double *)CALLOC((size_t)*n * *m,sizeof(double));
  bt=0;ct=0;mgcv_pmmult(Rymu,R,y,&bt,&ct,m,n,m,nt);  
  /* compute the log likelihood */
  for (*ll=0.0,p=Rymu,p1=p + *n * *m;p<p1;p++) *ll += *p * *p;
  *ll = - *ll/2 + ldetR * *n;  
  
  /* now the grad vector */
  
  p = lb;
  /* first the derivatives w.r.t. the coeffs of the linear predictors */
  for (l=0;l<*m;l++) { /* work through dimensions */
    if (l==0) { Xl = X;pl = lpi[0];} /* Xl is lth model matrix with pl columns */ 
    else { Xl = X + *n * lpi[l-1];pl = lpi[l]-lpi[l-1];} 
    for (i=0;i<pl;i++) { /* work through coefs for this dimension */
      *p = 0.0; /* clear lb */
      for (j=0;j<*n;j++) {
        xx = Xl[j + *n * i];
        for (p1=R + l * *m,p2 = p1 + l,p3 = Rymu + *m *j;p1<=p2;p1++,p3++) *p += xx * *p1 * *p3; 
      }
      p++;
    }
  }
  /* now the derivatives w.r.t. the parameters, theta, of R */ 
  
  for (k=0,i=0;i<*m;i++) { /* i is row */ 
    /* get tr(R^{-1}R R_theta^k) */
    xx = deriv_theta[k]; /* the non-zero element of R_theta^i at i,i */;
    *p += *n;
    k++; /* increment the theta index */ 
    /* quadratic form involves only ith dimension */
    for (zz=0.0,l=0,p1 = Rymu+i,p2=y+i;l<*n;l++,p1 += *m,p2 += *m) zz += *p1 * *p2 * xx;
    *p += -zz; 
    p++;
    for (j=i+1;j<*m;j++) { /* j is col */
      k++; /* increment the theta index */ 
      for (zz=0.0,l=0,p1 = Rymu+i,p2=y+j;l<*n;l++,p1 += *m,p2 += *m) zz += *p1 * *p2;
      *p += -zz;
      p++;
    }
  }  

  /* the Hessian is needed next */
  /* create index vector of dimension to which each coef relates ...*/ 
  din =  (int *)CALLOC((size_t)ncoef,sizeof(int));
  for (k=0,i=0;i<ncoef;i++) { 
    if (i==lpi[k]) k++; 
    din[i] = k;
  } 
  /* first the mean coef blocks */
  for (i=0;i<ncoef;i++) for (j=0;j<=i;j++) {
     l=din[i];k=din[j]; /* note l>=k */
     /* inner product of col l and col k of R ... */
     for (p=R+l * *m,p1=R+k * *m,rip=0.0,p2=p1+k;p1<=p2;p++,p1++) rip += *p * *p1;
     lbb[i + nb * j] = lbb[j + nb * i] = -XX[i + ncoef * j]*rip; /* -xx*rip; */ 
  }
  /* now the mixed blocks */
  for (i=0;i<ncoef;i++) for (j=0;j<ntheta;j++) {
     ri = rri[j]; /* row index of theta[j] in R */
     rj = rci[j]; /* col index of theta[j] in R */
     l = din[i]; /* which dimension does beta[i] relate to */
     xx = 0.0;
     zz = deriv_theta[j];  /* the non-zero derivative of R w.r.t. theta */
     /* term \bar x_i^{lT} R_\theta^{jT} R(y-\mu) is only non zero if l=rj... */
     if (l==rj) for (p = X + i* *n,p1=Rymu+ri,p2=p + *n;p<p2;p++,p1 += *m) xx += *p * *p1;
     xx *= zz;
     /* next term is inner product of ith col of X with row rj of y-mu (in y), multiplied by a constant*/
     if (ri<=l) { 
       for (yy=0.0,p = X + i* *n,p1=y+rj,p2=p + *n;p<p2;p++,p1 += *m) yy += *p * *p1;
       xx += yy * R[ri + *m * l]*zz;
     }
     lbb[i + nb * (j+ncoef)] = lbb[j + ncoef + nb * i] = xx;
  }
  /* the theta block completes the Hessian... */
  for (k=0;k<ntheta;k++) for (l=0;l<=k;l++) {
      xx=0.0; /* to accumulate term */
      if (k==l) {
        ri=rri[k];rj=rci[k];
        if (ri==rj) { 
          /* compute (y-\mu)'R'R_\theta^k_\theta^l(y-\mu) */
          for (zz=0.0,i=0,p=Rymu+ri,p2=y+ri;i<*n;i++,p += *m,p2+= *m) zz += *p * *p2;
          xx -= zz *  deriv_theta[k];
        }
      } 
      
      ri=rri[k];rj=rci[k];
      zz = deriv_theta[k]; 
      /* compute (y-\mu)'R_\theta^l'R_\theta^k_\theta^l(y-\mu) */
      ril=rri[l];rjl=rci[l];
      rik=rri[k];rjk=rci[k];
      if (ril==rik) { /* then term is non-zero */
	for (yy=0.0,i=0,p=y+rjl,p1=y+rjk;i<*n;i++,p+= *m, p1+= *m) yy += *p * *p1;
        yy *= zz;
        if (ril==rjl) yy *= deriv_theta[l];
        xx -= yy;
      }
      lbb[k + ncoef + nb * (l+ncoef)] = lbb[l + ncoef + nb * (k+ncoef)] = xx;
  }

  /* Now the derivatives of the Hessian, given the derivatives of the coefficients,
     wrt the smoothing parameters */
  if (*deriv) {
    yX = (double *)CALLOC((size_t)*m * ncoef,sizeof(double)); /* need (y-mu)X - m by ncoef */
    bt=0;ct=0;mgcv_pmmult(yX,y,X,&bt,&ct,m,&ncoef,n,nt); /* rows, dim, cols coef */   
    yRX = (double *)CALLOC((size_t)*m * ncoef,sizeof(double)); /* need R(y-mu)X - m by ncoef*/
    bt=0;ct=0;mgcv_pmmult(yRX,Rymu,X,&bt,&ct,m,&ncoef,n,nt); /* rows, dim, cols coef */  
    yty = (double *)CALLOC((size_t)*m * *m,sizeof(double)); /* need (y-mu)(y-mu)' - m by m */
    bt=0;ct=1;mgcv_pmmult(yty,y,y,&bt,&ct,m,m,n,nt); /* rows, cols dim */  
    for (r=0;r< *nsp;r++) { 
      db = dbeta + nb * r; /* d coefs / d rho_r */
      dtheta = db + ncoef; /* d theta / d rho_r */

      /* the derivatives of the hessian w.r.t. the smoothing parameters. First the portions 
         relating to the coefficients */
      for (i=0;i<ncoef;i++) for (j=0;j<=i;j++) {
	l = din[i];k = din[j]; /* dimensions for these elements */
	xx=0.0;p=R+l* *m;p1=R+k* *m;
        for (q=0;q<ntheta;q++) { /* sum over derivs of theta_q */
          ri=rri[q];rj=rci[q]; /* row and col of non zero element of deriv of R wrt theta_q */
          if (rj==l) xx += p1[ri]*deriv_theta[q]*dtheta[q];
          if (rj==k) xx += p[ri]*deriv_theta[q]*dtheta[q];
        }
        /* xx now contains the required element of (R_\theta^{qT} R + R^T R_\theta^q) d\theta_q/d \rho_r 
           which is multiplied by the inner product of columns i and j of X... */
        dH[i + nb * j] = dH[j + nb * i] = -xx * XX[i +  ncoef * j];
      } 
      /* now the mixed blocks */
      for (i=0;i<ncoef;i++) for (j=0;j<ntheta;j++) {  
	/* first the summation over the derivatives of beta */
        l=din[i];
	for (xx=0.0,q=0;q<ncoef;q++) {
          k=din[q];
          ri=rri[j];rj=rci[j]; /* row and col of non zero element of deriv of R wrt theta_j */
          zz = deriv_theta[j];/* deriv R w.r.t theta_l */
          if (rj==l) xx += -R[ri + *m * k]*zz*XX[i + ncoef * q] * db[q];
          if (rj==k) xx += -R[ri + *m * l]*zz*XX[i + ncoef * q] * db[q];
	}
        /* now the summation over the derivatives of theta */
        rij=rri[j];rjj=rci[j]; /* row and col of non zero element of deriv of R wrt theta_j */
        for (k=0;k<ntheta;k++) {
          zz=0.0;
          rik=rri[k];rjk=rci[k]; /* row and col of non zero element of deriv of R wrt theta_k */
          if (rij==rik&&(l==rjj||l==rjk)) { /* then term is non-zero */
            if (l==rjj) zz +=  yX[rjk + i * *m] * deriv_theta[j] * deriv_theta[k];
            if (l==rjk) zz +=  yX[rjj + i * *m] * deriv_theta[j] * deriv_theta[k];
            if (k==j&&rik==rjk) { /* then second deriv of R is non-zero */
               zz +=  deriv_theta[k] * yRX[rjj + *m * i];   
            }
            xx += zz * dtheta[k];
          }
          if (k==j&&rik==rjk) xx += dtheta[k]* deriv_theta[k] * R[rjj + *m * l] * yX[rjj + *m * i];/* x_i^l'R'R_tt^jk(y-mu) */
        }
        dH[i + (j+ncoef) * nb] = dH[j+ncoef + i * nb] = xx;   
      } /* mixed block loop */
      
      /* finally the theta block... */
      for (j=0;j<ntheta;j++) for (k=j;k<ntheta;k++) {
        rij=rri[j];rjj=rci[j];rik=rri[k];rjk=rci[k];
	/* first sum over the derivatives of beta... */
	xx = 0.0;
        for (i=0;i<ncoef;i++) {
	  zz=0.0;l=din[i];
          if (rij==rik&&(l==rjj||l==rjk)) { /* then term is non-zero */
            if (l==rjj) zz +=  yX[rjk + i * *m] * deriv_theta[j] * deriv_theta[k];
            if (l==rjk) zz +=  yX[rjj + i * *m] * deriv_theta[j] * deriv_theta[k];
            if (k==j&&rik==rjk) { /* then second deriv of R is non-zero */
              // following is suspicious... commenting it out improves the dodgy element from 7 out to 1 out, but 
              // spoils final element to 1 out.
	      // zz +=  deriv_theta[k] * R[rjj + *m * l] * yX[rjj + *m * i]; /* x_i^l'R'R_tt^jk(y-mu) */
              if (l==rjj) zz +=  deriv_theta[k] * yRX[rjj + *m * i];  /* x_i^l'R_tt^jk R(y-mu) */
            }
            xx += zz * db[i];
          }
          if (k==j&&rij==rjj) xx +=  db[i]*deriv_theta[k] * R[rjj + *m * l] * yX[rjj + *m * i]; /* x_i^l'R'R_tt^jk(y-mu) */
        }
        for (i=0;i<ntheta;i++) {
	  ri = rri[i];rj=rci[i];zz=0.0;
          if (j==k&&ri==rij&&rjk==rik) zz += deriv_theta[j]*deriv_theta[i]*yty[rj * *m + rjj];  /* row rjj, col rj */ 
          if (i==k&&rik==rij&&rj==ri) zz += deriv_theta[j]*deriv_theta[i]*yty[rjk * *m + rjj];  /* row rjj, col rjk */ 
          if (i==j&&rik==rij&&rj==ri) zz += deriv_theta[k]*deriv_theta[i]*yty[rjk * *m + rj];  /* row rjk, col rj */ 

          if (i==j&&j==k&&ri==rj) { /* pure derivative on diagonal of R: (y-mu)'R'R_ttt^iii(y-mu)*/
            for (yy=0.0,p=Rymu+ri,p1=y+ri,q=0;q<*n;p+= *m,p1+= *m,q++) yy += *p * *p1;           
            zz += deriv_theta[k]*yy;
          }
          xx += -zz*dtheta[i];
        }
	dH[k + ncoef + (j+ncoef) * nb] = dH[j+ncoef + (k+ncoef) * nb] = xx;
      }

      dH += nb * nb; /* move on to next Hessian */
    } /* smoothing parameter loop */ 
    FREE(yX);FREE(yRX);FREE(yty);
  } /* if (*deriv) */
  

  FREE(din); FREE(rri); FREE(rci);
  
  FREE(R);FREE(Rymu);FREE(deriv_theta);
} /* mvn_ll */



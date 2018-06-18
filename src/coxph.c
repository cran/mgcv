/* The Cox Proportional Hazard model for survival data, 
   for mgcv. (c) Simon N. Wood 2013-14 
   
   See for example, Hastie and Tibshirani (1990) for the log partial 
   likelihood (Peto's, 1972, approximation for ties). 

   For details of hazard estimation see...
 
   Klein, J.P and Moeschberger, M.L. (2003) Survival Analysis: Techniques for
   Censored and Truncated Data (2nd ed.) Springer

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include "mgcv.h"


void coxpred(double *X,double *t,double *beta,double *off,double *Vb,double *a,double *h,double *q,
             double *tr,int *n,int *p, int *nt,double *s,double *se) {
/* Function to predict the survivor function for the new data in 
   X (n by p), t, given fit results in a, h, q, Vb, and original event times 
   tr (length nt). 
   The new data are in descending order on entry, as is tr. 
   On exit n - vectors s and se contain the estimated survival function and its se.
*/
  double eta,*p1,*p2,*p3,*v,*pv,*pa,x,vVv,hi,exp_eta; 
  int ir=0,i=0;
  v = (double *)CALLOC((size_t)*p,sizeof(double)); 
  for (i=0;i<*n;i++) { /* loop through new data */
    while (ir < *nt && t[i]<tr[ir]) { /* find current interval */
      ir++;
      a += *p; /* moving the a pointer to a vector for this interval */
    } 
    if (ir == *nt) { /* before start of fit data */
      se[i] = 0; 
      s[i] = 1;
    } else { /* in the range */
      hi = h[ir]; /* cumulative hazard for this point */
      pv = v;pa = a;
      for (eta=0,p1=X,p2=beta + *p,p3=beta;p3<p2;pa++,p3++,pv++,p1+= *n) { 
        eta += *p1 * *p3; /* X beta */
        *pv = *pa - *p1 * hi; /* v = a - x * h */
      }
      exp_eta = exp(eta+off[i]);
      s[i] = exp(-hi*exp_eta); /* estimated survivor function */
      /* now get the s.e. for this... */
      p1 = Vb;pv = v;p2 = pv + *p;
      for (vVv=0;pv<p2;pv++) {
        for (x=0.0,p3 = v;p3<p2;p3++,p1++) x += *p3 * *p1;
        vVv += x * *pv; /* v'Vbv */
      }
      se[i] = exp_eta*s[i]*sqrt(q[ir] + vVv); /* standard error on survivor function */
    }
    X++; /* next prediction */
  } /* data loop */
  FREE(v);
} /* coxpred */

void coxpp(double *eta,double *X,int *r, int *d,double *h,double *q,double *km, 
            int *n,int *p, int *nt) {
/* Cox PH post-processing code computing 
   1. Baseline hazard + variance
   2. The a vectors used to compute the survival function variance.

   On entry 'eta' is X%*%beta - the linear predictor. rows of 'X' and 'eta'
   are arranged in reverse time order. There are 'nt' unique times. 
   r[i] is the index of the unique time corresponding to row i of 'X'.
   The latest times have the lowest indices. Notionally tr[r[i]] is the 
   time corresponding to row i, although this function does not use 'tr'.
   'X' is 'n' by 'p'.

   On exit:
   *  X is over written with the 'a' vectors. Each is length 'p' and all
      'nt' are stored one after the other. 
   * h is the cumulative hazard (h[i] at tr[i]) - an nt vector.
   * km is the basic Kaplan Meier hazard estimate
   * q is the variance of the hazard - an nt vector.

   - note that in R terms the log survivor function for the fit data is 
     -h[r+1]*eta ## r+1 to convert C indices to R indices.

   These ingredients are to be supplied to 'coxpred' to obtain the predicted 
   survivor function for individuals. 
*/
  double *b,*gamma_p,*gamma,*gamma_np,*bj,*bj1,*p1,*p2,gamma_i,*Xp,*aj,*aj1,x,y;
  int *dc,i,j;
  b = (double *)CALLOC((size_t) *nt * *p,sizeof(double)); /* storage for the b vectors */
  gamma_p = (double *)CALLOC((size_t) *nt,sizeof(double)); 
  gamma_np = (double *)CALLOC((size_t) *nt,sizeof(double));
  dc = (int *)CALLOC((size_t) *nt,sizeof(int)); /* storage for event counts at each time*/
  gamma = (double *)CALLOC((size_t)*n,sizeof(double)); 
  if (*p>0) for (i=0;i<*n;i++) gamma[i] = exp(eta[i]);
  else for (p1=gamma,p2=p1 + *n;p1<p2;p1++) *p1 = 1.0;

  bj1 = bj = b;
  for (i=0,j=0;j<*nt;j++) { /* work back in time */
    if (j>0) {
      gamma_p[j] = gamma_p[j-1]; gamma_np[j] = gamma_np[j-1];
      /* copy b^+_{j-1}, bj1, into b^+_j, bj */
      for (p1=bj,p2=p1 + *p;p1<p2;p1++,bj1++) *p1 = *bj1;
    }
    while (i < *n && r[i]==j+1) { /* accumulating this event's information */
      gamma_i = gamma[i];
      gamma_p[j] +=  gamma_i; gamma_np[j] += 1.0;
      dc[j] += d[i]; /* count the events */
      /* accumulate gamma[i]*X[i,] into bj */
      for (p1=bj,p2=p1 + *p,Xp = X + i;p1<p2;p1++,Xp += *n) *p1 += *Xp * gamma_i; 
      i++; /* increase the data counter */
    }
    bj += *p; /* move on to next b^+ vector */
  } /* back in time loop done */
  
  /* with gamma_p, dc and b computed, we can now do time forward accumulations 
     of h, q and a... */
  j = *nt - 1;
  x =  dc[j]/gamma_p[j];h[j] = x;km[j] = dc[j]/gamma_np[j];
  x /= gamma_p[j];q[j] = x;
  i = j * *p;
  for (aj=X+i,p1=aj+ *p,p2=b+i;aj<p1;p2++,aj++) *aj = *p2 * x;
  for (j--;j>=0;j--) { /* back recursion, forwards in time */
    y = dc[j];
    x = y/gamma_p[j];
    y/=gamma_np[j];
    h[j] = h[j+1] + x;
    km[j] = km[j+1] + y; /* kaplan meier hazard estimate */
    x /= gamma_p[j];
    q[j] = q[j+1] + x;
    /* now accumulate the a vectors into X for return */
    i = j * *p;
    //for (aj=X+i,aj1=p1=aj+ *p,p2=b+i;aj<p1;p2++,aj++) *aj = *aj1 + *p1 * x; 
    for (aj=X+i,aj1=p1=aj+ *p,p2=b+i;aj<p1;p2++,aj++,aj1++) *aj = *aj1 + *p2 * x;
  }
  FREE(b);FREE(gamma);FREE(dc);
  FREE(gamma_p);FREE(gamma_np);
} /* coxpp */



void coxlpl(double *eta,double *X,int *r, int *d,double *tr, 
            int *n,int *p, int *nt,double *lp,double *g,double *H,
            double *d1beta,
            double *d1H,
            double *d2beta,
            double *d2H,
            int *n_sp,int *deriv)
/* rows of n by p model matrix X are arranged in decreasing order
   of time. The unique event times are in nt vector tr in time reverse order.
   The ith row of X corresponds to event time tr[r[i]]. If d[i] is 0 then the 
   event is censoring. 

   On output:
   lp is the log partial likelihood.
   g is the p vector of derivatives of lp w.r.t. beta.
   H is the p by p second derivative matrix of lp wrt beta
   
   The d1* are the derivatives of H and beta wrt rho=log(lambda), 
   the log smoothing parameters. In each case there are n_sp replicates
   of the same dimension as the original object stored end to end. 
   
   The d1* & d2* are unused unless deriv is non-zero. d1/2beta contains the derivatives
   of beta wrt the log smoothing parameters, on entry.

   *deriv controls which derivatives are returned. 
          * < 0 and only lp is returned 
          * 0 only lp, g and H are returned.
          * 1 the first derivative of the leading diagonal of H w.r.t. rho is 
              returned in d1H, using d1b. 
          * 2 the first derivative of H w.r.t. rho is returned in d1H, using 
              d1b. 
	  * 3 is for first and second derivatives of H are returned in d1H
              and d2H. d2H contains only the derivatives of the leading 
              diagonal of H. This uses d1b and d2b.

   Except for the d2H structure, the d2* contain the 
   second derivative structures, packed as follows...   

   * v2 will contain d^2v/d\rho_0d\rho_0, d^2v/d\rho_1d\rho_0,... but rows will not be
     stored if they duplicate an existing row (e.g. d^2v/d\rho_0d\rho_1 would not be 
     stored as it already exists and can be accessed by interchanging the sp indices).
     So to get d^2v_k/d\rho_id\rho_j: 
     i)   if i<j interchange the indices
     ii)  off = (j*m-(j+1)*j/2+i)*q (m is number of rho, q is dim(v))
     iii) v2[off+k] is the required derivative.       

   d2H contains second derivatives of the leading diagonal of H, only.

*/

{ int dr,i,j,tB=0,tC=0,k,l,m,off,nhh;
  double lpl=0.0,*gamma,gamma_p=0.0,
    eta_sum,
    *b_p=NULL,*A_p=NULL,*p1,*p2,*p3,*p4,
    *d1gamma=NULL,
    *d1gamma_p=NULL,*d1eta=NULL,xx,xx0,xx1,xx2,xx3,*d1b_p=NULL,*d1A_p=NULL,
    *d2eta=NULL,*d2gamma=NULL,*d2gamma_p=NULL,*d2b_p=NULL,
    *d2ldA_p=NULL;
  gamma = (double *)CALLOC((size_t)*n,sizeof(double)); 
  if (*deriv >=0) {
    b_p = (double *)CALLOC((size_t)*p,sizeof(double));
    A_p = (double *)CALLOC((size_t)(*p * *p),sizeof(double));
  }
  /* form exponential of l.p. */

  for (i=0;i<*n;i++) gamma[i] = exp(eta[i]);

  if (*deriv>0) { /* prepare for first derivatives */
    /* Get basic first derivatives given d1beta */
    d1eta = (double *)CALLOC((size_t)(*n * *n_sp),sizeof(double));
    mgcv_mmult(d1eta,X,d1beta,&tB,&tC,n,n_sp,p);
    p1=d1gamma = (double *)CALLOC((size_t)(*n * *n_sp),sizeof(double));
    p2=d1eta;
    for (j=0;j<*n_sp;j++) 
    for (i=0;i<*n;i++) {
	*p1 = *p2 * gamma[i]; p1++; p2++;
    }
    /* accumulation storage */
    d1gamma_p = (double *)CALLOC((size_t)*n_sp,sizeof(double));
    d1b_p = (double *)CALLOC((size_t)(*n_sp * *p),sizeof(double));   
  }

  if (*deriv>2) { /* prepare for second derivative calculations */
    /* Basic second derivative derived from d2beta */ 
    nhh = *n_sp * (*n_sp+1) / 2; /* elements in `half hessian' */
    d2eta  = (double *)CALLOC((size_t)(*n * nhh),sizeof(double));
     
    mgcv_mmult(d2eta,X,d2beta,&tB,&tC,n,&nhh,p);
   
    p1=d2gamma  = (double *)CALLOC((size_t)(*n * nhh),sizeof(double));
    p2=d2eta;
    for (j=0;j<*n_sp;j++) {  /* create d2gamma */
      for (k=j;k<*n_sp;k++) {
        p3 = d1eta + j * *n;
        p4 = d1eta + k * *n; 
        for (i=0;i<*n;i++) {
          *p1 = gamma[i] * (*p2 + *p3 * *p4);
	    p1++;p2++;p3++;p4++;
        }
      }
    } /* end of d2gamma loop */  
    /* accumulation storage */
    d2gamma_p = (double *)CALLOC((size_t) nhh,sizeof(double));
    d2b_p = (double *)CALLOC((size_t)( nhh * *p),sizeof(double));
  }

  if (*deriv>0) { /* Derivatives of H are required */
    /* create storage for accumulating derivatives */
    d1A_p = (double *)CALLOC((size_t)(*n_sp * *p * *p),sizeof(double));
    /* clear incoming storage */
    for (j = *n_sp * *p * *p,k=0;k<j;k++) d1H[k] = 0.0;
    /* note that only leading diagonal of d2H is obtained and stored */ 
    if (*deriv>2) {   
      d2ldA_p = (double *)CALLOC((size_t)(nhh * *p),sizeof(double));
      for (j = nhh * *p,k=0;k<j;k++) d2H[k] = 0.0; 
    }
  }

  /* now accumulate the log partial likelihood */
  lpl=0.0;
  for (k=0;k<*p;k++) g[k] =0.0; 
  for (k = 0;k < *p;k++) for (m = 0;m < *p ;m++)  H[k + *p * m] = 0.0;
  i=0; /* the row index */

  for (j=0;j<*nt;j++) { /* work back in time */
    eta_sum=0.0;
    dr=0;
    while (i < *n && r[i]==j+1) { /* accumulating this event's information */
      /* lpl part first */ 
      gamma_p += gamma[i];
      if (d[i]==1) { dr++;eta_sum+=eta[i];}
      /* now the first derivatives */
      if (*deriv >= 0) {
        for (k=0;k<*p;k++) b_p[k] += gamma[i]*X[i + *n * k];
        if (d[i]==1) for (k=0;k<*p;k++) g[k] += X[i + *n * k];
        /* and second derivatives */
        for (k = 0;k < *p;k++) for (m = k;m < *p ;m++)
	    A_p[k + *p *m] +=  gamma[i]*X[i + *n * k] * X[i + *n * m];
      }
      /* derivatives w.r.t. smoothing parameters */
      if (*deriv >0 ) { /* first derivative stuff only */
        for (k=0;k<*n_sp;k++) d1gamma_p[k] += d1gamma[i + *n * k];
        for (m=0;m<*n_sp;m++) {
          xx = d1gamma[i + *n * m];
          for (k=0;k<*p;k++) d1b_p[k + *p * m] += xx * X[i + *n * k];
        }
      } /* end of first derivative accumulation */
    
      if (*deriv>2) { /* second derivative accumulation */
         off = 0;         
         for (m=0;m<*n_sp;m++)   
         for (k=m;k<*n_sp;k++) { /* second derivates loop */
	      d2gamma_p[off] += d2gamma[i+ off * *n];
              for (l=0;l<*p;l++) 
              d2b_p[l + off * *p] +=  d2gamma[i+ off * *n] * X[i + *n * l];
              off++;
         } /* end k-loop */
      }

      if (*deriv>0) { /* H derivatives needed */	                
          for (m=0;m<*n_sp;m++) { /* First derivatives of A_p */
	      xx = d1gamma[i + *n * m];
              for (k = 0;k < *p;k++) for (l = k;l < *p ;l++) 
              d1A_p[k + *p * l + m * *p * *p] += xx * X[i + *n * k] * X[i + *n * l];             
          }
          if (*deriv>2) {      
            off = 0;         
            for (m=0;m<*n_sp;m++)   
	    for (k=m;k<*n_sp;k++) { /* second derivates of leading diagonal of A_p loop */
                for (l=0;l<*p;l++) 
                d2ldA_p[l + off * *p] +=  d2gamma[i+ off * *n] * X[i + *n * l] * X[ i + *n *l];
                off++;
            } /* end m/k -loop */
	  }        
      }

      i++;
    } /* finished getting this event's information */

    lpl += eta_sum - dr * log(gamma_p);
    if (*deriv>=0) {
      for (k=0;k<*p;k++) g[k] += - dr/gamma_p * b_p[k]; 
      for (k = 0;k < *p;k++) for (m = k;m < *p ;m++) 
      H[k + *p * m] += - dr * A_p[k + *p *m] /gamma_p +
        	        dr * b_p[k]*b_p[m]/(gamma_p*gamma_p); 
    }

    if (*deriv>0) { /* need derivatives of H */
        for (m=0;m<*n_sp;m++) { /* first derivatives of H */
	    xx0 =dr/gamma_p;
            xx = d1gamma_p[m]*xx0/gamma_p;        
            xx1 = xx0/gamma_p;
            xx2 = xx1*2*d1gamma_p[m]/gamma_p;
            for (k = 0;k < *p;k++) for (l = k;l < *p ;l++) {
		off = k + *p * l + m * *p * *p;
                d1H[off] += xx1 * (d1b_p[k + *p *m] * b_p[l] + b_p[k] * d1b_p[l + *p *m]) -
		    xx2 * b_p[k] * b_p[l] + xx * A_p[k + *p * l] - xx0 * d1A_p[off];
            }
        } /* m-loop end */
        if (*deriv>2) {
          xx = dr/gamma_p;
          xx0 = xx/gamma_p; /* dr/gamma_p^2 */
          xx1 = xx0/gamma_p; /* dr/gamma_p^3 */
          xx2 = xx1/gamma_p;
          off = 0;         
          for (m=0;m<*n_sp;m++) {
	    xx3 = -2*xx1*d1gamma_p[m];
	    for (k=m;k<*n_sp;k++) { /* second derivates of leading diagonal of H */
              for (l=0;l<*p;l++) {
		  d2H[l + off * *p] += xx3 * (A_p[l + *p *l] * d1gamma_p[k] + 
                                              2 * d1b_p[l + *p * k] * b_p[l]) + 
		                       xx0 * (d1A_p[l + l * *p + m * *p * *p] * d1gamma_p[k] 
                                              + A_p[l + *p * l] * d2gamma_p[off] + 
                                              d2b_p[l + off * *p] * b_p[l] + 
                                              2 * d1b_p[l + *p * k] * d1b_p[ l + *p * m] + 
		                              b_p[l] * d2b_p[l + off * *p]) +
                                       xx0 * d1gamma_p[m] * d1A_p[l + l * *p + k * *p * *p] -
                                       xx * d2ldA_p[l + off * *p] + 
                                       6 * xx2 * d1gamma_p[m] * b_p[l] * b_p[l] * d1gamma_p[k] -
		                       2 * xx1 * (2*d1b_p[l + *p * m] * b_p[l] * d1gamma_p[k] +
						  b_p[l]*b_p[l]*d2gamma_p[off]);
               
              }
              off++;
            } /* end k -loop */    
	  } /* end m - loop */
        } /* end if (*deriv>2) */
     } /* end of H derivatives */
  } /* end of j loop (work back in time) */

  for (k=0;k<*p;k++) for (m=0;m<k;m++) H[k + *p *m] = H[m + *p *k];
  if (*deriv>1) for (m=0;m<*n_sp;m++) {
    off = *p * *p * m;
    for (k = 0;k < *p;k++) for (l = 0;l < k ;l++) 
	d1H[k + *p * l + off] = d1H[l + *p * k + off];
  }
 
  if (*deriv>=0) { FREE(A_p);FREE(b_p);}
  FREE(gamma);

  if (*deriv > 0) { /* clear up first derivative storage */
    FREE(d1eta);FREE(d1gamma);
    FREE(d1gamma_p);FREE(d1b_p);
    FREE(d1A_p);
  }

  if (*deriv > 2) { /* clear up second derivative storage */
    FREE(d2eta);FREE(d2gamma);
    FREE(d2gamma_p);FREE(d2b_p);
    FREE(d2ldA_p);    
  }
  *lp = lpl;
} /* end coxlpl */




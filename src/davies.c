/* Simon N. Wood (Feb 2020)

   C implementation of Davies, R.B. (1980) "The Distribution of a Linear Combination 
   of \chi^2 Random Variables" J. R. Statist. Soc. C 29,323-333.
   For the basic method see Davies, R. B. (1973). "Numerical inversion of a characteristic function" 
   Biometrika, 60(2), 415-417. The 1980 paper provides the detail for the error bounds needed
   for the 1973 method, and provides Algol 60 code.

   Hand translated from the original Algol 60, but removing global variables and use of goto,
   doing a slightly more efficient sort in place of original 'order' (and simplifying the way
   this is called), calling R functions for ln1, and combining the original intl1 and intl2 into
   their sum intl, and ersm1 and ersm2 into their sum ersm, since only the sums are actually used.

   For a C++ translation (leaving the globals and gotos in place) see CompQuadForm, function 'davies' 
   (they also provide an 'imhof' function).
*/

#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "mgcv.h"

double ln1(double x,int first) {
  if (first) return(log1p(x)); // log(1+x)
  else return(log1pmx(x)); // log(1+x)-x 
}

/* procedure order is th = order(lb,decreasing=TRUE)
   R_orderVector1(th,r,Rf_lang(lb),1,1)
 */ 

int counter(int clear) {
  static int count=0;
  int a=0;
  if (clear) { a=count;count=0;} else count++;
  return(a);					
}  

double errbd(double u,double *cx,double sigsq,int r,int *n,double *lb,double *nc) {
/* finds bound on the tail probability, cutoff returned as cx */
  double sum1,lj,ncj,x,y,xy;
  int j,nj;
  j = counter(0);
  *cx = u * sigsq;
  sum1 = u * *cx;
  u = u*2;
  for (j=r-1;j>=0;j--) {
    nj = n[j];lj = lb[j];
    ncj = nc[j]; x = u * lj;
    y = 1.0 - x; *cx +=  lj*(ncj/y+nj)/y;
    xy = x/y;
    sum1 = sum1 + ncj * xy*xy + nj*(x*xy+ln1(-x,0)); 
  }
  x = exp(-0.5*sum1);
  return(x);
} /* errbd */

double ctff(double accx, double *upn,double mean,double lmin,double lmax,double sigsq,
	    int r,int *n,double *lb,double *nc) {
/* find ctff so that Pr(qf>ctff)<accx if upn>0, Pr(qf<ctff)<accx othrwise 
   modifies upn.
*/
  double u2,u1,u,c1,c2,rb,cst;
  u2 = *upn; u1 = 0.0; c1 = mean;
  if (u2>0) rb = 2*lmax; else rb = 2*lmin;
  
  while (1) {
    if (errbd(u2/(1+u2*rb),&c2,sigsq,r,n,lb,nc)<=accx) break;
    u1=u2;c1=c2;u2 = u2*2;
  }
  while (1) {
    if ((c1-mean)/(c2-mean) >= 0.9) break;
    u = (u1+u2)*0.5;
    if (errbd(u/(1+u*rb),&cst,sigsq,r,n,lb,nc)>accx) {
      u1=u;c1=cst;
    } else {
      u2=u;c2=cst;
    }
  }
  *upn = u2;
  return(c2);
} /* ctff */

double truncation(double u,double tausq,double sigsq,int r,int *n,double *lb,double *nc) {
/* bound the integration error due to truncation...*/
  int j,s,nj;
  double pi,sum1,sum2,prod1,prod2,prod3,lj,ncj,x,y,err1,err2;
  j=counter(0);
  pi = asin(1.0)*2;
  sum1 = prod2 = prod3 = 0.0;
  s = 0;sum2 = (sigsq + tausq)*u*u;
  prod1 = 2*sum2; u = 2*u;
  for (j=0;j<r;j++) {
    lj=lb[j];ncj=nc[j];nj=n[j];
    x=u*lj;x = x*x;
    sum1 += ncj*x/(1+x);
    if (x>1) {
      prod2 += nj*log(x);
      prod3 += nj*ln1(x,1);
      s = s + nj;
    } else prod1 +=  nj * ln1(x,1);
  } /* j loop */
  sum1 = sum1 * 0.5; prod2 += prod1;prod3 += prod1;
  x = exp(-sum1 - 0.25*prod2)/pi;
  y = exp(-sum1 - 0.25*prod3)/pi;
  if (s==0) err1 = 1.0; else err1 = 2.0*x/s;
  if (prod3>1) err2 = 2.5*y; else err2 = 1.0;
  if (err2<err1) err1 = err2;
  x = 0.5*sum2;
  if (x<=y) err2 = 1.0; else err2 = y/x;
  if (err1<err2) x = err1; else x=err2;
  return(x);
} /*truncation*/


double findu(double utx,double accx,double sigsq,int r,int *n,double *lb,double *nc) {
/* find u s.t. truncation(u,...) < accx and truncation(u/1.2,...) > accx */
  double ut,u,a[4]={2.0,1.4,1.2,1.1};
  int i;
  ut = utx; u = ut * 0.25;
  if (truncation(u,0,sigsq,r,n,lb,nc) > accx) {
    while (truncation(ut,0,sigsq,r,n,lb,nc) > accx) ut = ut * 4;
  } else {
    ut = u;u = u/4;
    while (truncation(u,0,sigsq,r,n,lb,nc) <= accx) { ut = u;u = u / 4;}
  }
  for (i=0;i<4;i++) {
    u = ut/a[i];
    if (truncation(u,0,sigsq,r,n,lb,nc) <= accx) ut = u;
  }
  return(ut);
} /* findu */

void integrate(int nterm,double interv,double tausq,int main,double c,double acc,
	       double *intl,double *ersm,double sigsq,int r,int *n,double *lb,double *nc) {
/* carry out integration with nterms terms at stepsize interv. If !main then
   multiply integrand by 1.0 - exp(-0.5*tausq*u^2).
   Updates intl and ersm. This is a slight modification of original code 
   which computes 2 versions of each of these, but then only ever uses the
   sum of the versions. 
 */
  double pi,inpi,u,sum1,sum2,sum3,x,y,z;
  int j,k,nj;
  pi = asin(1.0)*2;
  inpi = interv/pi;
  for (k=nterm;k>=0;k--) {
    u = (k+.5)*interv; sum1 = -2*u*c;
    sum2 = fabs(sum1);sum3= -0.5*sigsq*u*u;
    for (j=r-1;j>=0;j--) {
      nj = n[j];x = 2*lb[j]*u; y = x*x;
      sum3 = sum3 - 0.25*nj*ln1(y,1);
      y = nc[j]*x/(1+y); z = nj*atan(x)+y;
      sum1 +=  z;sum2 += fabs(z);
      sum3 += - 0.5*x*y;
    } /* j loop */
    x = inpi*exp(sum3)/u;
    if (!main) x = x * (1-exp(-0.5*tausq*u*u));
    sum1 = sin(0.5*sum1)*x; sum2 = 0.5*sum2*x;
    *intl += sum1; *ersm +=  sum2; 
  } /*k loop */
} /* integrate */


double cfe(double x,int *th,double ln28,int r,int *n,double *lb,double *nc,int *fail) {
/* coef of tausq in error when convergence factor of exp(-0.5*tausq*u^2) is used when df
   is evaluated at x */
  double pi,axl,sum1,lj,axl1,axl2;
  int sxl,j,k,t;
  j=counter(0);
  pi = asin(1.0)*2;
  axl = fabs(x);
  if (x<0) sxl = -1; else sxl = 1;
  sum1=0.0;
  for (j=r-1;j>=0;j--) {
    t = th[j];
    if (lb[t]*sxl>0.0) {
      lj = fabs(lb[t]); axl1 = axl - lj *(n[t]+nc[t]);
      axl2 = lj / ln28;
      if (axl1>axl2) axl=axl1; else {
        if (axl>axl2) axl=axl2;
	sum1 = (axl-axl1)/lj;
	for (k=j-1;k>=0;k--) sum1 += n[th[k]] + nc[th[k]];
	break;
      }
    } /* if lb[t] */
  } /* j loop */
  if (sum1>100.0) {
    *fail=1;return(1.0);
  } else {
    lj=pow(2.0,sum1*.25)/(pi*axl*axl);*fail=0;
    return(lj);
  } 
} /* cfe */


void davies(double *lb,double *nc,int *n,int *r,double *sigma,double *c,int *lim,
	double *acc,double *trace,int *ifault) {
/* Evaluate Pr(Q<c) where Q = \sum_j^r lb[j] * X_j + \sigma X_0
   X_j is a chi^2_n[j](nc[j]) [i.e. n[j] is DoF and nc[j] non-centrality param delta_j^2]
   X_0 ~ N(0,1). 'lim' is upper bound on terms in integral. 'acc' is error bound.
   
   trace is a 7 vector of algorithm information: 
   {abs val sum, integration terms, number of integrations, main integration interval,
    truncation point in initial integral, sd of cov factor term, cycles to locate integration params}
   ifault a fault indicator:
   0 - fine; 1 - desired accuracy not obtained; 2 round off error possibly significant;
   3 - invalid parameters; 4 - unable to locate integration parameters. 
   result returned in c
*/
  double ln28,pi,intl,ersm,acc1,*alb,sd,sigsq,lmin,lmax,mean=0,lj,ncj,almx,
    utx,up,un,intv,intv1,d1,d2,tausq,x;
  int i,j,fail,*th,nj,ok,nt,ntm;
  j=counter(1);
  ln28 = log(2.0)/8;pi = asin(1.0)*2;
  for (j=0;j<7;j++) trace[j] = 0.0;
  *ifault = 0;
  intl =  ersm = 0.0;
  acc1 = *acc;
  fail = 0;
  alb = (double *)CALLOC((size_t)*r,sizeof(double));
  th = (int *)CALLOC((size_t)*r,sizeof(int));
  for (j=0;j < *r;j++) {
    alb[j]=fabs(lb[j]);
    th[j]=j;
  }  
  revsort(alb,th,*r); /* order |lb| */
  FREE(alb);
  /* ... equivalent to R call th = order(abs(lb),decreasing=TRUE) */
  /* find mean sd max and min of lb check params ok... */

  sd = sigsq = *sigma * *sigma;
  /* it seems that lmin and lmax are not quite as commented,
     and should be initialized to 0 here and not lb[1] which
     comment would suggest (fails badly otherwise!)...*/
  lmax=lmin=0;
  mean=0.0;
  for (j=0;j < *r;j++) {
    nj = n[j]; lj = lb[j]; ncj = nc[j];
    if (nj<0||ncj<0) {
      FREE(th);
      *ifault = 3;
      return;
    }	
    sd +=  lj*lj*(2.0*nj+4.0*ncj);
    mean +=  lj*(nj+ncj);
    if (lmax<lj) lmax = lj; else if (lmin>lj) lmin=lj;
  } /*j loop */
  if (sd==0.0) {
    if (*c>0.0) *c = 1.0; else *c = 0.0;
    FREE(th);
    return;
  }
  if (lmin == 0.0 && lmax == 0.0 && *sigma == 0.0) {
    *ifault = 3;FREE(th);
    return;
  }    
  sd = sqrt(sd);
  if (lmax < -lmin) almx = -lmin; else almx = lmax;

  /* starting values for findu and ctff... */
  utx = 16.0/sd; up = 4.5/sd; un = -up;

  /* truncation point with no convergence factor...*/
  utx = findu(utx,0.5*acc1,sigsq,*r,n,lb,nc);

  /* does convergence factor help? */
  if (*c != 0.0 && almx > 0.07*sd) {
    tausq = 0.25 * acc1 / cfe(*c,th,ln28,*r,n,lb,nc,&fail);
    if (!fail) {
      if (truncation(utx,tausq,sigsq,*r,n,lb,nc)<0.2*acc1) {
        sigsq = sigsq + tausq;
	utx = findu(utx,0.25*acc1,sigsq,*r,n,lb,nc);
	trace[5] = sqrt(tausq);
      }
    }
  } /* if (c != 0.0 ...) */
  trace[4] = utx; acc1 = 0.5 * acc1;

  /* find 'range' of distribution, quit if outside this... */
  ok = 1;
  while (ok) {
    d1 = ctff(acc1,&up,mean,lmin,lmax,sigsq,*r,n,lb,nc) - *c;
    if (d1<0.0) {
      FREE(th);*c=1.0;trace[6] = counter(1);
      return;
    }	
    d2 = *c - ctff(acc1,&un,mean,lmin,lmax,sigsq,*r,n,lb,nc);
    if (d2<0.0) {
      FREE(th);*c=0.0;trace[6] = counter(1);
      return;
    }
    /* find integration interval */
    if (d1>d2) intv = 2*pi/d1; else intv = 2*pi/d2;
    /* calculate number of terms required for main and auxiliary integrations */
    x = utx/intv;nt = (int)floor(x);if(x-nt>.5) nt++; /*round(utx/intv)*/
    x = 3.0/sqrt(acc1);ntm = (int)floor(x);if(x-ntm>.5) ntm++; /*round(3.0/sqrt(acc1))*/
    if (nt>ntm*1.5) {
      /* parameters for auxiliary integration */
      intv1 = utx/ntm; x = 2*pi/intv1;
      if (x<=fabs(*c)) break;
      /* calculate convergence factor */
      tausq = 0.33 * acc1 / (1.1*(cfe(*c-x,th,ln28,*r,n,lb,nc,&fail)+
				  cfe(*c+x,th,ln28,*r,n,lb,nc,&fail)));
      if (fail) break;
      acc1 = acc1 * 0.67;
      if (ntm > *lim) {
	FREE(th);*c = -1.0;trace[6] = counter(1);
	return;
      }	  
      /* auxiliary integration */
      integrate(ntm,intv1,tausq,0,*c,*acc,&intl,&ersm,sigsq,*r,n,lb,nc);
      *lim = *lim - ntm; sigsq = sigsq + tausq;
      trace[2] += 1.0; trace[1] += ntm + 1;
      /* find truncation point with new convergence factor */
      utx = findu(utx,0.25*acc1,sigsq,*r,n,lb,nc);
      acc1 = 0.75 * acc1;
    } else ok = 0;
  } /* while (ok) */
  /* main integration... */
  trace[3] = intv;
  if (nt > *lim) {
    FREE(th);*c=-1.0;*ifault=1;trace[6] = counter(1);
    return;
  }    
  integrate(nt,intv,0.0,1,*c,*acc,&intl,&ersm,sigsq,*r,n,lb,nc);
  trace[2] += 1; trace[1] += nt + 1;
  *c = 0.5 - intl;
  trace[0] = ersm;
  /* test whether round off error could be significant */
  x = ersm + *acc/10.0;
  j = 1;
  for (i=0;i<4;i++) {
    if (j*x == j*ersm) *ifault = 2;
    j = j * 2;
  }
  FREE(th);
  trace[6] = counter(1);
} /* davies */



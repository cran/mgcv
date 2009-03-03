/* Copyright (C) 2008 Simon N. Wood  simon.wood@r-project.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define ANSI
#include "matrix.h"
#include "mgcv.h"


double *forward_buf(double *buf,int *jal,int update)
/* extend buffer forward 1000 */
{ double *buf2,*p,*p1,*p2;
  int n=1000;
  buf2 = (double *)calloc((size_t)*jal+n,sizeof(double));
  for (p=buf,p1=buf + *jal,p2=buf2;p<p1;p++,p2++) *p2 = *p;
  free(buf);
  if (update) *jal += n;
  return(buf2);
}


double *backward_buf(double *buf,int *jal,int *j0,int *j_lo,int *j_hi,int update)
/* extend buffer backwards by 1000 elements, or to j=1 */
{ int n=1000;
  double *buf2,*p,*p1,*p2;
  if (n > *j0-1) n = *j0 - 1; /* only extend back to j=1 */
  if (n==0) return(buf);
  buf2 = (double *)calloc((size_t)*jal+n,sizeof(double));
  for (p=buf,p1=buf + *jal,p2=buf2 + n;p<p1;p++,p2++) *p2 = *p;  
  if (update) {
    *jal += n;
    *j_lo += n;
    *j_hi += n;  
    *j0 -= n;
  }
  free(buf);
  return(buf2);
}



void tweedious(double *w,double *w1,double *w2,double *y,double *phi,double *p,double *eps,int *n)
/* Routine to perform tedious series summation needed for Tweedie distribution
   evaluation, following Dunn & Smyth (2005) Statistics and Computing 15:267-280.
   Notation as in that paper. 
   
   log W returned in w. 
   d logW / dphi in w1, 
   d2 logW / d phi2 in w2.
 
   *Only* suitable for 1<p<2 (inequalities strict). 

   NOTE: still some redundancy for readability 
 
*/
{ int j_max,i,j_lo,j_hi,jb,jal,j0,j,ok;
  double x,ymax,ymin,alpha,*alogy,*p1,*p2,*p3,*wb,*wb1,*wb2,w_base,
    wmax,w1max,w2max,wmin,w1min,w2min,wi,w1i,w2i,
    log_eps,wj,w1j,w2j,jalogy;
  
  log_eps = log(*eps);
  alpha = (2 - *p)/(1 - *p);
  w_base = alpha * log(*p-1) - (1-alpha)*log(*phi) - log(2 - *p);

  /* initially establish the min and max y values, and hence the initial buffer range,
     at the same time produce the alpha log(y) vector. */ 
  
  alogy = (double *)calloc((size_t)*n,sizeof(double));
  ymax = ymin = *y;
  *alogy = alpha * log(*y);
  for (p1=y+1,p2=y+ *n,p3=alogy+1;p1<p2;p1++,p3++) {
    *p3 = alpha * log(*p1); /* alogy[i] = alpha * log(y[i]) */
    if (*p1 > ymax) ymax = *p1; else if (*p1 < ymin) ymin = *p1;
  }
  
  x = pow(ymin,2 - *p)/(*phi * (2 - *p));
  j_lo = (int) floor(x);if (j_lo<1) j_lo = 1;
   
  x = pow(ymax,2 - *p)/(*phi * (2 - *p));
  j_hi = (int) ceil(x);if (j_hi<j_lo) j_hi = j_lo;
  
  j0 = j_lo - 1000;if (j0<1) j0=1;
  jal = j_hi + 1000;
  
  jal -= j0-1;
  j_lo -= j0;
  j_hi -= j0;

  wb = (double *)calloc((size_t)jal,sizeof(double)); /* add -j*alogy[i] to get logW_j, for y[i] */
  wb1 = (double *)calloc((size_t)jal,sizeof(double)); /* add -j*alogy[i] to get logW_j', for y[i] */
  wb2 = (double *)calloc((size_t)jal,sizeof(double)); /* add -j*alogy[i] to get logW_j'', for y[i] */
  
  /* ... note that in the above it's log of derivative, not derivative of log */ 

  for (jb=j_lo,j=j_lo+j0;jb <= j_hi;jb++,j++) { /* jb is in buffer index, j is true index */ 
    wb[jb] = j * w_base - lgamma((double)j+1) - lgamma(-j * alpha);
    x = j*(alpha-1)/ *phi;
    wb1[jb] = wb[jb] + log(-x); /* note this is for log(-W_j')) */
    wb2[jb] = wb[jb]  + log(x*(x-1/ *phi));
  }

  /* Now j0 is the true j corresponding to buffer position 0. j starts at 1.
     jal is the number of buffer locations allocated. locations in the buffer between 
     j_lo and j_hi contain data. */



  for (i=0;i<*n;i++) { /* loop through y */
    /* first find the location of the series maximum... */
    x = pow(y[i],2 - *p)/(*phi * (2 - *p));
    j_max = (int) floor(x);
    if (x - j_max  > .5||j_max<1) j_max++; 
    j_max -= j0; /* converted to buffer index */
    
    j = j_max+j0;
    jalogy = j*alogy[i];
    wi=w1i=w2i=1.0;
    wmax = wb[j_max] - jalogy;wmin = wmax + log_eps; 
    w1max = wb1[j_max] - jalogy;w1min = w1max + log_eps;    
    w2max = wb2[j_max] - jalogy;w2min = w2max + log_eps;


    /* start upsweep to convergence or end of available buffered values */
    ok = 0;
    for (j=j_max+1+j0,jb=j_max+1;jb<=j_hi;jb++,j++) { 
      jalogy = j*alogy[i];
      wj = wb[jb] - jalogy;
      w1j = wb1[jb]  - jalogy;
      w2j = wb2[jb]  - jalogy;
      wi += exp(wj-wmax);
      w1i += exp(w1j-w1max);
      w2i += exp(w2j-w2max);
      if ((wj < wmin)&&(w1j < w1min)&&(w2j < w2min)) { ok=1;break;} /* converged on upsweep */
    } /* end of upsweep to buffer end */ 

    while (!ok) { /* while upsweep unconverged need to fill in more buffer */
      for (;jb<jal;jb++,j++) { /* fill buffers and calculate w terms */
        wb[jb] = j * w_base - lgamma((double)j+1) - lgamma(-j*alpha);
        x = j*(alpha-1)/ *phi;
        wb1[jb] = wb[jb] + log(-x);
        wb2[jb] = wb[jb] + log(x*(x-1/ *phi));
        jalogy = j*alogy[i];
        wj = wb[jb] - jalogy;
        w1j = wb1[jb]  - jalogy;
        w2j = wb2[jb]  - jalogy;
        wi += exp(wj-wmax);
        w1i += exp(w1j-w1max);
        w2i += exp(w2j-w2max);
        if ((wj < wmin)&&(w1j < w1min)&&(w2j < w2min)) { ok=1;break;} /* converged on upsweep */
        
      } 
      j_hi = jb; if (j_hi > jal-1) j_hi = jal-1; /* set j_hi to last element filled */
      if (!ok) { /* need to expand buffer storage*/
        wb = forward_buf(wb,&jal,0);
        wb1 = forward_buf(wb1,&jal,0);
        wb2 = forward_buf(wb2,&jal,1);
      }
    } /* finished upsweep and any buffer expansion */

    /* start downsweep to convergence or start of available buffered values */
    ok=0;
    for (j=j_max-1+j0,jb=j_max-1;jb>=j_lo;jb--,j--) { 
      jalogy = j*alogy[i];
      wj = wb[jb] - jalogy;
      w1j = wb1[jb]  - jalogy;
      w2j = wb2[jb]  - jalogy;
      wi += exp(wj-wmax);
      w1i += exp(w1j-w1max);
      w2i += exp(w2j-w2max);
      if ((wj < wmin)&&(w1j < w1min)&&(w2j < w2min)) { ok=1;break;} /* converged on downsweep */
    } /* end of downsweep to buffer end */ 
   
    if (j<=1&&j_lo==0) ok=1; /* don't care about element size if reached base */

    while (!ok) { /* while downsweep unconverged need to fill in more buffer */
      for (jb=j_lo-1;jb>=0;jb--,j--) { /* fill buffers and calculate w terms */
        wb[jb] = j * w_base - lgamma((double)j+1) - lgamma(-j*alpha);
        x = j*(alpha-1)/ *phi;
        wb1[jb] = wb[jb] + log(-x);
        wb2[jb] = wb[jb] + log(x*(x-1/ *phi));
        jalogy = j*alogy[i];
        wj = wb[jb] - jalogy;
        w1j = wb1[jb]  - jalogy;
        w2j = wb2[jb]  - jalogy;
        wi += exp(wj-wmax);
        w1i += exp(w1j-w1max);
        w2i += exp(w2j-w2max);
        if ((wj < wmin)&&(w1j < w1min)&&(w2j < w2min)) { ok=1;break;} /* converged on downsweep */
      } 
      if (j<=1) ok=1; /* don't care about element size if reached base */

      j_lo = jb; if (j_lo<0) j_lo=0; /* set j_lo to first element filled */
      if (!ok) { /* need to expand buffer storage*/
        wb = backward_buf(wb,&jal,&j0,&j_lo,&j_hi,0);
        wb1 = backward_buf(wb1,&jal,&j0,&j_lo,&j_hi,0);
        wb2 = backward_buf(wb2,&jal,&j0,&j_lo,&j_hi,1);
      }

    } /* finished downsweep and any buffer expansion */

    /* Summation now complete: need to do final transformations */
    w[i] = wmax + log(wi);    /* contains log W */
    w1i = w1max + log(w1i); /* contains log dW/dphi */
    w1[i] = -exp(w1i-w[i]);    /* d logW / d phi */
    w2[i] = w2max + log(w2i); /* contains log d2W/dphi2 */   
    w2[i] = exp(w2[i]-w[i]) - exp(2*w1i-2*w[i]); /* d2 logW / dphi2 */
 
  } /* end of looping through y */
  free(alogy);free(wb);free(wb1);free(wb2);
}


/* test code for tweedious...
library(mgcv);library(tweedie)
phi <- 2
p <- 1.1
mu <- .001
y <- c(1,1,2,1,3,0,0,30,67)
eps <- 1e-6
l0 <- colSums(mgcv:::ldTweedie(y,mu=mu,p=p,phi=phi))
l1 <- colSums(mgcv:::ldTweedie(y,mu=mu,p=p,phi=phi+eps))
  (l1-l0)/eps;l0

log(dtweedie(y,power=p,mu=mu,phi=phi))

j <- 1:100
alpha <- (2-p)/(1-p)
w <- -j*alpha*log(y)+alpha*j*log(p-1)-j*(1-alpha)*log(phi)-j*log(2-p)-lgamma(j+1) - lgamma(-j*alpha)
theta <- mu^(1-p)
k.theta <- mu*theta/(2-p)
theta <- theta/(1-p)  
(y*theta-k.theta)/phi - log(y) +  log(sum(exp(w)))

n <- 20
mu <- rep(1,n)
ml <- mgcv:::ldTweedie(1:n,mu,p=1.5,phi=1);ml
dl <- log(dtweedie.series(1:n,power=1.5,mu,phi=1));dl
x <- seq(.05,100,by=.1)
mu <- 1+x*0
sum(dtweedie(x,power=1.5,mu,phi=1))*.1 + dtweedie(0,power=1.5,1,phi=1)

sum(exp(mgcv:::ldTweedie(x,mu,p=1.5,phi=1)))*.1 + exp(mgcv:::ldTweedie(0,1,p=1.5,phi=1))


x <- rtweedie(10000,power=1.5,mu=1,phi=1)
  system.time(d1 <- dtweedie(x,power=1.5,mu=1,phi=1))
  system.time(d2 <- mgcv:::ldTweedie(x,mu=1,p=1.5,phi=1))
  range(d2-log(d1))

*/ 

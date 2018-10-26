/* Copyright (C) 2008-2014 Simon N. Wood  simon.wood@r-project.org

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
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rconfig.h>
#include "mgcv.h"

void *R_chk_calloc1(size_t nmemb,size_t size) {
  /* checks for zero or negative memory allocation calls...*/
  if (nmemb<=0) {
    Rprintf("adjusting %d memory allocation\n",nmemb);
    nmemb++;
  }  
  return(R_chk_calloc(nmemb,size));
}

/* Compute reproducing kernel for spline on the sphere */

void rksos(double *x,int *n,double *eps) {
/* Function to compute reproducing kernel for spline on the sphere,
   based on Jim Wendelberger's (1981) thesis. 
   Returns evaluated kernel rk(x) in n vector x.  
*/
  double dl1,xi,rk,xk,xx;
  int i,k;
  dl1 = acos(0)*2;
  dl1 = dl1*dl1/6; /* dilog(1) = pi^2/6, dilog(0)=0 */
  for (i=0;i< *n;i++) {
    xi = x[i];
    if (xi <= 0) {
      if (xi < -1) xi = -1;
      rk = 1.0 - dl1;
      xk = xi = xi/2 + 0.5;
      for (k=1;k<1000;k++) {
        xx = xk/(k*k);
        rk += xx;
        xk *= xi;
	if (xx < *eps) break;
      }
    } else {
      if (xi>1) xi=1;
      if (xi/2>=.5) rk=1.0; else
      rk = 1 - log(.5+xi/2)*log(.5-xi/2);
      xk = xi = .5 - xi/2;
      for (k=1;k<1000;k++) {
        xx = xk/(k*k);
        rk += -xx;
        xk *= xi;
	if (xk < *eps) break;
      } 
    }
    x[i] = rk;
  }
}


/* inside polygon tester.... */

void in_out(double *bx, double *by, double *break_code, double *x,double *y,int *in, int *nb, int *n)
/* finds out whether points in arrays x,y are inside boundary or outside, by counting boundary 
   crossings. The boundaries nodes are defined by bx, by.  bx[i] and by[i] less than or equal to 
   break_code signals a break in the boundary (e.g. between island and external boundary.) Each 
   section of boundary is assumed to be a closed loop. nb is dimenion of bx and by; n is dimension 
   of x and y. `in' will contain a 1 for an interior point and a 0 otherwise, on exit. 
   Both bx[i] and by[i] or neither must be less than the break_code.
*/ 
{ double xx,yy,dum,x0,x1,y0,y1;
  int i,j,count,start,swap;
  for (i=0;i<*n;i++) { /* loop through all test points */
      xx=x[i];yy=y[i]; /* the current test point */
      start=0; /* start of current boundary section */
      for (count=0,j=0;j<*nb;j++) { /* loop through entire boundary */
	x0 = bx[j]; /* start node */
        if (x0 <= *break_code) start=j+1; /* next segment start */
        else { /* not a new section start */
          if (j==*nb-1) x1=bx[start]; else x1 = bx[j+1];  /* end node */
          if (x1 <= *break_code) x1 = bx[start]; /* must join up segment end */
          if (x0!=x1) { /* x0==x1 => segment immaterial to decision */
	    if (x1<x0) { dum=x0;x0=x1;x1=dum;swap=1;} else swap=0; /* ordered */
            if (x0<xx&&x1>=xx) { /* might have a crossing */
	      y0 = by[j]; /* start node y co-ord */
              if (j==*nb-1) y1=by[start]; else
              y1 = by[j+1]; /* end node y co-ord */
              if (y1 <= *break_code) y1=by[start]; /* must join up */
              if (y0<=yy&&y1<=yy) count++; /* definite crossing */
              else { /* more detail needed to determine crossing */
		if (!(y0>yy&&y1>yy)) { /* could still be one */
                  if (swap) {dum=y0;y0=y1;y1=dum;}
                  dum = (xx-x0)*(y1-y0)/(x1-x0)+y0; /* at what y does vertical cross segment */
                  if (yy>=dum) count++; /* it's a crossing */
		} /* end - could still be one */
              } /* end - more detail */
            } /* end - might be a crossing */
          } /* end - does seg matter */
        } /* end - not skipped because break */
      } /* end boundary loop */
     if (count%2) in[i]=1;else in[i]=0; /* record result */
  } /* end x,y test loop */
} /* end of in_out */



/******************************/
/* Tweedie distribution stuff */
/******************************/

void psum(double *y, double *x,int *index,int *n) {
  /* y is of length max(index). x and index are of the same length, n.
     This routine fills y[index[i]-1] so that it contains the sum of 
     the x[i]'s sharing index[i]. It is assumed that y is cleared to zero
     on entry.
  */
  int i; 
  for (i=0;i< *n;i++) {
    y[index[i]-1] += x[i];
  }
}

double *forward_buf(double *buf,int *jal,int update)
/* extend buffer forward 1000 */
{ double *buf2,*p,*p1,*p2;
  int n=1000;
  buf2 = (double *)CALLOC((size_t)(*jal+n),sizeof(double));
  for (p=buf,p1=buf + *jal,p2=buf2;p<p1;p++,p2++) *p2 = *p;
  FREE(buf);
  if (update) *jal += n;
  return(buf2);
}


double *backward_buf(double *buf,int *jal,int *j0,int *j_lo,int *j_hi,int update)
/* extend buffer backwards by 1000 elements, or to j=1 */
{ int n=1000;
  double *buf2,*p,*p1,*p2;
  if (n > *j0-1) n = *j0 - 1; /* only extend back to j=1 */
  if (n==0) return(buf);
  buf2 = (double *)CALLOC((size_t)(*jal+n),sizeof(double));
  for (p=buf,p1=buf + *jal,p2=buf2 + n;p<p1;p++,p2++) *p2 = *p;  
  if (update) {
    *jal += n; /* number of buffer locations allocated */
    *j_lo += n; /* start of initialized elements within buffer */
    *j_hi += n;  /* end of initialized elements within buffer */
    *j0 -= n; /* j0 is the true j corresponding to buffer element 0 */
  }
  FREE(buf);
  return(buf2);
} /* backward_buf */


void tweedious(double *w,double *w1,double *w2,double *w1p,double *w2p,
	       double *w2pp,double *y,double *eps,int *n,
               double *th,double *rho,double *a, double *b)
/* Routine to perform tedious series summation needed for Tweedie distribution
   evaluation, following Dunn & Smyth (2005) Statistics and Computing 15:267-280.
   Notation as in that paper. For 
   
   log W returned in w. (where W means sum_j W_j)
   d logW / drho in w1, 
   d2 logW / d rho2 in w2.
   d logW / dth in w1p
   d2 logW / dth2 in w2p
   d2 logW / dth drho in W2pp

   rho=log(phi), and th defines p = (a + b * exp(th))/(exp(th)+1). 
   note, 1<a<b<2 (all strict) 

   The somewhat involved approach is all about avoiding overflow or underflow. 
   Extensive use is made of 
        log { sum_j exp(x_j)} = log { sum_j exp(x_j-x_max) } + x_max
   digamma and trigamma functions are from Rmath.h

   NOTE: still some redundancy for readability 
 
*/
{ int j_max,i,j_lo,j_hi,jb,jal,j0,j,ok;
  double x,x1,x2,xx,ymax,ymin,alpha,*alogy,*p1,*p2,*p3,*p4,*p5,
    *wb,*wb1,//*wb2,
    *wp1,*wp2,*wpp,dpth1=0,dpth2=0,//xmax,x1max,x2max,
    w_base,wp_base,wp2_base,wp1j,wp2j,wppj,wj_scaled,wdlogwdp,
    wdW2d2W,dWpp,exp_th,
    wmax,wmin,
    wi,w1i,w2i,//drho_const,
    log_eps,wj,w1j,jalogy,onep,onep2,*logy1p2,*logy1p3,p,phi;
  
  /* do everything in terms of working parameters, rho, th */
  phi = exp(*rho); 
  /* compute p and its derivatives w.r.t. th */
  if (*th>0) { 
      exp_th =  exp(- *th);
      //drho_const = (1+exp_th)/(1 - *b + (1 - *a)*exp_th);
      x = 1 + exp_th;p = (*b + *a * exp_th)/x;
      x1 = x*x;dpth1 = exp_th*(*b - *a)/x1;
      dpth2 =  ((*a - *b)*exp_th+(*b - *a)*exp_th*exp_th)/(x1*x);
  } else {
      exp_th =  exp(*th);
      //drho_const = (1+exp_th)/((1 - *b)*exp_th + 1 - *a);
      x = exp_th+1;p = (*b * exp_th + *a)/x;
      x1 = x*x;dpth1 = exp_th*(*b - *a)/x1;
      dpth2 = ((*a - *b)*exp_th*exp_th+(*b - *a)*exp_th)/(x*x1);
  }
  
  log_eps = log(*eps);
  onep = 1 - p;onep2 = onep * onep;
  alpha = (2 - p)/onep;
  /* get terms that are repeated in logWj etc., but simply multiplied by j */
  w_base = alpha * log(p-1) + *rho/onep - log(2 - p);
  wp_base = (log(-onep) + *rho)/onep2 - alpha/onep + 1/(2 - p);
  wp2_base= 2*(log(-onep) + *rho)/(onep2*onep) - (3*alpha-2)/(onep2) + 1/((2 - p)*(2 - p));
 
  /* initially establish the min and max y values, and hence the initial buffer range,
     at the same time produce the alpha log(y) log(y)/(1-p)^2 and log(y)/(1-p)^3 vectors. */ 
  
  alogy = (double *)CALLOC((size_t)*n,sizeof(double));
  logy1p2 = (double *)CALLOC((size_t)*n,sizeof(double));
  logy1p3 = (double *)CALLOC((size_t)*n,sizeof(double));

  ymax = ymin = *y;
  *alogy = alpha * log(*y);
  *logy1p2 = log(*y)/(onep2);
  *logy1p3 = *logy1p2/onep;
 
  for (p1=y+1,p2=y+ *n,p3=alogy+1,p4=logy1p2+1,p5=logy1p3+1;p1<p2;p1++,p3++,p4++,p5++) {
    x = log(*p1); /* log(y) */
    *p3 = alpha * x; /* alogy[i] = alpha * log(y[i]) */
    *p4 = x/onep2; /* log(y[i])/(1-p)^2 */ 
    *p5 = *p4/onep; /* log(y[i])/(1-p)^3 */
    if (*p1 > ymax) ymax = *p1; else if (*p1 < ymin) ymin = *p1;
  }
  
  x = pow(ymin,2 - p)/(phi * (2 - p));
  j_lo = (int) floor(x);if (j_lo<1) j_lo = 1;
   
  x = pow(ymax,2 - p)/(phi * (2 - p));
  j_hi = (int) ceil(x);if (j_hi<j_lo) j_hi = j_lo;

  j0 = j_lo - 1000;if (j0<1) j0=1;
  jal = j_hi + 1000;
  
  jal -= j0-1;
  j_lo -= j0;
  j_hi -= j0;

  /* prepare, up front, everything needed to form log W_j etc. except for the part depending on y[i]... */

  /* evaluation... */
  wb = (double *)CALLOC((size_t)jal,sizeof(double)); /* add -j*alogy[i] to get logW_j, for y[i] */
  /* first deriv wrt phi... */
  wb1 = (double *)CALLOC((size_t)jal,sizeof(double)); /* add -j*alogy[i] to get logW_j', for y[i] */
  /* second deriv wrt phi... */
  // wb2 = (double *)CALLOC((size_t)jal,sizeof(double)); /* add -j*alogy[i] to get logW_j'', for y[i] */
  /* ... note that in the above it's log of derivative, not derivative of log, but in the 
    following it's the derivative of the log... */

  /* first deriv wrt p... */
  wp1 = (double *)CALLOC((size_t)jal,sizeof(double));
  /* second deriv wrt p... */
  wp2 = (double *)CALLOC((size_t)jal,sizeof(double));
  /* second deriv wrt p and phi... */
  wpp = (double *)CALLOC((size_t)jal,sizeof(double));

  for (jb=j_lo,j=j_lo+j0;jb <= j_hi;jb++,j++) { /* jb is in buffer index, j is true index */ 
    wb[jb] = j * w_base - lgamma((double)j+1) - lgamma(-j * alpha);
    wb1[jb] = -j/onep;
    xx = j/onep2;
    x = xx*digamma(-j*alpha);
    wp1[jb] = j * wp_base + x; /* base for d logW_j/dp */
    xx = trigamma(-j*alpha) * xx * xx;
    wp2[jb] = j * wp2_base + 2*x/onep - xx;
    wpp[jb] = j /onep2;
  }

  /* Now j0 is the true j corresponding to buffer position 0. j starts at 1.
     jal is the number of buffer locations allocated. locations in the buffer between 
     j_lo and j_hi contain data. */


  for (i=0;i<*n;i++) { /* loop through y */
    /* first find the location of the series maximum... */
    x = pow(y[i],2 - p)/(phi * (2 - p));
    j_max = (int) floor(x);
    if (x - j_max  > .5||j_max<1) j_max++; 
    j_max -= j0; /* converted to buffer index */
    
    j = j_max+j0;
    jalogy = j*alogy[i];
    wdW2d2W= wdlogwdp=dWpp=0.0;
    wi=w1i=w2i=0.0; // 1.0;
    wmax = wb[j_max] - jalogy;wmin = wmax + log_eps; 
    //  w1max = wb1[j_max] - jalogy;w1min = w1max + log_eps;    
    // w2max = wb2[j_max] - jalogy;w2min = w2max + log_eps;
 
    /* start upsweep to convergence or end of available buffered values */
    ok = 0;//xmax=x1max=x2max=0.0;
    for (j=j_max+j0,jb=j_max;jb<=j_hi;jb++,j++) { // note initially wi etc initialized to 1 and summation starts 1 later
      jalogy = j * alogy[i];
      wj = wb[jb] - jalogy;
      w1j = wb1[jb];
      wp1j = wp1[jb] - j * logy1p2[i]; /* d log W / dp */
      wp2j = wp2[jb] - 2 * j * logy1p3[i]; /* d^2 log W/ dp^2 */
      /* transform to working parameterization ... */
      wp2j = wp1j * dpth2 + wp2j * dpth1 * dpth1; /* d^2 log W/ dth^2 */
      wp1j *= dpth1; /* d log W / dth */
      wppj = wpp[jb] * dpth1; 

      wj_scaled =  exp(wj-wmax);
      wi += wj_scaled; /* sum of the scaled W_j */

      w1i += wj_scaled * w1j; /* sum W_j dlogW_j / d rho */
      w2i += wj_scaled * w1j*w1j; /*  sum W_j d^2logW_j / d rho^2 */ 
 
      x = wj_scaled*wp1j;
      wdlogwdp += x;    /* sum_j W_j dlogW_j/dp */
      x1 = wj_scaled*(wp1j*wp1j + wp2j);
      wdW2d2W += x1;  /* sum_j  (dlog W_j/dp)^2 + W_j d^2logW_j/dp^2 */     
      
      x2 = wj_scaled*(wp1j*j/onep + wppj);
      dWpp += x2; 
      //  x=fabs(x);x1=fabs(x1);x2=fabs(x2);
 
      // if (x>xmax) {xmax=x;wp1jmin=x * *eps;}
      //if (x1>x1max) {x1max=x1;wdW2min=x1 * *eps;}
      //if (x2>x2max) {x2max=x2;Wppmin=x2 * *eps;}
      if (wj < wmin) { ok=1;break;}
      //&&(w1j < w1min)&&(w2j < w2min)&&
      //  (x < wp1jmin)&&(x1 < wdW2min)&&(x2 < Wppmin)) { ok=1;break;} /* converged on upsweep */
    } /* end of upsweep to buffer end */ 

    while (!ok) { /* while upsweep unconverged need to fill in more buffer */
      for (;jb<jal;jb++,j++) { /* fill buffers and calculate w terms */
        wb[jb] = j * w_base - lgamma((double)j+1) - lgamma(-j * alpha);
        wb1[jb] = -j/onep;
        xx = j/onep2;
        x = xx*digamma(-j*alpha);
        wp1[jb] = j * wp_base + x; /* base for d logW_j/dp */
        xx = trigamma(-j*alpha) * xx * xx;
        wp2[jb] = j * wp2_base + 2*x/onep - xx;
        wpp[jb] = j /onep2;

        jalogy = j * alogy[i];
        wj = wb[jb] - jalogy;
        w1j = wb1[jb];
        wp1j = wp1[jb] - j * logy1p2[i]; /* d log W / dp */
        wp2j = wp2[jb] - 2 * j * logy1p3[i]; /* d^2 log W/ dp^2 */
        /* transform to working parameterization ... */
        wp2j = wp1j * dpth2 + wp2j * dpth1 * dpth1; /* d^2 log W/ dth^2 */
        wp1j *= dpth1; /* d log W / dth */
        wppj = wpp[jb] * dpth1; 
    
        wj_scaled =  exp(wj-wmax);
        wi += wj_scaled; /* sum of the scaled W_j */

        w1i += wj_scaled * w1j; /* sum W_j dlogW_j / d rho */
        w2i += wj_scaled * w1j*w1j; /*  sum W_j d^2logW_j / d rho^2 */ 
 
        x = wj_scaled*wp1j;
        wdlogwdp += x;    /* sum_j W_j dlogW_j/dp */
        x1 = wj_scaled*(wp1j*wp1j + wp2j);
        wdW2d2W += x1;  /* sum_j  (dlog W_j/dp)^2 + W_j d^2logW_j/dp^2 */     
      
        x2 = wj_scaled*(wp1j*j/onep + wppj);
        dWpp += x2; 
    
        if (wj < wmin) { ok=1;break;} /* converged on upsweep */
        
      } 
      j_hi = jb; if (j_hi > jal-1) j_hi = jal-1; /* set j_hi to last element filled */
      if (!ok) { /* need to expand buffer storage*/
        /*Rprintf("forward buffer expansion\n");*/
        wb = forward_buf(wb,&jal,0);
        wb1 = forward_buf(wb1,&jal,0);
	// wb2 = forward_buf(wb2,&jal,0);
        wp1 = forward_buf(wp1,&jal,0);
        wp2 = forward_buf(wp2,&jal,0);
        wpp = forward_buf(wpp,&jal,1);
      }
    } /* finished upsweep and any buffer expansion */
  
    /* start downsweep to convergence or start of available buffered values */
    ok=0;
    for (j=j_max-1+j0,jb=j_max-1;jb>=j_lo;jb--,j--) { 
        jalogy = j * alogy[i];
      wj = wb[jb] - jalogy;
      w1j = wb1[jb];
      wp1j = wp1[jb] - j * logy1p2[i]; /* d log W / dp */
      wp2j = wp2[jb] - 2 * j * logy1p3[i]; /* d^2 log W/ dp^2 */
      /* transform to working parameterization ... */
      wp2j = wp1j * dpth2 + wp2j * dpth1 * dpth1; /* d^2 log W/ dth^2 */
      wp1j *= dpth1; /* d log W / dth */
      wppj = wpp[jb] * dpth1; 

      wj_scaled =  exp(wj-wmax);
      wi += wj_scaled; /* sum of the scaled W_j */

      w1i += wj_scaled * w1j; /* sum W_j dlogW_j / d rho */
      w2i += wj_scaled * w1j*w1j; /*  sum W_j d^2logW_j / d rho^2 */ 
 
      x = wj_scaled*wp1j;
      wdlogwdp += x;    /* sum_j W_j dlogW_j/dp */
      x1 = wj_scaled*(wp1j*wp1j + wp2j);
      wdW2d2W += x1;  /* sum_j  (dlog W_j/dp)^2 + W_j d^2logW_j/dp^2 */     
      
      x2 = wj_scaled*(wp1j*j/onep + wppj);
      dWpp += x2; 

      if (wj < wmin) { ok=1;break;} /* converged on downsweep */
    } /* end of downsweep to buffer end */ 
   
    if (j<=1&&j_lo==0) ok=1; /* don't care about element size if reached base */

    while (!ok) { /* while downsweep unconverged need to fill in more buffer */
      for (jb=j_lo-1;jb>=0;jb--,j--) { /* fill buffers and calculate w terms */
        wb[jb] = j * w_base - lgamma((double)j+1) - lgamma(-j * alpha);
        wb1[jb] = -j/onep;
        xx = j/onep2;
        x = xx*digamma(-j*alpha);
        wp1[jb] = j * wp_base + x; /* base for d logW_j/dp */
        xx = trigamma(-j*alpha) * xx * xx;
        wp2[jb] = j * wp2_base + 2*x/onep - xx;
        wpp[jb] = j /onep2;

        jalogy = j * alogy[i];
        wj = wb[jb] - jalogy;
        w1j = wb1[jb];
        wp1j = wp1[jb] - j * logy1p2[i]; /* d log W / dp */
        wp2j = wp2[jb] - 2 * j * logy1p3[i]; /* d^2 log W/ dp^2 */
        /* transform to working parameterization ... */
        wp2j = wp1j * dpth2 + wp2j * dpth1 * dpth1; /* d^2 log W/ dth^2 */
        wp1j *= dpth1; /* d log W / dth */
        wppj = wpp[jb] * dpth1; 
    
        wj_scaled =  exp(wj-wmax);
        wi += wj_scaled; /* sum of the scaled W_j */

        w1i += wj_scaled * w1j; /* sum W_j dlogW_j / d rho */
        w2i += wj_scaled * w1j*w1j; /*  sum W_j d^2logW_j / d rho^2 */ 
 
        x = wj_scaled*wp1j;
        wdlogwdp += x;    /* sum_j W_j dlogW_j/dp */
        x1 = wj_scaled*(wp1j*wp1j + wp2j);
        wdW2d2W += x1;  /* sum_j  (dlog W_j/dp)^2 + W_j d^2logW_j/dp^2 */     
      
        x2 = wj_scaled*(wp1j*j/onep + wppj);
        dWpp += x2; 
    
        if (wj < wmin) { ok=1;break;} /* converged on upsweep */
      } 

      if (j<=1) ok=1; /* don't care about element size if reached base */

      j_lo = jb; if (j_lo<0) j_lo=0; /* set j_lo to first element filled */
      if (!ok) { /* need to expand buffer storage*/
        /*Rprintf("backward buffer expansion\n");*/
        wb = backward_buf(wb,&jal,&j0,&j_lo,&j_hi,0);
        wb1 = backward_buf(wb1,&jal,&j0,&j_lo,&j_hi,0);
	// wb2 = backward_buf(wb2,&jal,&j0,&j_lo,&j_hi,0);
        wp1 = backward_buf(wp1,&jal,&j0,&j_lo,&j_hi,0);
        wp2 = backward_buf(wp2,&jal,&j0,&j_lo,&j_hi,0);
        wpp = backward_buf(wpp,&jal,&j0,&j_lo,&j_hi,1); /* final '1' updates jal,j0 etc. */
      }

    } /* finished downsweep and any buffer expansion */
    /* Summation now complete: need to do final transformations */
    w[i] = wmax + log(wi);    /* contains log W */
    w2[i] = w2i/wi - (w1i/wi)*(w1i/wi);
    w2p[i] = wdW2d2W/wi - (wdlogwdp/wi)*(wdlogwdp/wi);
    w2pp[i] = (w1i/wi)*(wdlogwdp/wi) + dWpp/wi;
    w1[i] = -w1i/wi;
    w1p[i] =  wdlogwdp/wi;

  } /* end of looping through y */
  FREE(alogy);FREE(wb);FREE(wb1);//FREE(wb2);
  FREE(logy1p2);FREE(logy1p3);FREE(wp1);FREE(wp2);FREE(wpp);
} /* tweedious */


void tweedious2(double *w,double *w1,double *w2,double *w1p,double *w2p,
	       double *w2pp,double *y,double *eps,int *n,
               double *th,double *rho,double *a, double *b)
/* Routine to perform tedious series summation needed for Tweedie distribution
   evaluation, following Dunn & Smyth (2005) Statistics and Computing 15:267-280.
   Notation as in that paper. For 
   
   log W returned in w. (where W means sum_j W_j)
   d logW / drho in w1, 
   d2 logW / d rho2 in w2.
   d logW / dth in w1p
   d2 logW / dth2 in w2p
   d2 logW / dth drho in W2pp

   rho=log(phi), and th defines p = (a + b * exp(th))/(exp(th)+1). 
   note, 1<a<b<2 (all strict) 

   The somewhat involved approach is all about avoiding overflow or underflow. 
   Extensive use is made of 
        log { sum_j exp(x_j)} = log { sum_j exp(x_j-x_max) } + x_max
   digamma and trigamma functions are from Rmath.h

   NOTE: still some redundancy for readability 
 
*/
{ int j_max,i,j,ok,incr;
  double x,x1,x2,xx,alpha,alogy,lgammaj1,
    wbj,wb1j,wp1jb,wp2jb,wppjb,
    wp1j,wp2j,wppj,dpth1,dpth2,
    w_base,wp_base,wp2_base,wj_scaled,wdlogwdp,
    wdW2d2W,dWpp,exp_th,
    wmax,wmin,
    wi,w1i,w2i,
    log_eps,wj,w1j,jalogy,onep,onep2,twop,logy1p2,logy1p3,p,phi;

  log_eps = log(*eps);

  for (i=0;i<*n;i++) { /* loop through y */
    phi = exp(rho[i]);
    if (th[i]>0) { 
        exp_th =  exp(-th[i]);
        x = 1 + exp_th;p = (*b + *a * exp_th)/x;
        x1 = x*x;dpth1 = exp_th*(*b - *a)/x1;
        dpth2 =  ((*a - *b)*exp_th+(*b - *a)*exp_th*exp_th)/(x1*x);
    } else {
        exp_th =  exp(th[i]);
        x = exp_th+1;p = (*b * exp_th + *a)/x;
        x1 = x*x;dpth1 = exp_th*(*b - *a)/x1;
        dpth2 = ((*a - *b)*exp_th*exp_th+(*b - *a)*exp_th)/(x*x1);
    }
    /* first find the location of the series maximum... */
    x = pow(y[i],2 - p)/(phi * (2 - p));
    j_max = (int) floor(x);
    if (x - j_max  > .5||j_max<1) j_max++; 
    
    j = j_max; 
    onep = 1 - p;onep2 = onep * onep;
    twop = 2 - p;
    alpha = twop/onep;
    alogy = log(y[i]);
    logy1p2 = alogy/(onep2);
    logy1p3 = logy1p2/onep;
    alogy *= alpha; /* alpha * log(y[i]) */

    wdW2d2W= wdlogwdp=dWpp=0.0;
    wi=w1i=w2i=0.0;

    /* get terms that are repeated in logWj etc., but simply multiplied by j */
    w_base = alpha * log(-onep) + rho[i]/onep - log(twop);
    wp_base = (log(-onep) + rho[i])/onep2 - alpha/onep + 1/twop;
    wp2_base= 2*(log(-onep) + rho[i])/(onep2*onep) - (3*alpha-2)/(onep2) + 1/(twop*twop);

    wmax = j * w_base - lgamma((double)j+1) - lgamma(-j * alpha) - j*alogy;
    wmin = wmax + log_eps; 
   
    /* start upsweep/downsweep to convergence */
    ok = 0;//xmax=x1max=x2max=0.0;
    //for (j=j_max+j0,jb=j_max;jb<=j_hi;jb++,j++) { // note initially wi etc initialized to 1 and summation starts 1 later
    incr = 1;
    lgammaj1 = lgamma((double)j+1); // lgamma(j+1) to be computed by recursion
    while (!ok) {
      wbj = j * w_base - lgammaj1 - lgamma(-j * alpha);
      wb1j = -j/onep;
      xx = j/onep2;
      x = xx*digamma(-j*alpha);
      wp1jb = j * wp_base + x; /* base for d logW_j/dp */
     
      xx = trigamma(-j*alpha) * xx * xx;
      wp2jb = j * wp2_base + 2*x/onep - xx;
      wppjb = j /onep2;

      jalogy = j * alogy;
      wj = wbj - jalogy;
      w1j = wb1j;
      wp1j = wp1jb - j * logy1p2; /* d log W / dp */
      
      wp2j = wp2jb - 2 * j * logy1p3; /* d^2 log W/ dp^2 */
      
      /* transform to working parameterization ... */
      wp2j = wp1j * dpth2 + wp2j * dpth1 * dpth1; /* d^2 log W/ dth^2 */
      wp1j *= dpth1; /* d log W / dth */
      wppj = wppjb * dpth1; 

      wj_scaled =  exp(wj-wmax);
      wi += wj_scaled; /* sum of the scaled W_j */

      w1i += wj_scaled * w1j; /* sum W_j dlogW_j / d rho */
      w2i += wj_scaled * w1j*w1j; /*  sum W_j d^2logW_j / d rho^2 */ 
 
      x = wj_scaled*wp1j;
      wdlogwdp += x;    /* sum_j W_j dlogW_j/dp */
      //Rprintf("wdlogwdp=%g   wj_scaled=%g  wp1j=%g\n",wdlogwdp,wj_scaled,wp1j);
      x1 = wj_scaled*(wp1j*wp1j + wp2j);
      wdW2d2W += x1;  /* sum_j  (dlog W_j/dp)^2 + W_j d^2logW_j/dp^2 */     
      
      x2 = wj_scaled*(wp1j*j/onep + wppj);
      dWpp += x2; 
      j += incr;
      if (incr>0) { // upsweep
	lgammaj1 += log(j);
        if (wj < wmin) { j = j_max - 1;incr = -1;
	  if (j==0) ok=1; // finished
          lgammaj1 = lgamma((double)j+1);
	} // change to downsweep
      } else {
	lgammaj1 += -log(j+1);
        if (wj < wmin||j<1) ok=1; // finished
      }
    } /* end of upsweep/downsweep */ 
    //Rprintf("wdlogwdp = %g\n",wdlogwdp);
    /* Summation now complete: need to do final transformations */
    w[i] = wmax + log(wi);    /* contains log W */
    w2[i] = w2i/wi - (w1i/wi)*(w1i/wi);
    w2p[i] = wdW2d2W/wi - (wdlogwdp/wi)*(wdlogwdp/wi);
    w2pp[i] = (w1i/wi)*(wdlogwdp/wi) + dWpp/wi;
    w1[i] = -w1i/wi;
    w1p[i] =  wdlogwdp/wi;

  } /* end of looping through y */
} /* tweedious2 */




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


/*******************************************************/
/** Fast re-weighting routines                         */
/*******************************************************/

void rwMatrix(int *stop,int *row,double *w,double *X,int *n,int *p,int *trans,double *work) {
/* Function to recombine rows of n by p matrix X (column ordered).
   ith row of X' is made up of row[stop[i-1]+1...stop[i]], weighted by 
   w[stop[i-1]+1...stop[i]]. stop[-1]=-1 by convention.
   stop is an n vector.     
   
   If (trans==0) the operation on a column x is x'[i] += w[row[j]] * X[row[j]] over the 
   j from stop[i-1]+1 to stop[i]. Otherwise the tranposed operation 
   x'[row[j]] += w[row[j]] * x[i] is used with the same j range. x' zero at outset.

   work is same dimension as X

   See rwMatrix in bam.r for call from R. 
*/
  ptrdiff_t i,j,jump,start=0,end,off;
  double *X1p,*Xp,weight,*Xpe,*X1;
  /* create storage for output matrix, cleared to zero */
  X1 = work; 
  jump = *n;
  off = *n * (ptrdiff_t) *p; 
  for (X1p=X1,Xpe=X1p+off;X1p<Xpe;X1p++) *X1p = 0.0;
  for (i=0;i<*n;i++) { /* loop through rows of output X1 */
    end = stop[i]+1;
    for (j=start;j<end;j++) { /* loop through the input rows */
      if (*trans) {
        X1p = X1 + row[j];
        Xp = X + i;
      } else { 
        X1p = X1 + i;    /* pointer to start of row i of output */
        Xp = X + row[j]; /* pointer to start of source row */
      }
      weight = w[j];   
      for (Xpe=Xp+off;Xp<Xpe;Xp+=jump,X1p+=jump) *X1p += weight * *Xp;
    }
    start = end;
  }
  /* copy output to input for return...*/
  for (Xp=X,X1p=X1,Xpe=Xp+off;Xp<Xpe;Xp++,X1p++) *Xp = *X1p;
}

/* Example code for rwMatrix in R....
   n <- 10;p<-5
   X <- matrix(runif(n*p),n,p)
   ## create transform to take AR1(rho) to independence...
   stop <- c(1:(n-1)*2,2*n-1)
   row <- rep(1:n,rep(2,n))[-1]
   rho <- .7;ld <- 1/sqrt(1-rho^2);sd <- -rho*ld
   w <- c(rep(c(ld,sd),n-1),1)
   mgcv:::rwMatrix(stop,row,w,X)
   
*/


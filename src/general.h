#include <R.h> /* required for R specific stuff */


/* comment in the following and out R.h for normal use */
/*#define Rprintf printf
#define DOUBLE_EPS 2.3e-16 */

/* DOUBLE_EPS is needed because R defaults to using compiler options that use floating point registers 
   in a silly way that e.g. makes it impossible to write strictly repeatable code! */ 

#include <R.h> /* required for R specific stuff */

/* comment in the following and out R.h for non-R use */

/* 
#define Rprintf printf
#define DOUBLE_EPS 2.3e-16 
*/

/* #define MEM_CHECK and compile using the dmalloc library (dmalloc.com) to check for memory leaks 
   For R package compilation this means creating a file Makevars in the src directory containing the 
   line:
        PKG_LIBS = -ldmalloc
   At shell prompt type:
   dmalloc -l /where/to/put/logfile
   (where logfile is the file name). Then start R and use the routines that call compiled code.
   Entry point functions all set to call dmalloc_log_unfreed() on exit, to allow leak chasing.
 
*/

#ifdef MEM_CHECK
#include <dmalloc.h>
#endif

/* DOUBLE_EPS is needed because R defaults to using compiler options that use floating point registers 
   in a silly way that e.g. makes it impossible to write strictly repeatable code! */ 

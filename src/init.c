/* Symbol registration initialization: original provided by Brian Ripley.
   Anything called from R should be registered here (and declared in mgcv.h).
   (See also NAMESPACE:1)
 */ 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "mgcv.h"


R_CMethodDef CEntries[] = {
    {"update_beta", (DL_FUNC) &update_beta, 22},
    {"RMonoCon", (DL_FUNC) &RMonoCon, 7},
    {"RuniqueCombs", (DL_FUNC) &RuniqueCombs, 4},
    {"RPCLS", (DL_FUNC) &RPCLS, 14},
    {"mgcv", (DL_FUNC) &mgcv, 27},
    {"construct_tprs", (DL_FUNC) &construct_tprs, 13},
    {"construct_cr", (DL_FUNC) &construct_cr, 8},
    {"predict_tprs", (DL_FUNC) &predict_tprs, 12},
    {"MinimumSeparation", (DL_FUNC) &MinimumSeparation, 7},
    {"magic", (DL_FUNC) &magic, 18},
    {"mgcv_mmult", (DL_FUNC) &mgcv_mmult,8},
    {"gdi",(DL_FUNC) &gdi,36},
    {"R_cond",(DL_FUNC) &R_cond,5} ,
    {"pls_fit",(DL_FUNC)&pls_fit,10},
    {NULL, NULL, 0}
};

void R_init_mgcv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

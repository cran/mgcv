/* Symbol registration initialization: original provided by Brian Ripley.
   Anything called from R should be registered here. (See also NAMESPACE:1)
 */ 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "mgcv.h"


R_CMethodDef CEntries[] = {
    {"update_beta", (DL_FUNC) &update_beta, 22},
    {"RMonoCon", (DL_FUNC) &RMonoCon, 7},
    {"RuniqueCombs", (DL_FUNC) &RuniqueCombs, 3},
    {"RPCLS", (DL_FUNC) &RPCLS, 14},
    {"mgcv", (DL_FUNC) &mgcv, 27},
    {"construct_tprs", (DL_FUNC) &construct_tprs, 14},
    {"construct_cr", (DL_FUNC) &construct_cr, 8},
    {"predict_tprs", (DL_FUNC) &predict_tprs, 12},
    {"MinimumSeparation", (DL_FUNC) &MinimumSeparation, 7},
    {"magic", (DL_FUNC) &magic, 16},
    {NULL, NULL, 0}
};

void R_init_mgcv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

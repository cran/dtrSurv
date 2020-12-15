#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(gettree)(int *iTree, int *nr, int *nc, double *nodes, double *survFunc, double *mean, double *survProb);
extern void F77_NAME(predictsurvtree)(int *n, int *np, double *xt, int *nCat, int *nt, int *nNodes, double *tsurvFunc, double *mean, double *survProb, int *nCols, double *tnodes, double *predSurvFunc, double *predMean, double *predSurvProb);
extern void F77_NAME(setupbasics)(int *t_nt, double *t_dt, double *t_rs, int *t_ERT, int *t_uniformSplit, int *t_nodeSize, int *t_minEvent, int *t_rule, int *t_sIndex, double *t_sFraction, double *t_stratifiedSplit, int *t_replace);
extern void F77_NAME(setupinners)(int *t_n, int *t_np, double *t_x, double *t_pr, int *t_delta, int *t_mTry, int *t_nCat, int *t_sampleSize, int *t_nTree, int *t_nrNodes);
extern void F77_NAME(survtree)(double *tSurvFunc, double *mean, double *survProb);
extern void F77_NAME(treedim)(int *iTree, int *nr, int *nc);

static const R_FortranMethodDef FortranEntries[] = {
    {"gettree",         (DL_FUNC) &F77_NAME(gettree),          7},
    {"predictsurvtree", (DL_FUNC) &F77_NAME(predictsurvtree), 14},
    {"setupbasics",     (DL_FUNC) &F77_NAME(setupbasics),     12},
    {"setupinners",     (DL_FUNC) &F77_NAME(setupinners),     10},
    {"survtree",        (DL_FUNC) &F77_NAME(survtree),         3},
    {"treedim",         (DL_FUNC) &F77_NAME(treedim),          3},
    {NULL, NULL, 0}
};

void R_init_dtrSurv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


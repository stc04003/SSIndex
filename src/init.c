#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void DouglasAn(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void DouglasUn(void *, void *, void *, void *, void *, void *, void *, void *);
extern void drankSmooth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kappa(void *, void *, void *, void *, void *, void *, void *);
extern void kappa2(void *, void *, void *, void *);
extern void rank(void *, void *, void *, void *, void *, void *, void *);
extern void rankSmooth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void shapeEq(void *, void *, void *, void *);
extern void shapeEq2(void *, void *, void *, void *, void *);
extern void shapeEqSmooth(void *, void *, void *, void *, void *, void *, void *);
extern void shapeFun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void shapeFun2(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"DouglasAn",     (DL_FUNC) &DouglasAn,      9},
    {"DouglasUn",     (DL_FUNC) &DouglasUn,      8},
    {"drankSmooth",   (DL_FUNC) &drankSmooth,   10},
    {"kappa",         (DL_FUNC) &kappa,          7},
    {"kappa2",        (DL_FUNC) &kappa2,         4},
    {"rank",          (DL_FUNC) &rank,           7},
    {"rankSmooth",    (DL_FUNC) &rankSmooth,    10},
    {"shapeEq",       (DL_FUNC) &shapeEq,        4},
    {"shapeEq2",      (DL_FUNC) &shapeEq2,       5},
    {"shapeEqSmooth", (DL_FUNC) &shapeEqSmooth,  7},
    {"shapeFun",      (DL_FUNC) &shapeFun,      10},
    {"shapeFun2",     (DL_FUNC) &shapeFun2,      9},
    {NULL, NULL, 0}
};

void R_init_SSIndex(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

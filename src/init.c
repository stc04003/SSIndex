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
extern void kappa3(void *, void *, void *, void *, void *, void *, void *);
extern void kappa4(void *, void *, void *, void *, void *, void *, void *);
extern void k0Mat(void *, void *, void *, void *, void *, void *);
extern void k02Mat(void *, void *, void *, void *, void *, void *);
extern void givek0s(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rank(void *, void *, void *, void *, void *, void *, void *);
extern void rankSmooth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sizeEq(void *, void *, void *, void *);
extern void sizeEq2(void *, void *, void *, void *, void *);
extern void sizeEqSmooth(void *, void *, void *, void *, void *, void *, void *);
extern void shapeFun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void shapeFun2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void outerC(void *, void *, void *, void *, void *);
extern void outerCsign(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"DouglasAn",     (DL_FUNC) &DouglasAn,      9},
    {"DouglasUn",     (DL_FUNC) &DouglasUn,      8},
    {"drankSmooth",   (DL_FUNC) &drankSmooth,   10},
    {"outerC",        (DL_FUNC) &outerC,         5},
    {"outerCsign",    (DL_FUNC) &outerCsign,     5},
    {"kappa",         (DL_FUNC) &kappa,          7},
    {"kappa2",        (DL_FUNC) &kappa2,         4},
    {"kappa3",        (DL_FUNC) &kappa3,         7},
    {"kappa4",        (DL_FUNC) &kappa4,         7},
    {"givek0s",       (DL_FUNC) &givek0s,       10},
    {"k0Mat",         (DL_FUNC) &k0Mat,          6},
    {"k02Mat",        (DL_FUNC) &k02Mat,         6},
    {"rank",          (DL_FUNC) &rank,           7},
    {"rankSmooth",    (DL_FUNC) &rankSmooth,    10},
    {"sizeEq",        (DL_FUNC) &sizeEq,        4},
    {"sizeEq2",       (DL_FUNC) &sizeEq2,       5},
    {"sizeEqSmooth",  (DL_FUNC) &sizeEqSmooth,  7},
    {"shapeFun",      (DL_FUNC) &shapeFun,      10},
    {"shapeFun2",     (DL_FUNC) &shapeFun2,      9},
    {NULL, NULL, 0}
};

void R_init_SSIndex(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void DouglasAn(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void DouglasUn(void *, void *, void *, void *, void *, void *, void *, void *);
extern void kappa(void *, void *, void *, void *, void *, void *, void *);
extern void rank(void *, void *, void *, void *, void *, void *, void *);
extern void shapeEq(void *, void *, void *, void *);
extern void shapeFun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void shapeFun2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void shapeFun3(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kernal(void*);



static const R_CMethodDef CEntries[] = {
    {"DouglasAn", (DL_FUNC) &DouglasAn,  9},
    {"DouglasUn", (DL_FUNC) &DouglasUn,  8},
    {"kappa",     (DL_FUNC) &kappa,      7},
    {"rank",      (DL_FUNC) &rank,       7},
    {"shapeEq",   (DL_FUNC) &shapeEq,    4},
    {"shapeFun",  (DL_FUNC) &shapeFun,  10},
    {"shapeFun2", (DL_FUNC) &shapeFun2,  9},
    {"shapeFun3", (DL_FUNC) &shapeFun3,  11},
    {"kernal", (DL_FUNC) &kernal,  1},
    {NULL, NULL, 0}
};

void R_init_GSM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

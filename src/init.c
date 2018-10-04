#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void DouglasAn(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void DouglasUn(void *, void *, void *, void *, void *, void *, void *, void *);
extern void rank(void *, void *, void *, void *, void *, void *, void *);
extern void shapeEq(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"DouglasAn", (DL_FUNC) &DouglasAn, 9},
    {"DouglasUn", (DL_FUNC) &DouglasUn, 8},
    {"rank",      (DL_FUNC) &rank,      7},
    {"shapeEq",   (DL_FUNC) &shapeEq,   8},
    {NULL, NULL, 0}
};

void R_init_GSM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

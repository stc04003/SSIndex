#include <R.h>
#include <Rmath.h>
#include <math.h>

// C version of outer
// This function minics R's outer(t1, t2, "<=") but returns a vecter.
void outerC(double *t1, double *t2, int *n1, int*n2, double *result) {
  int i, j, ind;
  ind = 0;
  for (j = 0; j < *n2; j++) {
    for (i = 0; i < *n1; i++) {
      result[ind] = (t1[i] <= t2[j]);
      ind += 1;
    }
  }  
}

// C version of outer
// This function minics R's outer(t1, t2, ">=") but returns a vecter.
void outerC2(double *t1, double *t2, int *n1, int*n2, double *result) {
  int i, j, ind;
  ind = 0;
  for (j = 0; j < *n2; j++) {
    for (i = 0; i < *n1; i++) {
      result[ind] = (t1[i] >= t2[j]);
      ind += 1;
    }
  }  
}

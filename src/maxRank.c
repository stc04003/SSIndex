#include <R.h>
#include <Rmath.h>
#include <math.h>

void rank(int *n, int *m, int *midx,
	  double *tij, double *yi, double *xb, double *result) {
  int i, j, k, l; // id index
  double M = 0.0;
  for (i = 0; i < (*n - 1); i++) {
    for (j = (i + 1); j < *n; j++) {
      for (k = 0; k < m[i]; k++) {
	if (tij[midx[i] + k] <= fmin(yi[i], yi[j])) {
	  for (l = 0; l < m[j]; l++) {
	    if (tij[midx[j] + l] <= fmin(yi[i], yi[j])) {
	      M += 1;
	      if ((xb[i] - xb[j]) * (tij[midx[i] + k] - tij[midx[j] + l]) > 0) {
		result[0] += 1;
	      }
	    }
	  } // end l
	} // end if I(t_ik \le \min(Y_i, Y_j))
      } // end k
    }
  }
  result[0] = result[0] / M;
}

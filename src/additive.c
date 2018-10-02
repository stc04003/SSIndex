#include <R.h>
#include <Rmath.h>
#include <math.h>

// The estimating equations used in Douglas et al. (2006)
// Un equation
void Douglas.Un(int *n, int *p, int *m, int *midx,
	double *tij, double *yi, double *X,
	double *result) {
  int i, j, k, r; // id index
  double de = 0.0;
  double *nu = Calloc(*p, double);
  for (i = 0; i < *n; i++) {
    for (j = 0; j < m[i]; j++) {
      for (k = 0; k < *n; k++) {
	if (yi[k] > tij[midx[i] + j]) {
	  de += 1;
	  for (r = 0; r < *p; r++) {
	    nu[r] += X[*n * r + k];
	  }
	} // end if
      } // end for k
      for (r = 0; r < *p; r++) {
	if (de > 0) result[r] += X[*n * r + i] - nu[r] / de;
	if (de == 0) result[r] += X[*n * r + i];
	nu[r] = 0;
      }
      de = 0;
    } // end for j
  }
  Free(nu);
  result;
}

// An equation
void Douglas.An(int *n, int *p, int *m, int *midx,
	double *tij, double *yi, double *X, double *dY,
	double *result) {
  int i, j, k, r; // id index
  double de = 0.0;
  double *nu = Calloc(*p, double);
  for (i = 0; i < *n; i++) {
    for (k = 0; k < *n; k++) {
      if (yi[k] > yi[i]) {
	de += 1;
	for (r = 0; r < *p; r++) {
	  nu[r] += X[*n * r + k];
	}
      } // end if
    } // end for k
    for (r = 0; r < *p; r++) {
      for (j = 0; j < *p; j++) {
	if (de > 0) {
	  result[r * *p + j] += dY[i] * (X[*n * r + i] - nu[r] / de) * (X[*n * j + i] - nu[j] / de);
	}
	if (de == 0) result[r * *p + j] += dY[i] * X[*n * r + i] * X[*n * j + i];
      }
    }
    for (r = 0; r < *p; r++) {nu[r] = 0;}
    de = 0;
  }
  Free(nu);
  result;
}

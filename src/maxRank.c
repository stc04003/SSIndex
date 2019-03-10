#include <R.h>
#include <Rmath.h>
#include <math.h>

// This is the C version of the maximum rank estimating euqation.
// Equation (2.1) in the paper

void rank(int *n, int *m, int *midx,
	  double *tij, double *yi, double *xb, double *result) {
  int i, j, k, l; // id index
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      if (i != j) {
	for (k = 0; k < m[i]; k++) {
	  if (tij[midx[i] + k] <= yi[j]) {
	    for (l = 0; l < m[j]; l++) {
	      if (tij[midx[j] + l] <= yi[i]) {
		if (xb[i] > xb[j] && tij[midx[i] + k] > tij[midx[j] + l]) {
		  result[0] += 1;
		}
	      }
	    }
	  } // end l
	} // end if I(t_ik \le \min(Y_i, Y_j))
      } // end k
    }
  }
  result[0] = result[0] / *n / (*n - 1);
}

// used for hypothesis test: H_0: \beta_0 = 0; the formula is similar to that in `rank`
void kappa(int *n, int *m, int *midx,
	  double *tij, double *yi, double *xb, double *result) {
  int i, j, k, l; // id index
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      if (i != j) {
	for (k = 0; k < m[i]; k++) {
	  if (tij[midx[i] + k] <= yi[j]) {
	    for (l = 0; l < m[j]; l++) {
	      if (tij[midx[j] + l] <= yi[i]) {
		if (xb[i] > xb[j] && tij[midx[i] + k] > tij[midx[j] + l]) {
		  result[0] += 1;
		}
		if (xb[i] < xb[j] && tij[midx[i] + k] > tij[midx[j] + l]) {
		  result[0] -= 1;
		}
	      }
	    }
	  } // end l
	} // end if I(t_ik \le \min(Y_i, Y_j))
      } // end k
    }
  }
  result[0] = result[0] / n[0] / (n[0] - 1);
}


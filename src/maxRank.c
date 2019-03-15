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
	      if (tij[midx[j] + l] <= yi[i] && xb[i] > xb[j] && tij[midx[i] + k] > tij[midx[j] + l]) {
		  result[0] += 1;
	      }
	      if (tij[midx[j] + l] <= yi[i] && xb[i] == xb[j] && tij[midx[i] + k] > tij[midx[j] + l]) {
		result[0] += 0.5;
	      } // accomdiate bootstrap
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

// Smooth rank equation??
double get_rikjl(double *X, double *sigma,
		 int *N, int *p, int ik_idx, int jl_idx) {
  double *xdif = Calloc(*p, double);
  double rikjl = 0.0;
  int m = 0, q = 0;
  for (m = 0; m < *p; m++) {
    xdif[m] = 0.0;
    xdif[m] = X[ik_idx + m * *N] - X[jl_idx + m * *N];
  }
  for (m = 0; m < *p; m++) {
    for (q = 0; q < *p; q++) {
      rikjl += xdif[m] * sigma[m * *p + q] * xdif[q];
    }
  }
  rikjl = sqrt(rikjl);
  Free(xdif);
  return(rikjl);
}

void rankSmooth(int *n, int *p, int *m, int *midx,
		double *sigma, 
		double *tij, double *yi, double *xb, double *xx, double *result) {
  int i, j, k, l; // id index
  double rijkl;
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      if (i != j) {
	for (k = 0; k < m[i]; k++) {
	  if (tij[midx[i] + k] <= yi[j]) {
	    for (l = 0; l < m[j]; l++) {
	      if (tij[midx[j] + l] <= yi[i] && tij[midx[i] + k] > tij[midx[j] + l]) {
		rijkl = 0;
		rijkl = get_rikjl(xx, sigma, n, p, i, j);
		// result[0] += pnorm(sqrt(n[0]) * (xb[i] - xb[j]) / rijkl, 0.0, 1.0, 1, 0);
		if (rijkl != 0) {
		  result[0] += pnorm(sqrt(n[0]) * (xb[i] - xb[j]) / rijkl, 0.0, 1.0, 1, 0);
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

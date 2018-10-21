#include <R.h>
#include <Rmath.h>
#include <math.h>

// n is the # of sbjects
// m is the # of recurrent events for each subject
// midx is an integer id strings that indicates where t_ij starts
// tij is the jth recurrent event time for the ith subject
// yi is the censoring time for the ith subject
// xb is a vector with the ith element being \beta^\top Z_i

double kernal(double dx) {
  double out;
  if (dx <= 1 && dx >= -1) {
    // out = 15 * (1 - dx * dx) * (1 - dx * dx) / 16; // Quartic/biweight
    out = 3 * (1 - dx * dx) / 4; // Epanechnikov
  }
  return(out);
}

// This function implements \hat{F}(t, x, \hat\beta)
void shapeFun(int *n, int *m, int *midx, double *tij, double *yi, double *xb,
	      double *x, double *t, double *result) {
  int i, j, k, l;
  double nu = 0.0;
  double de = 0.0;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < m[i]; k++) {
      if (tij[midx[i] + k] > t[0]) {
  	de = 0.0;
  	nu = 0.0;
  	nu = kernal(x[0] - xb[i]);
  	for (j = 0; j < *n; j++) {
  	  for (l = 0; l < m[j]; l++) {
  	    if (tij[midx[i] + k] >= tij[midx[j] + l] && tij[midx[i] + k] <= yi[j])
  	      de += kernal(x[0] - xb[j]);
  	  }
  	}
	result[0] += nu / de;
      }
    }
  }
}

void shapeEq(int *n, double *xr, double *mFhat, double *result) {
  int i, j;
  for (i = 0; i < *n; i++) {
    // tmp = shapeFun(n, m, midx, tij, yi, xb, xb[i], yi[i]);    
    for (j = 0; j < *n; j++) {
      if (xr[i] > xr[j]) {
	result[0] += mFhat[i];
      }
    }
  }
}

/* void Sn(int *n, int *p, int *m, int *midx, */
/* 	double *tij, double *yi, double *xr, double *xb, double *X, */
/* 	double *result) { */
/*   int i, j, k, r; // id index */
/*   double tmp; */
/*   double de = 0.0; */
/*   double *nu = Calloc(*p, double); */
/*   for (i = 0; i < *n; i++) { */
/*     for (j = 0; j < m[i]; j++) { */
/*       for (k = 0; k < *n; k++) { */
/* 	if (yi[k] >= tij[midx[i] + j]) { */
/* 	  tmp = xb[k] - xb[i]; */
/* 	  de += xr[k] * kernal(tmp); */
/* 	  for (r = 0; r < *p; r++) { */
/* 	    nu[r] += X[*n * r + k] * xr[k] * kernal(tmp); */
/* 	  } */
/* 	} */
/*       } // end for k */
/*       for (r = 0; r < *p; r++) { */
/* 	if (de > 0) result[r] += X[*n * r + i] - nu[r] / de; */
/* 	if (de == 0) result[r] += X[*n * r + i]; */
/* 	nu[r] = 0; */
/*       } */
/*       de = 0; */
/*     } */
/*   } */
/*   Free(nu); */
/*   result; */
/* } */

// gives \tau_n(\theta) in the variance estimation (Section 2.3)
void tauN(int *n, int *m, int *midx,
	  double *tij, double *yi, double *xb, double *result) {
  int i, j, k, l; // id index
  double *e1 = Calloc(*n, double);
  double *e2 = Calloc(*n, double);
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      for (k = 0; k < m[i]; k++) {
	if (tij[midx[i] + k] <= fmin(yi[i], yi[j])) {
	  for (l = 0; l < m[j]; l++) {
	    if (tij[midx[j] + l] <= fmin(yi[i], yi[j])) {
	      if ((xb[i] - xb[j]) * (tij[midx[i] + k] - tij[midx[j] + l]) > 0) {
		e1[i] += 1;
		e2[j] += 1;
	      }
	    }
	  } // end l
	} // end if I(t_ik \le \min(Y_i, Y_j))
      } // end k
    }
  }
  for (i = 0; i < *n; i++) {
    result[0] = (e1[i] + e2[j]) / 2;
  }
  Free(e1);
  Free(e2);
}

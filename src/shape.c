#include <R.h>
#include <Rmath.h>
#include <math.h>

// n is the # of sbjects
// m is the # of recurrent events for each subject
// midx is an integer id strings that indicates where t_ij starts
// tij is the jth recurrent event time for the ith subject
// yi is the censoring time for the ith subject
// xb is a vector with the ith element being \beta^\top Z_i

/* double kernal(double dx) { */
/*   double out = 0.0; */
/*   out = exp(-1 * dx * dx / 2) / 2.506628; */
/*   return(out); */
/* } */

double kernal(double dx) {
  // I added "= 0.0", not sure if this matters or not
  double out=0.0;
  if ((dx <= 1.0) && (dx >= -1.0)) {
    // out = 15 * (1 - dx * dx) * (1 - dx * dx) / 16; // Quartic/biweight
    out = 3 * (1 - dx * dx) / 4; // Epanechnikov
  }
  return(out);
}

// This function implements \hat{F}(t, x, \hat\beta)
// returns one value at \hat{F}(t, x, \hat\beta)
void shapeFun(int *n, int *m, int *midx, double *tij, double *yi, double *xb,
	      double *x, double *t, double *h, double *result) {
  int i, j, k, l;
  double nu = 0.0;
  double de = 0.0;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < m[i]; k++) {
      if (tij[midx[i] + k] >= t[0]) {
  	de = 0.0;
  	nu = 0.0;
  	nu = kernal((x[0] - xb[i]) / h[0]);
  	for (j = 0; j < *n; j++) {
  	  for (l = 0; l < m[j]; l++) {
  	    if (tij[midx[i] + k] >= tij[midx[j] + l] && tij[midx[i] + k] <= yi[j])
  	      de += kernal((x[0] - xb[j]) / h[0]);
  	  }
  	}
  	if(de == 0) {
  	  result[0] += 0;
  	} else {
  	  result[0] += nu / de;
  	}
      }
    }
  }
}


// This function implements \hat{F}(t, x, \hat\beta)
// assumes fixed x
// returns a vector for
// \frac{K_h(\beta_0^\top Z_i - x)}{\sum_{j=1}^n\sum_{l=1}^{m_i}K_h(\beta_0^\top Z_j - x) I(T_{jl} \le T_{ik} \le C_j)}
// for unique(tij)
void shapeFun2(int *n, int *m, int *midx, double *tij, double *yi, double *xb,
               double *x, double *h, double *result) {
  int i, j, k, l, tind;
  double nu = 0.0;
  double de = 0.0;
  tind = 0;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < m[i]; k++) {
      de = 0.0;
      nu = 0.0;
      nu = kernal((x[0] - xb[i]) / h[0]);
      for (j = 0; j < *n; j++) {
        for (l = 0; l < m[j]; l++) {
          if (tij[midx[i] + k] >= tij[midx[j] + l] && tij[midx[i] + k] <= yi[j])
            de += kernal((x[0] - xb[j]) / h[0]);
        }
      }
      if (de == 0) {
	result[tind] += 0;
      } else { 
	result[tind] += nu / de;
      }
      tind += 1;
    }
  }
}

// return the estimating equation for \gamma_0
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

// return the estimating equation for \gamma_0
// This is `shapeEq` but with a small \tau_0
void shapeEq2(int *n, double *xr, double *mFhat, double *y, double *result) {
  int i, j;
  for (i = 0; i < *n; i++) {
    // tmp = shapeFun(n, m, midx, tij, yi, xb, xb[i], yi[i]);    
    for (j = 0; j < *n; j++) {
      if (xr[i] > xr[j] && y[i] >= 0.03 && y[j] >= 0.03) {
	result[0] += mFhat[i];
      }
    }
  }
}

// return \hat r(t,x,\beta) from SS1117 (paper #2)
// h2 is the bandwidth for smoothing on the time scale
// I think this function can be improved 
void shapeFun3(int *n, int *m, int *midx, double *tij, double *yi, double *xb,
               double *x, double *t, double *h, double *h2, double *result) {
  int i, j, k, l;
  double nu = 0.0;
  double de = 0.0;

  for (i = 0; i < *n; i++) {
    for (k = 0; k < m[i]; k++) {
      de = 0.0;
      nu = 0.0;
      nu = kernal((x[0] - xb[i]) / h[0]);
      // new addition:
      nu = nu * kernal((tij[midx[i]+k] - t[0]) / h2[0]);
      for (j = 0; j < *n; j++) {
        for (l = 0; l < m[j]; l++) {
          if (tij[midx[i] + k] >= tij[midx[j] + l] && tij[midx[i] + k] <= yi[j])
            de += kernal((x[0] - xb[j]) / h[0]);
        }
      }
      if(de == 0)
      {
        result[0] += 0;
      }else{
        result[0] += nu / de;
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

void kappa2(int *n, double *xb, double *mFhat, double *result) {
  int i, j;
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      if (i != j) {
	if (xb[i] > xb[j]) {
	  result[0] += mFhat[i] - mFhat[j];
	}
	if (xb[i] <= xb[j]) {
	  result[0] += mFhat[j] - mFhat[i];	  
	}
      }
    }
  }
  result[0] = result[0] / n[0] / (n[0] - 1);
}


void kappa3(int *n, int *m, int *midx, double *xb, double *tij, double *yi, double *result) {
  int i, j, k;
  double N1, N2;
  for (i = 0; i < (*n - 1); i++) {
    for (j = (i + 1); j < *n; j++) {
      N1 = 0.0;
      N2 = 0.0;
      for (k = 0; k < m[i]; k++) {
	if (tij[midx[i] + k] <= yi[i] & tij[midx[i] + k] <= yi[j]) {
	  N1 += 1;
	}
      }
      for (k = 0; k < m[j]; k++) {
	if (tij[midx[j] + k] <= yi[i] & tij[midx[j] + k] <= yi[j]) {
	  N2 += 1;
	}
      }
      if (xb[i] > xb[j]) {
	result[0] += N1 - N2;
      }
      if (xb[i] <= xb[j]) {
	result[0] += N2 - N1;
      }
    }
  }
  result[0] = 2 * result[0] / n[0] / (n[0] - 1);
}

// gives k0 and k02 at once
void kappa4(int *n, int *m, int *midx, double *xb, double *tij, double *yi, double *result) {
  int i, j, k, l;
  double N1, N2;
  for (i = 0; i < (*n - 1); i++) {
    for (j = (i + 1); j < *n; j++) {
      N1 = 0.0;
      N2 = 0.0;
      for (k = 0; k < m[i]; k++) {
	if (tij[midx[i] + k] <= yi[j]) {
	  N1 += 1;
	}
      }
      for (k = 0; k < m[j]; k++) {
	if (tij[midx[j] + k] <= yi[i]) {
	  N2 += 1;
	}
      }
      if (xb[i] > xb[j]) {
	result[0] += N1 - N2;
      }
      if (xb[i] < xb[j]) {
	result[0] += N2 - N1;
      }
      for (k = 0; k < m[i]; k++) {
	for (l = 0; l < m[j]; l++) {
	  if (tij[midx[i] + k] <= yi[j] &&
	      tij[midx[j] + l] <= yi[i] &&
	      xb[i] > xb[j] && 
	      tij[midx[i] + k] > tij[midx[j] + l]) {
	    result[1] += 1;
	  }
	  if (tij[midx[i] + k] <= yi[j] &&
	      tij[midx[j] + l] <= yi[i] &&
	      xb[i] < xb[j] && 
	      tij[midx[i] + k] > tij[midx[j] + l]) {
	    result[1] -= 1;
	  }
	  if (tij[midx[i] + k] <= yi[j] &&
	      tij[midx[j] + l] <= yi[i] &&
	      xb[i] > xb[j] && 
	      tij[midx[i] + k] < tij[midx[j] + l]) {
	    result[1] -= 1;
	  }
	  if (tij[midx[i] + k] <= yi[j] &&
	      tij[midx[j] + l] <= yi[i] &&
	      xb[i] < xb[j] && 
	      tij[midx[i] + k] < tij[midx[j] + l]) {
	    result[1] += 1;
	  }
	}
      }
    }
  }
  result[0] = 2 * result[0] / n[0] / (n[0] - 1);
  result[1] = 2 * result[1] / n[0] / (n[0] - 1);
}

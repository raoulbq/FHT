//******************************************
//NAME: Robert Taintor
//DATE: 8/3/07
//PURPOSE: Performs a fast polynomial transform on a group of orthogonal
//         polynomials which satisfy a 3 term recurrance
//         from uniformly sampled points.

#include <string.h>
#include <math.h>
#include <fftw3.h>

#include "lalgebra.h"
#include "fht.h"
#include "parameters.h"


// Defines the data points
double xk(int k) {
  return (double) (((double)(k - LITN) / (double)LITN) * BIGC);
}

// Use Zl to calculate Zl + 1
// Z0 = Zl
// Z1 = Zl + 1
// n is the length of Z0
void calculateFirstZ(double *Z0, double *Z1, int n)
{
  memset(Z1, 0, sizeof(double) * n);
  int i;
  // Normalize the Chebyshev result
  for(i=0; i<n; ++i)
    Z0[i] *= D0;
  // Now compute Z1
  for(i=1; i<(n-2); ++i) {
    Z1[i] = AL(0) * (1.0/ALPHA * Z0[i+1] - BETA/ALPHA * Z0[i] - GAMMA/ALPHA * Z0[i-1]) + BL(0)*Z0[i];
  }
}

// Define An(l) as in the paper
// Result is a 8*n sized vector with each 2n representing a circulant block
void createAn(int n, int l, double *result) {
  // Top left is all zeros (we'll zero out everything else while we're at it.
  memset(result, 0, sizeof(double) * 8*n);
  // Top right is I2n
  result[2*n] = 1;
  // Bottom left is cl*I2n
  result[4*n] = CL(l);
  // Bottom right is Cn(wl,vl,ul);
  result[6*n]   = VL(l);
  result[6*n+1] = WL(l);
  result[8*n-1] = UL(l);
}

// Needed columns of circulant matrices listed vertically as top left, top right, bottom left, bottom right.
void precomputeRnAndStore(int n, int l, double* result) {
  int i;
  double *temp = (double*) fftw_malloc(sizeof(double) * 8*n);
  double *temp2 = (double*) fftw_malloc(sizeof(double) * 8*n);
  createAn(n, n/2+l, result);
  for(i=l+n/2-1; i>l; --i) {
    createAn(n, i, temp);
    fourBcirculantSqMatrixMultiply(result, temp, 4*n, temp2);
    memcpy(result, temp2, sizeof(double) * 8*n);
  }
  fftw_free(temp);
  fftw_free(temp2);
}

// Perform a Chebyshev transform in the most naive way possible directly from the data points defined by xl()
void naiveChebyshev(double *data, double *results) {
  double Lminus1;
  double Lminus2;
  double curVal;
  int x, y;
  memset(results, 0, sizeof(double)*BIGN);
  // For each data point
  for(x=0; x<=2*LITN; ++x) {
    Lminus1 = xk(x) / BIGC;
    Lminus2 = 2.0*pow(Lminus1, 2) - 1.0;
    // Go through the n Chebyshev polynomials
    for(y=0; y<BIGN; ++y) {
      curVal = (ALPHA * xk(x) + BETA) * Lminus1 + GAMMA * Lminus2;
      results[y] += data[x] * curVal;
      Lminus2 = Lminus1;
      Lminus1 = curVal;
    }
  }
}

// A recursive function which performs a Hermite transform in O(n(logn)^2) time
// Z0 and Z1 must be precomputed as defined in the paper.
// The value of l should be first set to 1.
// You must precompute all the necessary Rns.
void performTransform(double* Z0, double* Z1, int n, int l, double* result) {
  result[l-1] = Z0[n-1];
  result[l]   = Z1[n-1];

  if (n < 3) return;

  // temp to store the new data
  double* temp = (double*) fftw_malloc(sizeof(double) * 4*n);

  // Combine Z0 and Z1 into Z to get ready for the matrix multiply
  double *Z;
  Z = (double*) fftw_malloc(sizeof(double) * 4*n);
  memcpy(Z, Z0, sizeof(double) * 2*n);
  memcpy(Z+n*2, Z1, sizeof(double) * 2*n);
  preFourBcirculantVcMatrixMultiply(n, l-1, Z, temp);

  fftw_free(Z);
  int nover2 = n / 2;
  performTransform(Z0+nover2, Z1+nover2, nover2, l, result);
  performTransform(temp+nover2, temp+5*nover2, nover2, l+nover2, result);
  fftw_free(temp);
  return;
}

// Performs a Hermite transform in the most naive way possible directly from the data points given in xl()
void naiveTransform(double *data, double *results) {
  memset(results, 0, sizeof(double) * BIGN);
  double Lminus1;
  double Lminus2;
  double curVal;
  int x, y;
  // For each data point
  for(x=0; x<=2*LITN; ++x) {
    Lminus1 = 0;
    curVal = D0;
    // Go through the Hermite polynomials
    for(y=0; y<BIGN; ++y) {
      results[y] += data[x] * curVal;
      Lminus2 = Lminus1;
      Lminus1 = curVal;
      curVal = (AL(y) * xk(x) + BL(y)) * Lminus1 + CL(y) * Lminus2;
    }
  }
}

void oneDTransform(double *data, double* result) {
  int i;
  double *Z0 = (double*) fftw_malloc(sizeof(double) * 2*BIGN);
  double *Z1 = (double*) fftw_malloc(sizeof(double) * 2*BIGN);

  // Do a Chebyshev Transform
  naiveChebyshev(data, Z0+BIGN-1);
  Z0[2*BIGN - 1] = 0;

  // Expand the data
  for(i=0; i<BIGN; ++i)
    Z0[i] = Z0[2*BIGN-i-2];

  // Find the next data point
  calculateFirstZ(Z0, Z1, 2*BIGN);

  // Do the second part
  performTransform(Z0, Z1, BIGN, 1, result);

  fftw_free(Z0);
  fftw_free(Z1);
}

void twoDTransform(double *data, int n, double* result) {
  int i;

  // Initialize the plans for the Fourier Transforms
  initFastFouriers(n);
  initRns(n);

  // Transform the rows
  for(i=0; i<=2*n; ++i) {
    oneDTransform(&data[n*i], &result[n*i]);
  }
  //transform the columns
  double* column = (double*) fftw_malloc(sizeof(double) * n);
  double* temp = (double*) fftw_malloc(sizeof(double) * n);
  for(i=0; i<n; ++i) {
    getColumnFromSq(result, column, i, n);
    oneDTransform(column, temp);
    setColumnToSq(temp, result, i, n);
  }
  fftw_free(column);
  fftw_free(temp);

  destroyFastFouriers(n);
}

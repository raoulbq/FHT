//******************************************
//NAME: Robert Taintor
//DATE: 8/3/07
//PURPOSE: Performs a fast polynomial transform on a group of orthogonal
//         polynomials which satisfy a 3 term recurrance
//         from uniformly sampled points.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>
#include <string.h>

#include "lalgebra.h"
#include "FHT.h"

// Recursion parameters
#define PI 3.1415926535897932385
#define BIGN 128
#define LITN (BIGN*5)

#define BIGC sqrt((double)2*BIGN + 1)

#define ALPHA (2.0 / BIGC)
#define BETA 0
#define GAMMA (-1.0)

#define D0 ((double)1.0 / pow(PI, (double)1.0/4.0))

// Define al, bl, cl, ul, vl and wl
#define AL(l) (double)sqrt((double)2.0 / (l+1.0))
#define BL(l) (double)0
#define CL(l) (double)((-1.0) * sqrt((double)l / (l+1.0)))

#define UL(l) (double)(AL(l) / ALPHA)
#define VL(l) (double)(BL(l) + (-1) * (UL(l) * BETA))
#define WL(l) (double)((-1) * GAMMA * UL(l))


// Defines the data points
fftw_complex xk(int k) {
  return (((double)(k - LITN) / (double)LITN) * BIGC);
}

// Scan filename which should be a text file of long floats seperated by carriage returns.
void oneDFileReader(char *filename, int n, double *result) {
  FILE *inputFile;
  int x;

  inputFile = fopen(filename, "r");
  if(inputFile != NULL) {
    for(x=0; x<n; ++x) {
      fscanf(inputFile, "%lf", &result[x]);
    }
  }
  else
    exit(0);
  fclose(inputFile);
}

// Scan filename which should be a text file of long floats seperated by tabs within a line and carriage returns between lines
void twoDFileReader(char *filename, int n, double *result) {
  FILE *inputFile;
  int i, j;

  inputFile = fopen(filename, "r");
  if(inputFile != NULL) {
    for(i=0; i<n; ++i) {
      for(j=0; j<n-1; ++j) {
        fscanf(inputFile, "%lf", &result[n*i+j]);
      }
      fscanf(inputFile, "%lf\n", &result[(n)*i+n-1]);
    }
  }
  else
    exit(0);
  fclose(inputFile);
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
void naiveChebyshev(double *data, fftw_complex *results) {
  memset(results, 0, sizeof(fftw_complex)*BIGN);
  fftw_complex Lminus1;
  fftw_complex Lminus2;
  fftw_complex curVal;
  int x, y;
  for(x=0; x<=2*LITN; ++x) { //for each data point
    Lminus1 = (double)xk(x) / BIGC;     // x
    Lminus2 = 2.0*pow(Lminus1, 2) - 1;  // 2*x^2 - 1
    for(y=0; y<BIGN; ++y) { //go through the n Chebyshev polynomials
      curVal = (ALPHA * xk(x) + BETA) * Lminus1 + GAMMA * Lminus2;
      results[y] += ((fftw_complex)data[x]) * curVal;
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
  fftw_complex Lminus1;
  fftw_complex Lminus2;
  fftw_complex curVal;
  int x, y;
  for(x=0; x<=2*LITN; ++x) {//for each data point
    Lminus1 = 0;
    curVal = D0;
    for(y=0; y<BIGN; ++y) {//go through the n hermites
      results[y] += (data[x]) * curVal;
      Lminus2 = Lminus1;
      Lminus1 = curVal;
      curVal = (AL(y) * (xk(x)) + BL(y)) * Lminus1 + CL(y) * Lminus2;
    }
  }
}

void oneDTransform(double *data, double* result) {
  int i;
  int n = BIGN;
  fftw_complex *Z0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*n);
  fftw_complex *Z1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*n);
  double *dblZ0 = (double*) fftw_malloc(sizeof(fftw_complex) * 2*n);
  double *dblZ1 = (double*) fftw_malloc(sizeof(fftw_complex) * 2*n);

  // Do a chebyshev Transform
  naiveChebyshev(data, Z0+n-1);
  Z0[2*n - 1] = 0;

  // We only want the real parts
  for(i=0; i<n; ++i) {
    dblZ0[n-1+i] = creal(Z0[n-1+i]);
  }

  // Expand the data
  for(i=0; i<n; ++i)
    dblZ0[i] = dblZ0[2*n-i-2];

  // Find the next data point
  calculateFirstZ(dblZ0, dblZ1, 2*n);

  // Do the second part
  performTransform(dblZ0, dblZ1, n, 1, result);

  fftw_free(Z0);
  fftw_free(dblZ0);
  fftw_free(Z1);
  fftw_free(dblZ1);
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


int main(int argc, char *argv[]) {
  clock_t fclock, nclock;
  int n = LITN;
  int N = BIGN;
  int i;
  double *data = (double*) fftw_malloc(sizeof(double) * (2*n+1));
  double *diag = (double*) fftw_malloc(sizeof(double) * (2*n+1) * (2*n+1));
  double *naiveResult = (double*) fftw_malloc(sizeof(double) * N);
  double *fancyResult = (double*) fftw_malloc(sizeof(double) * N);


  // Load input data
  oneDFileReader("doc/data.txt", (2*n+1), data);
  memset(diag, 0, sizeof(double) * (2*n+1) * (2*n+1));
  for(i=0; i<=2*n; ++i) {
    diag[(2*n + 1) * i + i] = exp(0 - pow(xk(i), 2) / 2);
    data[i] *= diag[(2*n + 1) * i + i];
  }


  // Perform transformations
  initFastFouriers(N);
  initRns(N);

  nclock = clock();
  naiveTransform(data, naiveResult);
  nclock -= clock();

  fclock = clock();
  oneDTransform(data, fancyResult);
  fclock -= clock();

  destroyFastFouriers(N);


  // Save output data
  FILE *resultsn;
  FILE *resultsf;
  resultsn = fopen("doc/resultsn.txt", "wt+");
  resultsf = fopen("doc/resultsf.txt", "wt+");

  for(i=0; i<N; ++i){
    fancyResult[i] *= (2*BIGC) / n;
    naiveResult[i] *= (2*BIGC) / n;
    fprintf(resultsn, "%.16lf\n", naiveResult[i]);
    fprintf(resultsf, "%.16lf\n", fancyResult[i]);
  }
  fclose(resultsn);
  fclose(resultsf);


  fprintf(stdout, "Naive transform took %li\n", -nclock);
  fprintf(stdout, "Fancy transform took %li\n", -fclock);
  fprintf(stdout, "Clocks per second: %li\n", CLOCKS_PER_SEC);

  return 0;
}

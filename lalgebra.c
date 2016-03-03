#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>

#include "lalgebra.h"
#include "FHT.h"


fftw_plan *daPlans;
double *daRns;
int daRnsSize;


void getColumnFromSq(double* input, double* output, int whichColumn, int n) {
  int i;
  for(i=0; i<n; ++i) {
    output[i] = input[n*i + whichColumn];
  }
}

void setColumnToSq(double* input, double* output, int whichColumn, int n) {
  int i;
  for(i=0; i<n; ++i) {
    output[n*i + whichColumn] = input[i];
  }
}

void initFastFouriers(int n) {
  int i;
  daPlans = (fftw_plan*) fftw_malloc(sizeof(fftw_plan) * 2*n+1);
  for(i=8; i<=2*n; i*=2) {
    double* in = (double*) fftw_malloc(sizeof(double) * i);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * i);
    daPlans[i] = fftw_plan_dft_r2c_1d(i, in, out, FFTW_ESTIMATE);
    fftw_free(in);
    fftw_free(out);
  }
}

void destroyFastFouriers(int n) {
  int i;
  for(i=8; i<=2*n; i*=2) {
    fftw_destroy_plan(daPlans[i]);
  }
  fftw_free(daPlans);
}

void initRns(int n) {
  int i, j;
  daRns	= (double*) fftw_malloc(sizeof(double) * (n+1) * n * 8*n);
  daRnsSize = n;
  // Precompute the necessary Rns
  for(i=n; i>=4; i/=2) {
    for(j=0; j<n; j+=i) {
      precomputeRnAndStore(i, j, &daRns[8 * n * (n*i + j)]);
    }
  }
}

// Multiply two 2x2 block square matrices
// A B  *  E F  = AE + BG  AF + BH
// C D     G H    CE + DG  CF + DH
void fourBcirculantSqMatrixMultiply(double* M1, double* M2, int n, double* result) {
  double *A, *B, *C, *D, *E, *F, *G, *H;
  double* temp1 = (double*) fftw_malloc(sizeof(double) * n/2);
  double* temp2 = (double*) fftw_malloc(sizeof(double) * n/2);
  int i;

  // Fill up the columns
  A = M1;
  E = M2;
  B = M1 + n/2;
  F = M2 + n/2;
  C = M1 + n;
  G = M2 + n;
  D = M1 + 3*n/2;
  H = M2 + 3*n/2;

  // Top left: A*E + B*G
  circulantVcMatrixMultiply(A, E, n/2, temp1);
  circulantVcMatrixMultiply(B, G, n/2, temp2);
  for(i=0; i<n/2; ++i)
    result[i] = temp1[i] + temp2[i];

  // Top right: A*F + B*H
  circulantVcMatrixMultiply(A, F, n/2, temp1);
  circulantVcMatrixMultiply(B, H, n/2, temp2);
  for(i=0; i<n/2; ++i)
    result[i+n/2] = temp1[i] + temp2[i];

  // Bottom left: C*E + D*G
  circulantVcMatrixMultiply(C, E, n/2, temp1);
  circulantVcMatrixMultiply(D, G, n/2, temp2);
  for(i=0; i<n/2; ++i)
    result[i+n] = temp1[i] + temp2[i];

  // Bottom right: C*F + D*H
  circulantVcMatrixMultiply(C, F, n/2, temp1);
  circulantVcMatrixMultiply(D, H, n/2, temp2);
  for(i=0; i<n/2; ++i)
    result[i+3*n/2] = temp1[i] + temp2[i];

  fftw_free(temp1);
  fftw_free(temp2);
}

// Give me the first column of a circulant matrix M.
void circulantVcMatrixMultiply(double* c, double* VecCpy, int n, double* result) {
  int i;
  fftw_complex* fftc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  fftw_complex* fftVec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

  fftw_execute_dft_r2c(daPlans[n], c, fftc);
  fftw_execute_dft_r2c(daPlans[n], VecCpy, fftVec);

  fftw_complex* multiply = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  for(i=0; i<n; ++i) {
    multiply[i] = fftc[i] * fftVec[i];
  }

  fftw_plan Finalplan  = fftw_plan_dft_c2r_1d(n, multiply, result, FFTW_ESTIMATE);
  fftw_execute(Finalplan);

  for(i=0; i<n; ++i) {
    result[i] /= n;
  }

  fftw_destroy_plan(Finalplan);
  fftw_free(fftVec);
  fftw_free(multiply);
  fftw_free(fftc);

  return;
}

// Multiply Z by the Rn that was precomputed at n, l
void preFourBcirculantVcMatrixMultiply(int n, int l, double* Vec, double* result) {
  int m = 4*n;
  double* temp = (double*) fftw_malloc(sizeof(double) * m/2);
  double* temp2 = (double*) fftw_malloc(sizeof(double) * m/2);
  int x;

  // Top left: a column of the top left times first half of Vec
  circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*m/4+l)], Vec, m/2, result);

  // Top right: a column of the top right times second half of Vec
  circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*m/4+l)+m/2], Vec+m/2, m/2, temp);

  // Add top left and top right
  for(x=0; x<m/2; ++x)
    result[x] += temp[x];

  // Bottom left: a column of the bottom left times first half of Vec
  circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*m/4+l)+m], Vec, m/2, result+m/2);

  // Bottom right: a column of the bottom right times second half of Vec
  circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*m/4+l)+3*m/2], Vec+m/2, m/2, temp2);

  // Add bottom left and bottom right
  for(x=m/2; x<m; ++x)
    result[x] += temp2[x-m/2];

  fftw_free(temp);
  fftw_free(temp2);
}

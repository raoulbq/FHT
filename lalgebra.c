#include <stdio.h>
#include <stdlib.h>
#include <complex.h> //I is now defined as sqrt(-1)
#include <fftw3.h>
#include "testingstuff.h"
#include "lalgebra.h"
#include "FHT.h"

//****************************************************************************
//==========================LINEAR ALGEBRA TOOLS==============================
//****************************************************************************
FILE *tester;
fftw_plan *daPlans;
double  *daRns;
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
    output[n*i+whichColumn] = input[i];
  }
}

void initFastFouriers(int n) {
  daPlans = (fftw_plan*) fftw_malloc(sizeof(fftw_plan) * 2*n+1);

  int i;
  for(i=8; i<=2*n; i*=2) {
    double* in = (double*) fftw_malloc(sizeof(double) * i);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * i);
    daPlans[i] = fftw_plan_dft_r2c_1d(i,in,out,FFTW_ESTIMATE);
    fftw_free(in);
    fftw_free(out);
  }
}

void initRns(int n) {
  int i,j;
  daRns	= (double*) fftw_malloc(sizeof(double) * (n+1) * n * 8*n);
  daRnsSize = n;
  //precompute the necessary Rns
  for(i=n; i>=4;i/=2) {
    for(j=0; j<n; j+=i) {
      precomputeRnAndStore(i, j, &daRns[8*n*(n*i+j)]);
    }
  }
}

void destroyFastFouriers(int n) {
  int i;
  for(i=8; i<=2*n; i*=2) {
    fftw_destroy_plan(daPlans[i]);
  }
  fftw_free(daPlans);
}

// A B  *  E F  = AE + BG  AF + BH
// C D     G H    CE + DG  CF + DH
void fourBcirculantSqMatrixMultiply(double* M1, double* M2, int n, double* result) {
  printf("doing a multiplication of size %i\n",n);
  double *A,*B,*C,*D,*E,*F,*G,*H;
  double* temp1	  = (double *) fftw_malloc(sizeof(double) * n/2);
  double* temp2	  = (double *) fftw_malloc(sizeof(double) * n/2);
  int i;

  //fill up the columns
  A = M1;
  E = M2;
  B = M1+n/2;
  F = M2+n/2;
  C = M1+n;
  G = M2+n;
  D = M1+3*n/2;
  H = M2+3*n/2;

  //A*E+B*G top left
  circulantVcMatrixMultiply(A,E,n/2,temp1);
  circulantVcMatrixMultiply(B,G,n/2,temp2);
  //Add em up
  for(i=0;i<n/2;++i)
    result[i]=temp1[i]+temp2[i];

  //A*F+B*H top right
  circulantVcMatrixMultiply(A,F,n/2,temp1);
  circulantVcMatrixMultiply(B,H,n/2,temp2);
  //Add em up
  for(i=0;i<n/2;++i)
    result[i+n/2]=temp1[i]+temp2[i];

  //C*E+D*G bottom left
  circulantVcMatrixMultiply(C,E,n/2,temp1);
  circulantVcMatrixMultiply(D,G,n/2,temp2);
  //Add em up
  for(i=0;i<n/2;++i)
    result[i+n]=temp1[i]+temp2[i];

  //C*F+D*H bottom right
  circulantVcMatrixMultiply(C,F,n/2,temp1);
  circulantVcMatrixMultiply(D,H,n/2,temp2);
  //Add em up
  for(i=0;i<n/2;++i)
    result[i+3*n/2]=temp1[i]+temp2[i];

  fftw_free(temp1);
  fftw_free(temp2);
}

//give me the first column of a circulant matrix in M.
void circulantVcMatrixMultiply(double* c, double* VecCpy, int n, double* result) {
  int i;
  fftw_complex* fftc = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n);
  fftw_complex* fftVec = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n);

  //fast fourier c
  //     fftw_plan cplan = fftw_plan_dft_r2c_1d(n,c,fftc,FFTW_ESTIMATE);
  fftw_execute_dft_r2c(daPlans[n],c,fftc);

  //fast fourier Vec
  //     fftw_plan Vecplan  = fftw_plan_dft_r2c_1d(n,VecCpy,fftVec,FFTW_ESTIMATE);
  //     fftw_execute(Vecplan);
  fftw_execute_dft_r2c(daPlans[n],VecCpy,fftVec);

  fftw_complex* multiply = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n);
  for(i=0; i<n; ++i) {
    multiply[i] = fftc[i]*fftVec[i];
  }

  fftw_plan Finalplan  = fftw_plan_dft_c2r_1d(n,multiply,result,FFTW_ESTIMATE);
  fftw_execute(Finalplan);

  for(i=0; i<n; ++i) {
    result[i]/=n;
  }

  fftw_destroy_plan(Finalplan);
  //     fftw_destroy_plan(cplan);
  //     fftw_destroy_plan(Vecplan);
  fftw_free(fftVec);
  fftw_free(multiply);
  fftw_free(fftc);

  return;
}

//multiply Z by the Rn that was precomputed at n,l
void preFourBcirculantVcMatrixMultiply(int n, int l, double* Vec, double* result) {
  n*=4;
  double* temp = (double *) fftw_malloc(sizeof(double) * n/2);
  double* temp2 = (double *) fftw_malloc(sizeof(double) * n/2);

  int x;
  //top left (want a column of the top left times first half of Vec)
  circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*n/4+l)],Vec,n/2,result);

  //top right (want a column of the top right times second half of Vec)
  circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*n/4+l)+n/2],Vec+n/2,n/2,temp);

  //add top left and top right
  for(x=0; x<n/2; ++x)
    result[x] += temp[x];

  //bottom left (want a column of the bottom left times first half of Vec)
  circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*n/4+l)+n],Vec,n/2,result+n/2);

  //bottom right (want a column of the bottom right times second half of Vec)
  circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*n/4+l)+3*n/2],Vec+n/2,n/2,temp2);

  //add bottom left and bottom right
  for(x=n/2; x<n; ++x)
    result[x] += temp2[x-n/2];

  fftw_free(temp);
  fftw_free(temp2);
  n/=4;
}

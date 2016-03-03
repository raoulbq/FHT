#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>
#include <string.h>

#include "lalgebra.h"
#include "fht.h"
#include "parameters.h"


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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>
#include <string.h>

#include "lalgebra.h"
#include "fht.h"
#include "io.h"
#include "parameters.h"


int main(int argc, char *argv[]) {
  clock_t fclock, nclock;
  int n = LITN;
  int N = BIGN;
  int i;
  double *data = (double*) fftw_malloc(sizeof(double) * (2*n+1));
  double *naiveResult = (double*) fftw_malloc(sizeof(double) * N);
  double *fancyResult = (double*) fftw_malloc(sizeof(double) * N);


  // Load input data
  oneDFileReader("data/input.txt", (2*n+1), data);
  for(i=0; i<=2*n; ++i) {
    data[i] *= exp(0 - pow(xk(i), 2) / 2);
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
  resultsn = fopen("data/result_n.txt", "wt+");
  resultsf = fopen("data/result_f.txt", "wt+");

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

#include <stdlib.h>
#include <stdio.h>

#include "io.h"

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

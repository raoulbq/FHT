#ifndef FHT_H
#define FHT_H

double xk(int k);
void oneDFileReader(char *filename, int n, double *result);
void twoDFileReader(char *filename, int n, double *result);
void calculateFirstZ(double *Z0, double *Z1, int n);
void createAn(int n, int l, double *result);
void precomputeRnAndStore(int n, int l, double* result);
void naiveChebyshev(double *data, double *results);
void performTransform(double* Z0, double* Z1, int n, int l, double* result);
void naiveTransform(double *data, double *results);
void oneDTransform(double *data, double* result);
void twoDTransform(double *data, int n, double* result);

#endif

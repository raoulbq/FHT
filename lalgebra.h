#ifndef LALGEBRA_H
#define LALGEBRA_H

void getColumnFromSq(double* input, double* output, int whichColumn, int n);
void setColumnToSq(double* input, double* output, int whichColumn, int n);
void initFastFouriers(int n);
void destroyFastFouriers(int n);
void initRns(int n);
void fourBcirculantSqMatrixMultiply(double* M1, double* M2, int n, double* result);
void circulantVcMatrixMultiply(double* c, double* VecCpy, int n, double* result);
void preFourBcirculantVcMatrixMultiply(int n, int l, double* Vec, double* result);

#endif

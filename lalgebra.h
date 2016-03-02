//give me the first column of a circulant matrix in M.
extern void circulantVcMatrixMultiply(double* c, double* VecCpy, int n, double* result);

//multiply Z by the Rn that was precomputed at n,l
extern void preFourBcirculantVcMatrixMultiply(int n, int l, double* Vec, double* result);

extern void fourBcirculantSqMatrixMultiply(double* M1, double* M2, int n, double* result);

extern void getColumnFromSq(double* input, double* output, int whichColumn, int n);

extern void setColumnToSq(double* input, double* output, int whichColumn, int n);

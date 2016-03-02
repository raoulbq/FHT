#include <stdio.h>
#include <stdlib.h>
#include <testingstuff.h>

//****************************************************************************
//==========================TESTING TOOLS=====================================
//****************************************************************************
FILE *tester; //for test information
//prints a len by len matrix to tester
void printMatrix(char* name, double* mat,int len) {
	int i,j;
	fprintf(tester,name);
	fprintf(tester,"=[\n");
	for(i=0; i<len; ++i) {
		fprintf(tester,"\n");
		for(j=0; j<len; ++j)
		    fprintf(tester,"%.16Lf	",mat[len*i+j]);
	}
	fprintf(tester,"];\n\n");
}

void printNonSqMatrix(double* mat,int n1, int n2) {
   	FILE* flTJ = fopen("Tj.txt","wt+");
	int i,j;
	for(i=0; i<n1; ++i) {
	    if(i!=0)fprintf(flTJ,"\n");
		for(j=0; j<n2; ++j)
		    fprintf(flTJ,"%.16Lf ",mat[n2*i+j]);
	}
	fclose(flTJ);
}

//prints a len vector to tester
void printVector(char* name, double* Vec,int len) {
	int i,j;
	fprintf(tester,name);
	fprintf(tester,"=[\n");
	for(i=0; i<len; ++i) {
	    fprintf(tester,"%.16Lf	",Vec[i]);
	}
	fprintf(tester,"];\n\n");
}

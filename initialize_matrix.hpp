#ifndef HEADER1  
#define HEADER1
typedef double (*Function)(double);

void InitV (double* A, double* B, double* V, double* H, double* f, double tau, double h, double mui, double lambda, int N, int M, int n, Function p);

#endif
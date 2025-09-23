#ifndef INITIALIZE_HEADER 
#define INITIALIZE_HEADER

typedef double (*Function2d)(double, double);
typedef double (*Function)(double);

void InitV (double* A, double* B, double* V, double* H, Function2d f, double tau, double h, double mui, double lambda, int M, int n, Function p);
void InitH (double* A, double* B, double* V, double* upV, double* H, Function2d f0, double tau, double h, int M, int n);

#endif
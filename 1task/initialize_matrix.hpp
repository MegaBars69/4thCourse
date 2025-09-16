#ifndef INITIALIZE_HEADER 
#define INITIALIZE_HEADER

typedef double (*Function)(double);

void InitV (double* A, double* B, double* V, double* H, double* f, double tau, double h, double mui, double lambda, int M, Function p);
void InitH (double* A, double* B, double* V, double* upV, double* H, double* f, double tau, double h, double mui, double lambda, int M);

#endif
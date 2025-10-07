#ifndef INITIALIZE_HEADER 
#define INITIALIZE_HEADER

typedef double (*Function2d)(double, double);
typedef double (*Function)(double);

void InitV (double* A, double* B, double* V, double* H, double tau, double h, double mui, double lambda, int M, double X_a, double X_b, double t);
void InitH (double* A, double* B, double* V, double* upV, double* H,  double tau, double h, int M, double X_a, double X_b, double t);

#endif
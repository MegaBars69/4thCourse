#ifndef HEADER_ALGO 
#define HEADER_ALGO
typedef double (*Function)(double);

double FindLambda (double* H, double mui, int M);
void SolveScheme (double* f, double *f0, double mui, double T, double X, int N, int M, Function p, double* V, double* H);

#endif
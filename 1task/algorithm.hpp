#ifndef HEADER_ALGO 
#define HEADER_ALGO
typedef double (*Function)(double);
typedef double (*Function2d)(double, double);

double FindLambda (double* H, double mui, int M);
void SolveScheme (Function2d f, Function2d f0, double mui, double T, double X, int N, int M, Function p, double* V, double* H);

#endif
#ifndef HEADER_ALGO 
#define HEADER_ALGO
typedef double (*Function)(double);

double FindLambda (double* H, double mui, int M);
void SolveScheme (double* f, double mui, double T, double X, int N, int M, Function p);

#endif
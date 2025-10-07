#ifndef HEADER_ALGO 
#define HEADER_ALGO
typedef double (*Function)(double);
typedef double (*Function2d)(double, double);

double FindLambda (double* H, double mui, int M);
void SolveScheme (double mui, double T_a, double T_b, double X_a, double X_b, int N, int M, Function p, double* V, double* H);

#endif
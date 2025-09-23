#ifndef TEST_HEADER
#define TEST_HEADER

extern double Cp;
extern double mui;

void InitF (double* f, double* f0, double tau, double h, int N, int M);
void Test_Scheme (int N, int M);
double U (double t, double x);
double po (double t, double x);
double func_f0 (double t, double x);
double func_f (double t, double x);
#endif 
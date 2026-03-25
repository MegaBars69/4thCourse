#ifndef TEST_HEADER
#define TEST_HEADER

#include <cmath>

extern int Cp;
extern double mui;
extern bool liniar;  

void InitF (double* f, double* f0, double tau, double h, int N, int M);
void Test_Scheme (int N, int M);
double U (double x);
double po (double x);
double func_f0 (double t, double x);
double func_f (double t, double x);

double p (double x);

#endif
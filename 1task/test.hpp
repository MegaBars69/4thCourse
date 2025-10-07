#ifndef TEST_HEADER
#define TEST_HEADER

#include <cmath>

extern double Cp;
extern double mui;

void InitF (double* f, double* f0, double tau, double h, int N, int M);
void Test_Scheme (int N, int M);
double U (double t, double x);
double po (double t, double x);
double func_f0 (double t, double x);
double func_f (double t, double x);

inline double H_0 (double x)
{
    return std::cos (3. * M_PI * x) + 1.5;
}

inline double V_0 (double x)
{
    return std::sin (4. * M_PI * x);
}

#endif 
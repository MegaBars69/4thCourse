#include "test.hpp"
#include <iostream>
#include <cmath>

double Cp = 1;
double mui = 0.1;

double U (double t, double x) {return cos (2 * M_PI * t) * sin (4 * M_PI * x);}

double po (double t, double x) {return exp (t) * (cos (3 * M_PI * x) + 1.5);}

double func_f (double t, double x)
{
    double u = U(t,x);
    double p = po (t, x);
    double ut = -2 * M_PI * sin (2 * M_PI * t) * sin (4 * M_PI * x);
    double ux = 4 * M_PI * cos (4 * M_PI * x) * cos (2 * M_PI * t);
    double uxx = -(16 * M_PI * M_PI) * sin (4 * M_PI * x) * cos (2 * M_PI * t) ;
    double Pxp = -3 * M_PI * sin (3 * M_PI * x) / (1.5 + cos (3 * M_PI * x));
    return ut + u * ux + Cp * Pxp - mui * uxx / p;
}

double func_f0 (double t, double x)
{
    double u = U (t, x);
    double ux = 4 * M_PI * cos (4 * M_PI * x) * cos (2 * M_PI * t);
    double px = -3 * M_PI * exp (t) * sin (3 * M_PI * x);
    double p = po (t,x);
    return p * (1 + ux) + u * px;
}

void InitF (double* f, double* f0, double tau, double h, int N, int M)
{
    for (int n = 0; n <= N; n++)
    {
        for (int m = 0;  m <= M; m++)
        {
            f[n * (M + 1) + m] = func_f (n*tau, m*h);
            f0[n * (M + 1) + m] = func_f0 (n*tau, m*h);
        }
    }
}

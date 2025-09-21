#include "progonka_solver.hpp"
#include <iostream>
#include <cmath>
#define EPSILON 1e-15

bool SolveSystem (double* A, double* B, double* x, int n)
{
    double* a = A;
    double* b = A + n;
    double* c = A + 2*n;
    double* d = B;

    double denomenator = b[0];
    double* alpha = c;
    double* betha = b;
    if (fabs (denomenator) < EPSILON)
        return false;
    alpha[0] = -c[0] / denomenator;
    betha[0] = d[0] / denomenator;

    for (int i = 1; i < n; i++)
    {
        denomenator = b[i] + a[i]*alpha[i-1];
        if (fabs (denomenator) < EPSILON)
            return false;
        alpha[i] = -c[i]/denomenator;
        betha[i] = (d[i]-a[i]*betha[i-1])/denomenator;
    }
    x[n-1] = betha[n-1];
    for (int i = n-2; i >=0 ; i--)
    {
        x[i] = alpha[i]*x[i+1] + betha[i];
    }
    return true;
}

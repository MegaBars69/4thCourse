#include <iostream>
#include "test.hpp"
#include "algorithm.hpp"
#include "log.hpp"
#include <string.h>
#include <cmath>

double p (double x) 
{
    double Cp = 10;
    return Cp * x;
}


int main(int argc, char* argv[])
{
    double T, X, mui = 0.01;
    int N, M;
    if (!init_args (argc, argv, T, X, N, M)) return 1;

    double h = X / M;
    double *f = new double[M + 1];

    for (int m = 0; m <= M; m++)
    {
        f[m] = sin(m*h);
    }
    SolveScheme (f, mui, T, X, N, M, p);

    delete [] f;

    return 0;
}
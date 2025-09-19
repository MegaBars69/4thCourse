#include "norms.hpp"
#include <cmath>

double C_norm (double* x, double* y, int M)
{
    double norm = 0;
    double coord;
    for (int i = 0; i <= M; i++)
    {
        coord = fabs (x[i] - y[i]);
        if (coord > norm) norm = coord;
    }
    return norm;
}

double L_norm (double* x, double* y, int M, double h)
{
    double norm = 0;
    for (int i = 1; i < M; i++)
    {
        norm += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt (h * norm);
}

double W_norm (double* x, double* y, int M, double h)
{
    double norm = 0;
    for (int i = 1; i < M; i++)
    {
        norm += (x[i] - y[i]) * (x[i] - y[i]);
    }
    norm+=0.5 * ((x[0] - y[0]) * (x[0] - y[0]) + (x[M] - y[M]) * (x[M] - y[M]));
    return sqrt (h * norm);
}
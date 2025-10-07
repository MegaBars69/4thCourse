#include "initialize_matrix.hpp"
#include "test.hpp"
#include <iostream>

void InitV (double* A, double* B, double* V, double* H,  double tau, double h, double mui, double lambda, int M, double X_a, double /*X_b*/, double t, Function p)
{
    double *a = A, *b = A + M - 1, *c = A + 2 * (M - 1), *d = B;
    double th = tau/h;

    for (int m = 1; m <= M - 1; m++)
    {
        a[m - 1] = (-th * (V[m] + V[m - 1])/6 - th * lambda / h);
        c[m - 1] = ( th * (V[m + 1] + V[m])/6 - th * lambda / h);
        b[m - 1] = (1 + 2 * th * lambda / h);
        d[m - 1] = V[m] - 0.5 * th * (p (H[m + 1]) - p (H[m - 1])) / H[m] - th * (lambda - mui / H[m]) * (V[m + 1] - 2 * V[m] + V[m - 1]) / h + tau * func_f(t , X_a + m * h);
    }    
}

void InitH (double* A, double* B, double* V, double* upV, double* H, double tau, double h, int M, double X_a, double X_b,  double t)
{
    double *a = A, *b = A + M + 1, *c = A + 2*(M + 1), *d = B;
    double th = tau/h;

    b[0] = 1 - th * upV[0]/ 2;
    c[0] = th * upV[1]/ 2;
    d[0] = H[0] - 0.5 * th * (H[0] * upV[1] - H[0] * upV[0] - 2 * H[2]*V[2] + 2.5 * H[1] * V[1] - 2 * H[0] * V[0] + 0.5 * H [3] * V [3] - 2 * H[0] * V[2] + 2.5 * H[0] * V[1] + 0.5 * H[0] * V[3]) + tau * func_f0(t, X_a);

    for (int m = 1; m <= M - 1; m++)
    {
        a[m] = -th * (upV[m] + upV[m-1]) / 4;
        b[m] = 1.0;
        c[m] = th * (upV[m + 1] + upV[m]) / 4;
        d[m] = H[m] - H[m] * th * 0.25 * (upV[m + 1] - upV[m - 1]) + tau * func_f0(t, X_a + m * h);
    }

    a[M] = -th * upV[M - 1] / 2; 
    b[M] = 1 + th * upV[M] / 2;
    d[M] = H[M] - 0.5 * th * (H[M] * upV[M] - H[M] * upV[M - 1] + 2 * H[M] * V[M] - 2.5 * H[M] * V[M - 1] - 2.5 * H[M - 1] * V[M - 1] + 2 * H[M] * V[M - 2] + 2 * H[M - 2] * V[M - 2] - 0.5 * H[M - 3] * V[M - 3] - 0.5 * H[M] * V[M - 3]) + tau * func_f0(t, X_b); 
}
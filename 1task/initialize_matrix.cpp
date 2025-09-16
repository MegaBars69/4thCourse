#include "initialize_matrix.hpp"
#include <iostream>

void InitV (double* A, double* B, double* V, double* H, double* f, double tau, double h, double mui, double lambda, int M, Function p)
{
    double *a = A, *b = A + M - 1, *c = A + 2*(M - 1), *d = B;
    double th = tau/h;

    for (int m = 1; m <= M - 1; m++)
    {
        a[m - 1] = (m != 1 ? (-th * (V[m] + V[m-1])/6 - th * lambda / h) * H[m] : 0);
        c[m - 1] = (m != M - 1 ? (th * (V[m+1] + V[m])/6 - th * lambda) * H[m] : 0);
        b[m - 1] = (1 + 2 * th * lambda / h) * H[m];
        d[m - 1] = V[m] * H[m] - (p (H[m + 1]) - p (H[m - 1])) * th / 2 - th * (lambda * H[m] - mui) * (V[m + 1] - 2 * V[m] + V[m - 1]) / h + H[m] * tau * f[m];
    }    
}

void InitH (double* A, double* B, double* V, double* upV, double* H, double* /*f*/, double tau, double h, double /*mui*/, double /*lambda*/, int M)
{
    double *a = A, *b = A + M + 1, *c = A + 2*(M + 1), *d = B;
    double th = tau/h;

    b[0] = 1 - th * upV[0]/ 2;
    c[0] = th * upV[1]/ 2;
    d[0] = H[0] - th * (H[0] * upV[1] - H[0] * upV[0] - 2 * H[2]*V[2] + 2.5 * H[1] * V[1] - 2 * H[0] * V[0] + 0.5 * H [3] * V [3] - 2 * H[0] * V[2] + 2.5 * H[0] * V[1] + 0.5 * H[0] * V[3]) / 2;

    for (int m = 1; m <= M - 1; m++)
    {
        a[m] = -th * (upV[m] + upV[m-1]) / 4;
        b[m] = 1;
        c[m] = th * (upV[m + 1] + upV[m]) / 4;
        d[m] = H[m] * (1 + th * (upV[m - 1] - upV[m + 1])/ 2);
    }

    a[M] = -th * upV[M - 1]/ 2;
    b[M] = 1 + th * (upV[M])/ 2;
    d[M] = H[M] - th * (H[M] * upV[M] - H[M] * upV[M - 1] + 2 * H[M] * V[M] - 2.5 * H[M] * V[M - 1] - 2.5 * H[M - 1] * V[M - 1] + 2 * H[M] * V[M - 2] + 2 * H[M - 2] * V[M - 2] - 0.5 * H[M - 3] * V[M - 3] - 0.5 * H[M] * V[M - 3] ) / 2; 
}
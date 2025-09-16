#include "initialize_matrix.hpp"
#include <iostream>

void InitV (double* A, double* B, double* V, double* H, double* f, double tau, double h, double mui, double lambda, int N, int M, int n, Function p)
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

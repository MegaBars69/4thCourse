#include <iostream>
#include "test.hpp"
#include "algorithm.hpp"
#include "progonka_solver.hpp"
#include "log.hpp"
#include "norms.hpp"
#include <string.h>
#include <cmath>
#include <fenv.h>

double p (double x) 
{
    return Cp * x;
}


int main(int argc, char* argv[])
{
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    double T, X;
    int N, M;
    if (!init_args (argc, argv, T, X, N, M)) return 1;

    double h = X / M;
    double tau = T / N;
    double *f = new double[(N + 1) * (M + 1)];
    double *f0 = new double[(N + 1) * (M + 1)];
    double* V = new double [M + 1];
    double* H = new double [M + 1];
    double* u = new double [M + 1]; 
    double* NullMatrix = new double [M + 1]; 
    for (int i = 0; i <= M; i++)
    {
        u[i] = U(T, i*h);
        NullMatrix[i] = 0;
    }
    InitF (f, f0, tau, h, N, M);
 
    auto start = clock();
    SolveScheme (f, f0,  mui, T, X, N, M, p, V, H);
    auto end = clock ();

    auto t = static_cast<double> (end - start) / CLOCKS_PER_SEC;

    double w_norm = W_norm(V, u, M, h);
    double l_norm = L_norm(V, u, M, h);
    double c_norm = C_norm(V, u, M);
    double V_w_norm = W_norm(V, NullMatrix, M, h);
    double V_l_norm = L_norm(V, NullMatrix, M, h);
    double V_c_norm = C_norm(V, NullMatrix, M);

    printf ("Res1 = %e Res2 = %e Res3 = %e T = %.2f \n", c_norm / V_c_norm, l_norm / V_l_norm, w_norm / V_w_norm, t);

    delete[] f;
    delete[] f0;
    delete[] V;
    delete[] H;
    delete[] u;
    delete[] NullMatrix;

    return 0;
}
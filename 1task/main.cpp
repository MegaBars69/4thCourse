#include <iostream>
#include "test.hpp"
#include "algorithm.hpp"
#include "progonka_solver.hpp"
#include "log.hpp"
#include <string.h>
#include <cmath>
#include <fenv.h>

double p (double x) 
{
    double Cp = 10;
    return Cp * x;
}


int main(int argc, char* argv[])
{
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    double T, X, mui = 0.01;
    int N, M;
    if (!init_args (argc, argv, T, X, N, M)) return 1;

    double h = X / M;
    double tau = T / N;
    double *f = new double[(N + 1) * (M + 1)];
    double *f0 = new double[(N + 1) * (M + 1)];
    double* V = new double [M + 1];
    double* H = new double [M + 1];

    InitF (f, f0, tau, h, N, M);
    SolveScheme (f, f0,  mui, T, X, N, M, p, V, H);

    double norm = 0;
    for (int i = 0; i <= M; i++)
    {
        norm = (fabs (V[i] - U (T - tau, i * h)) > norm ? fabs (V[i] - U (T - tau, i * h)) : norm);
    }
    std::cout<<norm<<std::endl;

    delete[] f;
    delete[] f0;
    delete[] V;
    delete[] H;

    return 0;
}
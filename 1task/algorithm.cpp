#include "algorithm.hpp"
#include "initialize_matrix.hpp"
#include "progonka_solver.hpp"
#include "log.hpp"
#include "test.hpp"
#include <iostream>
#include <cmath>
#include <string.h>

double FindLambda (double* H, double mui, int M)
{
    double min = H[0];
    for (int i = 1; i <= M; i++)
    {
        min = H[i] > min ? min : H[i];
    }
    return mui/min;
}

void SolveScheme (double* f, double *f0, double mui, double T, double X, int N, int M, Function p, double* V, double* H)
{
    double tau = T / N, h = X / M, lambda = 0;

    double* upV = new double [M + 1];
    double* A = new double [3 * (M + 1)];
    double* b = new double [M + 1];

    for (int i = 0; i <= M; ++i) 
    {
        V[i] = U (0 , i*h);
        H[i] = po (0, i*h);
        upV[i] = 0;
    }
    V[0] = V[M] = 0; upV[0] = upV[M] = 0;

    for (int n = 0; n < N; n++)
    {
        lambda = FindLambda (H, mui, M);
        InitV (A, b, V, H, f, tau, h, mui, lambda, M, n, p);
        
        upV[0] = 0; upV[M] = 0;
        if (!SolveSystem (A, b, upV + 1, M - 1))
        {
            printf ("Something went wrong!\n");
            delete[] upV; 
            delete[] A; 
            delete[] b; 

            return;
        }

        InitH (A, b, V, upV, H, f0, tau, h, mui, lambda, M, n);
        if (!SolveSystem (A, b, H, M + 1))
        {
            printf ("Something went wrong!\n");
            delete[] upV; 
            delete[] A; 
            delete[] b; 
            return;
        }

        for (int i = 0; i <= M; i++)
            V[i] = upV[i];
    }
    printf ("PRINTING V:\n");
    PrintVector (V, M + 1);
    printf ("PRINTING H:\n");
    PrintVector (H, M + 1);

    delete[] upV; 
    delete[] A; 
    delete[] b; 
}

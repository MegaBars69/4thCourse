#include "algorithm.hpp"
#include "initialize_matrix.hpp"
#include "progonka_solver.hpp"
#include "log.hpp"
#include "test.hpp"
#include "norms.hpp"
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

void SolveScheme (Function2d f, Function2d f0, double mui, double T, double X, int N, int M, Function p, double* V, double* H)
{
    double tau = T / N, h = X / M, lambda = 0;

    double* upV = new double [M + 1];
    double* A = new double [3 * (M + 1)];
    double* b = new double [M + 1];
    double* u = new double [M + 1]; 
    double* rho = new double [M + 1];
    double* NullVector = new double [M + 1]; 
    for (int i = 0; i <= M; i++) 
    {
        V[i] = U (0, i*h);
        H[i] = po (0, i*h);
        upV[i] = 0;
    }
    
    auto something_went_wrong = [&]() 
    {
        printf("Something went wrong!\n");
        delete[] upV; 
        delete[] A; 
        delete[] b; 
    };


    V[0] = V[M] = 0; upV[0] = upV[M] = 0;

    for (int n = 0; n < N; n++)
    {
        lambda = FindLambda (H, mui, M);

        InitV (A, b, V, H, f, tau, h, mui, lambda, M, n, p);
        if (!SolveSystem (A, b, upV + 1, M - 1)) something_went_wrong ();
    
        InitH (A, b, V, upV, H, f0, tau, h, M, n);
        if (!SolveSystem (A, b, H, M + 1)) something_went_wrong ();

        memcpy (V, upV, (M + 1) * sizeof (double));
        for (int i = 0; i <= M; i++)
        {
            u[i] = U((n+1)*tau, i*h);
            rho[i] = po((n+1)*tau, i*h);
            NullVector[i] = 0;
        }
        double w_norm = W_norm(V, u, M, h);
        double l_norm = L_norm(V, u, M, h);
        double c_norm = C_norm(V, u, M);
        double V_w_norm = W_norm(V, NullVector, M, h);
        double V_l_norm = L_norm(V, NullVector, M, h);
        double V_c_norm = C_norm(V, NullVector, M);
        double H_w_norm = W_norm(H, NullVector, M, h);
        double H_l_norm = L_norm(H, NullVector, M, h);
        double H_c_norm = C_norm(H, NullVector, M);
        printf ("ResV1 = %e ResV2 = %e ResV3 = %e \n", c_norm / V_c_norm, l_norm / V_l_norm, w_norm / V_w_norm);

        w_norm = W_norm(H, rho, M, h);
        l_norm = L_norm(H, rho, M, h);
        c_norm = C_norm(H, rho, M);

        printf ("ResH1 = %e ResH2 = %e ResH3 = %e \n", c_norm / H_c_norm, l_norm / H_l_norm, w_norm / H_w_norm);
        
    }

    delete[] upV; 
    delete[] A; 
    delete[] b; 
}

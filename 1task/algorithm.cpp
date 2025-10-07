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

void SolveScheme (double mui, double T_a, double T_b, double X_a, double X_b, int N, int M, Function p, double* V, double* H)
{
    double tau = (T_b - T_a) / N, h = (X_b - X_a) / M, lambda = 0;
    double t = T_a;

    double* upV = new double [M + 1];
    double* A = new double [3 * (M + 1)];
    double* b = new double [M + 1];    
   
    for (int i = 0; i <= M; i++) 
    {
        V[i] = V_0 (X_a + i*h);
        H[i] = H_0 (X_a + i*h);
        upV[i] = 0;
    }
    
    auto something_went_wrong = [&](int n) 
    {
        std::cout<<"Something went wrong!\n"<<"Step:"<<n<<std::endl;
        delete[] upV; 
        delete[] A; 
        delete[] b; 
    };

    for (int n = 0; n < N; n++)
    {
        t = (n != N ? T_a + n * tau : T_b);
        lambda = FindLambda (H, mui, M);

        InitV (A, b, V, H, tau, h, mui, lambda, M, X_a, X_b, t, p);
        if (!SolveSystem (A, b, upV + 1, M - 1)) something_went_wrong (n);
        
        InitH (A, b, V, upV, H, tau, h, M, X_a, X_b, t);
        if (!SolveSystem (A, b, H, M + 1)) something_went_wrong (n);
              
        memcpy (V, upV, (M + 1) * sizeof (double));
    }

    delete[] upV; 
    delete[] A; 
    delete[] b; 
}

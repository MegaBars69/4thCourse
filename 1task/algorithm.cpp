#include "algorithm.hpp"
#include "initialize_matrix.hpp"
#include "progonka_solver.hpp"
#include "log.hpp"
#include "test.hpp"
#include "norms.hpp"
#include <iostream>
#include <cmath>
#include <string.h>

#define EPS 10e-4

double NormOfStable (double *V, double *H, const double &H_solution, int M)
{
    double norm = sqrt (pow (V[1], 2) + pow (H[1] - H_solution, 2));
    std::cout<<"H_solution: "<<H_solution<<" Norm: ";
    for (int i = 1; i < M; i++)
    {
        double cur = sqrt (pow (V[i], 2) + pow (H[i] - H_solution, 2));
        norm = (norm > cur ? norm : cur);
    }
    std::cout<<norm<<"\n";
    return norm;
}

bool SolutionIsStable (double *V, double *H, int M, int num)
{
    double H_solution = 0;
    for (int i = 1; i < M; i++){ H_solution += H[i];}
    H_solution /= num;
    return NormOfStable (V, H, H_solution, M) < EPS;   
}

double FindLambda (double* H, double mui, int M)
{
    double min = H[0];
    for (int i = 1; i <= M; i++)
    {
        min = H[i] > min ? min : H[i];
    }
    return mui/min;
}

void SolveScheme (double mui, double T_a, double T_b, double X_a, double X_b, int N, int M, double* V, double* H)
{
    double tau = (T_b - T_a) / N, h = (X_b - X_a) / M, lambda = 0;
    double t = T_a;
    int num = M;

    double* upV = new double [M + 1];
    double* A = new double [3 * (M + 1)];
    double* b = new double [M + 1];    
   
    for (int i = 0; i <= M; i++) 
    {
        V[i] = U (X_a + i*h);
        H[i] = po (X_a + i*h);
        upV[i] = 0;
    }
    
    auto something_went_wrong = [&](int n) 
    {
        std::cout<<"Something went wrong!\n"<<"Step:"<<n<<std::endl;
        delete[] upV; 
        delete[] A; 
        delete[] b; 
    };
    int n = 1;
    while (!SolutionIsStable (V, H, M, num))
    {    
        //num = n * M;
        lambda = FindLambda (H, mui, M);

        InitV (A, b, V, H, tau, h, mui, lambda, M, X_a, X_b, t);
        if (!SolveSystem (A, b, upV + 1, M - 1)) something_went_wrong (n);
        
        InitH (A, b, V, upV, H, tau, h, M, X_a, X_b, t);
        if (!SolveSystem (A, b, H, M + 1)) something_went_wrong (n);
        std::cout<<"V Norm: "<<C_norm (V, upV, M)<<"\n";      
        memcpy (V, upV, (M + 1) * sizeof (double));
        n++;
    }

    delete[] upV; 
    delete[] A; 
    delete[] b; 
}

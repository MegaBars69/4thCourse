#include "table_printer.hpp"
#include "test.hpp"
#include "algorithm.hpp"
#include "progonka_solver.hpp"
#include "norms.hpp"
#include "log.hpp"
#include <iostream>

void TestScheme ()
{
    double T_a = 0, T_b = 1, X_a = 0, X_b = 1;
    int MaxN = 100000, MaxM = 100000;
    double *f = new double[(MaxN + 1) * (MaxM + 1)];
    double *f0 = new double[(MaxN + 1) * (MaxM + 1)];
    double* V = new double [MaxM + 1];
    double* H = new double [MaxM + 1];
    double* u = new double [MaxM + 1]; 
    double* rho = new double [MaxM + 1];
    double* NullVector = new double [MaxM + 1]; 

    std::fill(NullVector, NullVector + MaxM + 1, 0.0);
    for (int lin = 1; lin >= 0; lin--)
    {
        liniar = (bool) lin;
        double mu = 0.1;
        for (int mu_count = 0; mu_count < 3; mu_count++, mu *= 0.1) 
        {
            mui = mu;
            for (int Cp_count = 0; Cp_count < 3; Cp_count++) 
            {
                Cp = pow (10, Cp_count);
                double res[4*4][4]; 
                if (liniar) { printf("\n\\begin{tabular}{ |l|l|l|l|l| }\n\\hline\n\\multicolumn{5}{|c|}{$\\mu = %g, p(\\rho)  =  %d \\rho$} \\\\\n\\hline\n$\\tau\\setminus h$ & $0.1$ & $0.01$ & $0.001$ & $0.0001$\\\\\n\\hline\n", mui, Cp); }
                else {printf("\n\\begin{tabular}{ |l|l|l|l|l| }\n\\hline\n\\multicolumn{5}{|c|}{$\\mu = %g, p(\\rho)  =  \\rho ^ {1.4}$} \\\\\n\\hline\n$\\tau\\setminus h$ & $0.1$ & $0.01$ & $0.001$ & $0.0001$\\\\\n\\hline\n", mui); }
       
                for (int row = 0; row < 4; row++) 
                {
                    double tau;
                    for (int col = 0; col < 4; col++) 
                    {
                        int N = pow (10, row + 1), M = pow (10, col + 1);
                        double h = (X_b - X_a) / M;
                        tau = (T_b - T_a) / N;
                        for (int i = 0; i <= M; i++)
                        {
                            u[i] = U(T_b , X_a + i*h);
                            rho[i] = po(T_b, X_a + i*h);
                        }

                        auto start = clock();
                        SolveScheme (mui, T_a, T_b, X_a, X_b, N, M, V, H);
                        auto end = clock ();

                        auto t = static_cast<double> (end - start) / CLOCKS_PER_SEC;

                        double w_norm = W_norm(V, u, M, h);
                        double l_norm = L_norm(V, u, M, h);
                        double c_norm = C_norm(V, u, M);
                        double V_w_norm = W_norm(V, NullVector, M, h);
                        double V_l_norm = L_norm(V, NullVector, M, h);
                        double V_c_norm = C_norm(V, NullVector, M);
                        res[row * 4 + col][0] = c_norm / V_c_norm;
                        res[row * 4 + col][1] = l_norm / V_l_norm;
                        res[row * 4 + col][2] = w_norm / V_w_norm;
                        res[row * 4 + col][3] = t;
                    }
                    
                    printf("$%g$ ", tau);
                    for (int cur_val = 0; cur_val < 4; cur_val++)
                    {
                        for (int cur_col = 0; cur_col < 4; cur_col++)
                        {
                            
                            printf ("& $%e$ ", res[row * 4 + cur_col][cur_val]);
                        }
                            printf("\\\\\n");
                    }
                    printf("\\hline\n");           
                }
                printf("\\end{tabular}\n\n");
            }
        }
    }

    delete[] f;
    delete[] f0;
    delete[] V;
    delete[] H;
    delete[] u;
    delete[] rho;
    delete[] NullVector;
}

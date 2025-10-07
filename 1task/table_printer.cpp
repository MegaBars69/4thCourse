#include "table_printer.hpp"
#include "test.hpp"
#include "algorithm.hpp"
#include "progonka_solver.hpp"
#include "norms.hpp"
#include "log.hpp"
#include <iostream>
#include <cmath>
#include <string.h>
#include <iomanip>
#include <ctime>
#include <sstream>

void RunConvergenceTest(double mu, int p_model, int max_power, double* results)
{
    double T_a = 0, T_b = 1, X_a = 0, X_b = 1;
    
    // Set pressure model using the extern variable
    liniar = (p_model == 0);
    
    int result_index = 0;
    
    for (int n_power = 1; n_power <= max_power; n_power++)
    {
        int N = (int)pow(10, n_power);
        
        for (int m_power = 1; m_power <= max_power; m_power++)
        {
            int M = (int)pow(10, m_power);
            
            double h = (X_b - X_a) / M;
            double* V = new double[M + 1];
            double* H = new double[M + 1];
            double* u_exact = new double[M + 1];
            double* rho_exact = new double[M + 1];
            double* NullVector = new double[M + 1];
            
            // Initialize exact solutions
            for (int i = 0; i <= M; i++)
            {
                u_exact[i] = U(T_b, X_a + i * h);
                rho_exact[i] = po(T_b, X_a + i * h);
                NullVector[i] = 0;
            }
            
            // Run the scheme
            auto start = clock();
            SolveScheme(mu, T_a, T_b, X_a, X_b, N, M, V, H);
            auto end = clock();
            double computation_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
            
            // Calculate norms for V (velocity)
            double c_norm_v = C_norm(V, u_exact, M);
            double l_norm_v = L_norm(V, u_exact, M, h);
            double w_norm_v = W_norm(V, u_exact, M, h);
            
            // Calculate norms for H (density)
            double c_norm_h = C_norm(H, rho_exact, M);
            double l_norm_h = L_norm(H, rho_exact, M, h);
            double w_norm_h = W_norm(H, rho_exact, M, h);
            
            // Store results: [C_v, L_v, W_v, C_h, L_h, W_h, time]
            results[result_index * 7 + 0] = c_norm_v;
            results[result_index * 7 + 1] = l_norm_v;
            results[result_index * 7 + 2] = w_norm_v;
            results[result_index * 7 + 3] = c_norm_h;
            results[result_index * 7 + 4] = l_norm_h;
            results[result_index * 7 + 5] = w_norm_h;
            results[result_index * 7 + 6] = computation_time;
            
            result_index++;
            
            // Cleanup
            delete[] V;
            delete[] H;
            delete[] u_exact;
            delete[] rho_exact;
            delete[] NullVector;
        }
    }
}

void PrintLatexTable(double mu, int p_model, double* results, int num_n, int num_m)
{
    std::cout << "\\begin{tabular}{ |l|l|l|l|l| }" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "\\multicolumn{5}{|c|}{$\\mu = " << std::fixed << std::setprecision(1) << mu << ", p = " << p_model << "$} \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$\\tau\\setminus h$ & $0.1$ & $0.01$ & $0.001$ & $0.0001$\\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    
    for (int n_power = 1; n_power <= num_n; n_power++)
    {
        int N = (int)pow(10, n_power);
        double tau = 1.0 / N;
        
        // Print tau value for this row
        std::cout << "$" << std::fixed << std::setprecision(n_power) << tau << "$ ";
        
        // Print C-norm for velocity (first data row)
        for (int m_power = 1; m_power <= num_m; m_power++)
        {
            int index = (n_power - 1) * num_m + (m_power - 1);
            std::cout << "& $" << std::scientific << std::setprecision(6) << results[index * 7 + 0] << "$ ";
        }
        std::cout << "\\\\" << std::endl;
        
        // Print L-norm for velocity (second data row)
        std::cout << " & ";
        for (int m_power = 1; m_power <= num_m; m_power++)
        {
            int index = (n_power - 1) * num_m + (m_power - 1);
            std::cout << "& $" << std::scientific << std::setprecision(6) << results[index * 7 + 1] << "$ ";
        }
        std::cout << "\\\\" << std::endl;
        
        // Print W-norm for velocity (third data row)
        std::cout << " & ";
        for (int m_power = 1; m_power <= num_m; m_power++)
        {
            int index = (n_power - 1) * num_m + (m_power - 1);
            std::cout << "& $" << std::scientific << std::setprecision(6) << results[index * 7 + 2] << "$ ";
        }
        std::cout << "\\\\" << std::endl;
        
        // Print computation time (fourth data row)
        std::cout << " & ";
        for (int m_power = 1; m_power <= num_m; m_power++)
        {
            int index = (n_power - 1) * num_m + (m_power - 1);
            std::cout << "& $" << std::scientific << std::setprecision(6) << results[index * 7 + 6] << "$ ";
        }
        std::cout << "\\\\" << std::endl;
        std::cout << "\\hline" << std::endl;
    }
    
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << std::endl;
}
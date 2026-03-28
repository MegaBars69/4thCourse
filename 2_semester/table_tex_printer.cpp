#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "matrix.h"

// Вспомогательная функция для печати одной таблицы для заданной переменной (mode = g, v1, v2)
void PrintTableForVariable(int var_mode,                      // g, v1 или v2
                           const std::string& var_name,       // название переменной для заголовка
                           double mu, double Cp, int mode,    // параметры газа
                           const std::vector<int>& N_points,
                           const std::vector<int>& T_steps_list,
                           double X, double Y, double T,
                           double gamma)
{
    int n_rows = T_steps_list.size();
    int n_cols = N_points.size();

    // Массив для хранения трёх норм: [row][col][norm]  norm: 0=C, 1=L2, 2=W1
    std::vector<std::vector<std::vector<double>>> norms(
        n_rows, std::vector<std::vector<double>>(n_cols, std::vector<double>(3, 0.0)));

    // Заполнение массива
    for (int row = 0; row < n_rows; ++row) {
        int T_steps = T_steps_list[row];
        double tau = T / T_steps;

        for (int col = 0; col < n_cols; ++col) {
            int Nx = N_points[col];
            int Ny = N_points[col];   // квадратная сетка
            double h = X / (Nx - 1);

            // Создание объектов
            P_gas gas(T, X, Y, Cp, gamma, mu, mode);
            Mesh mesh(X, Y, T, Nx, Ny, T_steps);
            Matrix matrix(gas, mesh);

            // Расчёт по всем временным слоям
            bool success = true;
            for (int step = 1; step <= T_steps; ++step) {
                matrix.step = step;
                if (matrix.init_and_solve_G() == -1) {
                    success = false;
                    break;
                }
                if (matrix.init_and_solve_V() == -1) {
                    success = false;
                    break;
                }
            }

            double errC, errL2, errW1;
            if (!success) {
                // Если ошибка, помечаем невязки большим числом
                errC = errL2 = errW1 = 1e10;
            } else {
                // Вычисление невязок на последнем слое
                errC = matrix.calc_res_C1(var_mode);
                errL2 = matrix.calc_res_L2(var_mode);
                errW1 = matrix.calc_res_W1(var_mode);
            }

            norms[row][col][0] = errC;
            norms[row][col][1] = errL2;
            norms[row][col][2] = errW1;
        }
    }

    // Печать заголовка таблицы
    std::cout << "\n\\begin{tabular}{|c|c|c|c|c|}\n";
    std::cout << "\\hline\n";
    std::cout << "\\multicolumn{5}{|c|}{"
              << "$\\mu = " << mu << ",\\; "
              << (mode == 1 ? "p(\\rho) = " + std::to_string(Cp) + "\\rho" : "p(\\rho) = \\rho^{1.4}")
              << "$, невязка для $" << var_name << "$} \\\\\n";
    std::cout << "\\hline\n";
    std::cout << "$\\tau \\setminus h$ ";
    for (int col = 0; col < n_cols; ++col) {
        double h = X / (N_points[col] - 1);
        std::cout << " & $" << h << "$ ";
    }
    std::cout << "\\\\\n";
    std::cout << "\\hline\n";

    // Печать строк таблицы (по одной для каждой нормы)
    for (int row = 0; row < n_rows; ++row) {
        double tau = T / T_steps_list[row];
        for (int norm = 0; norm < 3; ++norm) {
            // Первая колонка: значение tau и название нормы
            if (norm == 0)
                std::cout << "$" << tau << "$ ";
            else
                std::cout << " ";

            // Значения для каждого столбца h
            for (int col = 0; col < n_cols; ++col) {
                std::cout << " & $" << std::scientific << std::setprecision(3)
                          << norms[row][col][norm] << "$ ";
            }
            std::cout << "\\\\\n";
        }
        std::cout << "\\hline\n";
    }
    std::cout << "\\end{tabular}\n\n";
}

// Основная тестовая функция
void test_scheme() {
    const double X = 1.0, Y = 1.0, T = 1.0;
    const double gamma = 1.4;

    std::vector<double> mu_values = {0.1, 0.01/*, 0.001*/};
    std::vector<double> Cp_values = {1.0, 10.0/*, 100.0*/};
    std::vector<int> modes = {0, 1};

    std::vector<int> N_points = {40, 60, 120, 160};
    std::vector<int> T_steps_list = {40, 60, 120, 160};

    struct VarInfo { int id; std::string name; };
    std::vector<VarInfo> vars = {{g, "G"}, {v1, "V_1"}, {v2, "V_2"}};
    
    for (int mode : modes) {
        for (double mu : mu_values) {
            if (mode == 1) {
                for (double Cp : Cp_values) {
                    for (const auto& var : vars) {
                        PrintTableForVariable(var.id, var.name, mu, Cp, mode,
                                              N_points, T_steps_list, X, Y, T, gamma);
                    }
                }
            } else {
                double Cp_dummy = 0.0;
                for (const auto& var : vars) {
                    PrintTableForVariable(var.id, var.name, mu, Cp_dummy, mode,
                                          N_points, T_steps_list, X, Y, T, gamma);
                }
            }
        }
    }
    std::cerr << "\nAll done!\n";
}

#include "matrix.h"
#include "funcs.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

void nested_grid_test() {
    const double X = 1.0, Y = 1.0, T = 1.0;
    const double mu = 0.1;                     // viscosity
    const double Cp = 10.0;                   // coefficient for p = Cp * ρ (linear case)
    const int gas_mode = 1;                   // 1 = linear equation of state
    const double gamma = 1.4;                 // not used when gas_mode=1

    // Base grids: number of intervals in x, y and time steps
    std::vector<int> base_Nx = {40, 60, 80};
    std::vector<int> base_Ny = {40, 60, 80};
    std::vector<int> base_Nt = {40, 60, 80};

    const int max_refinement = 2;              // refine 2 times (k=1,2)

    // Structure to hold solution data for a grid
    struct GridResult {
        std::vector<double> G, V1, V2;
        Mesh mesh;
        std::vector<block_status> status;
        double hx, hy;
    };

    // ----- 1. Compute base solutions -----
    std::vector<GridResult> base_results;
    for (size_t idx = 0; idx < base_Nx.size(); ++idx) {
        int Nx = base_Nx[idx];
        int Ny = base_Ny[idx];
        int Nt = base_Nt[idx];
        std::cout << "Computing base grid " << Nx << "x" << Ny << "..." << std::endl;

        P_gas gas(T, X, Y, Cp, gamma, mu, gas_mode);
        Mesh mesh(X, Y, T, Nx, Ny, Nt);
        Matrix matrix(gas, mesh);

        // Time stepping
        for (int step = 1; step <= Nt; ++step) {
            matrix.step = step;
            if (matrix.init_and_solve_G() == -1) {
                std::cerr << "Error solving for G on base grid" << std::endl;
                return;
            }
            if (matrix.init_and_solve_V() == -1) {
                std::cerr << "Error solving for V on base grid" << std::endl;
                return;
            }
        }

        GridResult res;
        res.G = matrix.solution_G;
        res.V1 = matrix.solution_V1;
        res.V2 = matrix.solution_V2;
        res.mesh = mesh;
        res.hx = mesh.h_x;
        res.hy = mesh.h_y;
        for (int i = 0; i < mesh.Dim; ++i) {
            res.status.push_back(mesh.mesh_points[i].status);
        }
        base_results.push_back(res);
    }

    // ----- 2. Compute refined solutions and differences -----
    // Norms: [base_idx][ref_level][norm], norm: 0=C, 1=L2
    std::vector<std::vector<std::vector<double>>> diff_G_norms(
        base_results.size(),
        std::vector<std::vector<double>>(max_refinement, std::vector<double>(2, 0.0)));
    std::vector<std::vector<std::vector<double>>> diff_V1_norms(
        base_results.size(),
        std::vector<std::vector<double>>(max_refinement, std::vector<double>(2, 0.0)));
    std::vector<std::vector<std::vector<double>>> diff_V2_norms(
        base_results.size(),
        std::vector<std::vector<double>>(max_refinement, std::vector<double>(2, 0.0)));

    for (size_t base_idx = 0; base_idx < base_results.size(); ++base_idx) {
        const auto& base = base_results[base_idx];
        int Nx_base = base.mesh.N;
        int Ny_base = base.mesh.M;
        int Nt_base = base.mesh.T_segm;

        for (int k = 1; k <= max_refinement; ++k) {
            int factor = 1 << k;   // 2^k
            int Nx_ref = Nx_base * factor;
            int Ny_ref = Ny_base * factor;
            int Nt_ref = Nt_base * factor;
            std::cout << "Computing refined grid for base " << Nx_base << "x" << Ny_base
                      << ", refinement " << k << " (factor " << factor << ")" << std::endl;

            P_gas gas(T, X, Y, Cp, gamma, mu, gas_mode);
            Mesh mesh_ref(X, Y, T, Nx_ref, Ny_ref, Nt_ref);
            Matrix matrix_ref(gas, mesh_ref);

            for (int step = 1; step <= Nt_ref; ++step) {
                matrix_ref.step = step;
                if (matrix_ref.init_and_solve_G() == -1) {
                    std::cerr << "Error solving for G on refined grid" << std::endl;
                    return;
                }
                if (matrix_ref.init_and_solve_V() == -1) {
                    std::cerr << "Error solving for V on refined grid" << std::endl;
                    return;
                }
            }

            // Compute differences at coarse grid nodes
            std::vector<double> G_diff(base.G.size(), 0.0);
            std::vector<double> V1_diff(base.V1.size(), 0.0);
            std::vector<double> V2_diff(base.V2.size(), 0.0);
            for (int i = 0; i < base.G.size(); ++i) {
                int base_i = base.mesh.mesh_points[i].i;
                int base_j = base.mesh.mesh_points[i].j;
                int refined_i = base_i * factor;
                int refined_j = base_j * factor;
                int refined_idx = refined_j * Nx_ref + refined_i;
                G_diff[i] = base.G[i] - matrix_ref.solution_G[refined_idx];
                V1_diff[i] = base.V1[i] - matrix_ref.solution_V1[refined_idx];
                V2_diff[i] = base.V2[i] - matrix_ref.solution_V2[refined_idx];
            }

            // Compute C and L2 norms (only on inner nodes)
            double max_G = 0.0, max_V1 = 0.0, max_V2 = 0.0;
            double l2_G = 0.0, l2_V1 = 0.0, l2_V2 = 0.0;
            for (int i = 0; i < base.G.size(); ++i) {
                if (base.status[i] == block_status::inner) {
                    double diffG = std::abs(G_diff[i]);
                    double diffV1 = std::abs(V1_diff[i]);
                    double diffV2 = std::abs(V2_diff[i]);
                    max_G = std::max(max_G, diffG);
                    max_V1 = std::max(max_V1, diffV1);
                    max_V2 = std::max(max_V2, diffV2);
                    l2_G += diffG * diffG;
                    l2_V1 += diffV1 * diffV1;
                    l2_V2 += diffV2 * diffV2;
                }
            }
            l2_G = std::sqrt(base.hx * base.hy * l2_G);
            l2_V1 = std::sqrt(base.hx * base.hy * l2_V1);
            l2_V2 = std::sqrt(base.hx * base.hy * l2_V2);

            diff_G_norms[base_idx][k-1][0] = max_G;
            diff_G_norms[base_idx][k-1][1] = l2_G;
            diff_V1_norms[base_idx][k-1][0] = max_V1;
            diff_V1_norms[base_idx][k-1][1] = l2_V1;
            diff_V2_norms[base_idx][k-1][0] = max_V2;
            diff_V2_norms[base_idx][k-1][1] = l2_V2;
        }
    }

    // ----- 3. Compute exact errors for base grids (for the last row) -----
    std::vector<std::vector<double>> exact_G_norms(base_results.size(), std::vector<double>(2, 0.0));
    std::vector<std::vector<double>> exact_V1_norms(base_results.size(), std::vector<double>(2, 0.0));
    std::vector<std::vector<double>> exact_V2_norms(base_results.size(), std::vector<double>(2, 0.0));
    for (size_t base_idx = 0; base_idx < base_results.size(); ++base_idx) {
        const auto& base = base_results[base_idx];
        double t = T;
        double max_G = 0.0, max_V1 = 0.0, max_V2 = 0.0;
        double l2_G = 0.0, l2_V1 = 0.0, l2_V2 = 0.0;
        for (int i = 0; i < base.G.size(); ++i) {
            if (base.status[i] == block_status::inner) {
                double x = base.mesh.mesh_points[i].x;
                double y = base.mesh.mesh_points[i].y;
                double exact_G = std::log(rho(t, x, y));
                double exact_V1 = u1(t, x, y);
                double exact_V2 = u2(t, x, y);
                double diffG = std::abs(base.G[i] - exact_G);
                double diffV1 = std::abs(base.V1[i] - exact_V1);
                double diffV2 = std::abs(base.V2[i] - exact_V2);
                max_G = std::max(max_G, diffG);
                max_V1 = std::max(max_V1, diffV1);
                max_V2 = std::max(max_V2, diffV2);
                l2_G += diffG * diffG;
                l2_V1 += diffV1 * diffV1;
                l2_V2 += diffV2 * diffV2;
            }
        }
        l2_G = std::sqrt(base.hx * base.hy * l2_G);
        l2_V1 = std::sqrt(base.hx * base.hy * l2_V1);
        l2_V2 = std::sqrt(base.hx * base.hy * l2_V2);
        exact_G_norms[base_idx][0] = max_G;
        exact_G_norms[base_idx][1] = l2_G;
        exact_V1_norms[base_idx][0] = max_V1;
        exact_V1_norms[base_idx][1] = l2_V1;
        exact_V2_norms[base_idx][0] = max_V2;
        exact_V2_norms[base_idx][1] = l2_V2;
    }

    // ----- 4. Print LaTeX tables for C norms (can be extended to L2) -----
    // Table for G (C norm)
    std::cout << "\n\\begin{table}[h]\n\\centering\n";
    std::cout << "\\caption{Сравнение решений на вложенных сетках для $G = \\ln\\rho$ ($C_h$-норма)}\n";
    std::cout << "\\begin{tabular}{|c|" << std::string(base_results.size(), 'c') << "|}\n\\hline\n";
    std::cout << "\\multicolumn{1}{|c|}{$k$} & ";
    for (size_t i = 0; i < base_results.size(); ++i) {
        std::cout << "$" << base_results[i].mesh.N-1 << "\\times" << base_results[i].mesh.M-1 << "$";
        if (i != base_results.size()-1) std::cout << " & ";
    }
    std::cout << " \\\\\n\\hline\n";
    for (int k = 1; k <= max_refinement; ++k) {
        std::cout << "$" << k << "$";
        for (size_t base_idx = 0; base_idx < base_results.size(); ++base_idx) {
            std::cout << " & $" << std::scientific << std::setprecision(3)
                      << diff_G_norms[base_idx][k-1][0] << "$";
        }
        std::cout << " \\\\\n";
    }
    std::cout << "\\hline\n";
    std::cout << "Exact &";
    for (size_t base_idx = 0; base_idx < base_results.size(); ++base_idx) {
        std::cout << " $" << std::scientific << std::setprecision(3) << exact_G_norms[base_idx][0] << "$";
        if (base_idx != base_results.size()-1) std::cout << " &";
    }
    std::cout << " \\\\\n\\hline\n";
    std::cout << "\\end{tabular}\n\\end{table}\n\n";

    // Similarly, tables for V1 and V2 (C norms)
    // ... (repeat pattern for V1 and V2)
    // For brevity, only G table is shown here; you can copy and adapt for V1, V2.

    std::cout << "\\begin{table}[h]\n\\centering\n";
    std::cout << "\\caption{Сравнение решений на вложенных сетках для $V_1$ ($C_h$-норма)}\n";
    std::cout << "\\begin{tabular}{|c|" << std::string(base_results.size(), 'c') << "|}\n\\hline\n";
    std::cout << "\\multicolumn{1}{|c|}{$k$} & ";
    for (size_t i = 0; i < base_results.size(); ++i) {
        std::cout << "$" << base_results[i].mesh.N-1 << "\\times" << base_results[i].mesh.M-1 << "$";
        if (i != base_results.size()-1) std::cout << " & ";
    }
    std::cout << " \\\\\n\\hline\n";
    for (int k = 1; k <= max_refinement; ++k) {
        std::cout << "$" << k << "$";
        for (size_t base_idx = 0; base_idx < base_results.size(); ++base_idx) {
            std::cout << " & $" << std::scientific << std::setprecision(3)
                      << diff_V1_norms[base_idx][k-1][0] << "$";
        }
        std::cout << " \\\\\n";
    }
    std::cout << "\\hline\n";
    std::cout << "Exact &";
    for (size_t base_idx = 0; base_idx < base_results.size(); ++base_idx) {
        std::cout << " $" << std::scientific << std::setprecision(3) << exact_V1_norms[base_idx][0] << "$";
        if (base_idx != base_results.size()-1) std::cout << " &";
    }
    std::cout << " \\\\\n\\hline\n";
    std::cout << "\\end{tabular}\n\\end{table}\n\n";

    std::cout << "\\begin{table}[h]\n\\centering\n";
    std::cout << "\\caption{Сравнение решений на вложенных сетках для $V_2$ ($C_h$-норма)}\n";
    std::cout << "\\begin{tabular}{|c|" << std::string(base_results.size(), 'c') << "|}\n\\hline\n";
    std::cout << "\\multicolumn{1}{|c|}{$k$} & ";
    for (size_t i = 0; i < base_results.size(); ++i) {
        std::cout << "$" << base_results[i].mesh.N-1 << "\\times" << base_results[i].mesh.M-1 << "$";
        if (i != base_results.size()-1) std::cout << " & ";
    }
    std::cout << " \\\\\n\\hline\n";
    for (int k = 1; k <= max_refinement; ++k) {
        std::cout << "$" << k << "$";
        for (size_t base_idx = 0; base_idx < base_results.size(); ++base_idx) {
            std::cout << " & $" << std::scientific << std::setprecision(3)
                      << diff_V2_norms[base_idx][k-1][0] << "$";
        }
        std::cout << " \\\\\n";
    }
    std::cout << "\\hline\n";
    std::cout << "Exact &";
    for (size_t base_idx = 0; base_idx < base_results.size(); ++base_idx) {
        std::cout << " $" << std::scientific << std::setprecision(3) << exact_V2_norms[base_idx][0] << "$";
        if (base_idx != base_results.size()-1) std::cout << " &";
    }
    std::cout << " \\\\\n\\hline\n";
    std::cout << "\\end{tabular}\n\\end{table}\n\n";
}
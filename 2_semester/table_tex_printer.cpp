#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "matrix.h"

/*
// Вспомогательная функция для печати одной таблицы для заданной переменной (mode = g, v1, v2)
void PrintTableForVariable(int var_mode,                      // g, v1 или v2
                           const std::string& var_name,       // название переменной для заголовка
                           double mu, double Cp, int mode,    // параметры газа
                           const std::vector<int>& N_points,
                           const std::vector<int>& T_steps_list,
                           double X, double Y, double T,
                           double gamma, double w, double rho_in)
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
            double h = X / Nx;

            // Создание объектов
            P_gas gas(T, X, Y, Cp, gamma, mu, mode);
            Mesh mesh(X, Y, T, Nx, Ny, T_steps);
            Matrix matrix(gas, mesh);
            mesh.mesh_points.clear();
            mesh.mesh_points.resize(mesh.Dim);
            fill_mesh_domain12(mesh.mesh_points, mesh.N, mesh.M, X, Y, w, rho_in);

            Matrix matrix(gas, mesh, w, rho_in, solver_type::library);
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
        double h = X / N_points[col] ;
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
*/
#include <iomanip>
#include <algorithm>
#include "funcs.h"

void PrintNestedGridTable(double mu, double Cp, int mode,
                          const std::vector<int>& base_N,
                          const std::vector<int>& base_T,
                          double X, double Y, double T,
                          double gamma, double omega, double rho_in)
{
    const double final_time = T;
    std::vector<int> factors = {2, 4, 0};   // 0 means exact
    const int n_factors = factors.size();
    const int n_base = base_N.size();
    const int n_vars = 3;   // 0:G, 1:V1, 2:V2

    // Store relative errors: [var][factor][base][norm]
    std::vector<std::vector<std::vector<std::vector<double>>>> rel_errors(
        n_vars, std::vector<std::vector<std::vector<double>>>(
            n_factors, std::vector<std::vector<double>>(
                n_base, std::vector<double>(3, 0.0))));

    // Loop over each base grid
    for (int b = 0; b < n_base; ++b) {
        int N_base = base_N[b];
        int T_base = base_T[b];
        double h = X / (N_base - 1);
        double tau = T / T_base;

        // Coarse grid (base)
        Mesh coarse_mesh(X, Y, T, N_base, N_base, T_base);
        coarse_mesh.mesh_points.clear();
        coarse_mesh.mesh_points.resize(coarse_mesh.Dim);
        fill_mesh_domain12(coarse_mesh.mesh_points, coarse_mesh.N, coarse_mesh.M, X, Y, omega, rho_in);
        P_gas coarse_gas(T, X, Y, Cp, gamma, mu, mode);
        Matrix coarse_mat(coarse_gas, coarse_mesh, omega, rho_in, solver_type::library);

        // Solve on coarse grid
        bool coarse_failed = false;
        for (int step = 1; step <= T_base; ++step) {
            coarse_mat.step = step;
            if (coarse_mat.init_and_solve_G() == -1 ||
                coarse_mat.init_and_solve_V() == -1) {
                coarse_failed = true;
                break;
            }
        }
        if (coarse_failed) {
            // Mark all entries as failed
            for (int v = 0; v < n_vars; ++v)
                for (int f = 0; f < n_factors; ++f)
                    for (int norm = 0; norm < 3; ++norm)
                        rel_errors[v][f][b][norm] = 1e10;
            continue;
        }

        // Norms of the coarse solution for each variable
        std::vector<double> zero(coarse_mat.Dim, 0.0);
        std::vector<double> coarse_normC(3), coarse_normL2(3), coarse_normW1(3);
        coarse_normC[0] = coarse_mat.calc_res_C1(g, zero);
        coarse_normL2[0] = coarse_mat.calc_res_L2(g, zero);
        coarse_normW1[0] = coarse_mat.calc_res_W1(g, zero);
        coarse_normC[1] = coarse_mat.calc_res_C1(v1, zero);
        coarse_normL2[1] = coarse_mat.calc_res_L2(v1, zero);
        coarse_normW1[1] = coarse_mat.calc_res_W1(v1, zero);
        coarse_normC[2] = coarse_mat.calc_res_C1(v2, zero);
        coarse_normL2[2] = coarse_mat.calc_res_L2(v2, zero);
        coarse_normW1[2] = coarse_mat.calc_res_W1(v2, zero);

        // For each refinement factor
        for (int f = 0; f < n_factors; ++f) {
            int factor = factors[f];
            // Get comparison values on coarse grid points
            std::vector<double> comp_G(coarse_mat.Dim);
            std::vector<double> comp_V1(coarse_mat.Dim);
            std::vector<double> comp_V2(coarse_mat.Dim);

            if (factor == 0) {
                // Exact solution
                for (int idx = 0; idx < coarse_mat.Dim; ++idx) {
                    Point p = coarse_mesh.mesh_points[idx];
                    comp_G[idx] = log(rho(final_time, p.x, p.y));
                    comp_V1[idx] = u1(final_time, p.x, p.y);
                    comp_V2[idx] = u2(final_time, p.x, p.y);
                }
            } else {
                // Refined grid
                int N_fine = (N_base - 1) * factor + 1;
                int T_fine = T_base * factor;
                Mesh fine_mesh(X, Y, T, N_fine, N_fine, T_fine);
                fine_mesh.mesh_points.clear();
                fine_mesh.mesh_points.resize(fine_mesh.Dim);
                fill_mesh_domain12(fine_mesh.mesh_points, fine_mesh.N, fine_mesh.M, X, Y, omega, rho_in);
                P_gas fine_gas(T, X, Y, Cp, gamma, mu, mode);
                Matrix fine_mat(coarse_gas, coarse_mesh, omega, rho_in, solver_type::library);

                bool fine_failed = false;
                for (int step = 1; step <= T_fine; ++step) {
                    fine_mat.step = step;
                    if (fine_mat.init_and_solve_G() == -1 ||
                        fine_mat.init_and_solve_V() == -1) {
                        fine_failed = true;
                        break;
                    }
                }
                if (fine_failed) {
                    for (int v = 0; v < n_vars; ++v)
                        for (int norm = 0; norm < 3; ++norm)
                            rel_errors[v][f][b][norm] = 1e10;
                    continue;
                }

                // Restrict fine solution to coarse grid
                for (int j = 0; j < N_base; ++j) {
                    int j_fine = j * factor;
                    for (int i = 0; i < N_base; ++i) {
                        int i_fine = i * factor;
                        int idx_coarse = j * N_base + i;
                        int idx_fine = j_fine * N_fine + i_fine;
                        comp_G[idx_coarse] = fine_mat.solution_G[idx_fine];
                        comp_V1[idx_coarse] = fine_mat.solution_V1[idx_fine];
                        comp_V2[idx_coarse] = fine_mat.solution_V2[idx_fine];
                    }
                }
            }

            // Compute differences and relative errors
            double errC_G = coarse_mat.calc_res_C1(g, comp_G);
            double errL2_G = coarse_mat.calc_res_L2(g, comp_G);
            double errW1_G = coarse_mat.calc_res_W1(g, comp_G);
            double errC_V1 = coarse_mat.calc_res_C1(v1, comp_V1);
            double errL2_V1 = coarse_mat.calc_res_L2(v1, comp_V1);
            double errW1_V1 = coarse_mat.calc_res_W1(v1, comp_V1);
            double errC_V2 = coarse_mat.calc_res_C1(v2, comp_V2);
            double errL2_V2 = coarse_mat.calc_res_L2(v2, comp_V2);
            double errW1_V2 = coarse_mat.calc_res_W1(v2, comp_V2);

            // Store relative errors (divide by coarse norm)
            rel_errors[0][f][b][0] = errC_G / coarse_normC[0];
            rel_errors[0][f][b][1] = errL2_G / coarse_normL2[0];
            rel_errors[0][f][b][2] = errW1_G / coarse_normW1[0];
            rel_errors[1][f][b][0] = errC_V1 / coarse_normC[1];
            rel_errors[1][f][b][1] = errL2_V1 / coarse_normL2[1];
            rel_errors[1][f][b][2] = errW1_V1 / coarse_normW1[1];
            rel_errors[2][f][b][0] = errC_V2 / coarse_normC[2];
            rel_errors[2][f][b][1] = errL2_V2 / coarse_normL2[2];
            rel_errors[2][f][b][2] = errW1_V2 / coarse_normW1[2];
        }
    }

    // Now print tables for each variable
    const char* var_names[] = {"G", "V_1", "V_2"};
    const int var_ids[] = {g, v1, v2};

    for (int v = 0; v < n_vars; ++v) {
        // Header
        std::cout << "\n\\begin{tabular}{|c|c|c|c|c|}\n";
        std::cout << "\\hline\n";
        if (mode == 1)
            std::cout << "\\multicolumn{5}{|c|}{$\\mu = " << mu << ",\\; p(\\rho) = " << Cp << "\\rho$, невязка для $" << var_names[v] << "$} \\\\\n";
        else
            std::cout << "\\multicolumn{5}{|c|}{$\\mu = " << mu << ",\\; p(\\rho) = \\rho^{1.4}$, невязка для $" << var_names[v] << "$} \\\\\n";
        std::cout << "\\hline\n";
        std::cout << "$\\tau \\setminus h$ ";
        for (int b = 0; b < n_base; ++b) {
            double h = X / (base_N[b] - 1);
            std::cout << " & $" << std::scientific << std::setprecision(3) << h << "$ ";
        }
        std::cout << "\\\\\n\\hline\n";

        // Rows: for each factor (refinement level)
        for (int f = 0; f < n_factors; ++f) {
            double tau = T / base_T[0];  
            if (f < n_factors - 1) {
                std::cout << "$v - v^{" << (f+1) << "}$ ";
            } else {
                std::cout << "$v - u$ ";
            }
            // Print the three norms for this row (C, L2, W1) across columns
            for (int norm = 0; norm < 3; ++norm) {
                // We need to output three lines per row (one for each norm)
                if (norm == 0)
                    std::cout << " & $" << std::scientific << std::setprecision(3) << rel_errors[v][f][0][norm] << "$ ";
                else
                    std::cout << " & $" << rel_errors[v][f][0][norm] << "$ ";
                for (int b = 1; b < n_base; ++b) {
                    std::cout << " & $" << rel_errors[v][f][b][norm] << "$ ";
                }
                std::cout << "\\\\\n";
                if (norm == 2 && f != n_factors-1) std::cout << "\\hline\n";
            }
        }
        std::cout << "\\hline\n\\end{tabular}\n\n";
    }
}

void RunNestedGridTests()
{
    double X = 3.0, Y = 3.0, T = 3.0;
    double gamma = 1.4;
    double omega = 0.1, rho = 1.0;
    std::vector<double> mu_values = {0.1, 0.01};
    std::vector<double> Cp_values = {1.0, 10.0};
    std::vector<int> modes = {0, 1};

    std::vector<int> base_N = {40, 60, 80};
    std::vector<int> base_T = {40, 60, 80};

    for (int mode : modes) {
        for (double mu : mu_values) {
            if (mode == 1) {
                for (double Cp : Cp_values) {
                    PrintNestedGridTable(mu, Cp, mode, base_N, base_T, X, Y, T, gamma, omega, rho);
                }
            } else {
                double Cp_dummy = 0.0;
                PrintNestedGridTable(mu, Cp_dummy, mode, base_N, base_T, X, Y, T, gamma, omega, rho);
            }
        }
    }
}

/*
// Основная тестовая функция
void test_scheme() {
    const double X = 1.0, Y = 1.0, T = 1.0;
    const double gamma = 1.4;

    std::vector<double> mu_values = {0.1, 0.01};
    std::vector<double> Cp_values = {1.0, 10.0};
    std::vector<int> modes = {0, 1};

    std::vector<int> N_points = {40, 80, 160, 320};
    std::vector<int> T_steps_list = {40, 80, 160, 320};

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
*/
/*
void test_speed ()
{
                            
    double mu = 0.1; double Cp = 1; int mode = 1;  
    std::vector<int> N_points = {100, 200, 300};
    std::vector<int> T_steps_list = {100, 200, 30};
    double X = 1; double Y = 1; double T = 1;
    double gamma = 1.4;
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
            Matrix matrix1(gas, mesh, solver_type::own);
            Matrix matrix2(gas, mesh, solver_type::library);
            // Расчёт по всем временным слоям

            bool success = true;
            auto own_time = matrix1.solve_scheme ();
            auto lib_time = matrix2.solve_scheme ();

            norms[row][col][0] = own_time;
            norms[row][col][1] = lib_time;
        }
    }

    // Печать заголовка таблицы
    std::cout << "\n\\begin{tabular}{|c|c|c|c|c|}\n";
    std::cout << "\\hline\n";
    std::cout << "\\multicolumn{5}{|c|}{"
              << "$\\mu = " << mu << ",\\; "
              << (mode == 1 ? "p(\\rho) = " + std::to_string(Cp) + "\\rho" : "p(\\rho) = \\rho^{1.4}")
              <<  "$} \\\\\n";
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
        for (int norm = 0; norm < 2; ++norm) {
            // Первая колонка: значение tau и название нормы
            if (norm == 0)
                std::cout << "$" << tau << "$ ";
            else
                std::cout << " ";

            // Значения для каждого столбца h
            for (int col = 0; col < n_cols; ++col) {
                std::cout << " & $"  << norms[row][col][norm] << "$ ";
            }
            std::cout << "\\\\\n";
        }
        std::cout << "\\hline\n";
    }
    std::cout << "\\end{tabular}\n\n";
}
    */
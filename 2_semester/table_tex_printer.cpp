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
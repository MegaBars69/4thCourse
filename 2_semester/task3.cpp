#include "task3.h"
#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

// Вспомогательная функция: расчёт времени стабилизации для фиксированных параметров
double compute_steady_time(double mu, double Cp, int mode,
                           double w, double rho_in,
                           int ramp_steps,
                           int X_segm, int Y_segm, int T_steps_max,
                           double X, double Y, double T,
                           double tol_steady)
{
    P_gas gas(T, X, Y, Cp, 1.4, mu, mode);
    Mesh mesh(X, Y, T, X_segm, Y_segm, T_steps_max);
    mesh.mesh_points.clear();
    mesh.mesh_points.resize(mesh.Dim);
    fill_mesh_domain12(mesh.mesh_points, mesh.N, mesh.M, X, Y, w, rho_in);

    Matrix matrix(gas, mesh, w, rho_in, solver_type::library);
    matrix.ramp_steps = ramp_steps;

    double t_steady = -1.0;
    bool converged = matrix.solve_to_steady(T_steps_max, tol_steady, t_steady);
    return converged ? t_steady : -1.0;
}

// ---------- Пункт 1: таблица времени стабилизации ----------
void task3_steady_time_table()
{
    const double X = 3.0, Y = 3.0, T = 10.0;   // большое конечное время
    const double mu = 0.1, eta = 0.1, Cp = 1.0;
    const int mode = 1;   // p = Cp * rho
    const int X_segm = 50, Y_segm = 50;   // грубая сетка для быстроты
    const int T_steps_max = 500;          // максимум шагов
    const double tol_steady = 1e-2;

    // Наборы входных параметров из отчёта (пример)
    std::vector<double> w_list = {0.1, 0.5, 1, 2};      // скорости w
    std::vector<double> rho_list = {1, 2, 3, 4, 5};    // плотности на входе
    std::vector<int> ramp_list = {0, 10, 100};               // способы запуска

    // Печать заголовка таблицы LaTeX
    std::cout << "\\begin{tabular}{|c|c|c|c|c|}\n\\hline\n";
    std::cout << "\\multicolumn{5}{|c|}{Время установления, $\\mu=" << mu
              << ", \\; \\eta=" << eta << ", \\; C_\\rho=" << Cp
              << ", \\; \\Phi = e^G, \\; M1=M2=" << X_segm << "$}\\\\\n\\hline\n";
    std::cout << "$\\omega_0$ & $\\rho_0$ & $t$ (сразу) & $t$ (10 шагов) & $t$ (100 шагов)\\\\\n\\hline\n";

    for (size_t i = 0; i < w_list.size(); ++i) {
        double w = w_list[i];
        double rho = rho_list[i];
        std::cout << w << " & " << rho;
        for (int r : ramp_list) {
            double t = compute_steady_time(mu, Cp, mode, w, rho, r,
                                          X_segm, Y_segm, T_steps_max,
                                          X, Y, T, tol_steady);
            if (t < 0) std::cout << " & –";
            else std::cout << " & $" << std::fixed << std::setprecision(3) << t << "$";
        }
        std::cout << "\\\\\n";
    }
    std::cout << "\\hline\n\\end{tabular}\n";
}

// ---------- Пункт 2: визуализация (данные для gnuplot) ----------
void task3_visualization(double w, double rho_in, int ramp_steps)
{
    // ===== Параметры расчёта =====
    const double T = 10.0;               // конечное время моделирования
    const double X = 3.0, Y = 2.0;      // область
    const double mu = 0.1, Cp = 1.0;
    const int mode = 1;                 // p = Cp * ρ
    const int X_segm = 30, Y_segm = 20; // сетка (можно менять)
    const int T_steps = 300;            // число шагов по времени (желательно кратно 4)

    const double tau = T / T_steps;

    // Моменты времени, в которые нужно сохранить кадры
    const std::vector<double> snapshots = { T/4., T/2., 3.*T/4., T };
    std::vector<int> snapshot_steps;
    for (double t : snapshots) {
        int step = static_cast<int>(round(t / tau));
        if (step < 1) step = 1;
        if (step > T_steps) step = T_steps;
        snapshot_steps.push_back(step);
    }

    // ===== Инициализация задачи =====
    P_gas gas(T, X, Y, Cp, 1.4, mu, mode);
    Mesh mesh(X, Y, T, X_segm, Y_segm, T_steps);
    mesh.mesh_points.clear();
    mesh.mesh_points.resize(mesh.Dim);
    fill_mesh_domain12(mesh.mesh_points, mesh.N, mesh.M, X, Y, w, rho_in);

    Matrix matrix(gas, mesh, w, rho_in, solver_type::library);
    matrix.ramp_steps = ramp_steps;

    // ===== Временной цикл с сохранением кадров =====
    int next_snap = 0;
    for (int step = 1; step <= T_steps; ++step) {
        matrix.step = step;
        if (matrix.init_and_solve_G() == -1) {
            std::cerr << "G solver failed at step " << step << std::endl;
            return;
        }
        if (matrix.init_and_solve_V() == -1) {
            std::cerr << "V solver failed at step " << step << std::endl;
            return;
        }

        if (next_snap < (int)snapshot_steps.size() && step == snapshot_steps[next_snap]) {
            double current_t = step * tau;
            char fname[100];
            sprintf(fname, "result_T%.2f.dat", current_t);
            std::ofstream fout(fname);
            for (int j = 0; j < mesh.M; ++j) {
                for (int i = 0; i < mesh.N; ++i) {
                    int idx = j * mesh.N + i;
                    Point p = mesh.mesh_points[idx];
                    if (p.status != block_status::outer) {
                        fout << p.x << " " << p.y << " "
                             << exp(matrix.solution_G[idx]) << " "
                             << matrix.solution_V1[idx] << " "
                             << matrix.solution_V2[idx] << "\n";
                    } else {
                        fout << p.x << " " << p.y << " 0 0 0\n";
                    }
                }
                fout << "\n";
            }
            fout.close();
            std::cout << "Сохранён кадр t = " << current_t << " в " << fname << std::endl;
            next_snap++;
        }
    }

    // ===== Генерация скрипта gnuplot с ПРАВИЛЬНЫМ совмещением =====
    std::ofstream gp("plot_instants.gp");
    gp << "# Корректное совмещение плотности и векторов скорости\n";
    gp << "set terminal pngcairo size 800,600 enhanced font 'Arial,12'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set size ratio -1\n";
    gp << "set xrange [0:" << X << "]\n";
    gp << "set yrange [0:" << Y << "]\n";
    gp << "set view map\n";                // вид сверху для 3D-графика
    gp << "set palette defined (0 \"blue\", 1 \"red\")\n";
    gp << "\n";

    // Для каждого момента создаём отдельное изображение
    for (size_t i = 0; i < snapshots.size(); ++i) {
        double t = snapshots[i];
        char fname[100];
        sprintf(fname, "result_T%.2f.dat", t);
        gp << "set output 'frame_T" << i+1 << ".png'\n";
        gp << "set title 'Момент времени T/" << (i+1 == 4 ? "1" : (i==0 ? "4" : (i==1 ? "2" : "3*4"))) << "'\n";
        // Ключевая команда: единый splot для pm3d (плотность) и vectors (скорость)
        gp << "splot '" << fname << "' using 1:2:3 with pm3d notitle, \\\n"
           << "      '" << fname << "' using 1:2:(0.1):4:5:(0.0) with vectors head filled lc rgb 'black' notitle\n\n";
    }
    gp.close();

    std::cout << "\nДанные сохранены. Для построения графиков выполните:\n";
    std::cout << "  gnuplot plot_instants.gp\n";
}

// ---------- Пункт 3: таблица сходимости (зависимость времени установления) ----------
std::array<double, 3> compute_diff_norms(
    Matrix& coarse,
    Mesh& fine_mesh, const std::vector<double>& fine_G,
    std::vector<double>& fine_V1, const std::vector<double>& fine_V2,
    int var_mode,  // g, v1, v2
    int factor)    // кратность измельчения
{
    std::vector<double> comp(coarse.Dim);
    int Nx_coarse = coarse.mesh.N, Ny_coarse = coarse.mesh.M;
    for (int j = 0; j < Ny_coarse; ++j) {
        int j_f = j * factor;
        for (int i = 0; i < Nx_coarse; ++i) {
            int i_f = i * factor;
            int idx_coarse = j * Nx_coarse + i;
            int idx_fine = j_f * fine_mesh.N + i_f;
            if (var_mode == g) comp[idx_coarse] = fine_G[idx_fine];
            else if (var_mode == v1) comp[idx_coarse] = fine_V1[idx_fine];
            else comp[idx_coarse] = fine_V2[idx_fine];
        }
    }

    double eC, eL2, eW1;
    if (var_mode == g) {
        eC  = coarse.calc_res_C1(g, comp);
        eL2 = coarse.calc_res_L2(g, comp);
        eW1 = coarse.calc_res_W1(g, comp);
    } else if (var_mode == v1) {
        eC  = coarse.calc_res_C1(v1, comp);
        eL2 = coarse.calc_res_L2(v1, comp);
        eW1 = coarse.calc_res_W1(v1, comp);
    } else {
        eC  = coarse.calc_res_C1(v2, comp);
        eL2 = coarse.calc_res_L2(v2, comp);
        eW1 = coarse.calc_res_W1(v2, comp);
    }
    return {eC, eL2, eW1};
}

void task3_convergence_table()
{
    // Параметры задачи (как в отчёте)
    const double T = 1.0, X = 3.0, Y = 2.0;
    const double mu = 0.1, eta = 0.1, C_rho = 1.0;
    const int mode = 1;    // p = C_rho * rho
    const double w = 3.0, rho_in = 6.0;
    const int ramp = 0;    // граничные условия подаются сразу

    // Базовые шаги сетки по пространству и времени
    const std::vector<double> h_base = {0.05, 0.025, 0.0125, 0.00625};
    const std::vector<double> tau_base = h_base;  // τ = h (квадратная сетка во времени)

    // Множители измельчения для сравнения
    const std::vector<int> factors = {2, 4, 0};  // 0 означает сравнение с точным решением

    // Для печати таблиц: имена переменных и идентификаторы
    const std::array<std::string, 3> var_names = {"G - g", "V_1 - u_1", "V_2 - u_2"};
    const std::array<int, 3> var_ids = {g, v1, v2};

    // Цикл по переменным (три таблицы)
    for (size_t v = 0; v < 3; ++v) {
        std::cout << "\n\\begin{tabular}{|c|c|c|c|c|}\n\\hline\n";
        std::cout << "\\multicolumn{5}{|c|}{$\\mu=" << mu
                  << ",\\; \\eta=" << eta
                  << ",\\; C_\\rho=" << C_rho
                  << ",\\; \\Phi=e^G,\\; p(\\rho)=" << C_rho << "\\rho$, невязка для $"
                  << var_names[v] << "$}\\\\\n\\hline\n";
        std::cout << "$\\tau \\setminus h$";
        for (double h : h_base) std::cout << " & $" << h << "$ ";
        std::cout << "\\\\\n\\hline\n";

        // Строки таблицы соответствуют τ = h (по строчкам)
        for (size_t row = 0; row < tau_base.size(); ++row) {
            double tau = tau_base[row];
            double h = h_base[row];
            // Целое число шагов по пространству и времени
            int Nx = static_cast<int>(X / h + 0.5) + 1;
            int Ny = static_cast<int>(Y / h + 0.5) + 1;
            int T_steps = static_cast<int>(T / tau + 0.5);

            // Расчёт на крупной (базовой) сетке
            P_gas gas(T, X, Y, C_rho, 1.4, mu, mode);
            Mesh mesh(X, Y, T, Nx, Ny, T_steps);
            mesh.mesh_points.clear();
            mesh.mesh_points.resize(mesh.Dim);
            fill_mesh_domain12(mesh.mesh_points, mesh.N, mesh.M, X, Y, w, rho_in);
            Matrix coarse(gas, mesh, w, rho_in, solver_type::library);
            coarse.ramp_steps = ramp;
            bool coarse_ok = true;
            for (int step = 1; step <= T_steps; ++step) {
                coarse.step = step;
                if (coarse.init_and_solve_G() == -1 || coarse.init_and_solve_V() == -1) {
                    coarse_ok = false;
                    break;
                }
            }
            if (!coarse_ok) {
                std::cout << "$" << tau << "$ & error ...\\\\\n\\hline\n";
                continue;
            }

            // Вычисляем три нормы разности для каждого фактора
            for (size_t f = 0; f < factors.size(); ++f) {
                int factor = factors[f];
                std::array<double, 3> norms;
                if (factor == 0) {
                    // Сравнение с точным решением
                    if (v == 0) {
                        norms[0] = coarse.calc_res_C1(g);
                        norms[1] = coarse.calc_res_L2(g);
                        norms[2] = coarse.calc_res_W1(g);
                    } else if (v == 1) {
                        norms[0] = coarse.calc_res_C1(v1);
                        norms[1] = coarse.calc_res_L2(v1);
                        norms[2] = coarse.calc_res_W1(v1);
                    } else {
                        norms[0] = coarse.calc_res_C1(v2);
                        norms[1] = coarse.calc_res_L2(v2);
                        norms[2] = coarse.calc_res_W1(v2);
                    }
                } else {
                    // Измельчённая сетка
                    int Nx_fine = (Nx-1)*factor + 1;
                    int Ny_fine = (Ny-1)*factor + 1;
                    int T_fine = T_steps * factor;
                    Mesh mesh_fine(X, Y, T, Nx_fine, Ny_fine, T_fine);
                    mesh_fine.mesh_points.clear();
                    mesh_fine.mesh_points.resize(mesh_fine.Dim);
                    fill_mesh_domain12(mesh_fine.mesh_points, mesh_fine.N, mesh_fine.M, X, Y, w, rho_in);
                    Matrix fine(gas, mesh_fine, w, rho_in, solver_type::library);
                    fine.ramp_steps = ramp;
                    bool fine_ok = true;
                    for (int step = 1; step <= T_fine; ++step) {
                        fine.step = step;
                        if (fine.init_and_solve_G() == -1 || fine.init_and_solve_V() == -1) {
                            fine_ok = false;
                            break;
                        }
                    }
                    if (!fine_ok) {
                        norms = {1e10, 1e10, 1e10};
                    } else {
                        norms = compute_diff_norms(coarse, mesh_fine,
                            fine.solution_G, fine.solution_V1, fine.solution_V2,
                            var_ids[v], factor);
                    }
                }

                // Печатаем строку с тремя нормами
                if (f == 0)
                    std::cout << "$v - u$";     // сравнение с точным
                else
                    std::cout << "$v - v^{(" << factor << ")}$";
                for (int k = 0; k < 3; ++k) {
                    if (norms[k] >= 1e9) std::cout << " & $–$ ";
                    else std::cout << " & $" << std::scientific << std::setprecision(3) << norms[k] << "$ ";
                }
                std::cout << "\\\\\n";
            }
            std::cout << "\\hline\n";
        }
        std::cout << "\\end{tabular}\n\n";
    }
}
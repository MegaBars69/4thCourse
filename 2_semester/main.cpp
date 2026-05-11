
#include <iostream>
#include <sstream>       // для std::istringstream
#include <vector>
#include <string>
#include <iomanip>       // для std::setprecision
#include "matrix.h"
#include "table_tex_printer.h"


#include "task3.h"

int main()
{
    //task3_steady_time_table();
    task3_visualization(10, 10.0, 0);   // пример: w=3, ρ₀=6, сразу полный поток
    //task3_convergence_table();

    return 0;
}
/*
// Вспомогательная функция для разбора строки с числами через запятую
std::vector<double> parse_csv(const std::string& str) {
    std::vector<double> vals;
    std::istringstream ss(str);
    std::string token;
    while (std::getline(ss, token, ',')) {
        vals.push_back(std::stod(token));
    }
    return vals;
}

int main(int argc, char const *argv[]) 
{
    
        if (argc != 13) {   // исправлено: 1(программа) + 1(ключ) + 11 параметров = 13
            std::cout << "argc = " << argc << "\n";
            printf("Usage: %s --steady-table mu Cp X Y T X_segm Y_segm T_steps mode \"w1,w2,...\" \"rho1,rho2,...\"\n", argv[0]);
            return -1;
        }
        double mu, Cp, X, Y, T;
        int X_segm, Y_segm, T_steps, mode;
        if (sscanf(argv[2], "%lf", &mu) != 1 ||
            sscanf(argv[3], "%lf", &Cp) != 1 ||
            sscanf(argv[4], "%lf", &X) != 1 ||
            sscanf(argv[5], "%lf", &Y) != 1 ||
            sscanf(argv[6], "%lf", &T) != 1 ||
            sscanf(argv[7], "%d", &X_segm) != 1 ||
            sscanf(argv[8], "%d", &Y_segm) != 1 ||
            sscanf(argv[9], "%d", &T_steps) != 1 ||
            sscanf(argv[10], "%d", &mode) != 1) {
            printf("Can't read base parameters!\n");
            return -2;
        }
        std::vector<double> w_list = parse_csv(argv[11]);      // список скоростей w
        std::vector<double> rho_in_list = parse_csv(argv[12]); // список плотностей rho_in

        const double tol_steady = 1e-2;
        const int max_steps = T_steps;   // используем заданное число шагов как максимум

        // Таблица: строки – w, столбцы – rho_in
        std::cout << "\\begin{tabular}{|c|";
        for (size_t i = 0; i < rho_in_list.size(); ++i) std::cout << "c|";
        std::cout << "}\n\\hline\n";
        std::cout << "\\multicolumn{" << rho_in_list.size()+1 << "}{|c|}{"
                  << "Время установления, $\\mu=" << mu << ", $\\h=" << X/X_segm << "$, $\\tau=" << T/T_steps << "$,$p(\\rho)=";
        if (mode == 1) std::cout << Cp << "\\rho";
        else std::cout << "\\rho^{1.4}";
        std::cout << "$}\\\\\n\\hline\n";
        std::cout << "$w \\setminus \\rho_{\\text{in}}$";
        for (double rho : rho_in_list) std::cout << " & $" << rho << "$";
        std::cout << "\\\\\n\\hline\n";

        for (double w : w_list) {
            std::cout << "$" << w << "$";
            for (double rho_in : rho_in_list) {
                // Создаём объекты для текущей пары (w, rho_in)
                P_gas gas(T, X, Y, Cp, 1.4, mu, mode);
                Mesh mesh(X, Y, T, X_segm, Y_segm, T_steps);
                mesh.mesh_points.clear();
                mesh.mesh_points.resize(mesh.Dim);
                fill_mesh_domain12(mesh.mesh_points, mesh.N, mesh.M, X, Y, w, rho_in);

                Matrix matrix(gas, mesh, w, rho_in, solver_type::library);
                double t_steady = -1.0;
                bool converged = matrix.solve_to_steady(max_steps, tol_steady, t_steady);
                std::cout << " & $" << t_steady << "$";

            }
            std::cout << "\\\\\n\\hline\n";
        }
        std::cout << "\\end{tabular}\n";
        return 0;
}

int main(int argc, char const *argv[]) {
    double mu;
    double C_rho;
    double X = 3.0, Y = 3.0, T = 10.0;   // время увеличено для установления
    int X_segm;
    int Y_segm;
    int T_steps;
    int mode;

    double w;
    double rho_in;

    if (argc != 12) {
        printf("Usage: %s mu C_rho X Y T X_segm Y_segm T_steps mode w rho_in\n", argv[0]);
        return -1;
    }

    if ((sscanf(argv[1], "%lf", &mu) != 1) ||
        (sscanf(argv[2], "%lf", &C_rho) != 1) ||
        (sscanf(argv[3], "%lf", &X) != 1) ||
        (sscanf(argv[4], "%lf", &Y) != 1) ||
        (sscanf(argv[5], "%lf", &T) != 1) ||
        (sscanf(argv[6], "%d", &X_segm) != 1) ||
        (sscanf(argv[7], "%d", &Y_segm) != 1) ||
        (sscanf(argv[8], "%d", &T_steps) != 1) ||
        (sscanf(argv[9], "%d", &mode) != 1) ||
        (sscanf(argv[10], "%lf", &w) != 1) ||
        (sscanf(argv[11], "%lf", &rho_in) != 1)) {
        printf("Can't read initial values!\n");
        return -2;
    }

    P_gas gas(T, X, Y, C_rho, 1.4, mu, mode);
    Mesh mesh(X, Y, T, X_segm, Y_segm, T_steps);

    // Построение сетки области 12
    mesh.mesh_points.clear();
    mesh.mesh_points.resize(mesh.Dim);
    fill_mesh_domain12(mesh.mesh_points, mesh.N, mesh.M, X, Y, w, rho_in);

    Matrix matrix(gas, mesh, w, rho_in, solver_type::library);

    // Для отслеживания сходимости к стационару
    double prev_norm_G = 0.0;
    double tol_steady = 1e-6;

    clock_t start1 = clock();
    for (int step = 1; step <= T_steps; step++) {
        matrix.step = step;

        // Сохраняем старое G до решения
        std::vector<double> old_G(matrix.Dim);
        for (int i = 0; i < matrix.Dim; ++i)
            old_G[i] = matrix.func_points[i].G;

        if (matrix.init_and_solve_G() == -1) {
            printf("G solver failed at step %d\n", step);
            return -1;
        }
        if (matrix.init_and_solve_V() == -1) {
            printf("V solver failed at step %d\n", step);
            return -1;
        }

        // Вычисляем норму изменения G (используем новые значения func_points.G)
        double norm_G = 0.0;
        for (int i = 0; i < matrix.Dim; ++i) {
            double diff = matrix.func_points[i].G - old_G[i];
            norm_G += diff * diff;
        }
        norm_G = sqrt(norm_G / matrix.Dim);

        printf("Step %3d: dG = %.3e | G_max = %.3e | V1_max = %.3e | V2_max = %.3e\n",
               step, norm_G,
               *std::max_element(matrix.solution_G.begin(), matrix.solution_G.end()),
               *std::max_element(matrix.solution_V1.begin(), matrix.solution_V1.end()),
               *std::max_element(matrix.solution_V2.begin(), matrix.solution_V2.end()));

        if (step > 10 && fabs(norm_G - prev_norm_G) < tol_steady * prev_norm_G) {
            printf("Steady state reached at step %d\n", step);
            break;
        }
        prev_norm_G = norm_G;
    }

    clock_t end1 = clock();
    auto t1 = static_cast<double>(end1 - start1) / CLOCKS_PER_SEC;
    std::cout << "Time: " << t1 << " s\n";

    FILE* f = fopen("result.dat", "w");
    for (int j = 0; j < mesh.M; ++j) {
        for (int i = 0; i < mesh.N; ++i) {
            int idx = j * mesh.N + i;
            Point p = mesh.mesh_points[idx];
            // Для всех узлов записываем координаты, а значения только для не-внешних
            if (p.status != block_status::outer) {
                fprintf(f, "%f %f %f %f %f\n", p.x, p.y,
                        exp(matrix.solution_G[idx]),
                        matrix.solution_V1[idx],
                        matrix.solution_V2[idx]);
            } else {
                // Для внешних узлов записываем только координаты (значения не пишем)
                // Gnuplot пропустит строки с недостаточным количеством столбцов при using 1:2:3
                fprintf(f, "%f %f\n", p.x, p.y);
            }
        }
        fprintf(f, "\n");   // пустая строка между рядами сетки обязательна!
    }
    fclose(f);

    return 0;
}


int main ()
{
    RunNestedGridTests ();
    return 0;
}
*/
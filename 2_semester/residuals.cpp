#include "matrix.h"

double calc_res_point (int mode, double t, double x, double y, double val)
{
    if (mode == g)
        return  fabs (val - log (rho (t, x, y)));
    if (mode == v1)
        return fabs (val - u1 (t, x, y));
    if (mode == v2)
        return fabs (val - u2 (t, x, y));

    return 0;

}

double Matrix::calc_res_C1 (int mode)
{
    double max = -1;
    int max_i = -1;
    for (int i = 0; i < Dim; i++ )
    {
        Point check_point = mesh.mesh_points[i];
        if (check_point.status != block_status::outer)
        {
            double res {};
            if (mode == g)
              res = fabs (solution_G[i] - log (rho (step * mesh.tau, check_point.x, check_point.y)));
            if (mode == v1)
                res = fabs (solution_V1[i] - u1 (step * mesh.tau, check_point.x, check_point.y));
            if (mode == v2)
                res = fabs (solution_V2[i] - u2 (step * mesh.tau, check_point.x, check_point.y));

            if (res > max)
            {
                max = res;
                max_i = i;
            }
        }
    }
    return max;
}

double Matrix::calc_res_L2 (int mode)
{
    double sum = 0;
    for (int i = 0; i < Dim; i++ )
    {
        Point check_point = mesh.mesh_points[i];
        if (check_point.status == block_status::inner)
        {
            double res {};
            if (mode == g)
                res = fabs (solution_G[i] - log (rho (step * mesh.tau, check_point.x, check_point.y)));
            if (mode == v1)
                res = fabs (solution_V1[i] - u1 (step * mesh.tau, check_point.x, check_point.y));
            if (mode == v2)
                res = fabs (solution_V2[i] - u2 (step * mesh.tau, check_point.x, check_point.y));

            sum += res * res;
        }
    }
    return sqrt(mesh.h_x * mesh.h_y * sum);
}

double Matrix::calc_res_W1 (int mode)
{
    double sum = 0;
    double tau = mesh.tau * step;
    double prev_res {};
    if (mode == g)
        prev_res = fabs (solution_G[0] - log (rho (tau, 0, 0)));
    if (mode == v1)
        prev_res = fabs (solution_V1[0] - u1 (tau, 0, 0));
    if (mode == v2)
        prev_res = fabs (solution_V2[0] - u2 (tau, 0, 0));
    for (int i = 1; i < Dim; i++)
    {
        double res {};
        Point check_point = mesh.mesh_points[i];
        if (mode == g)
            res = fabs (solution_G[i] - log (rho (step * mesh.tau, check_point.x, check_point.y)));
        if (mode == v1)
            res = fabs (solution_V1[i] - u1 (step * mesh.tau, check_point.x, check_point.y));
        if (mode == v2)
            res = fabs (solution_V2[i] - u2 (step * mesh.tau, check_point.x, check_point.y));

        sum += (res - prev_res) *
               (res - prev_res);
        prev_res = res;
    }
    double res_l2 = calc_res_L2(mode);
    sum = res_l2 * res_l2 + sum / mesh.h_x;
    return sqrt (sum);
}
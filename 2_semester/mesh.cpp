#include "mesh.h"
#include <cmath>

static bool is_inside_domain12(double x, double y) {
    // Область 11: пять квадратов
    // Ω00: x∈[0,1], y∈[0,1]
    // Ω01: x∈[0,1], y∈[1,2]
    // Ω11: x∈[1,2], y∈[1,2]
    // Ω20: x∈[2,3], y∈[0,1]
    // Ω21: x∈[2,3], y∈[1,2]
    if (x < 0.0 || x > 3.0 || y < 0.0 || y > 3.0) return false;

    if ((x >= 0.0 && x <= 1.0 && y >= 0.0 && y <= 1.0) || // Ω00
        (x >= 0.0 && x <= 1.0 && y >= 1.0 && y <= 2.0) || // Ω01
        (x >= 1.0 && x <= 2.0 && y >= 1.0 && y <= 2.0) || // Ω11
        (x >= 2.0 && x <= 3.0 && y >= 0.0 && y <= 1.0) || // Ω20
        (x >= 2.0 && x <= 3.0 && y >= 1.0 && y <= 2.0))   // Ω21
        return true;
    return false;
}

void fill_mesh_domain12(std::vector<Point>& mesh_points, int N, int M,
                        double X, double Y,
                        double /*w*/, double /*rho_inlet*/) {
    double hx = X / (N - 1);
    double hy = Y / (M - 1);

    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < N; ++i) {
            double x = i * hx;
            double y = j * hy;
            Point pt(x, y, i, j, block_status::outer);

            if (!is_inside_domain12(x, y)) {
                pt.status = block_status::outer;
                mesh_points[j * N + i] = pt;
                continue;
            }

            bool left_out   = (i == 0)   || !is_inside_domain12(x - hx, y);
            bool right_out  = (i == N-1) || !is_inside_domain12(x + hx, y);
            bool bottom_out = (j == 0)   || !is_inside_domain12(x, y - hy);
            bool top_out    = (j == M-1) || !is_inside_domain12(x, y + hy);

            bool is_boundary = left_out || right_out || bottom_out || top_out;
            if (!is_boundary) {
                pt.status = block_status::inner;
                mesh_points[j * N + i] = pt;
                continue;
            }

            // Вход: левая граница при x = 0, y ∈ [0,1] (Ω00)
            bool inlet_cond = (fabs(x) < 1e-12) && (y >= 0.0 - 1e-12 && y <= 1.0 + 1e-12);
            // Выход снизу: y = 0, x ∈ [2,3] (Ω20)
            bool outlet_bottom_cond = (fabs(y) < 1e-12) && (x >= 2.0 - 1e-12 && x <= 3.0 + 1e-12);
            // Выход сверху в области 11 отсутствует → не задаём

            if (inlet_cond && left_out) {
                pt.status = block_status::inlet_left;
            } else if (outlet_bottom_cond && bottom_out) {
                pt.status = block_status::outlet_bottom;
            } else {
                // Все остальные границы – твёрдая стенка
                if (left_out && !right_out) pt.status = block_status::x1c_r;
                else if (right_out && !left_out) pt.status = block_status::x1c_l;
                else if (bottom_out && !top_out) pt.status = block_status::x2c_u;
                else if (top_out && !bottom_out) pt.status = block_status::x2c_d;
                else {
                    // Углы
                    if (left_out && bottom_out) pt.status = block_status::corner_l2u;
                    else if (left_out && top_out) pt.status = block_status::corner_u2r;
                    else if (right_out && bottom_out) pt.status = block_status::corner_d2l;
                    else if (right_out && top_out) pt.status = block_status::corner_r2d;
                    else pt.status = block_status::error;
                }
            }
            mesh_points[j * N + i] = pt;
        }
    }
}

int Mesh::get_point_neigb_glob_num(Point point, int dir) {
    int i = point.i;
    int j = point.j;
    if (dir == m_00) return j * N + i;
    if (dir == m_R0) return j * N + i + 1;
    if (dir == m_L0) return j * N + i - 1;
    if (dir == m_0R) return (j + 1) * N + i;
    if (dir == m_0L) return (j - 1) * N + i;
    if (dir == m_RR) return (j + 1) * N + i + 1;
    if (dir == m_LR) return (j + 1) * N + i - 1;
    if (dir == m_RL) return (j - 1) * N + i + 1;
    if (dir == m_LL) return (j - 1) * N + i - 1;
    if (dir == m_0RR) return (j + 2) * N + i;
    if (dir == m_0RRR) return (j + 3) * N + i;
    if (dir == m_0LL) return (j - 2) * N + i;
    if (dir == m_0LLL) return (j - 3) * N + i;
    if (dir == m_RR0) return j * N + i + 2;
    if (dir == m_RRR0) return j * N + i + 3;
    if (dir == m_LL0) return j * N + i - 2;
    if (dir == m_LLL0) return j * N + i - 3;
    return -1;
}

Point Mesh::get_neighbour(Point point, int dir) {
    int dest_num = get_point_neigb_glob_num(point, dir);
    if (dest_num < 0 || dest_num >= Dim)
        return Point();
    return mesh_points[dest_num];
}

void Mesh::print_mesh() {
    for (int j = M - 1; j >= 0; --j) {
        for (int i = 0; i < N; ++i) {
            Point point = get_point(i, j);
            switch (point.status) {
                case block_status::inner:       printf("x"); break;
                case block_status::x1c_r:
                case block_status::x1c_l:       printf("|"); break;
                case block_status::x2c_u:
                case block_status::x2c_d:       printf("_"); break;
                case block_status::inlet_left:  printf("I"); break;
                case block_status::outlet_bottom:
                case block_status::outlet_top:  printf("O"); break;
                case block_status::outer:       printf(" "); break;
                default:                        printf("?"); break;
            }
        }
        printf("\n");
    }
}
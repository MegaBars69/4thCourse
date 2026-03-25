#pragma once
#include <vector>
#include <cstdio>

enum class block_status {inner, x1c_r, x1c_l, x2c_u, x2c_d, corner_l2u, corner_u2r, corner_r2d, corner_d2l, corner_l2d, outer, error};
#define m_00 0
#define m_R0 1
#define m_L0 2
#define m_0R 3
#define m_0L 4
#define m_RR 5
#define m_LR 6
#define m_RL 7
#define m_LL 8
#define m_0RR 9
#define m_0RRR 10
#define m_0LL 11
#define m_0LLL 12
#define m_RR0 13
#define m_RRR0 14
#define m_LL0 15
#define m_LLL0 16

class Point
{
  public :
    double x{}, y{};
    int i{},j{};
    block_status status;

  Point (double x, double y, int i, int j, block_status status) :
    x(x), y(y), i(i), j(j), status(status)
  {}

  Point ()
  {}
};

void fill_mesh_blocks_vector (std::vector<Point> & mesh_points, int N, int M, double X, double Y);

class Mesh
{
public :
  double X{}, Y{}, T {};
  int N {}, M {}, Dim{};
  int T_segm;

  std::vector<Point> mesh_points {};
  double h_x, h_y, tau;


public :
Mesh (double X, double Y, double T, int N, int M, int T_segm) :
X(X), Y(Y), T(T), N(N), M(M), T_segm(T_segm)
{
  Dim = N*M;
  mesh_points.resize (Dim);
  h_x = X / (N - 1);
  h_y = Y / (M - 1);
  tau = T / (T_segm);
  fill_mesh_blocks_vector (mesh_points, N, M, X, Y);
}

Point get_neighbour (Point point, int dir);
Point get_point (int i, int j)
{
  return mesh_points[j*N + i];
}

int get_point_neigb_glob_num (Point point, int dir);


void print_mesh();

};

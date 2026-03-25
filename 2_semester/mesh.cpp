#include "mesh.h"

/*
block_status get_block_status (int i, int j, int bound_x, int bound_y, int N, int M)
{
  if (i == 0 && j == 0)
    return block_status::corner_l2u;

  if (i == 0 && j == M-1)
    return block_status::corner_u2r;

  if (j == 0 && i == bound_x)
    return block_status::corner_d2l;

  if (j == bound_y && i == bound_x)
    return block_status::corner_l2d;

  if (i == N - 1  && j == bound_y)
    return block_status::corner_d2l;

  if (i == N - 1 && j == M - 1)
    return block_status::corner_r2d;

  if (i == 0)
    return block_status::x1c_r;

  if (j == M - 1)
    return block_status::x2c_d;

  if (j==0 && i < bound_x)
    return block_status::x2c_u;

  if (i==bound_x && j < bound_y)
    return block_status::x1c_l;

  if (i > bound_x && j == bound_y)
    return block_status::x2c_u;

  if (i == N - 1 && j > bound_y)
    return block_status::x1c_l;

  if (i <= bound_x)
    return block_status::inner;

  if (i > bound_x)
    return j > bound_y ? block_status::inner : block_status::outer;

  printf ("INCORRECT CASE %d %d\n", i, j);
  return block_status::error;
}

void fill_mesh_blocks_vector (std::vector<Point> & mesh_points, int N, int M, double X, double Y)
{

  double square_size = X / 3;
  double h_x = X / (N - 1);
  double h_y = Y / (M - 1);
  int square_x = square_size / h_x;
  int square_y = square_size / h_y;
  printf ("square_x = %d sq_y = %d, x = %lf, y = %lf\n", square_x, square_y, square_x * h_x, square_y*h_y);

  for (int j = 0; j < M; j++)
  {
    for (int i = 0; i < N; i++)
    {
      block_status status = get_block_status (i, j, square_x, square_y, N, M);
      double x = i * h_x;
      double y = j * h_y;
      Point point (x, y, i, j, status);
      mesh_points [j * N + i] = point;
    }
  }
}*/

block_status get_block_status (int i, int j, int N, int M)
{
    // Углы
    if (i == 0 && j == 0)
        return block_status::corner_l2u;
    if (i == 0 && j == M-1)
        return block_status::corner_u2r;
    if (i == N-1 && j == 0)
        return block_status::corner_d2l;
    if (i == N-1 && j == M-1)
        return block_status::corner_r2d;

    // Грани (без углов)
    if (i == 0)
        return block_status::x1c_r;
    if (i == N-1)
        return block_status::x1c_l;
    if (j == 0)
        return block_status::x2c_u;
    if (j == M-1)
        return block_status::x2c_d;

    // Все остальные – внутренние
    return block_status::inner;
}

void fill_mesh_blocks_vector (std::vector<Point> & mesh_points, int N, int M, double X, double Y)
{
    double h_x = X / (N - 1);
    double h_y = Y / (M - 1);

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            block_status status = get_block_status (i, j, N, M);
            double x = i * h_x;
            double y = j * h_y;
            Point point (x, y, i, j, status);
            mesh_points [j * N + i] = point;
        }
    }
}

int Mesh::get_point_neigb_glob_num (Point point, int dir)
{
  int i = point.i;
  int j = point.j;

  if (dir == m_00)
    return j * N + i;

  if (dir == m_R0)
      return j * N + i + 1;

  if (dir == m_L0)
      return j * N + i - 1;

  if (dir == m_0R)
      return (j + 1) * N + i;

  if (dir == m_0L)
      return (j - 1) * N + i;

  if (dir == m_RR)
      return (j + 1) * N + i + 1;

  if (dir == m_LR)
      return (j + 1) * N + i - 1;

  if (dir == m_RL)
      return (j - 1) * N + i + 1;

  if (dir == m_LL)
      return (j - 1) * N + i - 1;

  if (dir == m_0RR)
      return (j+2) * N + i;

  if (dir == m_0RRR)
      return (j+3) * N + i;

  if (dir == m_0LL)
      return (j - 2) * N + i;

  if (dir == m_0LLL)
      return (j - 3) * N + i;

  if (dir == m_RR0)
      return (j) * N + i + 2;

  if (dir == m_RRR0)
      return (j) * N + i + 3;

  if (dir == m_LL0)
      return (j) * N + i - 2;

  if (dir == m_LLL0)
      return (j) * N + i - 3;

  return -1;
}

Point Mesh::get_neighbour (Point point, int dir)
{
    int dest_num = get_point_neigb_glob_num(point, dir);
    if (dest_num < 0 || dest_num >= Dim)
        return Point ();
    return mesh_points[dest_num];
}


void Mesh::print_mesh ()
{
  for (int j = M - 1; j >= 0; j--)
  {
    for (int i = 0; i < N; i++)
    {
      Point point = get_point (i, j);
      if (point.status == block_status::inner)
        {
          printf("x");
        }
          if (point.status == block_status::x1c_l ||
              point.status == block_status::x1c_r)
              printf("|" );
          if (point.status == block_status::x2c_u ||
              point.status == block_status::x2c_d )
              printf("_");
          if (point.status == block_status::corner_d2l)
            printf("_|");
          if (point.status == block_status::corner_l2d)
            printf("|-");
          if (point.status == block_status::corner_l2u)
            printf("|_");
          if (point.status == block_status::corner_r2d)
            printf("-|");
          if (point.status == block_status::corner_u2r)
            printf("|--");
          if (point.status == block_status::error)
            printf("error!");
          if  (point.status == block_status::outer)
          printf("0");

    }
    printf("\n");
  }
}

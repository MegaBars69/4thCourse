#pragma once
#include <cmath>
#include "P_scheme.h"
#include "P_gas.h"
#include "mesh.h"
#include <vector>
#include "funcs.h"
#include "Eigen/Sparse"

#define g 1
#define v1 2
#define v2 3

using eigen_matrix_t = Eigen::SparseMatrix<double>;
using eigen_triplet_t = Eigen::Triplet<double>;
using eigen_vector_t = Eigen::VectorXd;
//using eigen_precond_t = Eigen::LeastSquareDiagonalPreconditioner<double>;
using eigen_solver_t = Eigen::BiCGSTAB<eigen_matrix_t>;

class func_point
{
public:
  int i, j, global_num;
  double G, V1, V2, H;
  double GV1, GV2;
  double HV1, HV2;

  func_point (int i, int j, int N, double G, double V1, double V2) :
  i(i), j(j), G(G), V1(V1), V2(V2)
  {
      global_num = j * N + i;

      GV1 = G * V1;
      GV2 = G * V2;
      H = exp (G);
      HV1 = H * V1;
      HV2 = H*V2;
  }

  func_point (int glob_num, int N, double G, double V1, double V2) :
  i(glob_num % N), j (glob_num / N), G(G), V1 (V1), V2 (V2)
  {
      global_num =  glob_num;
      GV1 = G * V1;
      GV2 = G * V2;
      H = exp (G);
      HV1 = H * V1;
      HV2 = H*V2;
  }

    func_point()
    {
    }

    void update_G (double new_G)
    {
      G = new_G;
      GV1 = G*V1;
      GV2 = G*V2;
      H = exp (G);
      HV1 = H * V1;
      HV2 = H*V2;
    }

    void update_V1 (double new_V1)
    {
        V1 = new_V1;
        GV1 = G*V1;
        HV1 = H * V1;
    }

    void update_V2 (double new_V2)
    {
        V2 = new_V2;
        GV2 = G*V2;
        HV2 = H*V2;
    }
};

class Matrix
{
public:

  P_gas gas;
  Mesh mesh;

  int Dim;
  int step = 0;

  std::vector<func_point> func_points;

  eigen_matrix_t matrix;
  eigen_vector_t rhs;
  std::vector<eigen_triplet_t> triplets;

  std::vector<double> solution_G;
  std::vector<double> solution_V1;
  std::vector<double> solution_V2;

  func_point get_func_point_check (Point point, int dir);
  func_point init_func_point (Point point);

  double get_mu ();
  void init_matrix_G(std::vector<eigen_triplet_t> & triplets, eigen_vector_t & rhs);
  void fill_matrix_G (Point point, std::vector<eigen_triplet_t> & triplets, eigen_vector_t & rhs);

  void init_matrices_V(std::vector<eigen_triplet_t> & triplets1, eigen_vector_t & rhs1,
                       std::vector<eigen_triplet_t> & triplets2, eigen_vector_t & rhs2);
  void fill_matrix_V1 (Point point, std::vector<eigen_triplet_t> & triplets, eigen_vector_t & rhs);
  void fill_matrix_V2 (Point point, std::vector<eigen_triplet_t> & triplets, eigen_vector_t & rhs);

  int init_and_solve_G();
  int init_and_solve_V();

  void update_func_points (int mode);

  double calc_res_C1 (int mode);
  double calc_res_L2 (int mode);
  double calc_res_W1 (int mode);

public:
    Matrix (P_gas gas,  Mesh mesh) :
    gas(gas), mesh (mesh)
    {
      Dim = mesh.Dim;

      solution_G.resize(Dim);
      solution_V1.resize(Dim);
      solution_V2.resize(Dim);

      func_points.resize(Dim);

      for (int i = 0; i < Dim; i++)
        {
          Point point = mesh.mesh_points[i];
          func_point f_point = init_func_point (point);
          func_points[i] = f_point;

        }
      step = 0;
    }
};

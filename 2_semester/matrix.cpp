#include <iostream>
#include "matrix.h"
#include "bicstab.h"


func_point Matrix::init_func_point (Point point)
{
  double x = point.x;
  double y = point.y;
  double G{}, V1{}, V2{};

  if (point.status == block_status::inner)
  {
      double H = rho (0, x, y);
      G = log(rho (0, x, y));
      V1 = u1 (0, x, y);
      V2 = u2  (0, x, y);
  }
  else if (point.status == block_status::outer)
  {G = 0; V1 = 0; V2 = 0; }
  else
  {
      double H = rho (0, x, y);
      G = log(rho (0, x, y));
      V1 = 0;
      V2 = 0;
  }

    return func_point(point.i, point.j, mesh.N, G, V1, V2);

}

func_point Matrix::get_func_point_check (Point point, int dir)
{
  int global_num = mesh.get_point_neigb_glob_num(point, dir);
  if (global_num >=0 && global_num < Dim)
      return func_points[global_num];
  return func_point ();
}

void Matrix::fill_matrix_G(Point point, std::vector<eigen_triplet_t> & triplets, eigen_vector_t & rhs)
{
  func_point m00, mR0, mL0, m0R, m0L;

  m00 = get_func_point_check(point, m_00);
  mR0 = get_func_point_check(point, m_R0);
  mL0 = get_func_point_check(point, m_L0);
  m0R = get_func_point_check(point, m_0R);
  m0L = get_func_point_check(point, m_0L);


  double tau = mesh.tau;
  double h_x = mesh.h_x, h_y = mesh.h_y;
  int matr_str_num = m00.global_num;


  // if corner 270 - we consider it like inner point with one neigb that has null values
  if (point.status == block_status::inner )
  {
      // global_num = номер столбца, нужен номер элемента, который равен глоб.номеру точки point * Dim + номер столбца

      triplets.emplace_back(matr_str_num, m00.global_num, 4);

      if (fabs((tau / h_x) * (mR0.V1 + m00.V1)) > 1e-10)
        triplets.emplace_back(matr_str_num, mR0.global_num,(tau / h_x) * (mR0.V1 + m00.V1));

      if (fabs(-(tau / h_x) * (mL0.V1 + m00.V1)) > 1e-10)
        triplets.emplace_back(matr_str_num, mL0.global_num, -(tau / h_x) * (mL0.V1 + m00.V1));

      if (fabs((tau / h_y) * (m0R.V2 + m00.V2)) > 1e-10)
        triplets.emplace_back(matr_str_num, m0R.global_num, (tau / h_y) * (m0R.V2 + m00.V2));

      if (fabs(-(tau / h_y) * (m0L.V2 + m00.V2)) > 1e-10)
        triplets.emplace_back(matr_str_num, m0L.global_num, -(tau / h_y) * (m0L.V2 + m00.V2));

      rhs[m00.global_num] =  4 * m00.G - (tau / h_x) * (2 - m00.G) * (mR0.V1 - mL0.V1)
                                  - (tau / h_y) * (2 - m00.G) * (m0R.V2 - m0L.V2)
                                  + 4 * tau * Func_0(step * tau, point.x, point.y);
  }

  //outer point - all values should be zero
  if (point.status == block_status::outer)
  {
//      matrix_G[matr_str_num * Dim + m00.global_num] = 1;
      triplets.emplace_back(matr_str_num, m00.global_num, 1);
      rhs[matr_str_num] = 0;
  }

  if (point.status == block_status::x1c_r || point.status == block_status::corner_l2u || point.status == block_status::corner_u2r)
  {
      func_point mRR0 = get_func_point_check(point, m_RR0);
      func_point mRRR0 = get_func_point_check(point, m_RRR0);

      if (fabs(2 - (tau / h_x) * (m00.V1)) > 1e-10)
        triplets.emplace_back(matr_str_num, m00.global_num, 2 - (tau / h_x) * (m00.V1));

      if (fabs((tau / h_x) * (mR0.V1)) > 1e-10)
        triplets.emplace_back(matr_str_num, mR0.global_num, (tau / h_x) * (mR0.V1));

      rhs[matr_str_num] = 2 * m00.G - (tau/h_x) * (2 - m00.G) * (mR0.V1 - m00.V1)
                                   + (tau / h_x) * (-0.5 * mRRR0.GV1 + 2 * mRR0.GV1 - 2.5 * mR0.GV1 + m00.GV1)
                                   + (tau / h_x) * (2 - m00.G) *(-0.5 * mRRR0.V1 + 2 * mRR0.V1 - 2.5 * mR0.V1 + m00.V1)
                                   + 2 * tau * Func_0(step * tau, point.x, point.y);

  }

  if (point.status == block_status::x2c_u)
  {
      func_point m0RR = get_func_point_check(point, m_0RR);
      func_point m0RRR = get_func_point_check(point, m_0RRR);

      if (fabs(2 -  (tau / h_y) * m00.V2) > 1e-10)
        triplets.emplace_back(matr_str_num, m00.global_num, 2 -  (tau / h_y) * m00.V2);

      if (fabs((tau / h_y) * (m0R.V2)) > 1e-10)
        triplets.emplace_back(matr_str_num, m0R.global_num, (tau / h_y) * (m0R.V2));

      rhs[matr_str_num] = 2 * m00.G - (tau / h_y) * (2 - m00.G) * (m0R.V2 - m00.V2)
                                   + (tau / h_y) * (-0.5 * m0RRR.GV2 + 2*m0RR.GV2 - 2.5*m0R.GV2 + m00.GV2)
                                   + (tau / h_y) * (2 - m00.G) * (-0.5 * m0RRR.V2 + 2*m0RR.V2 - 2.5*m0R.V2 + m00.V2)
                                   + 2 * tau * Func_0(step * tau, point.x, point.y);
  }

  if (point.status == block_status::x1c_l || point.status == block_status::corner_d2l || point.status == block_status::corner_r2d || point.status == block_status::corner_l2d)
  {
      func_point mLL0 = get_func_point_check(point, m_LL0);
      func_point mLLL0 = get_func_point_check(point, m_LLL0);

      if (fabs(2 + (tau / h_x) * m00.V1) > 1e-10)
        triplets.emplace_back(matr_str_num, m00.global_num, 2 + (tau / h_x) * m00.V1);

      if (fabs(-(tau / h_x) * (mL0.V1)) > 1e-10)
        triplets.emplace_back(matr_str_num, mL0.global_num, -(tau / h_x) * (mL0.V1));

      rhs[matr_str_num] = 2 * m00.G - (tau / h_x) * (2 - m00.G) * (m00.V1 - mL0.V1)
                                   - (tau / h_x) * (m00.GV1 - 2.5*mL0.GV1 + 2 * mLL0.GV1 - 0.5 * mLLL0.GV1)
                                   - (tau / h_x) * (2 - m00.G) * (m00.V1 - 2.5*mL0.V1 + 2 * mLL0.V1 - 0.5 * mLLL0.V1)
                                   + 2 * tau * Func_0(step * tau, point.x, point.y);/*+ 2 * tau * Func_0(t,x,y)*/
  }

  if (point.status == block_status::x2c_d)
  {
      func_point m0LL = get_func_point_check(point, m_0LL);
      func_point m0LLL = get_func_point_check(point, m_0LLL);

      if (fabs(2 + (tau / h_y) * m00.V2) > 1e-10)
        triplets.emplace_back(matr_str_num, m00.global_num, 2 + (tau / h_y) * m00.V2) ;

      if (fabs(-(tau / h_y) * (m0L.V2)) > 1e-10)
        triplets.emplace_back(matr_str_num, m0L.global_num,  -(tau / h_y) * (m0L.V2));

      rhs[matr_str_num] = 2 * m00.G - (tau / h_y) * (2 - m00.G) * (m00.V2 - m0L.V2)
                                   - (tau / h_y) * (m00.GV2 - 2.5 * m0L.GV2 + 2 * m0LL.GV2 - 0.5 * m0LLL.GV2)
                                   - (tau / h_y) * (2 - m00.G) * (m00.V2 - 2.5 * m0L.V2 + 2 * m0LL.V2 - 0.5 * m0LLL.V2)
                                   + 2 * tau * Func_0(step * tau, point.x, point.y); /*+ 2 * tau * Func_0(t,x,y)*/
  }

}


// каждый узел дает одно уравнение в системе, то есть матрица будет Dim * Dim в каждой строчке есть Dim столбцов.
void Matrix::init_matrix_G(std::vector<eigen_triplet_t> & triplets, eigen_vector_t & rhs)
{
    //вообще надо бы делать так, чтобы матрица поступала пустая сюда
    // Проход по строкам матрицы.

    for (int i = 0; i < Dim; i++)
    {
        Point point = mesh.mesh_points[i];
        fill_matrix_G ( point, triplets, rhs);
    }
}


void Matrix::fill_matrix_V1 (Point point, std::vector<eigen_triplet_t> & triplets, eigen_vector_t & rhs)
{
    func_point m00  = get_func_point_check (point, m_00);
    int matr_str_num = m00.global_num;

    if (point.status != block_status::inner)
    {
        triplets.emplace_back(matr_str_num, m00.global_num, 1);
        rhs[matr_str_num] = 0;
    }

    if (point.status == block_status::inner)
    {
        func_point mR0, mL0, m0R, m0L, mRR, mLR, mRL, mLL;
        double tau = mesh.tau, mu = gas.mu;
        double h_x = mesh.h_x, h_y = mesh.h_y;

        mR0 = get_func_point_check(point, m_R0);
        mL0 = get_func_point_check(point, m_L0);
        m0R = get_func_point_check(point, m_0R);
        m0L = get_func_point_check(point, m_0L);
        mRR = get_func_point_check(point, m_RR);
        mLR = get_func_point_check(point, m_LR);
        mRL = get_func_point_check(point, m_RL);
        mLL = get_func_point_check(point, m_LL);

        if (fabs(m00.H / tau + (8 * mu) /(3 * h_x * h_x) + (2 * mu) /(h_y * h_y)) > 1e-10)
          triplets.emplace_back(matr_str_num, m00.global_num, m00.H / tau + (8 * mu) /(3 * h_x * h_x) + (2 * mu) /(h_y * h_y));

        if (fabs((1 / (6 * h_x)) * ( m00.HV1 + mR0.HV1) - (4 * mu) / (3 * h_x * h_x)) > 1e-10)
          triplets.emplace_back(matr_str_num, mR0.global_num, (1 / (6 * h_x)) * ( m00.HV1 + mR0.HV1) - (4 * mu) / (3 * h_x * h_x));

        if (fabs(-(1 / (6 * h_x)) * ( m00.HV1 + mL0.HV1) - (4 * mu) / (3 * h_x * h_x)) > 1e-10)
          triplets.emplace_back(matr_str_num, mL0.global_num, -(1 / (6 * h_x)) * ( m00.HV1 + mL0.HV1) - (4 * mu) / (3 * h_x * h_x));

        if (fabs((1 / (4 * h_y)) * (m00.HV2 + m0R.HV2) - mu / (h_y * h_y)) > 1e-10)
          triplets.emplace_back(matr_str_num, m0R.global_num, (1 / (4 * h_y)) * (m00.HV2 + m0R.HV2) - mu / (h_y * h_y));

        if (fabs(-(1 / (4 * h_y)) * (m00.HV2 + m0L.HV2) - mu / (h_y * h_y)) > 1e-10)
          triplets.emplace_back(matr_str_num, m0L.global_num, -(1 / (4 * h_y)) * (m00.HV2 + m0L.HV2) - mu / (h_y * h_y));

        rhs[matr_str_num] = m00.HV1 / tau + (1 / (6 * h_x)) * (m00.V1 * m00.V1) * (mR0.H - mL0.H)
                              + (1 / (4 * h_y)) * (m00.V1 * (m0R.HV2 - m0L.HV2)) - (gas.P(mR0.H) - gas.P(mL0.H)) * (1 / (2 * h_x))
                              + (mu / (12 * h_x * h_y)) * (mRR.V2 - mRL.V2 - mLR.V2 + mLL.V2)
                              + m00.H * Func_1(step * tau, point.x, point.y, gas);  // + m00.H * f1 (point) ;



    }

}

void Matrix::fill_matrix_V2 (Point point, std::vector<eigen_triplet_t> & triplets, eigen_vector_t & rhs)
{
    func_point m00  = get_func_point_check (point, m_00);
    int matr_str_num = m00.global_num;

    if (point.status != block_status::inner)
    {
        triplets.emplace_back(matr_str_num, m00.global_num, 1);
        rhs[matr_str_num] = 0;
    }

    if (point.status == block_status::inner)
    {
        func_point mR0, mL0, m0R, m0L, mRR, mLR, mRL, mLL;
        double tau = mesh.tau, mu = gas.mu;
        double h_x = mesh.h_x, h_y = mesh.h_y;

        mR0 = get_func_point_check(point, m_R0);
        mL0 = get_func_point_check(point, m_L0);
        m0R = get_func_point_check(point, m_0R);
        m0L = get_func_point_check(point, m_0L);
        mRR = get_func_point_check(point, m_RR);
        mLR = get_func_point_check(point, m_LR);
        mRL = get_func_point_check(point, m_RL);
        mLL = get_func_point_check(point, m_LL);

        if (fabs(m00.H / tau + (8 * mu) /(3 * h_y * h_y) + (2 * mu) /(h_x * h_x)) > 1e-10)
          triplets.emplace_back(matr_str_num, m00.global_num, m00.H / tau + (8 * mu) /(3 * h_y * h_y) + (2 * mu) /(h_x * h_x));

        if (fabs((1 / (6 * h_y)) * ( m00.HV2 + m0R.HV2) - (4 * mu) / (3 * h_y * h_y)) > 1e-10)
          triplets.emplace_back(matr_str_num, m0R.global_num, (1 / (6 * h_y)) * ( m00.HV2 + m0R.HV2) - (4 * mu) / (3 * h_y * h_y));

        if (fabs(-(1 / (6 * h_y)) * ( m00.HV2 + m0L.HV2) - (4 * mu) / (3 * h_y * h_y)) > 1e-10)
          triplets.emplace_back(matr_str_num, m0L.global_num, -(1 / (6 * h_y)) * ( m00.HV2 + m0L.HV2) - (4 * mu) / (3 * h_y * h_y));

        if (fabs((1 / (4 * h_x)) * (m00.HV1 + mR0.HV1) - mu / (h_x * h_x)) > 1e-10)
          triplets.emplace_back(matr_str_num, mR0.global_num, (1 / (4 * h_x)) * (m00.HV1 + mR0.HV1) - mu / (h_x * h_x));

        if (fabs(-(1 / (4 * h_x)) * (m00.HV1 + mL0.HV1) - mu / (h_x * h_x)) > 1e-10)
          triplets.emplace_back(matr_str_num, mL0.global_num, -(1 / (4 * h_x)) * (m00.HV1 + mL0.HV1) - mu / (h_x * h_x));

        rhs[matr_str_num] = m00.HV2 / tau + (1 / (6 * h_y)) * (m00.V2 * m00.V2) * (m0R.H - m0L.H)
                                 + (1 / (4 * h_x)) * (m00.V2 * (mR0.HV1 - mL0.HV1)) - (gas.P(m0R.H) - gas.P(m0L.H)) * (1 / (2 * h_y))
                                 + (mu / (12 * h_x * h_y)) * (mRR.V1 - mRL.V1 - mLR.V1 + mLL.V1)
                                 + m00.H * Func_2(step * tau, point.x, point.y, gas);  // + m00.H * f2 (point) ;

    }

}



void Matrix::init_matrices_V(std::vector<eigen_triplet_t> & triplets1, eigen_vector_t & rhs1,
                             std::vector<eigen_triplet_t> & triplets2, eigen_vector_t & rhs2)
{
    for (int i = 0; i < Dim; i++)
    {
        Point point = mesh.mesh_points[i];
        fill_matrix_V1 (point, triplets1,  rhs1 );
        fill_matrix_V2 (point, triplets2, rhs2);
    }
}

int Matrix::init_and_solve_G()
{
    eigen_matrix_t matrix (Dim, Dim);
    std::vector<eigen_triplet_t> triplets;
    eigen_vector_t rhs (Dim);
    eigen_vector_t guess (Dim);

    for (int i = 0; i < Dim; i++)
        guess[i] = func_points[i].G;

    init_matrix_G(triplets, rhs);
    matrix.setFromTriplets(triplets.begin(), triplets.end());

    eigen_solver_t solver;
    solver.compute(matrix);
    eigen_vector_t res = solver.solveWithGuess(rhs, guess);

    for (int i = 0; i < Dim; i++)
      solution_G[i] = res[i];

//    if (step % 6 == 0)
//    {
//        for (int i = 0; i < Dim; i++)
//        {
//            auto point = mesh.mesh_points[i];
//            solution_G[i] = log(rho(step*mesh.tau, point.x, point.y ));
//        }
//    }

    update_func_points(g);
    return 0;
}

int Matrix::init_and_solve_V ()
{
    eigen_matrix_t matrix1 (Dim, Dim);
    std::vector<eigen_triplet_t> triplets1;
    eigen_vector_t rhs1 (Dim);

    eigen_matrix_t matrix2 (Dim, Dim);
    std::vector<eigen_triplet_t> triplets2;
    eigen_vector_t rhs2 (Dim);

    eigen_vector_t guess1 (Dim);
    eigen_vector_t guess2 (Dim);

    for (int i = 0; i < Dim; i++)
    {
        guess1[i] = func_points[i].V1;
        guess2[i] = func_points[i].V2;
    }

    init_matrices_V(triplets1, rhs1, triplets2, rhs2 );
    matrix1.setFromTriplets(triplets1.begin(), triplets1.end());
    matrix2.setFromTriplets(triplets2.begin(), triplets2.end());

    eigen_solver_t solver;
    solver.compute(matrix1);
    eigen_vector_t res1 = solver.solveWithGuess(rhs1, guess1);

    solver.compute(matrix2);
    eigen_vector_t res2 = solver.solveWithGuess(rhs2, guess2);

    for (int i = 0; i < Dim; i++)
    {
        solution_V1[i] = res1[i];
        solution_V2[i] = res2[i];
    }

//    if (step % 6 == 0)
//    {
//        for (int i = 0; i < Dim; i++) {
//            Point point = mesh.mesh_points[i];
//            solution_V1[i] = u1(step * mesh.tau, point.x, point.y);
//            solution_V2[i] = u2(step * mesh.tau, point.x, point.y);
//        }
//    }

    update_func_points(v1);
    update_func_points(v2);
    return 0;
}

/*
int Matrix::init_and_solve_G()
{
    eigen_matrix_t matrix(Dim, Dim);
    std::vector<eigen_triplet_t> triplets;
    eigen_vector_t rhs(Dim);
    eigen_vector_t guess(Dim);

    for (int i = 0; i < Dim; i++)
        guess[i] = func_points[i].G;

    init_matrix_G(triplets, rhs);
    matrix.setFromTriplets(triplets.begin(), triplets.end());

    BiCGSTABSolver solver(matrix);
    eigen_vector_t res = guess;   // начальное приближение
    bool converged = solver.solve(res, rhs, 1000, 1e-10);
    if (!converged) {
        std::cerr << "BiCGSTAB не сошёлся для G на шаге " << step << std::endl;
        return -1;
    }

    for (int i = 0; i < Dim; i++)
        solution_G[i] = res[i];

    update_func_points(g);
    return 0;
}

int Matrix::init_and_solve_V()
{
    eigen_matrix_t matrix1(Dim, Dim);
    std::vector<eigen_triplet_t> triplets1;
    eigen_vector_t rhs1(Dim);

    eigen_matrix_t matrix2(Dim, Dim);
    std::vector<eigen_triplet_t> triplets2;
    eigen_vector_t rhs2(Dim);

    eigen_vector_t guess1(Dim);
    eigen_vector_t guess2(Dim);

    for (int i = 0; i < Dim; i++) {
        guess1[i] = func_points[i].V1;
        guess2[i] = func_points[i].V2;
    }

    init_matrices_V(triplets1, rhs1, triplets2, rhs2);
    matrix1.setFromTriplets(triplets1.begin(), triplets1.end());
    matrix2.setFromTriplets(triplets2.begin(), triplets2.end());

    // Решение для V1
    BiCGSTABSolver solver1(matrix1);
    eigen_vector_t res1 = guess1;
    bool conv1 = solver1.solve(res1, rhs1, 1000, 1e-10);
    if (!conv1) {
        std::cerr << "BiCGSTAB не сошёлся для V1 на шаге " << step << std::endl;
        return -1;
    }

    // Решение для V2
    BiCGSTABSolver solver2(matrix2);
    eigen_vector_t res2 = guess2;
    bool conv2 = solver2.solve(res2, rhs2, 1000, 1e-10);
    if (!conv2) {
        std::cerr << "BiCGSTAB не сошёлся для V2 на шаге " << step << std::endl;
        return -1;
    }

    for (int i = 0; i < Dim; i++) {
        solution_V1[i] = res1[i];
        solution_V2[i] = res2[i];
    }

    update_func_points(v1);
    update_func_points(v2);
    return 0;
}
*/

void Matrix::update_func_points(int mode)
{
    for (int i = 0; i < Dim; i++)
    {
        if (mode == g)
          func_points[i].update_G(solution_G[i]);

        if (mode == v1)
          func_points[i].update_V1(solution_V1[i]);

        if (mode == v2)
          func_points[i].update_V2(solution_V2[i]);
    }
}

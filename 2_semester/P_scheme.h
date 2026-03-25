#pragma once
// параметры сетки
class P_scheme
{
  public :
  int M_x; //число разбиений по оси Х
  int N; // по оси Т
  int Dim; // общее число узлов по оси Х (М_х + 1)
  double h_x; // шаг по оси Х = Segm_X / M_x
  double tau; // шаг по оси Т = Segm_T / N
  double eta; //  параметр пропорциональности искуственной вязкости

  P_scheme (int M_x, int N, double h_x, double tau) :
  M_x(M_x), N(N), h_x(h_x), tau(tau)
  {
    Dim = M_x + 1;
    eta = 0;
  }

  double get_point_x_by_i (int i)
  {
    return h_x * i;
  }
  double get_point_y_by_i (int j)
  {
    return tau * j;
  }
};

#include <iostream>
#include "matrix.h"

int main(int argc, char const *argv[])
{
  double mu;
  double C_rho;
  double tau;
  double h_x;
  double X = 1., Y = 1., T = 1.;
  int X_segm;
  int Y_segm;
  int T_steps;
  int mode;
  if (argc != 10 )
  {
      printf ("NULL\n");
      return -1;
  }

  if ((sscanf (argv[1], "%lf", &mu) != 1) ||
      (sscanf (argv[2], "%lf", &C_rho) != 1) ||
      (sscanf (argv[3], "%lf", &X) != 1) ||
      (sscanf (argv[4], "%lf", &Y) != 1) ||
      (sscanf (argv[5], "%lf", &T) != 1 ) ||
      (sscanf (argv[6], "%d", &X_segm) != 1) ||
      (sscanf (argv[7], "%d", &Y_segm) != 1) ||
      (sscanf (argv[8], "%d", &T_steps) != 1) ||
      (sscanf (argv[9], "%d", &mode) != 1)
      )
     {
       printf ("Can`t read initial values!\n");
       return -2;
     }


     P_gas gas(T, X, Y, C_rho, 1.4, mu, mode);
     Mesh mesh (X, Y, T, X_segm, Y_segm, T_steps);
     Matrix matrix (gas, mesh);
//     mesh.print_mesh();

//    matrix.step = 1;
    for (int step = 1; step <= T_steps; step ++)
    {
        matrix.step = step;
        matrix.init_and_solve_G();
        matrix.init_and_solve_V();

        double gl2 = matrix.calc_res_L2(g);
        double gc = matrix.calc_res_C1(g);
        double gw = matrix.calc_res_W1(g);

        double v1l2 = matrix.calc_res_L2(v1);
        double v1c = matrix.calc_res_C1(v1);
        double v1w = matrix.calc_res_W1(v1);

        double v2l2 = matrix.calc_res_L2(v2);
        double v2c = matrix.calc_res_C1(v2);
        double v2w = matrix.calc_res_W1(v2);

        printf ("G : %.3le %.3le %.3le\n", gc, gl2, gw);
        printf ("V1 : %.3le %.3le %.3le\n", v1c, v1l2, v1w);
        printf ("V2 : %.3le %.3le %.3le\n", v2c, v2l2, v2w);
        printf("---------------------------------------\n");


    }
    printf("\n");



  return 0;
}

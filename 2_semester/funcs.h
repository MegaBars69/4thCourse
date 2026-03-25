#pragma once
#include <cmath>
#include "P_gas.h"

double rho (double t, double x, double y);

double u1 (double t, double x, double y);

double u2 (double t, double x, double y);

double g_t (double t, double x, double y);

double g_x (double t, double x, double y);

double g_y (double t, double x, double y);

double u1_t (double t, double x, double y);

double u2_t (double t, double x, double y);

double u1_x (double t, double x, double y);

double u1_y (double t, double x, double y);

double uu1_x (double t, double x, double y);

double u1_xy (double t, double x, double y);

double u1_xx (double t, double x, double y);

double u2_x (double t, double x, double y);

double u2_y (double t, double x, double y);

double uu2_y (double t, double x, double y);

double u2_xy (double t, double x, double y);

double u2_yy (double t, double x, double y);

double Func_0 (double t, double x, double y);

double Func_1 (double t, double x, double y, P_gas gas);

double Func_2 (double t, double x, double y, P_gas gas);
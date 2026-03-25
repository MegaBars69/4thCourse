#include "funcs.h"

#define X t,x,y

double rho (double t, double x, double y)
{
    return (cos (2 * M_PI * x) + 1.5)
           * (sin (2 * M_PI * y) + 1.5) *
           exp (t);
}

double u1 (double t, double x, double y)
{
    return sin (2 * M_PI * x) * sin (2 * M_PI * y) * exp (t);
}

double u2 (double t, double x, double y)
{
    return sin (2 * M_PI * x) * sin (2 * M_PI * y) * exp (-t);
}

double g_t (double t, double x, double y)
{
    return 1;
}

double g_x (double t, double x, double y)
{
    return  - (2 * M_PI * sin (2 * M_PI * x)) / (cos (2 * M_PI * x) + 1.5);
}

double g_y (double t, double x, double y)
{
    return   (2 * M_PI * cos (2 * M_PI * y)) / (sin (2 * M_PI * y) + 1.5);
}

double u1_t (double t, double x, double y)
{
    return u1 (t,x,y);
}

double u2_t (double t, double x, double y)
{
    return -u2 (t,x,y);
}

double u1_x (double t, double x, double y)
{
    return 2 * M_PI * exp (t) * cos (2 * M_PI * x) * sin (2 * M_PI * y);
}

double u1_y (double t, double x, double y)
{
    return 2 * M_PI * exp (t) * cos (2 * M_PI * y) * sin (2 * M_PI * x);
}

double uu1_x (double t, double x, double y)
{
    return 4 * M_PI * exp(2 * t) * sin(2 * M_PI * x)*cos (2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * y);
}

double u1_xy (double t, double x, double y)
{
    return 4 * M_PI * M_PI * exp (t) * cos (2 * M_PI * x) * cos (2 * M_PI * y);
}

double u1_xx (double t, double x, double y)
{
    return -4 * M_PI * M_PI * exp (t) * sin (2 * M_PI * x) * sin (2 * M_PI * y);
}

double u1_yy (double t, double x, double y)
{
    return -4 * M_PI * M_PI * exp (t) * sin (2 * M_PI * x) * sin (2 * M_PI * y);
}

double u2_x (double t, double x, double y)
{
    return 2 * M_PI * exp (-t) * cos (2 * M_PI * x) * sin (2 * M_PI * y);
}

double u2_y (double t, double x, double y)
{
    return 2 * M_PI * exp (-t) * sin (2 * M_PI * x) * cos (2 * M_PI * y);
}

double u2_xy (double t, double x, double y)
{
    return 4 * M_PI * M_PI * exp (-t) * cos (2 * M_PI * x) * cos (2 * M_PI * y);
}

double uu2_y (double t, double x, double y)
{
    return 4 * M_PI * exp(-2 * t) * sin(2 * M_PI * y)*cos (2 * M_PI * y) * sin(2 * M_PI * x) * sin(2 * M_PI * x);
}

double u2_yy (double t, double x, double y)
{
    return -4 * M_PI * M_PI * exp (-t) * sin (2 * M_PI * x) * sin (2 * M_PI * y);
}

double u2_xx (double t, double x, double y)
{
    return -4 * M_PI * M_PI * exp (-t) * sin (2 * M_PI * x) * sin (2 * M_PI * y);
}

double Func_0(double t, double x, double y)
{
    return g_t (t,x,y) + u1(t,x,y) * g_x (t,x,y) + u1_x (t,x,y) + u2(t,x,y) * g_y (t,x,y) + u2_y (t,x,y);
}


double P_gas::P (double value)
{
    if (mode == 1) // C * rho
    {
        return p_ro * value;
    }
    else
        return pow(value, p_gamma);
}

double P_gas::dP(double value)
{
    return mode ? p_ro : p_gamma * pow(value, p_gamma - 1);
}



double Func_1 (double t, double x, double y, P_gas gas)
{
  return u1_t (X) +  (u1(X)* u1_x(X) + uu1_x(X)) / 3 + u2(X) * u1_y(X) + gas.dP(rho(X))* g_x(X)
         - (gas.mu / rho(X)) * ( (4./3.)*u1_xx (X) + (u1_yy(X) + u2_xy(X) / 3));
}

double Func_2 (double t, double x, double y, P_gas gas)
{
    return u2_t(X) + (u2(X) * u2_y(X) + uu2_y(X)) / 3 + u1(X) * u2_x(X) + gas.dP(rho(X)) * g_y(X)
           - (gas.mu / rho(X)) * ((4./3.)*u2_yy(X) + (u2_xx(X) + (u1_xy(X)) / 3));
}



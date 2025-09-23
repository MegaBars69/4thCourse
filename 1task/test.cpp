#include "test.hpp"
#include <iostream>
#include <cmath>

double Cp = 10;
double mui = 0.1;
/*
double U (double t, double x) {return cos (2 * M_PI * t) * sin (4 * M_PI * x);}

double po (double t, double x) {return exp (t) * (cos (3 * M_PI * x) + 1.5);}

double func_f0 (double t, double x)
{
    double u = U (t, x);
    double ux =  4 * M_PI * cos (4 * M_PI * x) * cos (2 * M_PI * t);
    double px = -3 * M_PI * exp (t) * sin (3 * M_PI * x);
    double p = po (t,x);
    return p * (1 + ux) + u * px;
}

double func_f (double t, double x)
{
    double u = U(t,x);
    double f0 = func_f0 (t, x);
    double p = po (t, x);
    double ut = -2 * M_PI * sin (2 * M_PI * t) * sin (4 * M_PI * x);
    double ux =  4 * M_PI * cos (4 * M_PI * x) * cos (2 * M_PI * t);
    double uxx = -(16 * M_PI * M_PI) * u ;
    double Pxp = -3 * M_PI * sin (3 * M_PI * x) / (1.5 + cos (3 * M_PI * x));
    return ut + u * ux + Cp * Pxp + (f0 * u  - mui * uxx) / p;
}*/

double U (double t, double x) {return x * x - x;}

double po (double t, double x) {return 10;}

double func_f0 (double t, double x)
{
    double u = U (t, x);
    double ux =  2*x - 1;
    double px = 0;
    double pt = 0;
    double p = po (t,x);
    return p * ux + pt + u * px;
}

double func_f (double t, double x)
{
    double u = U(t,x);
    double f0 = func_f0 (t, x);
    double p = po (t, x);
    double ut = 0;
    double ux = 2*x - 1;
    double uxx = 2;
    double Pxp = 0;
    return ut + u * ux + Cp * Pxp + (f0 * u  - mui * uxx) / p;
}
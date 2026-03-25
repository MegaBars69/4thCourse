#include "test.hpp"
#include <iostream>
#include <cmath>

int Cp = 10;
double Gamma = 1.4;
double mui = 0.1;
bool liniar = true;
int num_of_func = 2;

double p (double x) { return (liniar ? Cp * x : pow (x, Gamma)); }
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
    double p = po (t, x);
    double ut = -2 * M_PI * sin (2 * M_PI * t) * sin (4 * M_PI * x);
    double ux =  4 * M_PI * cos (4 * M_PI * x) * cos (2 * M_PI * t);
    double uxx = -(16 * M_PI * M_PI) * u;
    double Pxp = -3 * M_PI * sin (3 * M_PI * x) / (1.5 + cos (3 * M_PI * x));

    if (liniar) { return ut + u * ux + Cp * Pxp - mui * uxx / p; }
    else { return ut + u * ux + Gamma * Pxp * pow (p, 0.4) - mui * uxx / p; } 
}*/

double po1(double x)
{
    return (x > 5.5 || x < 4.5 ? 1 : 2);
}

double U1( double /*x*/)
{
    return 0;
}

double po2 (double /*x*/)
{
    return 1;
}

double U2 (double x)
{
    return (x > 5.5 || x < 4.5 ? 0 : 1);
}

double U (double x)
{
    return (num_of_func == 1 ? U1 (x) : U2 (x));
}

double po (double x)
{
    return (num_of_func == 1 ? po1 (x) : po2 (x));
}

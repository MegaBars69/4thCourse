#include <iostream>
#include "initialize_matrix.hpp"
#include "progonka_solver.hpp"
#include "log.hpp"
#include <string.h>

int main(int argc, char* argv[])
{
    double T, X;
    int N, M;
    
    if (!init_args (argc, argv, T, X, N, M))
        return 1;

    return 0;
}
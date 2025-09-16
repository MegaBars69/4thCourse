#include <iostream>
#include "log.hpp"

void PrintMatrix (double *A, int n)
{
    std::cout<<"a: ";
    for (int i = 1; i < n; i++)
    {
        std::cout<<A[i]<<" ";
    }
    printf ("\n");

    std::cout<<"b: ";
    for (int i = n; i < 2*n; i++)
    {
        std::cout<<A[i]<<" ";
    }
    printf ("\n");

    std::cout<<"c: ";
    for (int i = 2*n; i < 3*n-1; i++)
    {
        std::cout<<A[i]<<" ";
    }
    printf ("\n");
}

void PrintVector (double* x, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout<<x[i]<<" ";
    }
    printf ("\n");
}

bool init_args(int argc, char* argv[], double &T, double &X, int &N, int &M)
{
    if (argc != 5)
    {
        printf("Usage: ./a.out T X N M\n");
        printf("Where:\n");
        printf("  T - maximum time value\n");
        printf("  X - maximum spatial coordinate value\n");
        printf("  N - number of time steps\n");
        printf("  M - number of spatial steps\n");
        return false;
    }
    
    try {
        T = std::stod(argv[1]);
        X = std::stod(argv[2]);
        N = std::stoi(argv[3]);
        M = std::stoi(argv[4]);
        
        // Validate input values
        if (T <= 0 || X <= 0 || N <= 0 || M <= 0) {
            printf("Error: All parameters must be positive values\n");
            return false;
        }
        
        return true;
    }
    catch (const std::invalid_argument& e) {
        printf("Error: Invalid argument format. Please provide numeric values.\n");
        return false;
    }
    catch (const std::out_of_range& e) {
        printf("Error: Argument value out of range.\n");
        return false;
    }
}

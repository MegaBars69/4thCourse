#include <iostream>
#include "log.hpp"

void PrintTriDiagMatrix (double *A, int n)
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

void PrintMatrix (double* A, int N, int M)
{
    for (int i = 0; i <= N; i++)
    {
        for (int j = 0; j <= M; j++)
        {
            std::cout<<A[i * (M + 1) + j]<<" ";
        }
        printf ("\n");
    }
}

void PrintVector (double* x, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout<<x[i]<<" ";
    }
    printf ("\n");
}

bool init_args(int argc, char* argv[], int &N, int &M)
{
    if (argc != 3)
    {
        printf("Usage: ./a.out N M\n");
        printf("Where:\n");
        printf("  N - number of time steps\n");
        printf("  M - number of spatial steps\n");
        return false;
    }
    
    try {
        N = std::stoi(argv[1]);
        M = std::stoi(argv[2]);
        
        if (N <= 0 || M <= 0) {
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

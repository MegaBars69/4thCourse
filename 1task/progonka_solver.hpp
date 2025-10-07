#ifndef PROGONKA_HEADER 
#define PROGONKA_HEADER

bool SolveSystem (double* A, double* b, double* x, int n);
void PrintMatrix (double *A, int n);
void PrintVector (double* x, int n);
void solve_tree_diag (double *up_diag, double *diag, double *low_diag, double *rhs, double *H_solution, int n);

#endif 
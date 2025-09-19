#ifndef LOG_HEADER 
#define LOG_HEADER

void PrintTriDiagMatrix (double *A, int n);
void PrintMatrix (double* A, int N, int M);
void PrintVector (double* x, int n);
bool init_args(int argc, char* argv[], double &T, double &X, int &N, int &M);

#endif 
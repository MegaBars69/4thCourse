#ifndef TABLE_PRINTER_HEADER 
#define TABLE_PRINTER_HEADER

#include <string>

void PrintLatexTable(double mu, int p_model, double* results, int num_n, int num_m);
void RunConvergenceTest(double mu, int p_model, int max_power, double* results);

#endif
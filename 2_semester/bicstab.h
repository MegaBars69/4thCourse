#pragma once

#include <Eigen/Sparse>
#include <vector>
#include <cmath>

class BiCGSTABSolver {
public:
    using Vector = Eigen::VectorXd;

    // Конструктор принимает ссылку на матрицу в формате Eigen::SparseMatrix
    // Матрица должна оставаться неизменной во время работы решателя
    BiCGSTABSolver(const Eigen::SparseMatrix<double>& matrix);

    // Решение системы A * x = b с начальным приближением x
    // Параметры: максимальное число итераций и допустимая норма относительной невязки
    bool solve(Vector& x, const Vector& b, int max_iter = 1000, double tol = 1e-10);

private:
    int n;                          // размерность системы
    const int* outerStarts;          // указатель на начало строк (CSR)
    const int* innerIndices;         // индексы столбцов
    const double* values;            // значения элементов
    Vector diag;                     // диагональ матрицы (для предобуславливателя)

    // Умножение матрицы на вектор (y = A * x)
    void matvec(const Vector& x, Vector& y) const;

    // Применение диагонального предобуславливателя: y = M^{-1} * x (M = diag)
    void apply_precond(const Vector& x, Vector& y) const;
};
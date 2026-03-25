#include "bicstab.h"
#include <iostream>

BiCGSTABSolver::BiCGSTABSolver(const Eigen::SparseMatrix<double>& matrix)
    : n(matrix.rows())
    , outerStarts(matrix.outerIndexPtr())
    , innerIndices(matrix.innerIndexPtr())
    , values(matrix.valuePtr())
    , diag(n)
{
    // Извлечение диагонали (предполагается, что диагональные элементы присутствуют)
    for (int i = 0; i < n; ++i) {
        diag[i] = 0.0;
        for (int j = outerStarts[i]; j < outerStarts[i+1]; ++j) {
            if (innerIndices[j] == i) {
                diag[i] = values[j];
                break;
            }
        }
        // Если диагональ отсутствует (редко), считаем её равной 1 (защита)
        if (diag[i] == 0.0) diag[i] = 1.0;
    }
}

void BiCGSTABSolver::matvec(const Vector& x, Vector& y) const
{
    y.setZero();
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = outerStarts[i]; j < outerStarts[i+1]; ++j) {
            sum += values[j] * x[innerIndices[j]];
        }
        y[i] = sum;
    }
}

void BiCGSTABSolver::apply_precond(const Vector& x, Vector& y) const
{
    for (int i = 0; i < n; ++i) {
        y[i] = x[i] / diag[i];
    }
}

bool BiCGSTABSolver::solve(Vector& x, const Vector& b, int max_iter, double tol)
{
    // Алгоритм BiCGSTAB с левым диагональным предобуславливателем
    Vector r(n), r_tilde(n), p(n), v(n), s(n), t(n), y(n), z(n);
    double rho, rho_prev = 1.0, alpha = 1.0, omega = 1.0, beta;
    double norm_b = b.norm();
    if (norm_b == 0.0) norm_b = 1.0;

    // Начальная невязка
    matvec(x, r);
    r = b - r;
    r_tilde = r;          // произвольный вектор, обычно берётся равным начальной невязке

    p.setZero();
    v.setZero();

    for (int iter = 0; iter < max_iter; ++iter) {
        rho = r_tilde.dot(r);
        if (std::abs(rho) < 1e-15) {
            // Вырожденный случай
            return false;
        }

        if (iter == 0) {
            p = r;
        } else {
            beta = (rho / rho_prev) * (alpha / omega);
            p = r + beta * (p - omega * v);
        }

        // Предобуславливание: y = M^{-1} p
        apply_precond(p, y);
        matvec(y, v);          // v = A * y

        alpha = rho / r_tilde.dot(v);

        s = r - alpha * v;

        // Предобуславливание: z = M^{-1} s
        apply_precond(s, z);
        matvec(z, t);          // t = A * z

        omega = t.dot(s) / t.dot(t);
        if (std::abs(omega) < 1e-15) {
            return false;
        }

        x += alpha * y + omega * z;
        r = s - omega * t;

        double norm_r = r.norm();
        if (norm_r < tol * norm_b) {
            return true;
        }

        rho_prev = rho;
    }
    return false;  // достигнут max_iter без сходимости
}
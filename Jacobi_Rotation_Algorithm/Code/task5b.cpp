#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace arma;

// Function for task 3.
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
    if (A.n_rows != A.n_cols) {
        throw std::invalid_argument("Jacobi eigensolver requires a square matrix");
    }
    if (A.n_rows < 2) {
        k = -1;
        l = -1;
        return 0.0;
    }

    double max_val = std::abs(A(0, 1));
    k = 0;
    l = 1;

    for (arma::uword i = A.n_rows - 2; i < A.n_rows; --i) {
        for (arma::uword j = i + 1; j < A.n_cols; ++j) {
            double current_abs = std::abs(A(i, j));

            if (current_abs > max_val) {
                max_val = current_abs;
                k = static_cast<int>(i);
                l = static_cast<int>(j);
            }
        }
    }
    return max_val;
}

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l)
{
    const int N = static_cast<int>(A.n_rows);
    if (A.n_rows != A.n_cols || R.n_rows != A.n_rows || R.n_cols != A.n_cols) {
        throw std::invalid_argument("Jacobi rotation received incompatible matrix sizes");
    }
    if (k < 0 || l < 0 || k >= N || l >= N || k == l) {
        throw std::invalid_argument("Jacobi rotation received invalid pivot indices");
    }

    const mat B = A;
    const mat Q = R;
    const double a_kl = B(k, l);

    if (a_kl == 0.0) {
        return;
    }

    const double tau = (B(l, l) - B(k, k)) / (2.0 * a_kl);
    const double t = (tau >= 0.0) ? 1.0 / (tau + std::sqrt(1.0 + tau * tau))
                                  : -1.0 / (-tau + std::sqrt(1.0 + tau * tau));
    const double c = 1.0 / std::sqrt(1.0 + t * t);
    const double s = t * c;

    A(k, k) = c * c * B(k, k) - 2.0 * c * s * B(k, l) + s * s * B(l, l);
    A(l, l) = s * s * B(k, k) + 2.0 * c * s * B(k, l) + c * c * B(l, l);
    A(k, l) = 0.0;
    A(l, k) = 0.0;

    for (int i = 0; i < N; ++i) {
        if (i != k && i != l) {
            A(i, k) = A(k, i) = c * B(i, k) - s * B(i, l);
            A(i, l) = A(l, i) = c * B(i, l) + s * B(i, k);
        }
        R(i, k) = c * Q(i, k) - s * Q(i, l);
        R(i, l) = s * Q(i, k) + c * Q(i, l);
    }
}

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues,
    arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged)
{
    if (A.n_rows != A.n_cols) {
        throw std::invalid_argument("Jacobi eigensolver requires a square matrix");
    }
    if (eps < 0.0 || maxiter < 0) {
        throw std::invalid_argument(
            "Jacobi eigensolver received an invalid tolerance or iteration limit");
    }

    mat C = A;
    const int N = static_cast<int>(A.n_rows);
    eigenvectors.eye(N, N);
    iterations = 0;
    converged = false;

    while (iterations < maxiter) {
        int k = -1;
        int l = -1;
        const double max_off = max_offdiag_symmetric(C, k, l);
        if (max_off <= eps) {
            converged = true;
            break;
        }
        jacobi_rotate(C, eigenvectors, k, l);
        ++iterations;
    }

    eigenvalues = C.diag();
}

// initialize matrix
int main()
{
    std::ofstream output("Task5b_iterations_random.txt");
    if (!output) {
        std::cerr << "Could not open Task5b_iterations_random.txt for writing." << std::endl;
        return 1;
    }

    output << "N iterations converged\n";

    for (int N = 10; N <= 100; N += 5) {
        // Generate random N*N matrix
        arma::mat A = arma::mat(N, N).randn();

        // Symmetrize the matrix by reflecting the upper triangle to lower triangle
        A = arma::symmatu(A);

        vec eigval;
        mat eigvec(N, N, fill::zeros);
        eigvec.diag(0).fill(1.0);

        bool converged = false;
        int maxiter = 100000;
        int iteration = 0;
        double eps = 1e-15;

        jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iteration, converged);

        std::cout << "N = " << N << " | Jacobi solver converged: " << std::boolalpha << converged
                  << " after " << iteration << " iterations." << std::endl;

        // Write to file
        output << N << ' ' << iteration << ' ' << converged << '\n';
    }

    output.close();
    return 0;
}
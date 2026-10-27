#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

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

    const arma::mat B = A;
    const arma::mat Q = R;
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

    arma::mat C = A;
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

inline void generate_matrices(const int n, arma::mat& A, arma::mat& eigvec)
{
    const double h = 1. / (double)n;
    const double h2 = h * h;
    const arma::uword N = n - 1;
    A = { N, N, arma::fill::zeros };
    A.diag().fill(2.0 / h2);
    A.diag(1).fill(-1.0 / h2);
    A.diag(-1).fill(-1.0 / h2);

    eigvec = { N, N, arma::fill::zeros };
}

void save_lowest_3(
    const arma::vec& eigval, const arma::mat& eigvec, const int n, const bool analytical)
{
    arma::uvec indices = arma::sort_index(eigval);
    arma::uvec lowest_3_idx = indices.head(3);

    arma::vec lowest_3_vals = eigval(lowest_3_idx);
    arma::mat lowest_3_vecs = eigvec.cols(lowest_3_idx);

    std::stringstream filename;
    filename << "task6_" << n << (analytical ? "_an" : "") << ".txt";

    std::ofstream outfile(filename.str());
    outfile << std::setprecision(16) << std::scientific;
    outfile << "Eigenvalues:\n";
    lowest_3_vals.raw_print(outfile);
    outfile << "Eigenvectors:\n";
    outfile << "0 0 0\n";
    lowest_3_vecs.raw_print(outfile);
    outfile << "0 0 0\n";
}

void solve(const int n)
{
    arma::mat A, eigvec;
    generate_matrices(n, A, eigvec);

    arma::vec eigval;
    bool converged = false;
    int maxiter = 25000;
    int iterations = 0;
    double eps = 1e-20;

    jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);
    save_lowest_3(eigval, eigvec, n, false);
}

void solve_analytical(const int n)
{
    const arma::uword N = n - 1;
    arma::vec eigval(N);
    arma::mat eigvec(N, N);

    const double h = 1. / (double)n;
    const double a = -1.0 / (h * h);
    const double d = 2.0 / (h * h);

    for (arma::uword j = 1; j <= N; ++j) {
        double lam = d + 2.0 * a * std::cos((j * arma::datum::pi) / (N + 1.0));
        eigval(j - 1) = lam;

        for (arma::uword i = 1; i <= N; ++i) {
            eigvec(i - 1, j - 1) = std::sin((i * j * arma::datum::pi) / (N + 1.0));
        }
    }

    for (arma::uword col = 0; col < N; ++col) {
        eigvec.col(col) = normalise(eigvec.col(col));
    }

    save_lowest_3(eigval, eigvec, n, true);
}

int main()
{
    solve(10);
    solve_analytical(10);

    solve(100);
    solve_analytical(100);
}
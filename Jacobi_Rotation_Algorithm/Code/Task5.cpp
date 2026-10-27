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
    std::ofstream output("Task5_iterations.txt");
    if (!output) {
        std::cerr << "Could not open Task5_iterations.txt for writing." << std::endl;
        return 1;
    }

    output << "N iterations converged\n";

    for (int N = 10; N <= 100; N += 5) {
        const double h = 1.0 / (N + 1.0);

        // 1. numerical calculation of eigenvectors
        mat A(N, N, fill::zeros);
        A.diag(0).fill(2.0);
        A.diag(-1).fill(-1.0);
        A.diag(1).fill(-1.0);
        A = (1.0 / (h * h)) * A;

        // 2. Analytical calculations of eigenvectors and values
        vec lambda(N);
        mat anaeigen(N, N);

        const double a = -1.0 / (h * h);
        const double d = 2.0 / (h * h);

        for (uword j = 1; j <= N; ++j) {
            double lam = d + 2.0 * a * std::cos((j * datum::pi) / (N + 1.0));
            lambda(j - 1) = lam;

            for (uword i = 1; i <= N; ++i) {
                anaeigen(i - 1, j - 1) = std::sin((i * j * datum::pi) / (N + 1.0));
            }
        }

        for (uword col = 0; col < N; ++col) {
            anaeigen.col(col) = normalise(anaeigen.col(col));
        }

        vec eigval;
        mat eigvec(N, N, fill::zeros);
        eigvec.diag(0).fill(1.0);
        bool converged = 0;
        int maxiter = 25000;
        int iteration = 0;
        double eps = 1e-20;
        jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iteration, converged);
        std::cout << "Jacobi solver converged: " << std::boolalpha << converged << " after "
                  << iteration << " iterations." << std::endl;

        // 4. Check of Eigenvectors
        vec out(N, fill::zeros); // 1.0 found, 0.0 not found
        const double tol = 1e-12; // Tolerance for this check

        for (uword col = 0; col < N; ++col) {
            bool contained = false;

            for (uword i = 0; i < N; ++i) {
                vec numi = eigvec.col(i);
                for (uword col = 0; col < N; ++col) {
                    vec anaeigi = anaeigen.col(col);
                    double diff_plus = norm(numi - anaeigi);
                    double diff_minus = norm(numi + anaeigi);
                    double diff_lambda = std::abs(eigval(i) - lambda(col));

                    if ((diff_plus < tol || diff_minus < tol) && diff_lambda < tol) {
                        contained = true;
                        break;
                    }
                }
                if (contained) {
                    break;
                }
            }
            out(col) = contained ? 1.0 : 0.0;
        }

        output << N << ' ' << iteration << ' ' << converged << '\n';
        // out.print("The found numerical Eigenvectors are also one of the analytical eigenvectors
        // in the following cases (All eigenvectors listed: 1= Agree, 0 = Different): ");
    }
    return 0;
}
#define _USE_MATH_DEFINES //unlocking mathematical constants 
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo> 

int main() {
    //create a 6*6 matrix
    int N = 6;
    int n = N + 1;
    double h = 1.0 / n; 
    double d = 2.0 / (h*h);
    double a = -1.0 / (h*h);

int main()
{
    const uword n = 7;
    double x_min = 0;
    double x_max = 1;

    double h = (x_max - x_min) / n;
    const uword N = n - 1;

    // 1. numerical calculation of eigenvectors
    mat A(N, N, fill::zeros);
    A.diag(0).fill(2.0);
    A.diag(-1).fill(-1.0);
    A.diag(1).fill(-1.0);
    A = (1.0 / (h * h)) * A;

    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);

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

    // 3. Normalisation of eigenvectors
    for (uword col = 0; col < N; ++col) {
        eigvec.col(col) = normalise(eigvec.col(col));
        anaeigen.col(col) = normalise(anaeigen.col(col));
    }

    // 4. Check of Eigenvectors
    vec out(N, fill::zeros); // 1.0 found, 0.0 not found
    const double tol = 1e-5; // Tolerance for this check

    for (uword col = 0; col < N; ++col) {
        vec anaeigi = anaeigen.col(col);
        bool contained = false;

        for (uword i = 0; i < N; ++i) {
            vec numi = eigvec.col(i);
            double numlam = eigval(i);
            double analam = lambda(i);

            // Check of both v and -v, because sign is unknown
            double diff_plus = norm(numi - anaeigi);
            double diff_minus = norm(numi + anaeigi);
            double diff_lambda = std::abs(numlam - analam);

            if ((diff_plus < tol || diff_minus < tol) && (diff_lambda < tol)) {
                contained = true;
                break; // Eigenvector and eigenvalue found
            }
        }
        out(col) = contained ? 1.0 : 0.0;
    }

    out.print("The found numerical Eigenvectors are also one of the analytical eigenvectors in the "
              "following cases (All eigenvectors listed: 1= Agree, 0 = Different): ");

    std::ofstream file("Task_2.txt");

    if (!file.is_open()) {
        std::cerr << "Could not open file\n";
        return 1;
    }

    file << "Numerical Eigenvectors\n";
    eigvec.save(file, arma::raw_ascii);
    file << "\n";

    file << "Analytical Eigenvectors\n";
    anaeigen.save(file, arma::raw_ascii);
    file << "\n";

    file << "Numerical Eigenvalues\n";
    eigval.save(file, arma::raw_ascii);
    file << "\n";

    file << "Analytical Eigenvalues\n";
    lambda.save(file, arma::raw_ascii);

    file.close();

    return 0;
}
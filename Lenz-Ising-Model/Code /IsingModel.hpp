#pragma once

#include <armadillo>
#include <vector>
#include <random>

class IsingModel {
    private:
    int L;
    int N;
    double T;
    double current_energy;
    double current_magnetization;

    //define the Matrix with armadillo

    arma::Mat<int> spin_lattice;
    std::vector<double> boltzman_factors; 

    //acumulator for variance calculations 

    double sum_energy = 0;
    double sum_energy2 = 0;
    double sum_m = 0;
    double sum_m2 = 0;
}
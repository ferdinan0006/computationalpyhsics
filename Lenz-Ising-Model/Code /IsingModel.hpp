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

    //stores five possible values for delta E so the code does not have to call the exp function every time 
    std::vector<double> boltzman_factors; 

    //acumulator for variance calculations 

    double sum_energy = 0;
    double sum_energy2 = 0;
    double sum_m = 0;
    double sum_m2 = 0;

    //using the pseudo random generator engine as rand() from <cstdlib> has weak statistical randomness
    std::mt19937 rng;

    //Metropolis acceptance test , a spin gets accepted with boltzmann probability 
    std::uniform_real_distribution<double> uniform_dist;

    //selects the spin at of the lattice random
    std::uniform_int_distribution<int> spin_dist;

    void initialize_Lattice(bool random_start);
    void precompute_boltzmann_factors();
    //need periodic boundary logic 
    inline int pbc(int idx);

public:
    IsingModel(int L_in, double T_in, bool random_start);
    void perform_mc_cycle();
    void void calculate_macroscopic_quantities(double& e, double& abs_m, double& cv, double& chi);
};
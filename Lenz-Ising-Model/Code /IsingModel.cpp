#include <IsingModel.hpp>
#include <cmath>

//Building the Constructor to setup the variables
IsingModel::IsingModel(int L_in, double T_in, bool random_start) {
    L = L_in;
    N = L*L;
    T = T_in;

    //allocate matrix size with aramdillo 
    spin_lattice = arma::Mat<int>(L,L);

    //random seed generation
    rng = std::mt19937(std::random_device{}());
    uniform_dist = std::uniform_real_distribution<double>(0.0, 1.0);
    spin_dist = std::uniform_int_distribution<int>(0,L-1);

    //internal setup functions call:
    initialize_Lattice(random_start);
    precompute_boltzmann_factors();

};

//building the functions
void IsingModel::initialize_Lattice(bool random_start) {

};

//precompute the factors for every possible dE (8,4,0,-4,-8)
void IsingModel::precompute_boltzmann_factors() {
    boltzman_factors.resize(9);
    for (int dE = -8; dE <= 8; dE += 4) {
        if (dE > 0) {
            boltzmann_factors[dE] = std::exp(-dE / T);
        }
        else {
            boltzmann_factors[dE] = 1.0;
        }
    }
};

//setting the boundary conditions 
inline int IsingModel::pbc(int idx) {
    return (idx + L) % L;
};

//doing one round of monte carlo simulations (n attempted spin flips)
void IsingModel::perform_mc_cycle() {
    for (int step = 0; step < N; ++ step) {
        //pick random site 
        int i = spin_dist(rng);
        int j = spin_dist(rng);

        // Compute dE using pbc(...) for neighbors
        // Evaluate Metropolis test
        // Update spin_lattice, current_energy, and current_magnetization if accepted

    }
    //calculating the sums 
    sum_e += current_energy / N;
    sum_e2 += (current_energy / N) * (current_energy / N);
    sum_m += std::abs(current_magnetization) / N;
    sum_m2 += (current_magnetization / N) * (current_magnetization / N);
    //increment mc_cylces
    mc_cycles++;
};

//calculating final expected values

void IsingModel::calculate_macroscopic_quantities(double& e, double& abs_m, double& cv, double& chi) {
    //create the average over all completed cycles 
};



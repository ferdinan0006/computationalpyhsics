#include "IsingModel.hpp"
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

//filling the matrix (the lattice)
void IsingModel::initialize_Lattice(bool random_start) {
    //fill the lattice with +1 or -1
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            if (random_start) {
                // 50% chance for -1, 50% chance for +1
                spin_lattice(i, j) = (uniform_dist(rng) < 0.5) ? -1 : 1;
            } else {
                spin_lattice(i, j) = 1; //ordered 
            }
        }
    }

    //calculate the initial energy and magnetization
    current_energy = 0.0;
    current_magnetization = 0.0;
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            int s = spin_lattice(i, j);
            current_magnetization += s;
            
            int right_neighbor = spin_lattice(i, pbc(j + 1));
            int bottom_neighbor = spin_lattice(pbc(i + 1), j);
            current_energy -= s * (right_neighbor + bottom_neighbor);
        }
    }
};

//precompute the factors for every possible dE (8,4,0,-4,-8)
void IsingModel::precompute_boltzmann_factors() {
    //resize to 9, fill all with 1.0 by default
    boltzman_factors.assign(9, 1.0);
    
    //calculate two positive cases 
    boltzman_factors[4] = std::exp(-4.0 / T);
    boltzman_factors[8] = std::exp(-8.0 / T);
}

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

        //compute the inner loop
        //find the current spin and the sum of its 4 periodic partners 
        int current_spin = spin_lattice(i,j);
        int sum_neighbors = spin_lattice(pbc(i+1), j) + spin_lattice(pbc(i-1), j) + spin_lattice(i, pbc(j+1)) + spin_lattice(i, pbc(j-1));
        //compute dE 
        int dE = 2* current_spin * sum_neighbors;
        //Metropolis test
        bool accept = false; 
        if (dE <= 0) {
            accept = true; //always accept energ lowering or similar energy
        }
        else {
            //draw random number between [0,1) for comparison with precomputed factor 
            if (uniform_dist(rng) <= boltzman_factors[dE]) {
                accept = true;
            }
        }
        //update spin lattice 
        if (accept) {
            spin_lattice(i,j) *= -1; //spin
            current_energy += dE; //add difference

            //the spin jut flipped, so the net change is exactly twice the net spin value 
            current_magnetization += 2 * spin_lattice(i, j);

        }

    }

    //calculating the sums 
    sum_energy  += current_energy / N;
    sum_energy2 += (current_energy / N) * (current_energy / N);
    sum_m += std::abs(current_magnetization) / N;
    sum_m2 += (current_magnetization / N) * (current_magnetization / N);
    //increment mc_cylces
    mc_cycles++;
};

//calculating final expected values

void IsingModel::calculate_macroscopic_quantities(double& e, double& abs_m, double& cv, double& chi) {
    if (mc_cycles == 0) return;
    //create the average over all completed cycles 
    //first order values after each spin
    e = sum_energy / mc_cycles; 
    abs_m = sum_m / mc_cycles;

    //second order values per spin 
    double avg_energy2 = sum_energy2 / mc_cycles;
    double avg_m2 = sum_m2 / mc_cycles; 

    //per spin variance 
    //per spin variance 
    double var_energy = avg_energy2 - (e * e);
    double var_m2 = avg_m2 - (abs_m * abs_m);

    //calculate Cv/N and chi/N (unitless kB = 1)
    cv = (N / (T * T)) * var_energy;
    chi = (N / T) * var_m2;

};



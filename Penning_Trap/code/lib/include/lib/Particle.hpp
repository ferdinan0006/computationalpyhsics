#pragma once
#include <armadillo>
#include <string>
#include <vector>

struct Particle {
    double q; // charge
    double m; // mass

    arma::vec3 position; // position vector r
    arma::vec3 velocity; // velocity vector v

    arma::vec3 r0; // initial vector positon r0
    arma::vec3 v0; // initial vector velocity v0
    std::vector<arma::vec> r_history; // trajectory history, storing position over time

    // Constructor
    Particle(double q, double m, arma::vec3 pos, arma::vec3 vel, double t0 = 0.0);

    // resetting particle information
    void reset(double t = 0.0);

    // save current state to r_history
    void save_state(double time);

    // write entire history of particle to txt file
    void save_history(const std::string& filename);
};
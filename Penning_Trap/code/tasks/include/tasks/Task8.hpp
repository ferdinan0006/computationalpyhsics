#pragma once
#include "lib/Solvers.hpp"
#include <armadillo>
#include <fstream>

class Task8 {
public:
    void single_particle_simulation();
    void two_particle_simulation();

private:
    void write_position(
        std::ofstream& file, const arma::vec3& r, double t, double RKdiff, double FEdiff);

    // General constant definitions for trap / particles (Given in task)
    const double q_in = 1;
    const double m_particle = 39.962590863; // Use isotope Ca40
    const double B0 = 96.5;
    const double V0 = 2.41e6;
    const double d0 = 500;
    const double total_time = 50.0;
};
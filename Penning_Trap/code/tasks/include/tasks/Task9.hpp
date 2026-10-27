#pragma once
#include "lib/Solvers.hpp"
#include <armadillo>

class Task9 {
public:
    void coarse_resonance_sweep();
    void fine_resonance_sweep(bool interactions);

private:
    void make_sweep(std::ofstream& file, double w0, double w_max, double dw, RK4Solver& rk4,
        double f, int N_particles, int iterations, double dt, bool interactions);

    // General constant definitions for trap / particles (Given in task)
    const double q_in = 1;
    const double m_particle = 39.962590863; // Use isotope Ca40
    const double B0 = 96.5;
    const double V0 = 2.41e6;
    const double d0 = 500;
    const double total_time = 500.0;
};
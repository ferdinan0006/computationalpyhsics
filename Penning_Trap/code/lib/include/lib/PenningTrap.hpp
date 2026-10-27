#pragma once

#include "Particle.hpp"
#include "lib/Solvers.hpp"
#include <armadillo>
#include <vector>

class PenningTrap {
public:
    PenningTrap(ISolver& solver, double B0_in, double V0_in, double d_in);

    void add_particle(Particle& p_in);
    void add_random_particles(double q, double m, int N);

    void step(double dt, bool interactions, bool record_history = true);

    // Helper function to save the position of all particles
    void save_positions(const std::string& filename);

    // Helper function to reset all particles to their initial state
    void reset_trap();

    // Compute the analytical solution for single particle specified by initial conditions given in
    // project
    arma::vec3 Ana_solution();

    int get_trapped_count();
    double get_time();
    Particle& operator[](int idx);

protected:
    virtual arma::vec3 electric_field(const arma::vec3& r, double t);

private:
    arma::vec3 magnetic_field(const arma::vec3& r);
    arma::vec3 force_particles(int i, int j);
    arma::vec3 total_force_external(int i, double t);
    std::vector<arma::vec3>& compute_forces(double t);

    double B0; // external B_field in z direction
    double V0; // initial potential
    double d; // definition is inside the report
    double t_total = 0;
    bool interactions = false;
    std::vector<Particle> particles;

    ISolver& solver;

    std::vector<arma::vec3> forces;
};

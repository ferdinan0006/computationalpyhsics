#pragma once
#include "lib/Particle.hpp"
#include <cstddef>

using compute_forces_f = std::function<std::vector<arma::vec3>&(double)>;

struct ISolver {
    virtual void step(
        double dt, std::vector<Particle>& particles, const compute_forces_f& compute_forces)
        = 0;
};

struct RK4Solver : ISolver {
    RK4Solver(std::size_t N);
    void step(double dt, std::vector<Particle>& particles, const compute_forces_f& compute_forces);

private:
    std::vector<arma::vec3> r0, v0;
    std::vector<arma::vec3> k1_r, k1_v;
    std::vector<arma::vec3> k2_r, k2_v;
    std::vector<arma::vec3> k3_r, k3_v;
    std::vector<arma::vec3> k4_r, k4_v;
};

struct EulerCromerSolver : ISolver {
    void step(double dt, std::vector<Particle>& particles, const compute_forces_f& compute_forces);
};

struct ForwardEulerSolver : ISolver {
    void step(double dt, std::vector<Particle>& particles, const compute_forces_f& compute_forces);
};

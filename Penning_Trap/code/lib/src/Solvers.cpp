#include "lib/Solvers.hpp"
#include "lib/Particle.hpp"
#include <cstddef>
#include <stdexcept>

RK4Solver::RK4Solver(std::size_t N)
    : k1_r(N)
    , k1_v(N)
    , k2_r(N)
    , k2_v(N)
    , k3_r(N)
    , k3_v(N)
    , k4_r(N)
    , k4_v(N)
    , r0(N)
    , v0(N)
{
}

void RK4Solver::step(
    double dt, std::vector<Particle>& particles, const compute_forces_f& compute_forces)
{
    /*
    Function to evlolve Penning trap system for one time step (dt) using the Runge-Kutta 4th order
    method. dt: time step size for evolution
    */

    const std::size_t N = particles.size();
    if (N != k1_r.size()) {
        throw std::runtime_error {
            "Particle vector size in RK4 Solver must match the size declared on construction"
        };
    }

    for (std::size_t i = 0; i < N; ++i) {
        r0[i] = particles[i].position;
        v0[i] = particles[i].velocity;
    }

    // k1
    const auto& forces1 = compute_forces(0);
    for (std::size_t i = 0; i < N; ++i) {
        k1_r[i] = dt * v0[i];
        k1_v[i] = dt * (forces1[i] / particles[i].m);
    }

    // k2 calculation with updated positions and velocities
    for (std::size_t i = 0; i < N; ++i) {
        particles[i].position = r0[i] + 0.5 * k1_r[i];
        particles[i].velocity = v0[i] + 0.5 * k1_v[i];
    }
    const auto& forces2 = compute_forces(dt * 0.5);
    for (std::size_t i = 0; i < N; ++i) {
        k2_r[i] = dt * particles[i].velocity;
        k2_v[i] = dt * (forces2[i] / particles[i].m);
    }

    // k3 calculation with updated positions and velocities
    for (std::size_t i = 0; i < N; ++i) {
        particles[i].position = r0[i] + 0.5 * k2_r[i];
        particles[i].velocity = v0[i] + 0.5 * k2_v[i];
    }
    const auto& forces3 = compute_forces(dt * 0.5);
    for (std::size_t i = 0; i < N; ++i) {
        k3_r[i] = dt * particles[i].velocity;
        k3_v[i] = dt * (forces3[i] / particles[i].m);
    }

    // k4 calculation with updated positions and velocities
    for (std::size_t i = 0; i < N; ++i) {
        particles[i].position = r0[i] + k3_r[i];
        particles[i].velocity = v0[i] + k3_v[i];
    }
    const auto& forces4 = compute_forces(dt);
    for (std::size_t i = 0; i < N; ++i) {
        k4_r[i] = dt * particles[i].velocity;
        k4_v[i] = dt * (forces4[i] / particles[i].m);
    }

    // final update of positions and velocities using the RK4 formula
    for (std::size_t i = 0; i < N; ++i) {
        particles[i].position = r0[i] + (k1_r[i] + 2.0 * k2_r[i] + 2.0 * k3_r[i] + k4_r[i]) / 6.0;
        particles[i].velocity = v0[i] + (k1_v[i] + 2.0 * k2_v[i] + 2.0 * k3_v[i] + k4_v[i]) / 6.0;
    }
}

// evolve system for one time step (dt) using Forward Euler
void ForwardEulerSolver::step(
    double dt, std::vector<Particle>& particles, const compute_forces_f& compute_forces)
{
    /*
    Function to evlolve Penning trap system for one time step (dt) using the Forward Euler method.
    dt: time step size for evolution
    */

    std::vector<arma::vec3> accelerations(particles.size());

    const auto& forces = compute_forces(0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        accelerations[i] = forces[i] / particles[i].m;
    }

    for (std::size_t i = 0; i < particles.size(); ++i) {
        particles[i].position += dt * particles[i].velocity;
        particles[i].velocity += dt * accelerations[i];
    }
}

// evolve system for one time step (dt) using Euler-Cromer (Not included in report)
void EulerCromerSolver::step(
    double dt, std::vector<Particle>& particles, const compute_forces_f& compute_forces)
{
    /*
    Function to evlolve Penning trap system for one time step (dt) using the Euler-Cromer method.
    dt: time step size for evolution
    */
    std::vector<arma::vec3> accelerations(particles.size());

    const auto& forces = compute_forces(0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        accelerations[i] = forces[i] / particles[i].m;
    }

    for (std::size_t i = 0; i < particles.size(); ++i) {
        particles[i].velocity += dt * accelerations[i];
        particles[i].position += dt * particles[i].velocity;
    }
}
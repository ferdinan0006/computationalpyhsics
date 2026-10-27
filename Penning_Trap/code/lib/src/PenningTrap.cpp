#include "lib/PenningTrap.hpp"
#include "lib/Particle.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace {
constexpr double k_e = 1.38935333e5;
}

/*
Constructor
B0_in: magnetic field strength
V0_in: electric potential
d_in: characteristic dimension of the trap
*/
PenningTrap::PenningTrap(ISolver& solver, double B0_in, double V0_in, double d_in)
    : solver(solver)
    , B0(B0_in)
    , V0(V0_in)
    , d(d_in)
{
}

void PenningTrap::add_particle(Particle& p_in)
{
    /*
    Adds a particle to the trap
    p_in: Particle object to be added, type: Particle
    */
    particles.push_back(p_in);
}

void PenningTrap::add_random_particles(double q, double m, int N)
{
    particles.reserve(N);
    for (int i = 0; i < N; ++i) {
        arma::vec3 r = arma::vec3().randn() * 0.1 * d;
        arma::vec3 v = arma::vec3().randn() * 0.1 * d;
        particles.push_back(Particle { q, m, r, v });
    }
}

// external electric field
arma::vec3 PenningTrap::electric_field(const arma::vec3& r, double t)
{
    /*
    Calculates Electric field at given position
    r: 3 dimensional arma vector containing position
    */
    arma::vec3 E(arma::fill::zeros);
    // only valid inside the trap
    // outside, we ignore the field
    if (arma::norm(r) > d) {
        return E; // returns {0, 0, 0} if escaped
    }

    // E = V0 / d^2 * (x, y, -2z)
    E(0) = r(0); // E_x = x
    E(1) = r(1); // E_y = y
    E(2) = -2.0 * r(2); // E_z = -2z

    E *= (V0 / (d * d)); // scale by V0 / d^2
    return E;
}

// external B-field
arma::vec3 PenningTrap::magnetic_field(const arma::vec3& r)
{
    /*
    Calculates Magnetic field at given position
    r: 3 dimensional arma vector containing position
    */
    arma::vec3 B(arma::fill::zeros);
    // similar to the E-field
    if (arma::norm(r) > d) {
        return B;
    }
    // only z-direction
    B(2) = B0;
    return B;
}

// Force on the particles
arma::vec3 PenningTrap::force_particles(int i, int j)
{
    const arma::vec3& r_i = particles[i].position;
    const arma::vec3& r_j = particles[j].position;
    const arma::vec3 r_diff = r_i - r_j; // pointing from j to i

    //|r_i - r_j|, calculating the distance
    double distance = arma::norm(r_diff);
    if (distance == 0.0) {
        throw std::domain_error("Coulomb force is undefined for overlapping particles");
    }

    // coulombs law for two particles i and j
    double charge_product = particles[i].q * particles[j].q;

    return k_e * charge_product * r_diff / (distance * distance * distance);
}

// total force on particle_i from the external fields
arma::vec3 PenningTrap::total_force_external(int i, double t)
{
    const arma::vec3& r = particles[i].position;
    const arma::vec3& v = particles[i].velocity;
    double q = particles[i].q;

    arma::vec3 E = electric_field(r, t);
    arma::vec3 B = magnetic_field(r);

    // arma:cross for the cross product of two 3D vectors
    arma::vec3 lorentz_force = q * (E + arma::cross(v, B));

    return lorentz_force;
}

// total force on particle_i from both external fields and other particles
std::vector<arma::vec3>& PenningTrap::compute_forces(double t)
{
    if (forces.size() != particles.size()) {
        forces.resize(particles.size());
    }

    const int N = particles.size();
    for (int i = 0; i < N; i++) {
        forces[i] = total_force_external(i, t);
    }

    if (particles.size() < 2 || !interactions) {
        return forces;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            auto f = force_particles(i, j);
            forces[i] += f;
            forces[j] -= f;
        }
    }

    return forces;
}

void PenningTrap::save_positions(const std::string& filename)
{
    // Save history of every particle to a separate .txt file
    for (std::size_t i = 0; i < particles.size(); ++i) {
        std::string filename_i = filename + "_" + std::to_string(i);
        particles[i].save_history(filename_i);
    }
}

void PenningTrap::reset_trap()
{
    // Save history of every particle to a separate .txt file
    t_total = 0.0;
    for (auto& particle : particles) {
        particle.reset();
    }
}

arma::vec3 PenningTrap::Ana_solution()
{
    // Computes the analytical solution for the timestep at the current time
    // Only valid for one particle and initial conditions as given in the report
    if (particles.size() > 1) {
        throw std::invalid_argument(
            "There is no analytical solution for more than one Particle (currently implemented).");
    }
    if (particles.size() == 0) {
        throw std::invalid_argument("There is no Particle in the trap.");
    }

    // Calculate constants
    double omega_z = std::sqrt(2.0 * particles[0].q * V0 / (particles[0].m * d * d));
    double omega_0 = particles[0].q * B0 / particles[0].m;

    // Calculate z
    double A = particles[0].r0[2];
    double B = particles[0].v0[2] / omega_z;
    double z = A * std::cos(omega_z * t_total) + B * std::sin(omega_z * t_total);

    // Calculate x and y
    // Note: Only applicable to the initial conditions listed in the paper
    double omega_minus
        = 0.5 * (omega_0 - std::pow(((omega_0 * omega_0) - 2 * omega_z * omega_z), 0.5));
    double omega_plus
        = 0.5 * (omega_0 + std::pow(((omega_0 * omega_0) - 2 * omega_z * omega_z), 0.5));
    double Aplus
        = (particles[0].v0[1] + omega_minus * particles[0].r0[0]) / (omega_plus - omega_minus);
    double Aminus
        = (particles[0].v0[1] + omega_plus * particles[0].r0[0]) / (omega_plus - omega_minus);

    double x = Aplus * std::cos(-omega_plus * t_total + std::numbers::pi)
        + Aminus * std::cos(-omega_minus * t_total);
    double y = Aplus * std::sin(-omega_plus * t_total + std::numbers::pi)
        + Aminus * std::sin(-omega_minus * t_total);

    arma::vec3 r_ana = { x, y, z };
    return r_ana;
}

int PenningTrap::get_trapped_count()
{
    int count = 0;
    for (auto& particle : particles) {
        if (arma::norm(particle.position) < d) {
            count++;
        }
    }
    return count;
}

double PenningTrap::get_time() { return t_total; }

Particle& PenningTrap::operator[](int idx) { return particles[idx]; }

void PenningTrap::step(double dt, bool interactions, bool record_history)
{
    this->interactions = interactions;
    solver.step(dt, particles,
        [this](double t_offset) -> auto& { return this->compute_forces(t_total + t_offset); });
    t_total = t_total + dt;
    if (!record_history)
        return;
    for (auto& particle : particles) {
        particle.save_state(t_total);
    }
}
#include "tasks/Task8.hpp"
#include "lib/Particle.hpp"
#include "lib/PenningTrap.hpp"
#include <armadillo>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

// Helper function to write a vector containing the positions to file
void Task8::write_position(
    std::ofstream& file, const arma::vec3& r, double t, double RKdiff, double FEdiff)
{
    file << std::scientific << std::setprecision(16) << r(0) << ' ' << r(1) << ' ' << r(2) << ' '
         << RKdiff << ' ' << FEdiff << ' ' << t << '\n';
};

void Task8::single_particle_simulation()
{
    // Initial conditions
    arma::vec3 pos_in1("20.0 0.0 20.0");
    arma::vec3 vel_in1("0.0 25.0 0.0");

    // Create Particles for single Particle simulation
    Particle Particle_1FE = Particle(q_in, m_particle, pos_in1, vel_in1);
    Particle Particle_1RK = Particle(q_in, m_particle, pos_in1, vel_in1);

    // Create Trap for single Particle simulation
    ForwardEulerSolver fe;
    RK4Solver rk4 { 1 };
    PenningTrap Penning1FE = PenningTrap(fe, B0, V0, d0);
    PenningTrap Penning1RK4 = PenningTrap(rk4, B0, V0, d0);
    Penning1FE.add_particle(Particle_1FE);
    Penning1RK4.add_particle(Particle_1RK);

    // Setup for first simulations
    std::vector<int> iterations_compare = { 4000, 8000, 16000, 32000 };

    // Create files containing the numerical and analytical solution with n=4000,8000,16000,32000
    // for both RK4 and FE
    std::filesystem::create_directories("data");
    for (int iterations : iterations_compare) {
        const double dt = total_time / iterations;
        // Reset traps for next run
        Penning1FE.reset_trap();
        Penning1RK4.reset_trap();

        // Create files for analytical solution
        std::ofstream ana("./data/analytical_single_" + std::to_string(iterations) + "_0.txt");
        if (!ana) {
            throw std::runtime_error("Could not open analytical output files");
        }

        // Compute first relative distance to analytical solution (computes a zero)
        double diff = 0;
        arma::vec3 diff_RK4 = Penning1RK4.Ana_solution() - Penning1RK4[0].position;
        arma::vec3 diff_FE = Penning1RK4.Ana_solution() - Penning1FE[0].position;
        double rel_RK4 = arma::norm(diff_RK4) / arma::norm(Penning1RK4.Ana_solution());
        double rel_FE = arma::norm(diff_FE) / arma::norm(Penning1RK4.Ana_solution());

        // Write initial first analytical position
        write_position(ana, Penning1RK4.Ana_solution(), Penning1FE.get_time(), rel_RK4, rel_FE);

        for (int i = 0; i < iterations; ++i) {
            Penning1FE.step(dt, false);
            Penning1RK4.step(dt, false);
            diff_RK4 = Penning1RK4.Ana_solution() - Penning1RK4[0].position;
            diff_FE = Penning1RK4.Ana_solution() - Penning1FE[0].position;
            rel_RK4 = arma::norm(diff_RK4) / arma::norm(Penning1RK4.Ana_solution());
            rel_FE = arma::norm(diff_FE) / arma::norm(Penning1RK4.Ana_solution());
            write_position(ana, Penning1RK4.Ana_solution(), Penning1FE.get_time(), rel_RK4, rel_FE);
        }

        // Save timesteps for numeric solutions. (Saves complete history with (x,y,z,vx,vy,vz,t))
        Penning1RK4.save_positions("./data/RK4_single_" + std::to_string(iterations));
        Penning1FE.save_positions("./data/FE_single_" + std::to_string(iterations));
    }
}

void Task8::two_particle_simulation()
{
    // Timestep settings used by the following simulations.
    int iterations = 4000;
    double dt = total_time / iterations;

    // Simulation of two particles without interactions
    // Generate Traps for both RK and FE
    ForwardEulerSolver fe;
    RK4Solver rk4 { 2 };
    PenningTrap Penning2FE = PenningTrap(fe, B0, V0, d0);
    PenningTrap Penning2RK4 = PenningTrap(rk4, B0, V0, d0);

    // Initial conditions for particle 1
    arma::vec3 pos_in2("20.0 0.0 20.0");
    arma::vec3 vel_in2("0.0 25.0 0.0");

    // Initial conditions for particle 2
    arma::vec3 pos_in3("25.0 25.0 0.0");
    arma::vec3 vel_in3("0.0 40.0 5.0");

    // Create Particle 1 for both FE and RK
    Particle Particle_2FE_1 = Particle(q_in, m_particle, pos_in2, vel_in2);
    Particle Particle_2RK_1 = Particle(q_in, m_particle, pos_in2, vel_in2);

    // Add particle 1 to traps
    Penning2FE.add_particle(Particle_2FE_1);
    Penning2RK4.add_particle(Particle_2RK_1);

    // Create Particle 2 for both traps
    Particle Particle_2FE_2 = Particle(q_in, m_particle, pos_in3, vel_in3);
    Particle Particle_2RK_2 = Particle(q_in, m_particle, pos_in3, vel_in3);

    // Add particle 2 to traps
    Penning2FE.add_particle(Particle_2FE_2);
    Penning2RK4.add_particle(Particle_2RK_2);

    // Evolve Penning trap objects
    for (int i = 0; i < iterations; ++i) {
        Penning2FE.step(dt, false);
        Penning2RK4.step(dt, false);
    }

    // Save the positions at all timesteps
    Penning2RK4.save_positions("./data/RK4_2_no_inter");
    Penning2FE.save_positions("./data/FE_2_no_inter");

    // Simulation of two particles with interactions
    // Reset the positions of both particles
    Penning2FE.reset_trap();
    Penning2RK4.reset_trap();

    // Evolve Penning trap objects
    for (int i = 0; i < iterations; ++i) {
        Penning2FE.step(dt, true);
        Penning2RK4.step(dt, true);
    }

    // Save the positions at all timesteps
    Penning2RK4.save_positions("./data/RK4_2_w_inter");
    Penning2FE.save_positions("./data/FE_2_w_inter");
}
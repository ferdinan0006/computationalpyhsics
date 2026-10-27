#include "tasks/Task9.hpp"
#include "lib/PenningTrapPerturbed.hpp"
#include <filesystem>

void Task9::coarse_resonance_sweep()
{
    const double fs[] = { 0.1, 0.4, 0.7 };
    const double dw = 0.02;
    const double w0 = 0.2;
    const double w_max = 2.5;
    const int N_particles = 100;
    const int iterations = 40000;
    const double dt = total_time / iterations;
    const double eps = 1e-6;
    RK4Solver rk4 { N_particles };

    std::filesystem::create_directories("data");
    std::ofstream file("./data/coarse_resonance_sweep.txt");
    if (!file) {
        throw std::runtime_error("Could not open output files");
    }
    for (auto f : fs) {
        make_sweep(file, w0, w_max, dw, rk4, f, N_particles, iterations, dt, false);
    }
    file.close();
}

void Task9::fine_resonance_sweep(bool interactions)
{
    const double fs[] = { 0.1, 0.4, 0.7 };
    const double dw = 0.002;
    const double w0 = 1.0;
    const double w_max = 1.8;
    const int N_particles = 100;
    const int iterations = 40000;
    const double dt = total_time / iterations;
    RK4Solver rk4 { N_particles };

    std::filesystem::create_directories("data");
    std::ofstream file(interactions ? "./data/fine_resonance_sweep_inter.txt"
                                    : "./data/fine_resonance_sweep_nointer.txt");
    if (!file) {
        throw std::runtime_error("Could not open output files");
    }

    for (auto f : fs) {
        make_sweep(file, w0, w_max, dw, rk4, f, N_particles, iterations, dt, interactions);
    }

    file.close();
}

void Task9::make_sweep(std::ofstream& file, double w0, double w_max, double dw, RK4Solver& rk4,
    double f, int N_particles, int iterations, double dt, bool interactions)
{
    const double eps = 1e-6;
    for (double w = w0; w <= w_max + eps; w += dw) {
        PenningTrapPerturbed trap { rk4, B0, V0, d0, f, w };
        trap.add_random_particles(q_in, m_particle, N_particles);

        for (int i = 0; i < iterations; ++i) {
            trap.step(dt, interactions, false);
        }

        auto count = trap.get_trapped_count();
        std::cout << std::format("f={:.1f}  w={:5.3f}  -  N={:3d}", f, w, count) << std::endl;
        file << count << " ";
    }
    file << std::endl;
}
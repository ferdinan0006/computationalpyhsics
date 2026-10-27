#include "tasks/Task8.hpp"
#include "tasks/Task9.hpp"

int main()
{
    Task8 task8;
    std::cout << "Performing single particle simulation..." << std::flush;
    task8.single_particle_simulation();
    std::cout << " Done." << std::endl;
    std::cout << "Performing two particle simulation..." << std::flush;
    task8.two_particle_simulation();
    std::cout << " Done." << std::endl;

    Task9 task9;
    std::cout << "Performing coarse resonance sweep..." << std::endl;
    task9.coarse_resonance_sweep();
    std::cout << "Performing fine resonance sweep with particle interactions..." << std::endl;
    task9.fine_resonance_sweep(true);
    std::cout << "Performing fine resonance sweep without particle interactions..." << std::endl;
    task9.fine_resonance_sweep(false);
}
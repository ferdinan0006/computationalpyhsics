#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

//time measurments and pre-sized vector
#include <cmath>
#include <chrono>
#include <vector>
#include <omp.h>

#include "IsingModel.hpp"

int main(int argc, char* argv[]) {
    //the line makes sure that exactly 4 user parameters were passed (L, T, mc_cycles, filename)
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <L> <MC_cycles> <output_file.csv>\n";
        return 1; //exit with err
    }

    //converting string to integer (e.g. L)
    int L = std::stoi(argv[1]);
    int mc_cycles = std::stoi(argv[2]);
    
    //parsing floating point values for Temp
    
    //would implement the temperature loop here for task 8 and 9 and use OpenMD parallel computing 
    std::string filename = argv[3];

    double T_min = 2.2;
    double T_max = 2.35;
    double dT = 0.005;

    struct TemperatureResult {
        double T;
        double e;
        double abs_m;
        double chi;
        double cv;
        double current_energy_per_spin;
        double simulation_seconds;
    };

    const int n_temperatures = static_cast<int>(std::lround((T_max - T_min) / dT)) +1;
    std::vector<TemperatureResult> results(n_temperatures);

    //time measuring 
    const auto start = std::chrono::steady_clock::now(); 

    std::cout << "Starting simulation..." << std::endl;

    //parallel computing with OpenMP
    int N = L * L;

   //print actual thread count
    #pragma omp parallel
    {
        #pragma omp single
        {
            std::cout << "OpenMP threads: "
                    << omp_get_num_threads()
                    << std::endl;
        }
    }

    //parallelize the temperature loop
    #pragma omp parallel for schedule(static)
    for (int step = 0; step < n_temperatures; ++step) {
        const auto temperature_start = std::chrono::steady_clock::now();
        const double T = T_min + step * dT;
        IsingModel model(L,T, true);

        for (int i = 0; i< mc_cycles; ++i) {
            model.perform_mc_cycle();
        }

        double e, abs_m, cv, chi; 
        model.calculate_macroscopic_quantities(e, abs_m, cv, chi);

        double current_energy_per_spin = model.get_current_energy() / (double)N;
        const auto temperature_end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> temperature_elapsed =
            temperature_end - temperature_start;

        results[step] = {
            T, e, abs_m, chi, cv, current_energy_per_spin,
            temperature_elapsed.count()
        };
    }

    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = end - start; 

    std::cout << "All temperatures completed!" << std::endl;

    std::cout <<"Simulation time: " << elapsed.count() << "s\n";

    //open outfile 
    std::ofstream outfile(filename, std::ios::app);

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file\n";
        return 1;
    }

    if (outfile.tellp() == 0) {
        //write outfile header 
        outfile << "L,T,mc_cycles,e,abs_m,cv,chi,current_energy_per_spin,simulation_seconds,total_simulation_second\n";
    }

    //write outputfile after all simulations are completed 
    for (const auto& result : results) {
        outfile << std::fixed << std::setprecision(8)
                << L << ","
                << result.T << ","
                << mc_cycles << ","
                << result.e << ","
                << result.abs_m << ","
                << result.cv << ","
                << result.chi << ","
                << result.current_energy_per_spin << ","
                << result.simulation_seconds << "," 
                << elapsed.count() << "\n";
    }

    return 0;
}
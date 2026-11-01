#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
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
    
    // I would implement the temperature loop here for task 8 and 9 and use OpenMD parallel computing 
    std::string filename = argv[3];

    double T_min = 1.0;
    double T_max = 2.4;
    double dT = 0.2;
    int n_steps = static_cast<int>((T_max - T_min) / dT) + 1;

    for (int step = 0; step <= n_steps; ++step) {
        double T = T_min + step * dT;

        //executing the class constructor from IsingModel
        IsingModel model(L, T, true);

        //run monte carlo cycles 
        for (int i = 0; i < mc_cycles; ++i) {
            model.perform_mc_cycle();
        }

        //calculate macroscopic quantities
        double e, abs_m, cv, chi;
        model.calculate_macroscopic_quantities(e, abs_m, cv, chi);

        //write results to output file
        std::ofstream outfile;
        outfile.open(filename, std::ios::app);
        if (!outfile.is_open()) {
            std::cerr << "Error: Could not open output file.\n";
            return 1; //report err when file could not be open
        }

        //writing header for empty files (AI suggestion)
        outfile.seekp(0, std::ios::end);
        if (outfile.tellp() == 0) {
            outfile << "L,T,MC_cycles,e,abs_m,cv,chi\n";
        }

        //appending the computed data with high precision!
        outfile << std::fixed << std::setprecision(8)
                << L << ","
                << T << ","
                << mc_cycles << ","
                << e << ","
                << abs_m << ","
                << cv << ","
                << chi << "\n";

        outfile.close();
    }

    return 0;
}
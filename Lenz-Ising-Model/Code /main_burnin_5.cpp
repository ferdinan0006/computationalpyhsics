#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "IsingModel.hpp"

int main(int argc, char* argv[]) {
    //the code expects five parameter, first checking if the request delivers all of them
    if (argc != 6) {
        std::cerr <<'Usage: ' <<argv[0] <<'<L>, <mc_cycles>, <T>, <Ordered Start (1 = Yes, 0 = No)>, <output_filename.csv>\n';
        return 1; 
    }

    int L = std::stoi(argv[1]);
    int mc_cycles = std::stoi(argv[2]);
    double T = std::stod(argv[3]);

    //if ordered state is true, random start has to be set to zero so every spin will be filled with + 1
    int ordered_input = std::stoi(argv[4]);
    bool random_state = (ordered_input == 0);

    std::string filename = argv[5];

    //initializing the model 
    IsingModel model(L, T, random_state);

    //if the output file is urrently open, return error since its a specif burn in trace 
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr <<"Error: Could not open the file\n";
        return 1;
    }
    //write the files header 
    outfile <<"cycle, cumulative_e, cumulative_abs\n";

    std::cout <<"Running Burn in simulation for T=" << T << (random_state ? " (Random Start)" : " (Ordered Start)") << "...\n";

    //runnning mc with metropolis 
    for (int cycle = 1; cycle <= mc_cycles; ++cycle) {
        model.perform_mc_cycle();

        double e,abs_m, cv, chi;
        //calculating the average using all values calculated so far 
        model.calculate_macroscopic_quantities(e, abs_m, cv, chi);
        //write cycle number and current epsilon to output file
        outfile << std::fixed << std::setprecision(8) << cycle << ","  << e << ","  << abs_m << "\n";

    }

    outfile.close();
    std::cout << "Done! Burn-in trace saved to " << filename << "\n";
    return 0;

};
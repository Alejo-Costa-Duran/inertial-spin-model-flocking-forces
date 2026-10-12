#include "simulation.h"
#include <iostream>
#include <vector>
#include <chrono>

int main(int arg_count, char* args[]) {
    if (arg_count > 1) {
        std::string folder = args[1];
        std::string data_file = args[3];
        std::cout << "Saving files to directory: " << folder << "\n";
        
        Simulation s(data_file);
        int counter = 0;
        std::vector<double> inertias = {0.0003, 0.005, 0.2, 1.0, 2.0, 5.0};
        std::vector<double> fields = {0.0, 0.01, 0.1, 0.5, 1.0, 2.0};
        
        for (auto chi : inertias) {
            for (auto f : fields) {
                int seed = rand();
                s.set_par("Inertia", chi);
                s.set_par("Field", f);
                s.set_par("Seed", seed);
                
                std::string run_id = std::to_string(counter);
                run_id.insert(0, 3 - run_id.length(), '0');
                counter++;
                
                auto start = std::chrono::high_resolution_clock::now();
                s.simMosquito(s.parameters, folder, run_id, "/Info");
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
                
                std::cout << "Finished run " << run_id << " in " << duration.count() << " seconds\n";
            }
        }
    } else {
        std::cout << "Usage: <executable> <output_folder> <dummy> <params_file>\n";
        return 1;
    }
    return 0;
}
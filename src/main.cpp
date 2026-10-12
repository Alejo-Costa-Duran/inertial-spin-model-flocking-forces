#include "simulation.h"
#include <iostream>
#include <vector>
#include <chrono>

int main(int arg_count, char* args[]) {
    if (arg_count > 3) {
        std::string folder = args[1];
        std::string data_file = args[3];
        std::cout << "Saving files to directory: " << folder << "\n";
        
        Simulation s(data_file);
        int counter = 0;
        std::vector<double> inertias = {0.0003, 0.005, 0.2, 1.0, 2.0, 5.0};
        std::vector<double> fields = {0.0, 0.01, 0.1, 0.5, 1.0, 2.0};
        
        for (auto chi : inertias) {
            for (auto f : fields) {
                int seed = 1337 + counter;
                s.set_par("Inertia", chi);
                s.set_par("QuarticK", f); // Maps your loop variable 'f' to our quartic confinement/field parameter
                s.set_par("Seed", seed);
                
                std::string run_id = std::to_string(counter);
                run_id.insert(0, 3 - run_id.length(), '0');
                counter++;
                
                std::string run_folder = folder + "/run_" + run_id;
                // Create subfolders or pass folder path accordingly
                
                auto start = std::chrono::high_resolution_clock::now();
                s.runSimulation(folder); // Calls the refactored simulation runner
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
                
                std::cout << "Finished run " << run_id << " in " << duration.count() << " seconds\n";
            }
        }
    } else {
        std::cout << "Usage: ism_batch <output_folder> <dummy> <params_file>\n";
        return 1;
    }
    return 0;
}
#include "simulation.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> vec;
    if (num <= 1) {
        vec.push_back(start);
        return vec;
    }
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        vec.push_back(start + i * step);
    }
    return vec;
}

int main(int arg_count, char* args[]) {
    if (arg_count > 2) {
        std::string base_folder = args[1];
        std::string data_file = args[2];
        std::cout << "Base output directory: " << base_folder << "\n";
        
        Simulation s(data_file);
        
        double chi_init = s.parameters.count("InertiaInit") ? s.parameters["InertiaInit"] : 0.0003;
        double chi_final = s.parameters.count("InertiaFinal") ? s.parameters["InertiaFinal"] : 5.0;
        int chi_steps = s.parameters.count("InertiaSteps") ? static_cast<int>(s.parameters["InertiaSteps"]) : 6;
        
        double temp_init = s.parameters.count("TempInit") ? s.parameters["TempInit"] : 0.0;
        double temp_final = s.parameters.count("TempFinal") ? s.parameters["TempFinal"] : 0.5;
        int temp_steps = s.parameters.count("TempSteps") ? static_cast<int>(s.parameters["TempSteps"]) : 6;
        
        std::vector<double> inertias = linspace(chi_init, chi_final, chi_steps);
        std::vector<double> temperatures = linspace(temp_init, temp_final, temp_steps);
        
        s.set_par("QuarticK", 0.5);   
        s.set_par("HarmonicK", 0.0);  
        s.set_par("RadialStr",  0.0);  
        s.set_par("SwirlStr",   0.0);   
        
        int counter = 0;
        int total_runs = inertias.size() * temperatures.size();
        
        std::cout << "Starting parameter sweep: " << inertias.size() << " Inertia values x " 
                  << temperatures.size() << " Temperature values (" << total_runs << " total runs).\n";
        
        for (auto chi : inertias) {
            for (auto temp : temperatures) {
                int seed = 1337 + counter;
                s.set_par("Inertia", chi);
                s.set_par("Temperature", temp);
                s.set_par("Seed", seed);
                
                // Format run ID with leading zeros (e.g., run_000, run_001)
                std::ostringstream oss;
                oss << std::setw(3) << std::setfill('0') << counter;
                std::string run_id = oss.str();
                
                // Construct a unique subfolder for this specific run
                std::string run_folder = base_folder + "/run_" + run_id;
                counter++;
                
                auto start = std::chrono::high_resolution_clock::now();
                s.runSimulation(run_folder); // Pass the unique run folder!
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
                
                std::cout << "[" << counter << "/" << total_runs << "] Finished run " << run_id 
                          << " (Inertia: " << chi 
                          << ", Temp: " << temp 
                          << ") -> Saved to " << run_folder << " in " << duration.count() << "s\n";
            }
        }
        std::cout << "All simulations completed successfully with isolated run folders!\n";
    } else {
        std::cout << "Usage: ism_batch <output_folder> <params_file>\n";
        return 1;
    }
    return 0;
}
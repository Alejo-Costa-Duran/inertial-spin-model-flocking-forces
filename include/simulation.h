#pragma once
#include <string>
#include <map>
#include <memory>
#include "flock.h"
#include "forcefield.h"

class Simulation {
public:
    std::string simulationConfigFile;
    std::map<std::string, double> parameters;

    Simulation() = default;
    explicit Simulation(const std::string& sim_file);
    ~Simulation() = default;

    void set_par(const std::string& par, double value);
    void readSimFile(const std::string& file);
    
    // Core simulation runner incorporating the modular ForceField framework
    void runSimulation(const std::string& outputFolder);
};
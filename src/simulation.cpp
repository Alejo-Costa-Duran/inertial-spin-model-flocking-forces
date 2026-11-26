#include "simulation.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

Simulation::Simulation(const std::string& sim_file) {
    simulationConfigFile = sim_file;
    readSimFile(simulationConfigFile);
}

void Simulation::set_par(const std::string& par, double value) {
    parameters[par] = value;
}

void Simulation::readSimFile(const std::string& file) {
    std::ifstream in(file);
    if (!in.is_open()) {
        std::cerr << "Error: Could not open configuration file: " << file << std::endl;
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '/') continue;

        std::stringstream ss(line);
        std::string key;
        double value;
        if (ss >> key >> value) {
            parameters[key] = value;
        }
    }
    in.close();
}

void Simulation::runSimulation(const std::string& outputFolder) {
    // 1. Extract configuration parameters with default fallbacks
    int seed = static_cast<int>(parameters.count("Seed") ? parameters["Seed"] : 42);
    int numBirds = static_cast<int>(parameters.count("NumBirds") ? parameters["NumBirds"] : 100);
    int totalSteps = static_cast<int>(parameters.count("TotalSteps") ? parameters["TotalSteps"] : 10000);
    int logSteps = static_cast<int>(parameters.count("LogSteps") ? parameters["LogSteps"] : 100);
    int saveSteps = static_cast<int>(parameters.count("SaveSteps") ? parameters["SaveSteps"] : 1000);

    
    double boxSize = parameters.count("BoxSize") ? parameters["BoxSize"] : 30.0;
    double timeStep = parameters.count("TimeStep") ? parameters["TimeStep"] : 0.001;
    double inertia = parameters.count("Inertia") ? parameters["Inertia"] : 1.0;
    double friction = parameters.count("Friction") ? parameters["Friction"] : 1.0;
    double temperature = parameters.count("Temperature") ? parameters["Temperature"] : 0.0;
    double coupling = parameters.count("Interaction") ? parameters["Interaction"] : 1.0;
    double maxSpeed = parameters.count("MaxSpeed") ? parameters["MaxSpeed"] : 1.0;
    double WellStrength = parameters.count("WellStrength") ? parameters["WellStrength"] : 0.0;
    int WellExponent = parameters.count("WellExponent") ? parameters["WellExponent"] : 12;
    double WellRadius = parameters.count("WellRadius") ? parameters["WellRadius"] : 10.0;
    // Force parameters
    double harmonicK = parameters.count("HarmonicK") ? parameters["HarmonicK"] : 0.0;
    double quarticK = parameters.count("QuarticK") ? parameters["QuarticK"] : 0.0;
    double radialStr = parameters.count("RadialStr") ? parameters["RadialStr"] : 0.0;
    double swirlStr = parameters.count("SwirlStr") ? parameters["SwirlStr"] : 0.0;
    bool periodicBoundary = parameters.count("PeriodicBoundary") &&
                            parameters["PeriodicBoundary"] != 0.0;
    bool randomInitialConditions = parameters.count("RandomInitialConditions") &&
                                   parameters["RandomInitialConditions"] != 0.0;

    // 2. Build the modular ForceField structure dynamically
    auto compositeForce = std::make_shared<CompositeForce>();
    bool hasForces = false;
    if (WellStrength > 0.0) {
        compositeForce->addForce(std::make_shared<GeneralizedPowerConfinement>(WellStrength, WellRadius, WellExponent));
        hasForces = true;
    }
    if (harmonicK > 0.0) {
        compositeForce->addForce(std::make_shared<HarmonicConfinement>(harmonicK));
        hasForces = true;
    }
    if (quarticK > 0.0) {
        compositeForce->addForce(std::make_shared<QuarticConfinement>(quarticK)); // Steep r^4 well
        hasForces = true;
    }
    if (radialStr != 0.0 || swirlStr != 0.0) {
        compositeForce->addForce(std::make_shared<CentralSwirlForce>(radialStr, swirlStr)); // Swirling / central driving force
        hasForces = true;
    }

    if (periodicBoundary && hasForces) {
        throw std::invalid_argument(
            "PeriodicBoundary=1 cannot be combined with active position-dependent force fields."
        );
    }

    const ForceField* activeField = hasForces ? compositeForce.get() : nullptr;

    std::filesystem::create_directories(outputFolder);

    // 3. Initialize Flock and spawn birds in a swirling/circular initial distribution
    Flock flock(seed, boxSize, timeStep, maxSpeed, inertia, friction, temperature, coupling, 1.0, 0.0, periodicBoundary);

    for (int idx = 0; idx < numBirds; ++idx) {
        Pvector pos(0.0, 0.0);
        Pvector vel(0.0, 0.0);
        Pvector acc(0.0, 0.0);
        if (randomInitialConditions) {
            const double initialRadius = std::min(WellRadius, boxSize);
            const double radialSample = initialRadius * std::sqrt(gsl_rng_uniform(flock.m_mt));
            const double positionAngle = 2.0 * M_PI * gsl_rng_uniform(flock.m_mt);
            pos.x = radialSample * std::cos(positionAngle);
            pos.y = radialSample * std::sin(positionAngle);
            const double heading = 2.0 * M_PI * gsl_rng_uniform(flock.m_mt);
            vel.x = maxSpeed * std::cos(heading);
            vel.y = maxSpeed * std::sin(heading);
        } else {
            const double r0 = 5.0 + (gsl_rng_uniform(flock.m_mt) * 5.0);
            const double theta0 = 2.0 * M_PI * gsl_rng_uniform(flock.m_mt);
            pos.x = r0 * std::cos(theta0);
            pos.y = r0 * std::sin(theta0);
            vel.x = -maxSpeed * std::sin(theta0);
            vel.y = maxSpeed * std::cos(theta0);
            acc.x = -maxSpeed * maxSpeed * std::cos(theta0) / r0;
            acc.y = -maxSpeed * maxSpeed * std::sin(theta0) / r0;
        }
        
        flock.bandada.push_back(Bird(false, idx, pos, vel, acc));
    }
    if (periodicBoundary) flock.boundary();
    flock.calculateForces(activeField);
    if (randomInitialConditions) flock.initializeAccelerations();

    // 4. Open CSV files for logging metrics and system states
    std::ofstream metricsFile(outputFolder + "/metrics.csv");
    if (!metricsFile.is_open()) {
        std::cerr << "Error: Could not open metrics file in " << outputFolder << "\n";
        return;
    }
    metricsFile << "Step,TotalSpin,SpinPerBird,AvgVelX,AvgVelY,TangentialSpeed,AvgRnorm,cmX,cmY,TotalEnergy\n";    std::cout << "Starting simulation run for " << numBirds << " birds over " << totalSteps << " steps...\n";

    for (int step = 0; step < totalSteps; ++step) {
        // Log metrics every N steps
        if (step % logSteps == 0) {
            double totalSpin = 0.0;
            double sumVelX = 0.0;
            double sumVelY = 0.0;
            double totalTangSpeed = 0.0;
            double totalRnorm = 0.0;
            // New variables for Center of Mass
            double sumPosX = 0.0;
            double sumPosY = 0.0;
            


            for (const auto& bird : flock.bandada) {
                totalSpin += bird.getSpin(inertia, maxSpeed);
                sumVelX += bird.velocity.x;
                sumVelY += bird.velocity.y;
                
                double rNorm = bird.position.getNorm();
                totalRnorm += rNorm;
                sumPosX += bird.position.x;
                sumPosY += bird.position.y;
                if (rNorm > 1e-12) {
                    // Tangential velocity relative to origin: (v x r) / |r|
                    totalTangSpeed += bird.velocity.crossProd(bird.position) / rNorm;
                }
            }

            double spinPerBird = totalSpin / numBirds;
            double avgVelX = sumVelX / numBirds;
            double avgVelY = sumVelY / numBirds;
            // Compute Center of Mass
            double cmX = sumPosX / numBirds;
            double cmY = sumPosY / numBirds;
            double avgTangSpeed = totalTangSpeed / numBirds;
            double totalEnergy = flock.totalEnergy();
            double avgRnorm = totalRnorm / numBirds;

            metricsFile << step << "," 
                        << totalSpin << "," 
                        << spinPerBird << "," 
                        << avgVelX << "," 
                        << avgVelY << "," 
                        << avgTangSpeed << ","
                        << avgRnorm << ","
                        << cmX << ","
                        << cmY << "," 
                        << totalEnergy << "\n";
        }

        // Save full system state snapshot every M steps
        if (step % saveSteps == 0) {
            std::string snapEnv = outputFolder + "/env_step_" + std::to_string(step) + ".csv";
            std::string snapBirds = outputFolder + "/birds_step_" + std::to_string(step) + ".csv";
            std::string rngFile = outputFolder + "/rng_step_" + std::to_string(step) + ".bin";
            flock.saveFile(const_cast<char*>(rngFile.c_str()), snapEnv, snapBirds, timeStep);
        }

        // Step the system forward, passing our modular ForceField pointer
        flock.updateFlock(timeStep, activeField);
    }

    metricsFile.close();
    std::cout << "Simulation completed successfully. Results saved to: " << outputFolder << "\n";
}
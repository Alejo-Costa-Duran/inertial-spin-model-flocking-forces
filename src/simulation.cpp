#include "simulation.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <vector>

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

    flock.saveEnvironment(outputFolder + "/environment.csv", timeStep);

    // 4. Open CSV files for logging metrics and system states
    std::ofstream metricsFile(outputFolder + "/metrics.csv");
    if (!metricsFile.is_open()) {
        std::cerr << "Error: Could not open metrics file in " << outputFolder << "\n";
        return;
    }
    metricsFile << "Step,TotalSpin,SpinPerBird,AvgVelX,AvgVelY,TangentialSpeed,AvgRnorm,cmX,cmY,TotalEnergy,"
                   "Polarization,RadiusGyration,MeanNeighborCount,NeighborCountVariance,SpinVariance,"
                   "EnergyVariance\n";
    std::cout << "Starting simulation run for " << numBirds << " birds over " << totalSteps << " steps...\n";

    for (int step = 0; step < totalSteps; ++step) {
        // Log metrics every N steps
        if (step % logSteps == 0) {
            double totalSpin = 0.0;
            double sumSpinSquared = 0.0;
            double sumVelX = 0.0;
            double sumVelY = 0.0;
            double totalTangSpeed = 0.0;
            double totalRnorm = 0.0;
            double sumPosX = 0.0;
            double sumPosY = 0.0;
            double sumCosX = 0.0;
            double sumSinX = 0.0;
            double sumCosY = 0.0;
            double sumSinY = 0.0;
            double totalEnergy = 0.0;
            double sumEnergySquared = 0.0;
            double totalNeighborCount = 0.0;
            double sumNeighborCountSquared = 0.0;
            const double birdCount = static_cast<double>(flock.bandada.size());
            const double interactionEnergyFactor =
                coupling / (2.0 * maxSpeed * maxSpeed);

            for (std::size_t idx = 0; idx < flock.bandada.size(); ++idx) {
                const Bird& bird = flock.bandada[idx];
                const double spin = bird.getSpin(inertia, maxSpeed);
                totalSpin += spin;
                sumSpinSquared += spin * spin;
                sumVelX += bird.velocity.x;
                sumVelY += bird.velocity.y;

                const double rNorm = bird.position.getNorm();
                totalRnorm += rNorm;
                sumPosX += bird.position.x;
                sumPosY += bird.position.y;
                if (periodicBoundary) {
                    const double xPhase = M_PI * bird.position.x / boxSize;
                    const double yPhase = M_PI * bird.position.y / boxSize;
                    sumCosX += std::cos(xPhase);
                    sumSinX += std::sin(xPhase);
                    sumCosY += std::cos(yPhase);
                    sumSinY += std::sin(yPhase);
                }
                if (rNorm > 1e-12) {
                    // Tangential velocity relative to origin: (v x r) / |r|
                    totalTangSpeed += bird.velocity.crossProd(bird.position) / rNorm;
                }

                const std::vector<Bird> neighbors = flock.flocking(static_cast<int>(idx));
                const double neighborCount = static_cast<double>(neighbors.size());
                totalNeighborCount += neighborCount;
                sumNeighborCountSquared += neighborCount * neighborCount;

                double birdEnergy = spin * spin / (2.0 * inertia);
                for (const Bird& neighbor : neighbors) {
                    birdEnergy -= interactionEnergyFactor *
                                  bird.velocity.dotProd(neighbor.velocity);
                }
                totalEnergy += birdEnergy;
                sumEnergySquared += birdEnergy * birdEnergy;
            }

            const double spinPerBird = totalSpin / birdCount;
            const double avgVelX = sumVelX / birdCount;
            const double avgVelY = sumVelY / birdCount;
            const double cmX = periodicBoundary
                ? (std::hypot(sumCosX, sumSinX) > 1e-12
                    ? std::atan2(sumSinX, sumCosX) * boxSize / M_PI
                    : 0.0)
                : sumPosX / birdCount;
            const double cmY = periodicBoundary
                ? (std::hypot(sumCosY, sumSinY) > 1e-12
                    ? std::atan2(sumSinY, sumCosY) * boxSize / M_PI
                    : 0.0)
                : sumPosY / birdCount;
            const double avgTangSpeed = totalTangSpeed / birdCount;
            const double avgRnorm = totalRnorm / birdCount;
            const double meanEnergy = totalEnergy / birdCount;
            const double meanNeighborCount = totalNeighborCount / birdCount;
            double squaredRadiusFromCm = 0.0;
            for (const Bird& bird : flock.bandada) {
                double dx = bird.position.x - cmX;
                double dy = bird.position.y - cmY;
                if (periodicBoundary) {
                    const double period = 2.0 * boxSize;
                    dx = std::remainder(dx, period);
                    dy = std::remainder(dy, period);
                }
                squaredRadiusFromCm += dx * dx + dy * dy;
            }
            const double radiusGyration =
                std::sqrt(squaredRadiusFromCm / birdCount);
            const double polarization =
                std::hypot(avgVelX, avgVelY) / maxSpeed;
            const double spinVariance = std::max(
                0.0, sumSpinSquared / birdCount - spinPerBird * spinPerBird
            );
            const double energyVariance = std::max(
                0.0, sumEnergySquared / birdCount - meanEnergy * meanEnergy
            );
            const double neighborCountVariance = std::max(
                0.0, sumNeighborCountSquared / birdCount -
                         meanNeighborCount * meanNeighborCount
            );

            metricsFile << step << "," 
                        << totalSpin << "," 
                        << spinPerBird << "," 
                        << avgVelX << "," 
                        << avgVelY << "," 
                        << avgTangSpeed << ","
                        << avgRnorm << ","
                        << cmX << ","
                        << cmY << "," 
                        << totalEnergy << ","
                        << polarization << ","
                        << radiusGyration << ","
                        << meanNeighborCount << ","
                        << neighborCountVariance << ","
                        << spinVariance << ","
                        << energyVariance << "\n";
        }

        // Save full system state snapshot every M steps
        if (step % saveSteps == 0) {
            std::string snapBirds = outputFolder + "/birds_step_" + std::to_string(step) + ".csv";
            std::string rngFile = outputFolder + "/rng_step_" + std::to_string(step) + ".bin";
            flock.saveFile(rngFile, snapBirds);
        }

        // Step the system forward, passing our modular ForceField pointer
        flock.updateFlock(timeStep, activeField);
    }

    metricsFile.close();
    std::cout << "Simulation completed successfully. Results saved to: " << outputFolder << "\n";
}
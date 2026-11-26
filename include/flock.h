#pragma once
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <unordered_map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bird.h"

// Forward declaration
class ForceField;

class Flock {
public:
    std::vector<Bird> bandada;
    double maxSpeed, inertia, friction, temperature, coupling;
    double numNeigh;
    double c0, c1, c2, exi;
    double siga, sigv, rho;
    double boxSize;
    bool periodicBoundary = false;
    unsigned long int seed;
    double fil;
    Pvector field;
    gsl_rng *m_mt;

    Flock() = default;
    Flock(unsigned long int seed, double boxSize, double dt, double maxs, double iner, double fric, double T, double J, double N, double f, bool periodic = false);
    ~Flock();

    void writeRng(char *c);
    void saveFile(char *c, std::string envs, std::string birds, double dt);
    void loadFile(FILE *rng, std::string envs, std::string birds);
    void randVecs(Bird &b);
    void boundary();
    void circularBoundary(double innerRad, double outerRad);
    std::vector<Bird> flocking(int index) const;
    void calculateForces(const ForceField* forceField);
    void initializeAccelerations();
    
    // >> UPDATE THIS SIGNATURE <<
    void updateFlock(double dt, const ForceField* forceField);
    
    double totalSpin() const;
    double totalEnergy() const;
    double totalTangSpeed() const;
    double avgRad() const;
    void setConstants(double dt);
    std::vector<double> measurements();
    std::vector<double> measurementsFerro();

private:
    mutable int cellsPerAxis = 0;
    mutable double cellWidth = 0.0;
    mutable bool cellListReady = false;
    mutable std::unordered_map<std::uint64_t, std::vector<std::size_t>> cellList;
    mutable std::vector<std::uint64_t> activeCellKeys;

    void buildCellList() const;
    void advanceFlockStep(double dt, const ForceField* forceField);
};
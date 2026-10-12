#include "forcefield.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

// Harmonic Confinement
HarmonicConfinement::HarmonicConfinement(double k) : k0(k) {}

Pvector HarmonicConfinement::calculateForce(const Pvector& pos) const {
    return Pvector(-k0 * pos.x, -k0 * pos.y);
}

// Quartic Confinement
QuarticConfinement::QuarticConfinement(double s) : strength(s) {}

Pvector QuarticConfinement::calculateForce(const Pvector& pos) const {
    double r2 = pos.dotProd(pos);
    double factor = -strength * r2;
    return Pvector(factor * pos.x, factor * pos.y);
}

// Central Swirl Force
CentralSwirlForce::CentralSwirlForce(double radS, double swS) 
    : radialStrength(radS), swirlStrength(swS) {}

Pvector CentralSwirlForce::calculateForce(const Pvector& pos) const {
    double r = pos.getNorm();
    if (r < 1e-12) return Pvector(0, 0);
    
    double fx = -radialStrength * pos.x;
    double fy = -radialStrength * pos.y;
    
    fx -= swirlStrength * (pos.y / r);
    fy += swirlStrength * (pos.x / r);
    
    return Pvector(fx, fy);
}

// Composite Force
void CompositeForce::addForce(std::shared_ptr<ForceField> f) {
    forces.push_back(f);
}

Pvector CompositeForce::calculateForce(const Pvector& pos) const {
    Pvector totalF(0, 0);
    for (const auto& f : forces) {
        Pvector fi = f->calculateForce(pos);
        totalF.addVector(fi);
    }
    return totalF;
}

// File Table Force
FileTableForce::FileTableForce(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open force table file: " + filepath);
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        double r, f_val;
        if (ss >> r >> f_val) {
            r_grid.push_back(r);
            f_grid.push_back(f_val);
        }
    }
}

Pvector FileTableForce::calculateForce(const Pvector& pos) const {
    double r = pos.getNorm();
    if (r < 1e-12 || r_grid.empty()) return Pvector(0, 0);

    double f_mag = 0.0;
    if (r <= r_grid.front()) {
        f_mag = f_grid.front();
    } else if (r >= r_grid.back()) {
        f_mag = f_grid.back();
    } else {
        for (size_t i = 0; i < r_grid.size() - 1; ++i) {
            if (r >= r_grid[i] && r <= r_grid[i + 1]) {
                double t = (r - r_grid[i]) / (r_grid[i + 1] - r_grid[i]);
                f_mag = f_grid[i] + t * (f_grid[i + 1] - f_grid[i]);
                break;
            }
        }
    }
    return Pvector(f_mag * (pos.x / r), f_mag * (pos.y / r));
}
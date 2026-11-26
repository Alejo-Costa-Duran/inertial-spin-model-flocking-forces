#include "flock.h"
#include "forcefield.h" // Include force field definitions
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace {
std::uint64_t cellKey(int x, int y) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(x)) << 32) |
           static_cast<std::uint32_t>(y);
}
}

Flock::Flock(unsigned long int seed, double boxSize, double dt, double maxSp, double iner, double fric, double T, double J, double N, double f, bool periodic)
    : maxSpeed(maxSp), inertia(iner), friction(fric), temperature(T), coupling(J),
      numNeigh(N), boxSize(boxSize), periodicBoundary(periodic), seed(seed), fil(f),
      field(0, 0), m_mt(nullptr)
{
    setConstants(dt);
    
    m_mt = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_mt, seed);
}

Flock::~Flock() {
    if (m_mt) {
        gsl_rng_free(m_mt);
    }
}

void Flock::loadFile(FILE *rng, std::string env, std::string birds) {
    gsl_rng_fread(rng, m_mt);
    std::ifstream ifs(env, std::ifstream::in);
}

void Flock::randVecs(Bird &b) {
    gsl_ran_bivariate_gaussian(m_mt, siga, sigv, rho, &b.randA.x, &b.randV.x);
    gsl_ran_bivariate_gaussian(m_mt, siga, sigv, rho, &b.randA.y, &b.randV.y);
}

void Flock::writeRng(char *c) {
    FILE *rngFile = fopen(c, "wb");
    if (rngFile) {
        gsl_rng_fwrite(rngFile, m_mt);
        fclose(rngFile);
    }
}

void Flock::saveFile(char *c, std::string envs, std::string birds, double dt) {
    writeRng(c);

    std::ofstream variables(envs);
    variables << std::fixed << std::setprecision(17) 
              << inertia << "," << seed << "," << friction << "," << temperature << "," 
              << coupling << "," << numNeigh << "," << boxSize << "," << maxSpeed << "," << dt
              << "," << (periodicBoundary ? 1 : 0);
    variables.close();

    std::ofstream data(birds);
    for (size_t idx = 0; idx < bandada.size(); ++idx) {
        const Bird &b = bandada[idx];
        data << std::setprecision(18) 
             << b.idx << "," << b.position.x << "," << b.position.y << "," 
             << b.velocity.x << "," << b.velocity.y << "," << b.acc.x << "," << b.acc.y << ","
             << b.lambda << "," << b.mu << "," << b.deltaV.x << "," << b.deltaV.y << "," <<
             b.force.x << "," << b.force.y << "," << b.randA.x << "," << b.randA.y << "," 
             << b.randV.x << "," << b.randV.y << "\n";
    }
    data.close();
}

void Flock::setConstants(double dt) {
    double xi = friction / inertia;
    double xidt = xi * dt;
    exi = std::exp(-xidt);
    
    if (xidt < 1e-3) {
        c0 = 1.0 - xidt + (xidt * xidt) / 2.0 - (xidt * xidt * xidt) / 6.0;
        c1 = 1.0 - xidt / 2.0 + (xidt * xidt) / 6.0 - (xidt * xidt * xidt) / 24.0;
        c2 = 0.5 - xidt / 6.0 + (xidt * xidt) / 24.0;
        rho = std::sqrt(3.0) * (0.5 - xidt / 16.0 - (17.0 / 1280.0) * xidt * xidt + (17.0 / 6144.0) * xidt * xidt * xidt);
    } else {
        c0 = exi;
        c1 = (1.0 - c0) / xidt;
        c2 = (1.0 - c1) / xidt;
        rho = (1.0 - exi) * (1.0 - exi) / std::sqrt((1.0 - exi * exi) * (2.0 * xidt - 3.0 + 4.0 * exi - exi * exi));
    }
    
    siga = (temperature / inertia) * (1.0 - exi * exi);
    siga = std::sqrt(siga);
    
    if (friction < 1e-3 || xidt < 1e-3) {
        sigv = temperature * dt * dt * dt * friction * (2.0 / 3.0 - 0.5 * xidt) / (inertia * inertia);
    } else {
        sigv = (temperature / friction) * (2.0 * dt - (3.0 - 4.0 * exi + exi * exi) / xi);
    }
    sigv = std::sqrt(sigv);
}

void Flock::boundary() {
    for (auto &b : bandada) {
        b.boundary(boxSize, boxSize);
    }
    cellListReady = false;
}

void Flock::circularBoundary(double innerRad, double outerRad) {
    for (auto &b : bandada) {
        double rad = b.position.getNorm();
        if (rad < innerRad || rad > outerRad) {
            double theta = 2.0 * M_PI * gsl_rng_uniform(m_mt);
            b.position.x = outerRad * std::cos(theta);
            b.position.y = outerRad * std::sin(theta);
            b.acc.x = 0.0;
            b.acc.y = 0.0;
        }
    }
    cellListReady = false;
}

double Flock::totalTangSpeed() const {
    double tang = 0.0;
    Pvector centerMass(0, 0);
    for (const auto &b : bandada) {
        centerMass.x += b.position.x;
        centerMass.y += b.position.y;
    }
    centerMass.mulScalar(1.0 / bandada.size());
    
    for (const auto &b : bandada) {
        Pvector dist(centerMass.x - b.position.x, centerMass.y - b.position.y);
        double distNorm = dist.getNorm();
        if (distNorm > 1e-12) {
            tang += b.velocity.crossProd(dist) / distNorm;
        }
    }
    return tang / bandada.size();
}

std::vector<Bird> Flock::flocking(int index) const {
    if (index < 0 || static_cast<std::size_t>(index) >= bandada.size() || numNeigh <= 0.0) {
        return {};
    }

    if (!cellListReady) {
        buildCellList();
    }

    const Bird &pajaro = bandada[index];
    std::vector<Bird> vecinos;
    if (cellWidth <= 0.0) return vecinos;

    int cellX = static_cast<int>(std::floor((pajaro.position.x + boxSize) / cellWidth));
    int cellY = static_cast<int>(std::floor((pajaro.position.y + boxSize) / cellWidth));
    if (periodicBoundary) {
        cellX = (cellX % cellsPerAxis + cellsPerAxis) % cellsPerAxis;
        cellY = (cellY % cellsPerAxis + cellsPerAxis) % cellsPerAxis;
    }

    int searchRange = std::max(1, static_cast<int>(std::ceil(numNeigh / cellWidth)));
    if (periodicBoundary) searchRange = std::min(searchRange, cellsPerAxis);

    auto addCellNeighbors = [&](int neighborCellX, int neighborCellY) {
        const auto cell = cellList.find(cellKey(neighborCellX, neighborCellY));
        if (cell == cellList.end()) return;

        for (std::size_t candidateIndex : cell->second) {
            if (candidateIndex == static_cast<std::size_t>(index)) continue;
            const Bird& candidate = bandada[candidateIndex];
            const double dist = pajaro.getDistance(candidate, boxSize, periodicBoundary);
            if (dist < numNeigh) vecinos.push_back(candidate);
        }
    };

    const bool visitEveryPeriodicCell =
        periodicBoundary && searchRange > (cellsPerAxis - 1) / 2;
    if (visitEveryPeriodicCell) {
        for (int x = 0; x < cellsPerAxis; ++x) {
            for (int y = 0; y < cellsPerAxis; ++y) {
                addCellNeighbors(x, y);
            }
        }
    } else {
        for (int dx = -searchRange; dx <= searchRange; ++dx) {
            for (int dy = -searchRange; dy <= searchRange; ++dy) {
                int neighborCellX = cellX + dx;
                int neighborCellY = cellY + dy;
                if (periodicBoundary) {
                    neighborCellX = (neighborCellX % cellsPerAxis + cellsPerAxis) % cellsPerAxis;
                    neighborCellY = (neighborCellY % cellsPerAxis + cellsPerAxis) % cellsPerAxis;
                }
                addCellNeighbors(neighborCellX, neighborCellY);
            }
        }
    }
    return vecinos;
}

void Flock::buildCellList() const {
    for (std::uint64_t key : activeCellKeys) {
        const auto cell = cellList.find(key);
        if (cell != cellList.end()) cell->second.clear();
    }
    activeCellKeys.clear();
    cellListReady = true;
    cellsPerAxis = 0;
    cellWidth = 0.0;

    if (boxSize <= 0.0 || numNeigh <= 0.0) return;

    const double boxWidth = 2.0 * boxSize;
    cellsPerAxis = std::max(1, static_cast<int>(std::floor(boxWidth / numNeigh)));
    cellWidth = boxWidth / cellsPerAxis;
    if (!periodicBoundary) cellWidth = std::max(cellWidth, numNeigh);

    const std::size_t maxRetainedCells = std::max<std::size_t>(4096, bandada.size() * 16);
    if (cellList.size() > maxRetainedCells) {
        for (auto cell = cellList.begin(); cell != cellList.end();) {
            if (cell->second.empty()) cell = cellList.erase(cell);
            else ++cell;
        }
    }

    for (std::size_t index = 0; index < bandada.size(); ++index) {
        const Bird& bird = bandada[index];
        int cellX = static_cast<int>(std::floor((bird.position.x + boxSize) / cellWidth));
        int cellY = static_cast<int>(std::floor((bird.position.y + boxSize) / cellWidth));
        if (periodicBoundary) {
            cellX = (cellX % cellsPerAxis + cellsPerAxis) % cellsPerAxis;
            cellY = (cellY % cellsPerAxis + cellsPerAxis) % cellsPerAxis;
        }
        const std::uint64_t key = cellKey(cellX, cellY);
        auto& cell = cellList[key];
        if (cell.empty()) activeCellKeys.push_back(key);
        cell.push_back(index);
    }
}

void Flock::calculateForces(const ForceField* forceField) {
    buildCellList();
    for (std::size_t index = 0; index < bandada.size(); ++index) {
        std::vector<Bird> vec = flocking(static_cast<int>(index));
        bandada[index].calcForce(inertia, coupling, vec, forceField);
    }
}

void Flock::initializeAccelerations() {
    const double speedSquared = maxSpeed * maxSpeed;
    for (auto& bird : bandada) {
        if (bird.fixed) continue;

        const double parallelForce =
            (bird.velocity.x * bird.force.x + bird.velocity.y * bird.force.y) / speedSquared;
        bird.acc.x = bird.force.x - parallelForce * bird.velocity.x;
        bird.acc.y = bird.force.y - parallelForce * bird.velocity.y;
    }
}

double Flock::avgRad() const {
    double rad = 0.0;
    for (const auto &b : bandada) {
        rad += std::sqrt(b.position.x * b.position.x + b.position.y * b.position.y);
    }
    return rad / bandada.size();
}

void Flock::advanceFlockStep(double dt, const ForceField* forceField) {
    cellListReady = false;

    // Solve every velocity constraint before mutating positions or velocities.
    for (auto &b : bandada) {
        if (!b.fixed) {
            randVecs(b);

            b.deltaV.x = c1 * dt * b.acc.x + c2 * dt * dt * b.force.x + b.randV.x;
            b.deltaV.y = c1 * dt * b.acc.y + c2 * dt * dt * b.force.y + b.randV.y;

            double b0 = 2.0 * (b.velocity.x * b.deltaV.x + b.deltaV.y * b.velocity.y);
            double a = maxSpeed * maxSpeed;
            if (a <= 0.0) {
                throw std::runtime_error("Maximum speed must be positive for the velocity constraint.");
            }
            double c = b.deltaV.x * b.deltaV.x + b.deltaV.y * b.deltaV.y - a;
            double w;
            double discriminant = b0 * b0 - 4.0 * a * c;
            const double discriminantScale = b0 * b0 + 4.0 * a * std::abs(c);
            const double discriminantTolerance =
                1e-12 * std::max(1.0, discriminantScale);
            if (!std::isfinite(discriminant) || discriminant < -discriminantTolerance) {
                throw std::runtime_error(
                    "Velocity constraint remained infeasible after adaptive step sizing for bird " +
                    std::to_string(b.idx) + "."
                );
            }
            discriminant = std::max(0.0, discriminant);
            
            if (std::abs(b0) < 1e-12) {
                w = std::sqrt(std::max(0.0, -c / a));
            } else {
                double sgnb = (b0 > 0) ? 1.0 : -1.0;
                double q = -0.5 * (b0 + sgnb * std::sqrt(discriminant));
                w = (b0 > 0) ? (c / q) : (q / a);
            }
            const double lambdaDenominator = c2 * dt * dt;
            if (!std::isfinite(w) || !std::isfinite(lambdaDenominator) ||
                lambdaDenominator == 0.0) {
                throw std::runtime_error("Adaptive integration produced an invalid constraint multiplier.");
            }
            b.lambda = (w - 1.0) / lambdaDenominator;
            if (!std::isfinite(b.lambda)) {
                throw std::runtime_error("Adaptive integration produced a non-finite constraint multiplier.");
            }
        }
    }

    // 1. Apply the accepted constrained update.
    for (auto &b : bandada) {
        if (!b.fixed) {
            b.position.x += b.velocity.x * dt;
            b.position.y += b.velocity.y * dt;
            b.acc.x = c0 * b.acc.x + (c1 - c2) * dt * (b.force.x + b.lambda * b.velocity.x) + b.randA.x;
            b.acc.y = c0 * b.acc.y + (c1 - c2) * dt * (b.force.y + b.lambda * b.velocity.y) + b.randA.y;

            b.velocity.x = (1.0 + c2 * dt * dt * b.lambda) * b.velocity.x + b.deltaV.x;
            b.velocity.y = (1.0 + c2 * dt * dt * b.lambda) * b.velocity.y + b.deltaV.y;
        }
    }

    if (periodicBoundary) boundary();

    // 2. Update interaction and positional forces using the modular ForceField pointer
    calculateForces(forceField);

    // 3. Final acceleration updates and second constraint multiplier (mu)
    for (auto &b : bandada) {
        if (!b.fixed) {
            b.acc.x += c2 * dt * b.force.x;
            b.acc.y += c2 * dt * b.force.y;

            b.mu = -b.velocity.dotProd(b.acc) / (maxSpeed * maxSpeed);

            b.acc.x += b.velocity.x * b.mu;
            b.acc.y += b.velocity.y * b.mu;
        }
    }
}

void Flock::updateFlock(double dt, const ForceField* forceField) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::invalid_argument("Time step must be finite and positive.");
    }
    if (!std::isfinite(maxSpeed) || maxSpeed <= 0.0) {
        throw std::invalid_argument("Maximum speed must be finite and positive.");
    }

    constexpr std::size_t maxSubsteps = 100000;
    const double safeVelocityIncrement = 0.5 * maxSpeed;
    const double minimumStep = std::numeric_limits<double>::epsilon() * dt;
    double remaining = dt;
    std::size_t acceptedSubsteps = 0;

    try {
        calculateForces(forceField);

        while (remaining > 0.0) {
            if (acceptedSubsteps >= maxSubsteps) {
                throw std::runtime_error(
                    "Adaptive integration exceeded 100000 substeps in one time step."
                );
            }

            double substep = remaining;
            while (true) {
                if (substep <= minimumStep) {
                    throw std::runtime_error(
                        "Adaptive integration reached the floating-point time-step limit."
                    );
                }

                setConstants(substep);
                double maximumPredictedIncrement = 0.0;
                for (const auto& bird : bandada) {
                    if (bird.fixed) continue;

                    const double accelerationIncrement =
                        std::abs(c1 * substep) * bird.acc.getNorm();
                    const double forceIncrement =
                        std::abs(c2 * substep * substep) * bird.force.getNorm();
                    const double noiseIncrement =
                        8.0 * std::sqrt(2.0) * std::abs(sigv);
                    maximumPredictedIncrement = std::max(
                        maximumPredictedIncrement,
                        accelerationIncrement + forceIncrement + noiseIncrement
                    );
                }
                if (!std::isfinite(maximumPredictedIncrement)) {
                    throw std::runtime_error(
                        "Non-finite force or acceleration encountered during adaptive integration."
                    );
                }

                if (maximumPredictedIncrement > safeVelocityIncrement) {
                    substep *= 0.5;
                    continue;
                }
                break;
            }

            advanceFlockStep(substep, forceField);
            remaining -= substep;
            ++acceptedSubsteps;
        }
    } catch (...) {
        setConstants(dt);
        throw;
    }

    setConstants(dt);
}

double Flock::totalSpin() const {
    double totSpin = 0.0;
    for (const auto &b : bandada) {
        totSpin += b.getSpin(inertia, maxSpeed);
    }
    return totSpin;
}

double Flock::totalEnergy() const {
    double energy = 0.0;
    for (size_t index = 0; index < bandada.size(); ++index) {
        energy += std::pow(bandada[index].getSpin(inertia, maxSpeed), 2) / (2.0 * inertia);
        std::vector<Bird> vec = flocking(index);
        for (const auto &b : vec) {
            energy += -(coupling / (2.0 * maxSpeed * maxSpeed)) * bandada[index].velocity.dotProd(b.velocity);
        }
    }
    return energy;
}

std::vector<double> Flock::measurements() {
    int N = bandada.size();
    double spin = 0.0, spinSq = 0.0, energy = 0.0, energySq = 0.0;
    double tangSpeed = 0.0, tangSpeedSq = 0.0;
    double vecX = 0.0, vecSqX = 0.0, vecY = 0.0, vecSqY = 0.0;
    
    for (int idx = 0; idx < N; ++idx) {
        const Bird &pajaro = bandada[idx];
        double s = pajaro.getSpin(inertia, maxSpeed);
        spin += s;
        spinSq += s * s;
        
        double vx = pajaro.velocity.x;
        double vy = pajaro.velocity.y;
        vecX += vx;
        vecY += vy;
        vecSqX += vx * vx;
        vecSqY += vy * vy;
        
        double normPos = pajaro.position.getNorm();
        double tsp = (normPos > 1e-12) ? (pajaro.velocity.crossProd(pajaro.position) / normPos) : 0.0;
        tangSpeed += tsp;
        tangSpeedSq += tsp * tsp;
        
        double e = (s * s) / (2.0 * inertia);
        std::vector<Bird> vec = flocking(idx);
        for (const auto &b : vec) { 
            e += -(coupling / (2.0 * maxSpeed * maxSpeed)) * pajaro.velocity.dotProd(b.velocity); 
        }
        energy += e;
        energySq += e * e;
    }
    return {spin / N, spinSq / N, energy / N, energySq / N, tangSpeed / N, tangSpeedSq / N, vecX / N, vecSqX / N, vecY / N, vecSqY / N};
}
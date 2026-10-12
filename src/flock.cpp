#include "flock.h"
#include "forcefield.h" // Include force field definitions
#include <algorithm>
#include <cmath>
#include <iomanip>

Flock::Flock(unsigned long int seed, double boxSize, double dt, double maxSp, double iner, double fric, double T, double J, double N, double f)
    : seed(seed), maxSpeed(maxSp), inertia(iner), friction(fric), temperature(T), coupling(J), 
      numNeigh(N), boxSize(boxSize), fil(f), field(0, 0)
{
    cells = (boxSize * 2.0) / std::floor(2.0 * boxSize);
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
              << coupling << "," << numNeigh << "," << boxSize << "," << maxSpeed << "," << dt;
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
        rho = (1.0 - exi) * (1.0 - exi) / std::sqrt((1.0 - exi * exi) * (3.0 * xidt - 3.0 + 4.0 * exi - exi * exi));
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
    const Bird &pajaro = bandada[index];
    std::vector<Bird> vecinos;
    int box = static_cast<int>(2.0 * boxSize - 1.0);
    
    for (const auto &b : bandada) {
        if ((std::abs(b.boxX - pajaro.boxX) % box < 3) && (std::abs(b.boxY - pajaro.boxY) % box < 3)) {
            double dist = pajaro.getDistance(b, boxSize);
            if (dist < numNeigh && dist != 0.0) {
                vecinos.push_back(b);
            }
        }
    }
    return vecinos;
}

double Flock::avgRad() const {
    double rad = 0.0;
    for (const auto &b : bandada) {
        rad += std::sqrt(b.position.x * b.position.x + b.position.y * b.position.y);
    }
    return rad / bandada.size();
}

void Flock::updateFlock(double dt, const ForceField* forceField) {
    // 1. Partial updates, noise generation, and constraint multipliers (lambda)
    for (auto &b : bandada) {
        if (!b.fixed) {
            randVecs(b);

            b.position.x += b.velocity.x * dt;
            b.position.y += b.velocity.y * dt;

            b.boxX = static_cast<int>(std::floor((b.position.x + boxSize) / cells));
            b.boxY = static_cast<int>(std::floor((b.position.y + boxSize) / cells));

            b.deltaV.x = c1 * dt * b.acc.x + c2 * dt * dt * b.force.x + b.randV.x;
            b.deltaV.y = c1 * dt * b.acc.y + c2 * dt * dt * b.force.y + b.randV.y;

            double b0 = 2.0 * (b.velocity.x * b.deltaV.x + b.deltaV.y * b.velocity.y);
            double a = maxSpeed * maxSpeed;
            double c = b.deltaV.x * b.deltaV.x + b.deltaV.y * b.deltaV.y - a;
            double w;
            
            if (std::abs(b0) < 1e-12) {
                w = std::sqrt(-c / a);
            } else {
                double sgnb = (b0 > 0) ? 1.0 : -1.0;
                double q = -0.5 * (b0 + sgnb * std::sqrt(b0 * b0 - 4.0 * a * c));
                w = (b0 > 0) ? (c / q) : (q / a);
            }
            b.lambda = (w - 1.0) / (c2 * dt * dt);

            b.acc.x = c0 * b.acc.x + (c1 - c2) * dt * (b.force.x + b.lambda * b.velocity.x) + b.randA.x;
            b.acc.y = c0 * b.acc.y + (c1 - c2) * dt * (b.force.y + b.lambda * b.velocity.y) + b.randA.y;

            b.velocity.x = (1.0 + c2 * dt * dt * b.lambda) * b.velocity.x + b.deltaV.x;
            b.velocity.y = (1.0 + c2 * dt * dt * b.lambda) * b.velocity.y + b.deltaV.y;
        }
    }

    // 2. Update interaction and positional forces using the modular ForceField pointer
    for (auto &b : bandada) {
        std::vector<Bird> vec = flocking(b.idx);
        b.calcForce(inertia, coupling, vec, forceField);
    }

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
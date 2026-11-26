#include "flock.h"
#include "forcefield.h"

#include <cmath>
#include <exception>
#include <iostream>

int main() {
    constexpr double timeStep = 0.001;
    constexpr double maxSpeed = 1.0;

    Flock freeFlock(321, 20.0, timeStep, maxSpeed, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    freeFlock.bandada.emplace_back(
        false, 0, Pvector(0.0, 0.0), Pvector(maxSpeed, 0.0), Pvector(0.0, 0.0)
    );
    freeFlock.updateFlock(timeStep, nullptr);
    if (std::abs(freeFlock.bandada.front().position.x - maxSpeed * timeStep) > 1e-12) {
        std::cerr << "Adaptive integration did not advance the full requested time step.\n";
        return 1;
    }

    Flock flock(123, 20.0, timeStep, maxSpeed, 0.0001, 1.0, 0.0, 0.0, 1.0, 0.0);
    flock.bandada.emplace_back(
        false, 0, Pvector(8.0, 0.0), Pvector(0.0, maxSpeed), Pvector(0.0, 0.0)
    );

    GeneralizedPowerConfinement strongWell(100.0, 8.0, 6);
    flock.calculateForces(&strongWell);
    flock.initializeAccelerations();

    try {
        flock.updateFlock(timeStep, &strongWell);
    } catch (const std::exception& error) {
        std::cerr << "Adaptive step failed: " << error.what() << '\n';
        return 1;
    }

    const Bird& bird = flock.bandada.front();
    const double speed = bird.velocity.getNorm();
    if (!std::isfinite(bird.position.x) || !std::isfinite(bird.position.y) ||
        !std::isfinite(speed) || std::abs(speed - maxSpeed) > 1e-8) {
        std::cerr << "Adaptive step produced a non-finite state or violated fixed speed.\n";
        return 1;
    }

    return 0;
}

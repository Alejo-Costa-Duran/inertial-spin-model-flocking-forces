#pragma once
#include "pvector.h"
#include <vector>

// Forward declaration to avoid circular dependencies
class ForceField;

class Bird {
public:
    // ... (keep all other attributes and constructors as they are)
    Pvector position, velocity, acc, force, deltaV;
    int idx, boxX, boxY;
    double lambda, mu;
    Pvector randA, randV;
    bool fixed;

    Bird() = default;
    Bird(bool fix, int idx, Pvector position, Pvector velocity, Pvector acc, double inertia, double maxSpeed, double length);

    // ... (keep other methods unchanged)
    void updatePosition(double dt);
    void updateVelocity(double dt, double c1, double c2);
    void calcLambda(double c1, double c2, double dt, double maxSpeed);
    void calcMu(double c2, double dt, double maxSpeed);
    void boundary(double lengthX, double lengthY);
    double getDistance(const Bird& vecino, double length) const;
    double getSpin(double inertia, double maxSpeed) const;
    Pvector polarToScreen(double width, double height, double lengthX, double lengthY) const;

    // >> UPDATE THIS SIGNATURE <<
    // Replaced 'double fieldStrength' with 'const ForceField* forceField'
    void calcForce(double inertia, double coupling, const std::vector<Bird>& vecinos, const ForceField* forceField);
};
#include "bird.h"
#include "forcefield.h" // Include the new force field definitions
#include <cmath>
Bird::Bird(bool fix, int idx, Pvector position, Pvector velocity, Pvector acc, double inertia, double maxSpeed, double length)
    : fixed(fix), idx(idx), position(position), velocity(velocity), acc(acc), 
      boxX(0), boxY(0), lambda(0.0), mu(0.0) 
{
    boundary(length, length);
    boxX = static_cast<int>(std::floor((position.x + length) * std::floor(2 * length) / (2 * length)));
    boxY = static_cast<int>(std::floor((position.y + length) * std::floor(2 * length) / (2 * length)));
    deltaV = Pvector(0, 0);
    randA = Pvector(0, 0);
    randV = Pvector(0, 0);
}

Pvector Bird::polarToScreen(double width, double height, double lengthX, double lengthY) const {
    Pvector screenPosition;
    screenPosition.x = (width / (2 * lengthX)) * (lengthX + position.x);
    screenPosition.y = (height / (2 * lengthY)) * (lengthY - position.y);
    return screenPosition;
}

void Bird::calcLambda(double c1, double c2, double dt, double maxSpeed) {
    Pvector delV = Pvector(dt * c1 * acc.x + c2 * dt * dt * force.x + randV.x, 
                           dt * c1 * acc.y + c2 * dt * dt * force.y + randV.y);
    double b = 2.0 * (velocity.x * delV.x + delV.y * delV.y);
    double a = maxSpeed * maxSpeed;
    double c = delV.x * delV.x + delV.y * delV.y - a;
    double w;
    
    if (std::abs(b) < 1e-12) { 
        w = std::sqrt(-c / a); 
    } else {
        double sgnb = (b > 0) ? 1.0 : -1.0;
        double q = -0.5 * (b + sgnb * std::sqrt(b * b - 4 * a * c));
        w = (b > 0) ? (c / q) : (q / a);
    }
    lambda = (w - 1.0) / (c2 * dt * dt);
}

void Bird::updatePosition(double dt) {
    position.x += dt * velocity.x;
    position.y += dt * velocity.y;
}

void Bird::calcMu(double c2, double dt, double maxSpeed) {
    mu = -(velocity.x * acc.x + velocity.y * acc.y) / (maxSpeed * maxSpeed * dt * c2);
}

void Bird::updateVelocity(double dt, double c1, double c2) {
    velocity.x += c1 * dt * acc.x + dt * dt * c2 * force.x + dt * dt * c2 * lambda * velocity.x + randV.x;
    velocity.y += c1 * dt * acc.y + dt * dt * c2 * force.y + dt * dt * c2 * lambda * velocity.y + randV.y;
}

void Bird::calcForce(double inertia, double coupling, const std::vector<Bird>& vecinos, const ForceField* forceField) {
    double currentFx = 0.0;
    double currentFy = 0.0;
    
    // 1. Vicsek-style velocity alignment/imitation interaction
    for (const auto& bird : vecinos) {
        currentFx += bird.velocity.x;
        currentFy += bird.velocity.y;
    }

    // Base force from neighbor coupling (F_coupling / inertia)
    force.x = (currentFx * coupling) / inertia;
    force.y = (currentFy * coupling) / inertia;

    // 2. Add external position-dependent F forces if a force field is provided
    if (forceField != nullptr) {
        Pvector extForce = forceField->calculateForce(position);
        
        // F forces act on acceleration/velocity dynamics (divided by inertia)
        force.x += extForce.x / inertia;
        force.y += extForce.y / inertia;
    }
}

double Bird::getSpin(double inertia, double maxSpeed) const {
    return (inertia / (maxSpeed * maxSpeed)) * (velocity.x * acc.y - velocity.y * acc.x);
}

void Bird::boundary(double lengthX, double lengthY) {
    if (position.x > lengthX)  { position.x -= 2 * lengthX; }
    if (position.y > lengthY)  { position.y -= 2 * lengthY; }
    if (position.x < -lengthX) { position.x += 2 * lengthX; }
    if (position.y < -lengthY) { position.y += 2 * lengthY; }
}

double Bird::getDistance(const Bird& vecino, double length) const {
    double deltaX = std::abs(position.x - vecino.position.x);
    if (deltaX > length) { deltaX = deltaX - 2 * length; }
    double deltaY = std::abs(position.y - vecino.position.y);
    if (deltaY > length) { deltaY = deltaY - 2 * length; }
    return std::sqrt(deltaX * deltaX + deltaY * deltaY);
}
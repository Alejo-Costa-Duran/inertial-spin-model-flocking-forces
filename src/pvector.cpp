#include "pvector.h"

Pvector::Pvector() : x(0.0), y(0.0) {}

Pvector::Pvector(double x_val, double y_val) : x(x_val), y(y_val) {}

double Pvector::getNorm() const {
    return std::sqrt(x * x + y * y);
}

void Pvector::normalize(double targetNorm) {
    double currentNorm = getNorm();
    if (currentNorm > 1e-12) {
        x = targetNorm * (x / currentNorm);
        y = targetNorm * (y / currentNorm);
    }
}

void Pvector::addVector(const Pvector& v) {
    x += v.x;
    y += v.y;
}

void Pvector::subVector(const Pvector& v) {
    x -= v.x;
    y -= v.y;
}

void Pvector::mulScalar(double s) {
    x *= s;
    y *= s;
}

double Pvector::dotProd(const Pvector& v) const {
    return x * v.x + y * v.y;
}

double Pvector::crossProd(const Pvector& v) const {
    return x * v.y - y * v.x; // Correct 2D scalar cross product magnitude
}
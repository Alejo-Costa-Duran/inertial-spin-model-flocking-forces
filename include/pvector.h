#pragma once
#include <cmath>

class Pvector {
public:
    double x, y;

    // Constructors
    Pvector();
    Pvector(double x_val, double y_val);

    // Vector operations
    double getNorm() const;
    void normalize(double targetNorm = 1.0);
    void addVector(const Pvector& v);
    void subVector(const Pvector& v);
    void mulScalar(double s);
    double dotProd(const Pvector& v) const;
    double crossProd(const Pvector& v) const;
};
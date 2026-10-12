#pragma once
#include "pvector.h"
#include <vector>
#include <memory>
#include <string>

class ForceField {
public:
    virtual ~ForceField() = default;
    virtual Pvector calculateForce(const Pvector& pos) const = 0;
};

class HarmonicConfinement : public ForceField {
    double k0;
public:
    explicit HarmonicConfinement(double k);
    Pvector calculateForce(const Pvector& pos) const override;
};

class QuarticConfinement : public ForceField {
    double strength;
public:
    explicit QuarticConfinement(double s);
    Pvector calculateForce(const Pvector& pos) const override;
};

class CentralSwirlForce : public ForceField {
    double radialStrength;
    double swirlStrength;
public:
    CentralSwirlForce(double radS, double swS);
    Pvector calculateForce(const Pvector& pos) const override;
};

class CompositeForce : public ForceField {
    std::vector<std::shared_ptr<ForceField>> forces;
public:
    void addForce(std::shared_ptr<ForceField> f);
    Pvector calculateForce(const Pvector& pos) const override;
};

class FileTableForce : public ForceField {
    std::vector<double> r_grid;
    std::vector<double> f_grid;
public:
    explicit FileTableForce(const std::string& filepath);
    Pvector calculateForce(const Pvector& pos) const override;
};
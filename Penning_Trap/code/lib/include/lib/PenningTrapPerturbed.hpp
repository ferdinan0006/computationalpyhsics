#pragma once

#include "lib/PenningTrap.hpp"
#include <armadillo>

class PenningTrapPerturbed : public PenningTrap {
public:
    PenningTrapPerturbed(
        ISolver& solver, double B0_in, double V0_in, double d_in, double f, double w_V);

protected:
    virtual arma::vec3 electric_field(const arma::vec3& r, double t) override;

    double f;
    double w_V;
};

#include "lib/PenningTrapPerturbed.hpp"
#include "lib/PenningTrap.hpp"

PenningTrapPerturbed::PenningTrapPerturbed(
    ISolver& solver, double B0_in, double V0_in, double d_in, double f, double w_V)
    : PenningTrap(solver, B0_in, V0_in, d_in)
    , f(f)
    , w_V(w_V)
{
}

arma::vec3 PenningTrapPerturbed::electric_field(const arma::vec3& r, double t)
{
    return PenningTrap::electric_field(r, 0) * (1 + f * std::cos(w_V * t));
}

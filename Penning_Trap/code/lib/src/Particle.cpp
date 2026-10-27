#include "lib/Particle.hpp"
#include <armadillo>
#include <string>
#include <vector>

// Constructor
Particle::Particle(double q, double m, arma::vec3 pos, arma::vec3 vel, double t0)
    : q(q)
    , m(m)
    , position(pos)
    , r0(pos)
    , velocity(vel)
    , v0(vel)
{
    // storing initial position as first entry in trajectory history
    save_state(t0);
}

// resetting particle information
void Particle::reset(double t)
{
    position = r0;
    velocity = v0;
    r_history.clear();
    save_state(t);
}

// save current state to r_history
void Particle::save_state(double time)
{
    arma::vec timevec = { time };
    arma::vec pos_time = join_cols(position, velocity, timevec);
    r_history.push_back(pos_time);
}

// Save the full stored trajectory as a matrix: rows = time steps, cols = x,y,z
void Particle::save_history(const std::string& filename)
{
    if (r_history.empty()) {
        return;
    }

    arma::mat traj(r_history.size(), 7);
    for (std::size_t i = 0; i < r_history.size(); ++i) {
        traj.row(static_cast<arma::uword>(i)) = r_history[i].t();
    }

    traj.save(filename + ".txt", arma::raw_ascii);
}

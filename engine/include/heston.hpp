#pragma once

#include <vector>
#include <cstddef>

struct HestonParams {
    double S0    = 100.0;
    double v0    = 0.04;
    double kappa = 2.0;
    double theta = 0.04;
    double xi    = 0.3;
    double rho   = -0.7;
    double r     = 0.02;
    double q     = 0.0;
};

using Path     = std::vector<double>;
using Paths    = std::vector<Path>;
using TwoPaths = std::pair<Paths, Path>; 

class HestonSimulator {
    public: 
        explicit HestonSimulator(const HestonParams& p);

        TwoPaths simulate(std::size_t n_paths, std::size_t n_steps, double T = 1.0) const;
    
    private:
        HestonParams params_;
        double dt_;
};
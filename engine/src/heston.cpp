#include "../include/heston.hpp"
#include <random>
#include <cmath>

HestonSimulator::HestonSimulator(const HestonParams& p)
    : params_(p), dt_(0.0)
{
    TwoPaths HestonSimulator::simulate(std::size_t n_paths, std::size_t n_steps, double T = 1.0) const {
        
    }
}
<
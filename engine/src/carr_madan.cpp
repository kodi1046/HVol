#include "../include/carr_madan.hpp"
// #include <cmath>
#include <complex>

static constexpr std::complex<double> i(0, 1);

std::complex<double> CarrMadanModel::Pricer::d(double u) const 
{
    return std::sqrt(std::pow(rho_ * xi_ * i * u - kappa_, 2) + xi_ * xi_ * (i * u + u * u));
}


std::complex<double> CarrMadanModel::Pricer::g(double u) const 
{
    return (kappa_ - rho_ * xi_ * i * u - d(u)) / (kappa_ - rho_ * xi_ * i * u + d(u));
}

std::complex<double> CarrMadanModel::Pricer::C(double u) const {
    std::complex<double> t1 = r_ * i * u * tau_;
    std::complex<double> t2 = (kappa_ * theta_) / (xi_ * xi_);
    std::complex<double> t3 = (kappa_ - rho_ * xi_ * i * u) * tau_;
    std::complex<double> t4 = (std::complex<double>(1) - g(u) * std::exp(-d(u) * tau_)) / (std::complex<double>(1)  - g(u));
    return t1 + t2 * (t3 - std::complex<double>(2) * std::log(t4));
}

std::complex<double> CarrMadanModel::Pricer::D(double u) const {
    std::complex<double> t1 = (kappa_ - rho_ * xi_ * i * u - d(u)) / (xi_ * xi_);
    std::complex<double> t2 = (std::complex<double>(1) - std::exp(-d(u) * tau_)) / (std::complex<double>(1) - g(u) * std::exp(-d(u) * tau_));
    return t1 * t2;
}

std::complex<double> CarrMadanModel::Pricer::characteristic_function(double u) const {
    return std::exp(C(u) + D(u) * v0_ + i * u * std::log(S0_));
}

double CarrMadanModel::Pricer::c(double k) const {
    return 1.1;
}

double CarrMadanModel::Pricer::carr_madan_transform(double v) const {
    return 1.1;
}

double CarrMadanModel::Pricer::operator()() const {
    return 1.1;
}

double CarrMadanModel::price(
    const Option& option,
    double S0,
    std::optional<double> exercise_t,
    double alpha) const 
{
    exercise_t = option.maturity_;
    return Pricer{*this, option, S0, exercise_t, alpha}();
}




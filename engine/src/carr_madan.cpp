#include "../include/carr_madan.hpp"
#include "../include/fft.hpp"
// #include <cmath>
#include <complex>
#include <vector>

using complex_t = CarrMadanModel::complex_t;

static constexpr complex_t i(0, 1);

complex_t CarrMadanModel::Pricer::d(complex_t u) const 
{
    return std::sqrt(std::pow(rho_ * xi_ * i * u - kappa_, 2) + xi_ * xi_ * (i * u + u * u));
}

complex_t CarrMadanModel::Pricer::g(complex_t u) const 
{
    return (kappa_ - rho_ * xi_ * i * u - d(u)) / (kappa_ - rho_ * xi_ * i * u + d(u));
}

complex_t CarrMadanModel::Pricer::C(complex_t u) const {
    complex_t t1 = r_ * i * u * tau_;
    complex_t t2 = (kappa_ * theta_) / (xi_ * xi_);
    complex_t t3 = (kappa_ - rho_ * xi_ * i * u) * tau_;
    complex_t t4 = (complex_t(1) - g(u) * std::exp(-d(u) * tau_)) / (complex_t(1)  - g(u));
    return t1 + t2 * (t3 - complex_t(2) * std::log(t4));
}

complex_t CarrMadanModel::Pricer::D(complex_t u) const {
    complex_t t1 = (kappa_ - rho_ * xi_ * i * u - d(u)) / (xi_ * xi_);
    complex_t t2 = (complex_t(1) - std::exp(-d(u) * tau_)) / (complex_t(1) - g(u) * std::exp(-d(u) * tau_));
    return t1 * t2;
}

complex_t CarrMadanModel::Pricer::characteristic_function(complex_t u) const {
    return std::exp(C(u) + D(u) * v0_ + i * u * std::log(S0_));
}

complex_t CarrMadanModel::Pricer::c(complex_t k) const {
    return std::exp(alpha_ * k) * C(k);
}

complex_t CarrMadanModel::Pricer::carr_madan_transform(complex_t u) const {
    complex_t t1 = (std::exp(-r_ * tau_) * characteristic_function(complex_t(u) - (alpha_ + 1) * i));
    complex_t t2 = (alpha_ * alpha_ + alpha_ - u * u + i * (complex_t(2) * alpha_ + complex_t(1)) * u);
    return t1 / t2;
}

double CarrMadanModel::Pricer::carr_madan_inverse_transform(complex_t k) const {

}

double CarrMadanModel::Pricer::operator()() const {
    return 0;
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




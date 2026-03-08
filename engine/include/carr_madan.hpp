#include "base.hpp"
#include <complex>

namespace {
    inline constexpr double DEFAULT_ALPHA = 1.5;
}

class CarrMadanModel : public HestonModel {
    public:
        CarrMadanModel(
            double kappa,
            double theta,
            double xi,
            double rho,
            double v0,
            double r)
            : HestonModel(kappa, theta, xi, rho, v0, r)
        {}
        
        double price(
            const Option& option,
            double S0,
            std::optional<double> exercise_t = std::nullopt,
            double alpha = DEFAULT_ALPHA) const;
        
        using complex_t = std::complex<double>;
    
    private:
        struct Pricer {
            const CarrMadanModel& model_;
            const double& kappa_;
            const double& theta_;
            const double& xi_;
            const double& rho_;
            const double& v0_;
            const double& r_;
            
            const Option& option_;
            double alpha_;
            double S0_;
            double tau_;
            
            Pricer(
                const CarrMadanModel& model,
                const Option& option,
                double S0,
                std::optional<double> tau = std::nullopt,
                double alpha = DEFAULT_ALPHA)
                : model_(model), option_(option), S0_(S0), tau_(tau.value()), alpha_(alpha),
                  kappa_(model.kappa_), 
                  theta_(model.theta_), 
                  xi_(model.xi_), 
                  rho_(model.rho_), 
                  v0_(model.v0_), 
                  r_(model.r_)
            {}
        
            
            complex_t g(complex_t u) const;
            complex_t d(complex_t u) const;
            complex_t C(complex_t u) const;
            complex_t D(complex_t u) const;
            complex_t characteristic_function(complex_t u) const;
            complex_t c(complex_t k) const;
            complex_t carr_madan_transform(complex_t u) const; 
            double carr_madan_inverse_transform(complex_t k) const;
            double operator()() const;     
        };     
};
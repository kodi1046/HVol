#include "base.hpp"
#include <cmath>

class VanillaOption : public Option {
    public:
        VanillaOption(
            double strike,
            double maturity,
            OptionType option_type
        )
            : Option(strike, maturity, option_type)
        {}
        
        double payoff(double S0) const override {
            return (option_type_ == OptionType::Call) ? std::max(S0 - strike_, 0.0)
                                                      : std::max(strike_ - S0, 0.0);
        }
};

class EuropeanOption : public VanillaOption {
    public:
        using VanillaOption::VanillaOption;
}



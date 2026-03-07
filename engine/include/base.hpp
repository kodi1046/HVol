#pragma once
#include <optional>

class HestonModel;
class Option;

enum OptionType {
    Call,
    Put,
};  

class HestonModel {
    public: 
        virtual ~HestonModel() = default;
        
        HestonModel(
            double kappa,
            double theta,
            double xi,
            double rho,
            double v0,
            double r)
            : kappa_(kappa), theta_(theta), xi_(xi), rho_(rho), v0_(v0), r_(r) 
        {}
        
        virtual double price(
            const Option& option,
            double S0,
            std::optional<double> exercise_t = std::nullopt
        ) const = 0;
    
    protected:
        double kappa_, theta_, xi_, rho_, v0_, r_;
};

class Option {
    public:
        Option(
            double strike,
            double maturity,
            OptionType option_type
        )
            : strike_(strike), maturity_(maturity), option_type_(option_type)
        {};
 
        virtual ~Option() = default;
        
        virtual double payoff(double S0) const = 0;

    
    const double strike_, maturity_;
    const OptionType option_type_;
};



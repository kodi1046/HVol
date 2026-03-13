#include "base.hpp"
#include <vector>

namespace {
    inline constexpr double DEFAULT_WEIGHT = 1.0;
}

namespace calibration {
    
    struct CalibrationTarget {
        Option& option;
        double weight = DEFAULT_WEIGHT;
    };

    class HestonCalibration {
        public:
            virtual ~HestonCalibration() = default;
            
            HestonCalibration(
                HestonModel& model
            )
                : model_(model)
            {}
            
            void calibrate();
        
        private:
            HestonModel& model_;
            std::vector<CalibrationTarget> targets_;
            
            double cost_function(const std::vector<double>& params) const;
            
            bool is_valid(const std::vector<double>& params) const;
        };

}
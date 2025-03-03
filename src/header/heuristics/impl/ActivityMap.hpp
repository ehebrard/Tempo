//
// Created by tluchterha on 22/11/22.
//

#ifndef SCHEDCL_ACTIVITYMAP_HPP
#define SCHEDCL_ACTIVITYMAP_HPP

#include <concepts>

#include "util/traits.hpp"


namespace tempo::heuristics::impl {
/**
 * @brief Class that can be used to record activity of variables
 */
class ActivityMap : public std::vector<double> {
public:
    /**
     * CTor. Initializes activity of all literals with baseIncrement.
     */
    ActivityMap(const double d=.999) : decay(d) {}
    
    bool incrementActivity(const var_t x) noexcept {
        bool need_rescaling{false};
        this->operator[](x) += increment;
        need_rescaling = this->operator[](x) > maxActivity;
        return need_rescaling;
    }
    
    template<typename Iterable>
    void update(const Iterable& vars) noexcept {
      bool should_normalize = false;
      for (const auto x : vars) {
        should_normalize |= incrementActivity(x);
      }
        // protect against overflow
        if (should_normalize) {
          normalize();
        } else {
          applyDecay();
        }
    }
    void applyDecay() { increment /= decay; }

    void normalize() {
        
        

      // do not consider constants for normalization
      auto [lp, up] =
          std::ranges::minmax_element(this->begin() + 1, this->end());
      double l{std::numeric_limits<double>::max()};
      if (lp != this->end())
        l = *lp;

      double u{-std::numeric_limits<double>::max()};
      if (up != this->end())
        u = *up;
        
        
//        std::cout << "normalize weights [" << l << "," << u << "]\n";

      auto factor{baseGap / (u - l)};
      for (auto a{this->begin() + 1}; a != this->end(); ++a) {
          
//          std::cout << *a << " -> ";
        *a = baseIncrement + (*a - l) * factor;
//          std::cout << *a << "\n";
      }
      increment = baseIncrement;
        
        
//        auto [lpa, upa] =
//            std::ranges::minmax_element(this->begin() + 1, this->end());
//        l = std::numeric_limits<double>::max();
//        if (lpa != this->end())
//          l = *lpa;
//
//        u = -std::numeric_limits<double>::max();
//        if (upa != this->end())
//          u = *upa;
//        
//        std::cout << " ===> [" << l << "," << u << "]\n";
        
//        int i{1};
//        for (auto a{this->begin() + 1}; a != this->end(); ++a) {
//            if(*a <= 1e-10)
//                std::cout << " " << i ;
//            ++i;
//        }
//        std::cout << std::endl;
    }

    constexpr static const double baseIncrement{1e-12};
    constexpr static const double maxActivity{1e24};
    constexpr static const double baseGap{baseIncrement};
    
protected:

    const double decay{};
    double increment{baseIncrement};
};

}

#endif //SCHEDCL_EVENTACTIVITYMAP_HPP

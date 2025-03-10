//
// Created by Emmanuel Hebrard on 5/3/24.
//

#ifndef SCHEDCL_LRBMAP_HPP
#define SCHEDCL_LRBMAP_HPP

#include <concepts>

#include "util/traits.hpp"


namespace tempo::heuristics::impl {
/**
 * @brief Class that can be used to record activity of variables
 */
class LearningRateMap : public std::vector<double> {
public:
    /**
     * CTor. Initializes activity of all literals with baseIncrement.
     */
    LearningRateMap(const double a=.4) : alpha(a) {}
    
    
    // called after the set of variables 'vars' has been assigned at fail step 'stamp'
    // record the fail step and reset the activity counter
    template<typename Iterable>
    void resetActivity(const Iterable& vars, const long unsigned int stamp) noexcept {
        if(assigned_at.size() < this->size()) {
            assigned_at.resize(this->size(), 0);
            participated.resize(this->size(), 0);
        }
        for(auto x : vars) {
            assigned_at[x] = stamp;
            participated[x] = 0;
        }
    }
    
    // called after the set of variables 'vars' has been involved in a conflict
    // increment the activity counter
    template<typename Iterable>
    void incrementActivity(const Iterable& vars) noexcept {
      for (const auto x : vars) {
          ++participated[x];
      }
    }
   
//     called after the set of variables 'vars' has been unassigned
//     update the learning rate (alpha-skewed average over the period wher it was assigned)
    template<typename Iterable>
    void update(const Iterable& vars, const long unsigned int stamp) {
        if(assigned_at.size() < this->size()) {
            assigned_at.resize(this->size(), 0);
            participated.resize(this->size(), 0);
        }
        for(auto x : vars) {
            update(x, stamp);
        }
    }
    
    void update(const var_t x, const long unsigned int stamp) {
        this->operator[](x) *= (1.0 - alpha);
        this->operator[](x) +=
        (static_cast<double>(participated[x]) /
         static_cast<double>(stamp - assigned_at[x] + 1)) *
        alpha;
    }
    
    constexpr static const double defaultValue{0};
    
public:
    
    const double alpha{.4};
        
    // [for each variable] the number of times it participated to a conflict
    std::vector<long unsigned int> participated;
    
    // [for each variable] the number of conflicts when it was assigned
    std::vector<long unsigned int> assigned_at;
    
    // [for each variable] its current learning rate (*this)

};

}

#endif //SCHEDCL_EVENTLearningRateMap_HPP

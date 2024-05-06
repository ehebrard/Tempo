//
// Created by tluchterha on 21/11/22.
//

#ifndef TEMPO_BASEHEURISTIC_HPP
#define TEMPO_BASEHEURISTIC_HPP

//#define DEBUG_HEURISTICS

#include <concepts>

#include "util/crtp.hpp"

namespace tempo {
    template<typename T>
    class Scheduler;
}

namespace tempo::heuristics {
    /**
     * @brief Requirement for a class that derives from BaseHeuristic
     * @details @copybrief
     * Requires a member function named getCost with valid signature
     * @tparam Impl
     */
    template<typename Impl>
    concept HeuristicImplementation = requires(Impl instance, var x) {
        { instance.getCost(x) } -> std::convertible_to<double>;
    };

    /**
     * @brief CRTP base class used to build heuristics that select a choice point by minimizing a cost function
     * @tparam Impl type of concrete Implementation
     */
    template<typename Impl>
    class BaseHeuristic : public crtp<Impl, BaseHeuristic>{
    public:
        /**
         * @brief Static polymorphic interface for the caller of the heuristic. Selects a choice point from a scheduler by
         * minimizing a cost function.
         * @details @copybrief
         * Also removes ground instances from the index sequence of the scheduler
         * @tparam T type of scheduler
         * @param scheduler scheduler for which to generate a choice point
         * @return selected choice point or DistanceConstraint::none if no further choices can be made
         */
        template<typename T> requires(HeuristicImplementation<Impl>)
        auto nextChoicePoint(Scheduler<T> & scheduler) {
            var best_var{NoVar};
            auto &indexSequence = scheduler.getBranch();
            double minCost = std::numeric_limits<double>::infinity();
            
            assert(not indexSequence.empty());
            
            for(auto x : indexSequence) {
                const auto cost = this->getImpl().getCost(x);
                
#ifdef DEBUG_HEURISTICS_CHOICE
                std::cout << scheduler.getEdge(POS(x)) << "<>" << scheduler.getEdge(NEG(x)) << ": " << cost;
#endif
                
                if (cost < minCost) {
                    minCost = cost;
                    best_var = x;
                    
#ifdef DEBUG_HEURISTICS_CHOICE
                    std::cout << "*";
#endif
                }
                
#ifdef DEBUG_HEURISTICS_CHOICE
                std::cout << std::endl;
#endif
                
            }
            
            assert(best_var != NoVar);
            return best_var;
        }

    private:
        BaseHeuristic() = default;
        friend Impl;
    };
}

#endif //TEMPO_BASEHEURISTIC_HPP

/**
* @author Tim Luchterhand
* @date 06.05.24
* @brief Base value selection heuristic
*/

#ifndef TEMPO_BASEVALUEHEURISTIC_HPP
#define TEMPO_BASEVALUEHEURISTIC_HPP

#include <concepts>
#include <stdexcept>

#include "Global.hpp"
#include "util/crtp.hpp"

namespace tempo::heuristics {

    /**
     * @brief Contains all information necessary for instantiating value selection heuristics
     */
    struct ValueHeuristicConfig {
        double epsilon;
    };

    /**
     * Interface for value selection heuristic implementations that derive from BaseValueHeuristic
     * @tparam H heuristic type
     * @tparam Sched information provider (usually the scheduler)
     */
    template<typename H, typename Sched>
    concept value_heuristic_implementation = requires(H heuristic, var x, const Sched &scheduler) {
        {heuristic.choose(x, scheduler)} -> std::same_as<lit>;
    };

    /**
     * Interface for value selection heuristics
     * @tparam H heuristic type
     * @tparam Sched information provider (usually the scheduler)
     */
    template<typename H, typename Sched>
    concept value_heuristic = requires(H heuristic, var x, const Sched &scheduler) {
        {heuristic.choosePolarity(x, scheduler)} -> std::same_as<lit>;
    };

    /**
     * @brief CRTP base class for value selection heuristics.
     * @details @copybrief
     * Implements epsilon-greedy mechanism and implements the value selection heuristic interface
     * @tparam Impl
     */
    template<typename Impl>
    class BaseValueHeuristic : public crtp<Impl, BaseValueHeuristic> {
        static constexpr auto EpsScale = 10000ul;
        unsigned long epsilon;
    public:
        /**
         * Ctor
         * @param epsilon epsilon value for epsilon greedy selection (0: deterministic, 1: purely random)
         */
        explicit constexpr BaseValueHeuristic(double epsilon): epsilon(static_cast<unsigned long>(epsilon * EpsScale)) {
            if (epsilon > 1 or epsilon < 0) {
                throw std::runtime_error("invalid epsilon value");
            }
        }

        /**
         * Value selection heuristic interface
         * @tparam Sched class that provides additional information for the actual implementation
         * @param cp choice point
         * @param scheduler scheduler instance
         * @return either POS(cp) or NEG(cp)
         * @return polarity selected by the actual heuristic with probability 1-epsilon, random polarity otherwise
         */
        template<typename Sched> requires(value_heuristic_implementation<Impl, Sched>)
        constexpr lit choosePolarity(var cp, const Sched &scheduler) {
            auto rval = tempo::random();
            if (rval % EpsScale < epsilon) {
                return rval % 2 == 0 ? POS(cp) : NEG(cp);
            }

            return this->getImpl().choose(cp, scheduler);
        }
    };
}

#endif //TEMPO_BASEVALUEHEURISTIC_HPP

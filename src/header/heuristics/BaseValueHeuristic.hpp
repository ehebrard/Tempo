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
#include "util/traits.hpp"

namespace tempo {
    template<typename T>
    class Scheduler;
}

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
     * @tparam T timing type
     */
    template<typename H, typename T>
    concept value_heuristic_implementation = requires(H heuristic, var x, const Scheduler<T> scheduler) {
        {heuristic.choose(x, scheduler)} -> std::same_as<lit>;
    };

    /**
     * Interface for value selection heuristics
     * @tparam H heuristic type
     * @tparam T timing type
     */
    template<typename H, typename T>
    concept value_heuristic = requires(H heuristic, var x, const Scheduler<T> scheduler) {
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
         * @tparam T timing type
         * @param cp choice point
         * @param scheduler scheduler instance
         * @return either POS(cp) or NEG(cp)
         * @return polarity selected by the actual heuristic with probability 1-epsilon, random polarity otherwise
         */
        template<concepts::scalar T> requires(value_heuristic_implementation<Impl, T>)
        constexpr lit choosePolarity(var cp, const Scheduler<T> &scheduler) {
            if (tempo::random() % EpsScale < epsilon) {
                return tempo::random() % 2 == 0 ? POS(cp) : NEG(cp);
            }

            return this->getImpl().choose(cp, scheduler);
        }
    };
}

#endif //TEMPO_BASEVALUEHEURISTIC_HPP

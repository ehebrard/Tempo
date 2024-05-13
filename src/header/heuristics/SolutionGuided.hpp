/**
* @author Tim Luchterhand
* @date 06.05.24
* @brief Solution guided value selection
*/

#ifndef TEMPO_SOLUTIONGUIDED_HPP
#define TEMPO_SOLUTIONGUIDED_HPP

#include <vector>

#include "BaseValueHeuristic.hpp"
#include "TightestValue.hpp"
#include "util/SubscribableEvent.hpp"

namespace tempo::heuristics {

    namespace detail {
        template<typename Sched>
        concept solution_provider = requires(const Sched &s, var x) {
            { s.hasSolution() } -> std::convertible_to<bool>;
            { s.getSolution()[x] } -> std::convertible_to<bool>;
        };
    }

    /**
     * @brief Solution guided value selection heuristic.
     * @details @copybrief
     * Performs the first decent using tempo::heuristics::TightestValue. After that follows the most recent solution
     */
    class SolutionGuided : public BaseValueHeuristic<SolutionGuided> {
    public:
        /**
         * Ctor.
         * @param epsilon see tempo::heuristics::BaseValueHeuristic
         */
        explicit SolutionGuided(double epsilon): BaseValueHeuristic<SolutionGuided>(epsilon) {}

        /**
         * heuristic interface
         * @tparam Sched class that provides previously encountered solutions
         * @param cp choice point
         * @param scheduler scheduler instance
         * @return either POS(cp) or NEG(cp)
         */
        template<detail::solution_provider Sched>
        [[nodiscard]] lit choose(var cp, const Sched &scheduler) const {
            if (not scheduler.hasSolution()) {
                return TightestValue::choose(cp, scheduler);
            }

            return scheduler.getSolution()[cp] ? POS(cp) : NEG(cp);
        }
    };

    MAKE_FACTORY(SolutionGuided, const ValueHeuristicConfig &config) {
            return SolutionGuided(config.epsilon);
    }};
}

#endif //TEMPO_SOLUTIONGUIDED_HPP

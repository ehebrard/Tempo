/**
* @author Tim Luchterhand
* @date 06.05.24
* @brief
*/

#ifndef TEMPO_SOLUTIONGUIDED_HPP
#define TEMPO_SOLUTIONGUIDED_HPP

#include <vector>

#include "BaseValueHeuristic.hpp"
#include "TightestValue.hpp"
#include "util/traits.hpp"
#include "util/SubscribableEvent.hpp"

namespace tempo::heuristics {
    class SolutionGuided : public BaseValueHeuristic<SolutionGuided> {
    public:
        explicit SolutionGuided(double epsilon): BaseValueHeuristic<SolutionGuided>(epsilon) {}

        template<concepts::scalar T>
        [[nodiscard]] lit choose(var cp, const Scheduler<T> &scheduler) const {
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

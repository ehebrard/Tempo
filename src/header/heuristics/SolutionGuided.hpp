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
        std::vector<bool> polarityCache;
        SubscriberHandle handle;
    public:
        template<concepts::scalar T>
        SolutionGuided(double epsilon, const Scheduler<T> &scheduler): BaseValueHeuristic<SolutionGuided>(epsilon),
                handle(scheduler.SolutionFound.subscribe_handled([this](const auto &sched) {
                    polarityCache = sched.getSolution();
                })) {}

        template<concepts::scalar T>
        [[nodiscard]] lit choose(var cp, const Scheduler<T> &scheduler) const {
            if (polarityCache.empty()) {
                return TightestValue::choose(cp, scheduler);
            }

            return polarityCache[cp] ? POS(cp) : NEG(cp);
        }
    };

    MAKE_TEMPLATE_FACTORY(SolutionGuided, concepts::scalar T, const ValueHeuristicConfig<T> &config) {
            return SolutionGuided(config.epsilon, config.scheduler);
    }};
}

#endif //TEMPO_SOLUTIONGUIDED_HPP

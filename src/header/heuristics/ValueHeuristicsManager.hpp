/**
* @author Tim Luchterhand
* @date 06.05.24
* @brief
*/

#ifndef TEMPO_VALUEHEURISTICSMANAGER_HPP
#define TEMPO_VALUEHEURISTICSMANAGER_HPP
#include <string>
#include <optional>

#include "util/factory_pattern.hpp"
#include "util/Options.hpp"
#include "util/traits.hpp"
#include "TightestValue.hpp"
#include "SolutionGuided.hpp"
#include "RandomValue.hpp"

namespace tempo {
    template<typename T>
    class Scheduler;
}


namespace tempo::heuristics {

    auto valHeuristicTypeToString(Options::PolarityHeuristic type) -> std::string;

    MAKE_FACTORY_PATTERN(ValueHeuristic, ValueHeuristicConfig, TightestValue, SolutionGuided, RandomValue)

    class ValueHeuristicsManager {
        std::optional<ValueHeuristic> impl;

    public:
        template<concepts::scalar T>
        explicit ValueHeuristicsManager(const Scheduler<T> &scheduler) {
            impl.emplace(ValueHeuristicFactory::getInstance().create(
                    valHeuristicTypeToString(scheduler.getOptions().polarity_heuristic),
                    ValueHeuristicConfig{.epsilon = 0.1}));
        }

        template<concepts::scalar T>
        lit choosePolarity(var cp, const Scheduler<T> &scheduler) {
            return std::visit([cp, &scheduler](auto &h) { return h.choosePolarity(cp, scheduler); }, *impl);
        }
    };

}

#endif //TEMPO_VALUEHEURISTICSMANAGER_HPP

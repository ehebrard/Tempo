/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief
 */

#ifndef TEMPO_VALUEHEURISTICSMANAGER_HPP
#define TEMPO_VALUEHEURISTICSMANAGER_HPP
#include <optional>
#include <string>

#include "RandomValue.hpp"
#include "SolutionGuided.hpp"
#include "TightestValue.hpp"
#include "util/Options.hpp"
#include "util/factory_pattern.hpp"
#include "util/traits.hpp"

namespace tempo {
template <typename T> class Scheduler;
}

namespace tempo::heuristics {

auto valHeuristicTypeToString(Options::PolarityHeuristic type) -> std::string;

MAKE_FACTORY_PATTERN(ValueHeuristic, ValueHeuristicConfig, TightestValue,
                     SolutionGuided, RandomValue)

class ValueHeuristicsManager {
  std::optional<ValueHeuristic> impl;
  template <typename T> struct FactoryChecker {
    HOLDS_FOR_ALL(ValueHeuristic, value_heuristic, T)
  };

public:
  template <concepts::scalar T>
  explicit ValueHeuristicsManager(const Scheduler<T> &scheduler) {
    impl.emplace(ValueHeuristicFactory::getInstance().create(
        valHeuristicTypeToString(scheduler.getOptions().polarity_heuristic),
        ValueHeuristicConfig{.epsilon =
                                 scheduler.getOptions().polarity_epsilon}));
  }

  template <typename Sched> lit choosePolarity(var cp, const Sched &scheduler) {
    static_assert(FactoryChecker<Sched>::template __ValueHeuristic_tester__<
                      ValueHeuristic>::value,
                  "At least one heuristic has an invalid signature");
    return std::visit(
        [cp, &scheduler](auto &h) { return h.choosePolarity(cp, scheduler); },
        *impl);
  }
};

} // namespace tempo::heuristics

#endif // TEMPO_VALUEHEURISTICSMANAGER_HPP

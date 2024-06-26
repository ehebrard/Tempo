/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief
 */

#ifndef TEMPO_VALUEHEURISTICSMANAGER_HPP
#define TEMPO_VALUEHEURISTICSMANAGER_HPP
#include <optional>
#include <string>

#include "RandomBinaryValue.hpp"
#include "SolutionGuided.hpp"
#include "TightestValue.hpp"
#include "util/Options.hpp"
#include "util/factory_pattern.hpp"
#include "util/traits.hpp"

namespace tempo {
template <typename T> class Solver;
}

namespace tempo::heuristics {

auto valHeuristicTypeToString(Options::PolarityHeuristic type) -> std::string;

MAKE_FACTORY_PATTERN(ValueHeuristic, ValueHeuristicConfig, TightestValue,
                     SolutionGuided, RandomBinaryValue)

class ValueHeuristicsManager {
  std::optional<ValueHeuristic> impl;
  template <typename T> struct FactoryChecker {
    HOLDS_FOR_ALL(ValueHeuristic, value_heuristic, T)
  };

public:
  template <concepts::scalar T>
  explicit ValueHeuristicsManager(const Solver<T> &scheduler) {
    impl.emplace(ValueHeuristicFactory::getInstance().create(
        valHeuristicTypeToString(scheduler.getOptions().polarity_heuristic),
        ValueHeuristicConfig{.epsilon =
                                 scheduler.getOptions().polarity_epsilon}));
  }

  template <typename Sched>
  auto valueDecision(const VariableSelection &selection, const Sched &scheduler) {
    static_assert(FactoryChecker<Sched>::template __ValueHeuristic_tester__<
                      ValueHeuristic>::value,
                  "At least one heuristic has an invalid signature");
    return std::visit([&](auto &h) { return h.valueDecision(selection, scheduler); }, *impl);
  }
};

} // namespace tempo::heuristics

#endif // TEMPO_VALUEHEURISTICSMANAGER_HPP

/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief Solution guided value selection
 */

#ifndef TEMPO_SOLUTIONGUIDED_HPP
#define TEMPO_SOLUTIONGUIDED_HPP

#include <vector>

#include "BaseBooleanHeuristic.hpp"
#include "TightestValue.hpp"
#include "util/SubscribableEvent.hpp"

namespace tempo::heuristics {

namespace detail {
template <typename Sched>
concept solution_provider = requires(const Sched &s, var_t x) {
  { s.hasSolution() } -> std::convertible_to<bool>;
  { s.getSolution()[x] } -> std::convertible_to<bool>;
};
}

/**
 * @brief Solution guided value selection heuristic.
 * @details @copybrief
 * Performs the first decent using tempo::heuristics::TightestValue. After that
 * follows the most recent solution
 */
class SolutionGuided : public BaseBooleanHeuristic<SolutionGuided> {
public:
  /**
   * Ctor.
   * @param epsilon see tempo::heuristics::BaseValueHeuristic
   */
  explicit SolutionGuided(double epsilon)
      : BaseBooleanHeuristic<SolutionGuided>(epsilon) {}

  /**
   * heuristic interface
   * @tparam T timing type
   * @tparam S class that provides previously encountered solutions
   */
  template <typename S>
  [[nodiscard]] auto choose(var_t x, const S &solver) const {
      throw std::runtime_error("needs implementation");
      return makeBooleanLiteral<int>(true, 0, 0);
      // something like if solver has not solution yet => return tightest
      // else lookup corresponding edge in solver.boolean.bestSolution and return a polarity
  }
};

MAKE_FACTORY(SolutionGuided, const ValueHeuristicConfig &config) {
  return SolutionGuided(config.epsilon);
}
};
}

#endif // TEMPO_SOLUTIONGUIDED_HPP

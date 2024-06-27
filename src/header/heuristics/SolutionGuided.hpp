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

template <typename T> class Solver;

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
template <class BaseHeuristic>
class SolutionGuided
    : public BaseBooleanHeuristic<SolutionGuided<BaseHeuristic>> {
public:
  /**
   * Ctor.
   * @param epsilon see tempo::heuristics::BaseValueHeuristic
   */
  template <class... arguments>
  explicit SolutionGuided(double epsilon, arguments &&...args)
      : BaseBooleanHeuristic<SolutionGuided>(epsilon),
        h(std::forward<arguments>(args)...) {}

  /**
   * heuristic interface
   * @tparam T timing type
   * @tparam S class that provides previously encountered solutions
   */
  template <typename T>
  [[nodiscard]] auto choose(var_t, const Solver<T> &) const {

    //    if (solver.boolean.hasSolution()) {
    //      return solver.boolean.value(x);
    //    } else {
    //      if (not solver.boolean.hasSemantic(x)) {
    //        return true;
    //      }
    //
    //      auto edgePos = solver.boolean.getEdge(true, x);
    //      auto edgeNeg = solver.boolean.getEdge(false, x);
    //      auto gapPos =
    //          solver.numeric.upper(edgePos.from) -
    //          solver.numeric.lower(edgePos.to);
    //      auto gapNeg =
    //          solver.numeric.upper(edgeNeg.from) -
    //          solver.numeric.lower(edgeNeg.to);
    //      return makeBooleanLiteral<T>(gapPos <= gapNeg, x);

    return makeBooleanLiteral<int>(false, 0, 0);
  }

  BaseHeuristic h;
};

// MAKE_FACTORY(SolutionGuided<TightestValue>, const ValueHeuristicConfig
// &config) {
//   return SolutionGuided<TightestValue>(config.epsilon);
// }
// };
}

#endif // TEMPO_SOLUTIONGUIDED_HPP

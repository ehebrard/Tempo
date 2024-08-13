/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief Tightest value selection
 */

#ifndef TEMPO_TIGHTESTVALUE_HPP
#define TEMPO_TIGHTESTVALUE_HPP

#include "BaseBooleanHeuristic.hpp"
#include "Constant.hpp"
#include "DistanceConstraint.hpp"
#include "Global.hpp"
#include "Literal.hpp"
#include "util/traits.hpp"


namespace tempo::heuristics {

namespace detail {
template <typename Solver>
concept edge_distance_provider = concepts::distance_provider<Solver> and requires(const Solver s, var_t x) {
  { s.boolean.getEdge(true, x) } -> concepts::same_template<DistanceConstraint>;
};
}

/**
 * @brief Tightest value selection heuristic.
 * @details @copybrief
 * Chooses the polarity that would leave the most slack in the timing network
 */
class TightestValue : public BaseBooleanHeuristic<TightestValue> {
public:
  /**
   * Ctor
   * @param epsilon see tempo::heuristics::BaseValueHeuristic
   */
  explicit TightestValue(double epsilon)
      : BaseBooleanHeuristic<TightestValue>(epsilon) {}

  /**
   * heuristic interface
   * @tparam Solver class that provides distances between events and a mapping
   * from literals to edges
   * @param cp choice point
   * @param solver scheduler instance
   * @return either POS(cp) or NEG(cp)
   */
  template <detail::edge_distance_provider Solver>
  requires(boolean_info_provider<Solver>) static auto choose(
      var_t x, const Solver &solver) {
    // @TODO no gap info available -> what should I return?
    if (not solver.boolean.hasSemantic(x)) {
      return solver.boolean.getLiteral(true, x);
    }

    auto edgePos = solver.boolean.getEdge(true, x);
    auto edgeNeg = solver.boolean.getEdge(false, x);
      
    auto gapPos = (edgePos.isNull() ? 1000000 :
        solver.numeric.upper(edgePos.from) - solver.numeric.lower(edgePos.to));
    auto gapNeg = (edgeNeg.isNull() ? 1000000 :
        solver.numeric.upper(edgeNeg.from) - solver.numeric.lower(edgeNeg.to));
      
    return solver.boolean.getLiteral(gapPos >= gapNeg, x);
  }
};
}

#endif // TEMPO_TIGHTESTVALUE_HPP

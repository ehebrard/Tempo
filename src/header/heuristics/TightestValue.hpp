/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief Tightest value selection
 */

#ifndef TEMPO_TIGHTESTVALUE_HPP
#define TEMPO_TIGHTESTVALUE_HPP

#include "BaseBooleanHeuristic.hpp"
#include "Global.hpp"
#include "util/traits.hpp"
#include "util/edge_distance.hpp"


namespace tempo::heuristics {

namespace detail {
}

/**
 * @brief Tightest value selection heuristic.
 * @details @copybrief
 * Chooses the polarity that would leave the most slack in the timing network
 */
template<concepts::scalar T>
class TightestValue : public BaseBooleanHeuristic<TightestValue<T>, T> {
public:
  /**
   * Ctor
   * @param epsilon see tempo::heuristics::BaseValueHeuristic
   */
  explicit TightestValue(double epsilon)
      : BaseBooleanHeuristic<TightestValue, T>(epsilon) {}
    
    explicit TightestValue(const Solver<T> &solver)
        : BaseBooleanHeuristic<TightestValue, T>(solver.getOptions().polarity_epsilon) {}

  /**
   * heuristic interface
   * @tparam Solver class that provides distances between events and a mapping
   * from literals to edges
   * @param cp choice point
   * @param solver scheduler instance
   * @return either POS(cp) or NEG(cp)
   */
  template <edge_distance_provider Solver>
  requires(boolean_info_provider<Solver>) static auto choose(
      var_t x, const Solver &solver) {
      
      
//      std::cout << solver.pretty(solver.boolean.getLiteral(true,x)) << " <> "
//      << solver.pretty(solver.boolean.getLiteral(false,x)) << "\n";
//      
    // @TODO no gap info available -> what should I return?
    if (not solver.boolean.hasSemantic(x)) {
        return solver.boolean.getLiteral((tempo::random() % 2), x);
//        return solver.boolean.getLiteral(false, x);
    }

    auto gapPos = boundEstimation(true, x, solver);
    auto gapNeg = boundEstimation(false, x, solver);
//    assert(gapPos.has_value() and gapNeg.has_value());
      // case where only one side is an edge
      if(not gapPos.has_value()) {
          return solver.boolean.getLiteral(false, x);
      } else if(not gapNeg.has_value()) {
          return solver.boolean.getLiteral(true, x);
      } else {
          return solver.boolean.getLiteral(gapPos >= gapNeg, x);
      }
  }
};

    struct TightestValueFactory : MakeValueHeuristicFactory<TightestValueFactory> {
        TightestValueFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(const Solver<T> &solver) const -> ValueHeuristic<T> {
            return std::make_unique<TightestValue<T>>(solver);
        }
    };
}

#endif // TEMPO_TIGHTESTVALUE_HPP

//
// Created by tluchterha on 21/11/22.
//

#ifndef TEMPO_TIGHTEST_HPP
#define TEMPO_TIGHTEST_HPP

#include "BaseHeuristic.hpp"
//#include "Scheduler.hpp"

namespace tempo::heuristics {

/**
 * @brief Heuristic selects the choice point with minimum Distance between two
 * nodes
 * @tparam Distance Type of Distance
 */
template <typename T> class Tightest : public BaseHeuristic<Tightest<T>> {
public:
  /**
   * CTor
   * @tparam DistFun type of Distance function
   * @param distFun Distance function
   */
//  explicit Tightest(Scheduler<T> &s) : distance(s) {}
    explicit Tightest() {}

  /**
   * Calculates the cost for a choice point which is the maximum of the Distance
   * between the nodes in both directions
   * @tparam T type of DistanceConstraint
   * @param choicePoint choice point to evaluate
   * @return maximum of the Distances in both directions between the nodes
   */
  T getCost(var x, const Scheduler<T>& sched) const {

    // to - from <= d // e_i - s_j <= 0

    auto prec_a{sched.getEdge(POS(x))};
    auto prec_b{sched.getEdge(NEG(x))};
    T gap_a{0};
    if (prec_a != Constant::NoEdge<T>)
      gap_a = sched.upper(prec_a.from) - sched.lower(prec_a.to) +
              prec_a.distance;
    T gap_b{0};
    if (prec_b != Constant::NoEdge<T>)
      gap_b = sched.upper(prec_b.from) - sched.lower(prec_b.to) +
              prec_b.distance;
    return std::max(gap_a, gap_b);
  }
    
    T getCost(var_t x, const Solver<T>& solver) const {
        double dom{1};
        
        if(solver.boolean.hasSemantic(x)) {
            auto p{solver.boolean.getLiteral(true,x)};
            auto n{solver.boolean.getLiteral(false,x)};
            
            auto prec_a{solver.boolean.getEdge(p)};
            auto prec_b{solver.boolean.getEdge(n)};
            
            auto gap_a = solver.numeric.upper(prec_a.from) - solver.numeric.lower(prec_a.to);
            auto gap_b = solver.numeric.upper(prec_b.from) - solver.numeric.lower(prec_b.to);
            
            dom = static_cast<double>(std::max(gap_a, gap_b));
        }
        
      return dom;
    }

//private:
//  const Scheduler<T> &distance;
};
}


#endif //TEMPO_TIGHTEST_HPP

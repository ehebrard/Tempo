/**
* @author Tim Luchterhand
* @date 16.09.24
* @brief Variable selection heuristics that orders variables by their respective gap in the temporal network divided by
 * variable activity
*/

#ifndef TEMPO_GAPOVERACTIVITY_HPP
#define TEMPO_GAPOVERACTIVITY_HPP

#include <algorithm>

#include "Global.hpp"
#include "Constant.hpp"
#include "RankingHeuristic.hpp"
#ifdef OLDVSIDS
#include "heuristics/impl/DecayingEventActivityMap.hpp"
#else
#include "heuristics/impl/ActivityMap.hpp"
#endif
#include "util/SubscribableEvent.hpp"
#include "util/edge_distance.hpp"
#include "Solver.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {

    /**
     * @brief Variable selection heuristics that orders variables by their respective gap in the temporal network
     * divided by variable activity.
     * @details @copybrief
     * This is the base class for heuristics like VSIDS and WeightedDegree
     */
    template<concepts::scalar T>
    class GapOverActivity : public RankingHeuristic<GapOverActivity<T>>, public BaseVariableHeuristic<T> {
    public:

        /**
         * CTor. Initializes activity map.
         * @param solver target solver
         * @param handle subscriber handle of event handler that updates activity map
         */
      explicit GapOverActivity(Solver<T> &solver, SubscriberHandle handle)
          :
#ifdef OLDVSIDS
            activity(solver, solver.getOptions().vsids_decay),
#else
            boolean_activity(solver.getBooleanActivity()),
            numeric_activity(solver.getNumericActivity()),
#endif
            handlerToken(std::move(handle)) {
#ifndef OLDVSIDS
          boolean_activity.resize(solver.boolean.size(),
                                  impl::ActivityMap::baseIncrement);
          numeric_activity.resize(solver.numeric.size(),
                                  impl::ActivityMap::baseIncrement);
#endif
        }

      // Copy and move are disabled. Otherwise, a call to the subscribed event
      // handler will cause undefined behavior
      GapOverActivity(const GapOverActivity &) = delete;
      GapOverActivity(GapOverActivity &&) = delete;
      GapOverActivity &operator=(const GapOverActivity &) = delete;
      GapOverActivity &operator=(GapOverActivity &&) = delete;
      ~GapOverActivity() override = default;

      [[nodiscard]] double getCost(const var_t x, const Solver<T> &solver) {

#ifdef OLDVSIDS
        activity.resize(solver);
#endif

        //@TODO: there shoud be a normalization thingy and Boolean variables
        //without semantic should get the highest value
        double dom{1};
        if (solver.boolean.hasSemantic(x)) {
          auto gapA =
              boundEstimation(true, x, solver).value_or(Constant::Infinity<T>);
          auto gapB =
              boundEstimation(false, x, solver).value_or(Constant::Infinity<T>);
          if (gapA == Constant::Infinity<T>) {
            dom = static_cast<double>(gapA / 2 + gapB / 2);
          } else if (gapB == Constant::Infinity<T>) {
            dom = static_cast<double>(gapA / 2 + gapB / 2);
          } else {
            dom = static_cast<double>(std::max(gapA, gapB));
          }
        }
#ifdef OLDVSIDS
        auto act{activity.get(x, solver)};
#else
        auto act{getActivity(x, solver)};
#endif
        return dom / act;
      }

#ifndef OLDVSIDS
      constexpr double getActivity(const var_t x,
                                   const Solver<T> &solver) const noexcept {

        double a{boolean_activity[x]};
        if (solver.boolean.hasSemantic(x)) {
          auto pe{solver.boolean.getEdge(true, x)};
          if (pe != Constant::NoEdge<T>)
            a += (numeric_activity[pe.from] + numeric_activity[pe.to]);

          auto ne{solver.boolean.getEdge(false, x)};
          if (ne != Constant::NoEdge<T>)
            a += (numeric_activity[ne.from] + numeric_activity[ne.to]);
        }
        return a;
      }
#endif

      /**
       * Ranking heuristic interface. Out of two variables selects the one for
       * which gap_v / activity_v is lower. Herby gap_v is the maximum distance
       * on the associated distance constraint arcs and activity_v is the
       * variables activity.
       * @tparam T timing type
       * @param x variable x
       * @param y variable y
       * @param solver solver that provides distance information
       * @return variable according to criteria described above
       */
      [[nodiscard]] var_t chooseBest(var_t x, var_t y,
                                     const Solver<T> &solver) {
        return getCost(x, solver) <= getCost(y, solver) ? x : y;
      }

      /**
       * Heuristic interface.
       * @param solver
       * @todo currently only selects boolean variables
       */
      auto nextVariable(const Solver<T> &solver) -> VariableSelection override {
        return {this->bestVariable(solver.getBranch(), solver),
                VariableType::Boolean};
      }

    protected:
#ifdef OLDVSIDS
      impl::DecayingEventActivityMap activity;
#else
      impl::ActivityMap &boolean_activity;
      impl::ActivityMap &numeric_activity;
#endif
      SubscriberHandle handlerToken;
    };
}



#endif //TEMPO_GAPOVERACTIVITY_HPP

//
// Created by tluchterha on 22/11/22.
//

#ifndef TEMPO_WEIGHTEDDEGREE_HPP
#define TEMPO_WEIGHTEDDEGREE_HPP

#include <vector>

#include "heuristic_interface.hpp"
#include "Literal.hpp"
#include "GapOverActivity.hpp"
#include "Solver.hpp"

namespace tempo::heuristics {
    /**
     * @brief Weighted degree heuristic. Prioritized choice points according to
     * their recent activity
     * @details @copybrief
     * Registers the activity of literals. The activity is a decaying counter which
     * counts the number of times the literal has been involved in a conflict.
     * @note The difference to VSIDS is that the activity of a literal depends on
     * how often it is involved in a conflict and not how often it is contained in a
     * learned clause
     */
    template<concepts::scalar T>
    class WeightedDegree : public GapOverActivity<T> {
    public:

        /**
         * CTor. Infers parameters from given scheduler. Subscribes to events of interest
         * @tparam T timing type
         * @param solver target solver
         */
      explicit WeightedDegree(Solver<T> &solver)
          : GapOverActivity<T>(
                solver,
                solver.ConflictEncountered.subscribe_handled(
                    [this
#ifdef OLDVSIDS
                     ,
                     &solver
#endif
                     ,
                     clause = std::vector<Literal<T>>{}](auto &expl) mutable {
                      clause.clear();
                      expl.explain(Contradiction<T>, clause);
#ifdef OLDVSIDS
                      this->activity.update(clause, solver);
#else
                      auto num_normalize{false};
                      auto bool_normalize{false};
                      for (auto l : clause) {
                        if (l.isNumeric()) {
                          num_normalize |=
                              this->numeric_activity.incrementActivity(
                                  l.variable());
                        } else {
                          bool_normalize |=
                              this->boolean_activity.incrementActivity(
                                  l.variable());
                        }
                      }
                      if (num_normalize) {
                        this->numeric_activity.normalize();
                      }
                      if (bool_normalize) {
                        this->boolean_activity.normalize();
                      }
#endif
                    })) {
      }
    };


    struct WeightedDegreeFactory : MakeVariableHeuristicFactory<WeightedDegreeFactory> {
        WeightedDegreeFactory();

        template<concepts::scalar T>
        auto build_impl(Solver<T>& solver) const -> VariableHeuristic<T> {
            return std::make_unique<WeightedDegree<T>>(solver);
        }
    };
}

#endif //SCHEDCL_WEIGHTEDDEGREE_HPP

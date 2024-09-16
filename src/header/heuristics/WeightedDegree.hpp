//
// Created by tluchterha on 22/11/22.
//

#ifndef TEMPO_WEIGHTEDDEGREE_HPP
#define TEMPO_WEIGHTEDDEGREE_HPP

#include <vector>

#include "Literal.hpp"
#include "GapOverActivity.hpp"

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
    class WeightedDegree : public GapOverActivity {
    public:

        template<concepts::scalar T>
        explicit WeightedDegree(Solver<T> &solver):
                  GapOverActivity(solver, solver.ConflictEncountered.subscribe_handled(
                          [this, &solver, clause=std::vector<Literal<T>>{}](auto &expl) mutable {
                              clause.clear();
                              expl.explain(Solver<T>::Contradiction, clause);
                              this->activity.update(clause, solver);
                          })) {}

        WeightedDegree(const WeightedDegree &) = delete;
        WeightedDegree(WeightedDegree &&) = delete;
        WeightedDegree &operator=(const WeightedDegree &) = delete;
        WeightedDegree &operator=(WeightedDegree &&) = delete;
        ~WeightedDegree() = default;
    };
}

#endif //SCHEDCL_WEIGHTEDDEGREE_HPP

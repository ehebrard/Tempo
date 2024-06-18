//
// Created by tluchterha on 22/11/22.
//

#ifndef TEMPO_WEIGHTEDDEGREE_HPP
#define TEMPO_WEIGHTEDDEGREE_HPP

#include "util/Options.hpp"
#include "util/SubscribableEvent.hpp"
#include "heuristics/RankingHeuristic.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
//#include <array>

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
    template<typename T>
    class WeightedDegree : public RankingHeuristic<WeightedDegree<T>> {
    public:

        explicit WeightedDegree(Solver<T> &solver)
                : activity(solver, solver.getOptions().vsids_decay),
                  handlerToken(solver.ConflictEncountered.subscribe_handled(
                          [this, &solver](auto &expl) {
                              this->nclause.clear();
                              expl.explain(Solver<T>::Contradiction, this->nclause);
                              this->activity.update(this->nclause, solver);
                          })) {}

        WeightedDegree(const WeightedDegree &) = delete;
        WeightedDegree(WeightedDegree &&) = delete;
        WeightedDegree &operator=(const WeightedDegree &) = delete;
        WeightedDegree &operator=(WeightedDegree &&) = delete;
        ~WeightedDegree() = default;


        [[nodiscard]] double getCost(const var x, const Solver<T> &solver) const {
            double dom{1};

            if (solver.boolean.hasSemantic(x)) {
                auto p{solver.boolean.getLiteral(true, x)};
                auto n{solver.boolean.getLiteral(false, x)};

                auto prec_a{solver.boolean.getEdge(p)};
                auto prec_b{solver.boolean.getEdge(n)};

                auto gap_a = solver.numeric.upper(prec_a.from) - solver.numeric.lower(prec_a.to);
                auto gap_b = solver.numeric.upper(prec_b.from) - solver.numeric.lower(prec_b.to);

                dom = static_cast<double>(std::max(gap_a, gap_b));
            }

            return dom / activity.get(x, solver);
        }

        /**
         * @param solver
         * @todo currently only selects boolean variables
         */
        auto nextVariable(const Solver<T> &solver) const -> std::pair<var_t, VariableType> {
            return {this->bestVariable(solver.getBranch(), solver), VariableType::Boolean};
        }

    private:
        impl::DecayingEventActivityMap<T> activity;
        SubscriberHandle handlerToken;
        std::vector<genlit> clause;
        std::vector<Literal<T>> nclause;
    };
}

#endif //SCHEDCL_WEIGHTEDDEGREE_HPP

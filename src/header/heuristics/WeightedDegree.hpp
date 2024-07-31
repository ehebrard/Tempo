//
// Created by tluchterha on 22/11/22.
//

#ifndef TEMPO_WEIGHTEDDEGREE_HPP
#define TEMPO_WEIGHTEDDEGREE_HPP

#include "Global.hpp"
#include "util/Options.hpp"
#include "util/SubscribableEvent.hpp"
#include "heuristics/RankingHeuristic.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"

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
    class WeightedDegree : public RankingHeuristic<WeightedDegree> {
    public:

        template<concepts::scalar T>
        explicit WeightedDegree(Solver<T> &solver)
                : activity(solver, solver.getOptions().vsids_decay),
                  handlerToken(solver.ConflictEncountered.subscribe_handled(
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


        template<concepts::scalar T>
        [[nodiscard]] double getCost(const var_t x, const Solver<T> &solver) const {
            double dom{1};

            if (solver.boolean.hasSemantic(x)) {
                auto p{solver.boolean.getLiteral(true, x)};
                auto n{solver.boolean.getLiteral(false, x)};

                auto prec_a{solver.boolean.getEdge(p)};
                auto prec_b{solver.boolean.getEdge(n)};

                auto gap_a = (prec_a == Constant::NoEdge<T> ? Constant::Infinity<T> : solver.numeric.upper(prec_a.from) - solver.numeric.lower(prec_a.to));
                auto gap_b = (prec_b == Constant::NoEdge<T> ? Constant::Infinity<T> : solver.numeric.upper(prec_b.from) - solver.numeric.lower(prec_b.to));
                
                if(gap_a == Constant::Infinity<T>) {
                  //                    assert(gap_b != Constant::Infinity<T>);
                  dom = static_cast<double>(gap_a / 2 + gap_b / 2);
                } else if(gap_b == Constant::Infinity<T>) {
                  //                    assert(gap_a != Constant::Infinity<T>);
                  dom = static_cast<double>(gap_a / 2 + gap_b / 2);
                } else {
                    dom = static_cast<double>(std::max(gap_a, gap_b));
                }
            }

            return dom / activity.get(x, solver);
        }

        /**
         * @param solver
         * @todo currently only selects boolean variables
         */
        template<concepts::scalar T>
        auto nextVariable(const Solver<T> &solver) const -> VariableSelection {
            return {this->bestVariable(solver.getBranch(), solver), VariableType::Boolean};
        }

private:
  impl::DecayingEventActivityMap activity;
  SubscriberHandle handlerToken;
};
}

#endif //SCHEDCL_WEIGHTEDDEGREE_HPP

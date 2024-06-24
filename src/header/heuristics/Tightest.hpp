//
// Created by tluchterha on 21/11/22.
//

#ifndef TEMPO_TIGHTEST_HPP
#define TEMPO_TIGHTEST_HPP

#include "RankingHeuristic.hpp"
#include "Global.hpp"

namespace tempo::heuristics {

    class Tightest : public RankingHeuristic<Tightest> {
    public:
        /**
         * @tparam T
         * @param x
         * @param solver
         */
        template<concepts::scalar T>
        [[nodiscard]] T getCost(var_t x, const Solver<T> &solver) const {
            T dom{1};
            if (solver.boolean.hasSemantic(x)) {
                auto p{solver.boolean.getLiteral(true, x)};
                auto n{solver.boolean.getLiteral(false, x)};

                auto prec_a{solver.boolean.getEdge(p)};
                auto prec_b{solver.boolean.getEdge(n)};

                auto gap_a = solver.numeric.upper(prec_a.from) - solver.numeric.lower(prec_a.to);
                auto gap_b = solver.numeric.upper(prec_b.from) - solver.numeric.lower(prec_b.to);

                dom = std::max(gap_a, gap_b);
            }

            return dom;
        }

        /**
         * @tparam T
         * @param solver
         * @todo currently only selects boolean variables
         */
        template<concepts::scalar T>
        [[nodiscard]] auto nextVariable(const Solver<T> &solver) const -> VariableSelection {
            return {this->bestVariable(solver.getBranch(), solver), VariableType::Boolean};
        }

    };
}


#endif //TEMPO_TIGHTEST_HPP

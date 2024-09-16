//
// Created by tluchterha on 21/11/22.
//

#ifndef TEMPO_TIGHTEST_HPP
#define TEMPO_TIGHTEST_HPP

#include "RankingHeuristic.hpp"
#include "Global.hpp"
#include "util/distance.hpp"

namespace tempo::heuristics {

    class Tightest : public RankingHeuristic<Tightest> {
    public:
        Tightest() = default;

        /**
         * Calculates the cost for a boolean variable which is the maximum of the distance
         * between the nodes in both directions of the associated distance constraint
         * @tparam T timing type
         * @param x variable identifier
         * @param solver solver for which to select a variable
         * @return maximum of the Distances in both directions between the nodes
         */
        template<edge_distance_provider S>
        [[nodiscard]] auto getCost(var_t x, const S &solver) const {
            using T = decltype(boundEstimation(true, 0, solver));
            T dom{1};
            if (solver.boolean.hasSemantic(x)) {
                auto gapA = boundEstimation(true, x, solver);
                auto gapB = boundEstimation(false, x, solver);
                dom = std::max(gapA, gapB);
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

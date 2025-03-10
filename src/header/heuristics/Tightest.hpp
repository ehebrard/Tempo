//
// Created by tluchterha on 21/11/22.
//

#ifndef TEMPO_TIGHTEST_HPP
#define TEMPO_TIGHTEST_HPP

#include "RankingHeuristic.hpp"
#include "Global.hpp"
#include "util/edge_distance.hpp"
#include "util/traits.hpp"
#include "Solver.hpp"
#include "heuristic_interface.hpp"

namespace tempo::heuristics {

    template<concepts::scalar T>
    class Tightest : public RankingHeuristic<Tightest<T>>, public BaseVariableHeuristic<T>{
    public:
        Tightest() = default;

        /**
         * Calculates the cost for a boolean variable which is the maximum of the distance
         * between the nodes in both directions of the associated distance constraint
         * @param x variable identifier
         * @param solver solver for which to select a variable
         * @return maximum of the Distances in both directions between the nodes
         */
        template<edge_distance_provider S>
        [[nodiscard]] auto getCost(var_t x, const S &solver) const {
            using Time = std::remove_cvref_t<decltype(*boundEstimation(true, 0, solver))>;
            assert(x != Constant::NoVar);
            Time dom{1};
            if (solver.boolean.hasSemantic(x)) {
                auto gapA = boundEstimation(true, x, solver);
                auto gapB = boundEstimation(false, x, solver);
                dom = std::max(gapA, gapB).value_or(InfiniteDistance<T>);
            }

            return dom;
        }

        template<edge_distance_provider S>
        [[nodiscard]] var_t chooseBest(var_t x, var_t y, const S &solver) const {
            return getCost(x, solver) <= getCost(y, solver) ? x : y;
        }

        /**
         * @tparam T
         * @param solver
         * @todo currently only selects boolean variables
         */
        [[nodiscard]] auto nextVariable(const Solver<T> &solver) -> VariableSelection override {
            return {this->bestVariable(solver.getBranch(), solver), VariableType::Boolean};
        }

    };

    struct TightestFactory : MakeVariableHeuristicFactory<TightestFactory> {
        TightestFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(const Solver<T>&) const -> VariableHeuristic<T> {
            return std::make_unique<Tightest<T>>();
        }
    };
}


#endif //TEMPO_TIGHTEST_HPP

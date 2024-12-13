//
// Created by tluchterha on 21/11/22.
//

#ifndef TEMPO_LRB_HPP
#define TEMPO_LRB_HPP

#include "RankingHeuristic.hpp"
#include "Global.hpp"
#include "util/edge_distance.hpp"

namespace tempo::heuristics {

    class LRB : public RankingHeuristic<LRB> {
    public:
        LRB() = default;

        /**
         * Calculates the cost for a boolean variable which is the maximum of the distance
         * between the nodes in both directions of the associated distance constraint
         * @tparam T timing type
         * @param x variable identifier
         * @param solver solver for which to select a variable
         * @return maximum of the Distances in both directions between the nodes
         */
        template<edge_distance_provider S>
        [[nodiscard]] auto getWindowCost(var_t x, const S &solver) const {
            using T = std::remove_cvref_t<decltype(*boundEstimation(true, 0, solver))>;
            assert(x != Constant::NoVar);
            T dom{1};
            if (solver.boolean.hasSemantic(x)) {
                auto gapA = boundEstimation(true, x, solver);
                auto gapB = boundEstimation(false, x, solver);
                dom = std::max(gapA, gapB).value_or(InfiniteDistance<T>);
            }

            return dom;
        }

        /**
         * Calculates the learning cost for a boolean variabletaking into
         * account the learning rates of the involved numeric variables between
         * the nodes in both directions of the associated distance constraint
         * @tparam T timing type
         * @param x variable identifier
         * @param solver solver for which to select a variable
         * @return maximum of the Distances in both directions between the nodes
         */
        template <edge_distance_provider S>
        [[nodiscard]] auto getActivityCost(var_t x, const S &solver) const {
          using T =
              std::remove_cvref_t<decltype(*boundEstimation(true, 0, solver))>;

          //            std::cout << "cost of " <<
          //            solver.pretty(solver.boolean.getLiteral(true,x)) << " <>
          //            " << solver.pretty(solver.boolean.getLiteral(false,x))
          //            << ":";

          auto c{100 * solver.boolean.getLearningRate(x)};

          //            std::cout << " " << c ;

          if (solver.boolean.hasSemantic(x)) {
            auto e1{solver.boolean.getEdge(true, x)};

            if (e1 != Constant::NoEdge<T>) {

              //                    std::cout << " + c(x" << e1.from << ")=" <<
              //                    solver.numeric.getLearningRate(e1.from);
              //                    std::cout << " + c(x" << e1.to << ")=" <<
              //                    solver.numeric.getLearningRate(e1.to);

              c += solver.numeric.getLearningRate(e1.from);
              c += solver.numeric.getLearningRate(e1.to);
            }
          }

          //            std::cout << " = " << c << std::endl;

          return c;
        }

//        template<edge_distance_provider S>
//        [[nodiscard]] auto getLearningRate(var_t x, const S &solver) const {
//            double solver.getLearningRate(x);
//        }

        template<edge_distance_provider S>
        [[nodiscard]] var_t chooseBest(var_t x, var_t y, const S &solver) const {
            if(solver.num_fails < threshold)
                return getWindowCost(x, solver) <= getWindowCost(y, solver) ? x : y;
            else
              return getActivityCost(x, solver) >= getActivityCost(y, solver)
                         ? x
                         : y;
            //                return solver.boolean.getLearningRate(x) >=
            //                solver.boolean.getLearningRate(y) ? x : y;
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

    private:
        unsigned threshold{100};
        
    };
}


#endif //TEMPO_LRB_HPP

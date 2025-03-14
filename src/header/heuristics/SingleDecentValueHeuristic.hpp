/**
* @author Tim Luchterhand
* @date 13.03.25
* @file SingleDecentValueHeuristic.hpp
* @brief Value selection heuristic wrapper that uses two heuristics, one until the first solution is found, then
* the other for the rest of the time
*/

#ifndef SINGLEDECENTVALUEHEURISTIC_HPP
#define SINGLEDECENTVALUEHEURISTIC_HPP

#include "util/traits.hpp"
#include "heuristic_interface.hpp"
#include "Solver.hpp"

namespace tempo::heuristics {
    /**
     * @brief Value selection heuristic wrapper that uses two heuristics, one until the first solution is found, then
     * the other for the rest of the time.
     * @tparam InitHeuristic Value heuristic for the first decent (until a solution is found)
     * @tparam BaseHeuristic Base heuristic that is used the rest of the time
     */
    template<typename InitHeuristic, typename BaseHeuristic>
    class SingleDecentValueHeuristic {
        InitHeuristic hInit;
        BaseHeuristic hBase;
    public:
        /**
         * Ctor
         * @tparam IH Initial heuristic type
         * @tparam BH Base heuristic type
         * @param initHeuristic heuristic to use until the first solution is found
         * @param baseHeuristic heuristic to use for the rest of the time
         */
        template<typename IH, typename BH>
        SingleDecentValueHeuristic(IH &&initHeuristic, BH &&baseHeuristic): hInit(std::forward<IH>(initHeuristic)),
                                                                            hBase(std::forward<BH>(baseHeuristic)) {}

        /**
         * Value heuristic interface
         * @tparam T timing type
         * @param x variable selection
         * @param solver target solver
         * @return decided literal
         */
        template<concepts::scalar T>
        requires(value_heuristic<InitHeuristic, Solver<T>> and value_heuristic<BaseHeuristic, Solver<T>>)
        auto valueDecision(const VariableSelection &x, const Solver<T> &solver) -> Literal<T> {
            if (solver.boolean.hasSolution()) {
                return hBase.valueDecision(x, solver);
            }

            return hInit.valueDecision(x, solver);
        }
    };

    /**
     * Helper method to create SingleDecentValueHeuristics
     * @tparam IH Initial heuristic type
     * @tparam BH Base heuristic type
     * @param initHeuristic heuristic to use until the first solution is found
     * @param baseHeuristic heuristic to use for the rest of the time
     * @return SingleDecentValueHeuristic
     */
    template<typename IH, typename BH>
    auto make_single_decent_heuristic(IH &&initHeuristic, BH &&baseHeuristic) {
        return SingleDecentValueHeuristic<IH, BH>(std::forward<IH>(initHeuristic), std::forward<BH>(baseHeuristic));
    }
}

#endif //SINGLEDECENTVALUEHEURISTIC_HPP

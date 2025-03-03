/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief Solution guided value selection
 */

#ifndef TEMPO_SOLUTIONGUIDED_HPP
#define TEMPO_SOLUTIONGUIDED_HPP

#include <vector>

#include "BaseBooleanHeuristic.hpp"
#include "heuristics/heuristic_interface.hpp"

namespace tempo::heuristics {

    namespace detail {
        template<typename Sched>
        concept solution_provider = requires(const Sched &s, var_t x) {
            { s.boolean.hasSolution() } -> std::convertible_to<bool>;
            { s.boolean.bestSolution() } -> std::same_as<const std::vector<bool> &>;
        };
    }

    /**
     * @brief Solution guided value selection heuristic.
     * @details @copybrief
     * Performs the first decent using some base heuristic. After that follows the most recent solution
     */
    template<class BaseHeuristic>
    class SolutionGuided : public BaseBooleanHeuristic<SolutionGuided<BaseHeuristic>> {
        BaseHeuristic h;
    public:
        /**
         * Ctor.
         * @tparam Args
         * @param epsilon see tempo::heuristics::BaseValueHeuristic
         * @param args arguments to base heuristic
         */
        template<class... Args>
        explicit SolutionGuided(double epsilon, Args &&...args): BaseBooleanHeuristic<SolutionGuided>(epsilon),
                                                                 h(std::forward<Args>(args)...) {}

        /**
         * heuristic interface
         * @tparam S class that provides previously encountered solutions
         */
        template<detail::solution_provider S>
        requires(heuristics::value_heuristic<BaseHeuristic, S> and boolean_info_provider<S>)
        [[nodiscard]] auto choose(var_t x, const S &solver) {
            const auto &b = solver.boolean;
            if (not b.hasSolution()) {
                return h.valueDecision({x, VariableType::Boolean}, solver);
            }

            auto posLit = b.getLiteral(true, x);
            auto negLit = b.getLiteral(false, x);
            const auto &sol = b.bestSolution();
            
            if(sol.size() <= posLit)
                return h.valueDecision({x, VariableType::Boolean}, solver);
            
            assert(sol[posLit] != sol[negLit]);
            return sol[posLit] ? posLit : negLit;
        }
    };
}

#endif // TEMPO_SOLUTIONGUIDED_HPP

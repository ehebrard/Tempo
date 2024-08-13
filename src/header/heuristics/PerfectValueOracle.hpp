/**
* @author Tim Luchterhand
* @date 13.08.24
* @brief Contains the perfect oracle value heuristic
*/

#ifndef TEMPO_PERFECTVALUEORACLE_HPP
#define TEMPO_PERFECTVALUEORACLE_HPP

#include <ranges>
#include <vector>

#include "BaseBooleanHeuristic.hpp"
#include "util/serialization.hpp"
#include "util/traits.hpp"

namespace tempo::heuristics {
    class PerfectValueHeuristic: public BaseBooleanHeuristic<PerfectValueHeuristic> {
        std::vector<bool> polarities;
    public:
        /**
         * @brief The perfect value heuristic. Follows the path to a given (possibly optimal) solution
         * @tparam T timing type
         * @param epsilon error chance
         * @param solution solution to follow
         */
        template<concepts::scalar T>
        PerfectValueHeuristic(double epsilon, const serialization::Solution <T> &solution)
                : BaseBooleanHeuristic<PerfectValueHeuristic>(epsilon), polarities(
                std::ranges::max(solution.decisions | std::views::elements<0>)) {
            for (auto [var, val] : solution.decisions) {
                polarities[var] = val;
            }
        }

        /**
         * heuristic interface
         * @tparam Solver solver class that provides boolean literal information for variables
         * @param x variable to fix
         * @param solver solver providing variable information
         * @return literal for x as in the specified solution
         */
        template<boolean_info_provider Solver>
        auto choose(var_t x, const Solver &solver) const {
            return solver.boolean.getLiteral(polarities[x], x);
        }
    };
}

#endif //TEMPO_PERFECTVALUEORACLE_HPP

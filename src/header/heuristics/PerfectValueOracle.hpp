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
#include "util/random.hpp"

namespace tempo::heuristics {
    /**
     * @brief The perfect value heuristic. Follows the path to a given (possibly optimal) solution.
     * @details @copybrief
     */
    template<concepts::scalar T>
    class PerfectValueHeuristic: public BaseBooleanHeuristic<PerfectValueHeuristic<T>, T> {
        std::vector<bool> polarities;
        static constexpr auto EpsScale = 10000ul;
    public:
        /**
         * Ctor
         * @param epsilon error chance
         * @param solution solution to follow
         */
        PerfectValueHeuristic(double epsilon, const serialization::Solution<T> &solution)
                : BaseBooleanHeuristic<PerfectValueHeuristic, T>(0), polarities(
                std::ranges::max(solution.decisions | std::views::elements<0>) + 1) {
            if (epsilon < 0 or epsilon > 1) {
                throw std::runtime_error("epsilon must be between 0 and one");
            }

            for (auto [var, val] : solution.decisions) {
                polarities.at(var) = (random() % EpsScale >= static_cast<unsigned long>(epsilon * EpsScale)) == val;
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

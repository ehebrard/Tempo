/**
* @author Tim Luchterhand
* @date 12.09.24
* @brief
*/

#ifndef TEMPO_RANDOMVARIABLESELECTION_HPP
#define TEMPO_RANDOMVARIABLESELECTION_HPP

#include "util/traits.hpp"
#include "heuristic_interface.hpp"
#include "util/SparseSet.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {
    /**
     * @brief Random variable selection strategy
     * @details @copybrief Randomly chooses a variable from the remenaing search variables
     * @note Right now, only binary variables are selected
     */
    struct RandomVariableSelection {

        /**
         * Heuristic interface
         * @tparam T timing type
         * @param solver solver for which to select the variable
         * @return randomly selected variable
         */
        template<concepts::scalar T>
        auto nextVariable(const Solver<T> &solver) const noexcept -> VariableSelection {
            const concepts::same_template<SparseSet> auto &variables = solver.getBranch();
            assert(not variables.empty());
            return {variables.any(), VariableType::Boolean};
        }
    };
}

#endif //TEMPO_RANDOMVARIABLESELECTION_HPP

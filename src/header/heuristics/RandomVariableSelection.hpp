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
#include "Solver.hpp"

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
    template<concepts::scalar T>
    struct RandomVariableSelection : BaseVariableHeuristic<T>{

        /**
         * Heuristic interface
         * @tparam T timing type
         * @param solver solver for which to select the variable
         * @return randomly selected variable
         */
        auto nextVariable(const Solver<T> &solver) noexcept -> VariableSelection override {
            const concepts::same_template<SparseSet> auto &variables = solver.getBranch();
            assert(not variables.empty());
            return {variables.any(), VariableType::Boolean};
        }
    };


    struct RandomVariableSelectionFactory : MakeVariableHeuristicFactory<RandomVariableSelectionFactory> {
        RandomVariableSelectionFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(const Solver<T>&) const -> VariableHeuristic<T> {
            return std::make_unique<RandomVariableSelection<T>>();
        }
    };
}

#endif //TEMPO_RANDOMVARIABLESELECTION_HPP

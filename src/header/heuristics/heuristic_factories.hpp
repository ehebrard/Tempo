/**
* @author Tim Luchterhand
* @date 09.07.24
* @brief
*/

#ifndef TEMPO_HEURISTIC_FACTORIES_HPP
#define TEMPO_HEURISTIC_FACTORIES_HPP

#include <iostream>

#include "heuristic_interface.hpp"
#include "util/Options.hpp"
#include "util/traits.hpp"

namespace tempo::heuristics {

    /**
     * Infers variable selection heuristic from solver options
     * @tparam T timing type
     * @param solver solver for which to create variable selection heuristic
     * @return variable selection heuristic inferred from solver options
     */
    template<concepts::scalar T>
    auto make_variable_heuristic(Solver<T> &solver) -> VariableHeuristic<T> {
        const auto &options = solver.getOptions();
        auto hType = options.choice_point_heuristics;
        if ((hType == Options::ChoicePointHeuristics::VSIDS or
             hType == Options::ChoicePointHeuristics::VSIDSHeap) and not options.learning) {
            if (options.verbosity >= Options::QUIET) {
                std::cout << "-- no learning, cannot use " << hType << std::endl;
            }

            hType = Options::ChoicePointHeuristics::WeightedDegree;
        }
        if (options.verbosity >= Options::QUIET) {
            std::cout << "-- using variable selection strategy '" << hType << "'" << std::endl;
        }

        return VariableHeuristicFactory::get().build(solver, hType);
    }

    /**
     * Infers value selection heuristic from solver options
     * @tparam T timing type
     * @param solver solver for which to create value selection heuristic
     * @return value selection heuristic inferred from solver options
     */
    template<concepts::scalar T>
    auto make_value_heuristic(Solver<T> &solver) -> ValueHeuristic<T> {
        const auto &options = solver.getOptions();
        if (options.verbosity >= Options::QUIET) {
            std::cout << "-- using value selection strategy '" << options.polarity_heuristic << "'" << std::endl;
        }

        return ValueHeuristicFactory::get().build(solver, options.polarity_heuristic);
    }

    /**
     * Infers branching heuristics to be used from solver options
     * @tparam T timing type
     * @param solver for which to create branching heuristic
     * @return branching heuristic inferred from solver options
     */
    template<concepts::scalar T>
    auto make_heuristic(Solver<T> &solver) {
        return make_compound_heuristic(make_variable_heuristic(solver), make_value_heuristic(solver));
    }
}

#endif //TEMPO_HEURISTIC_FACTORIES_HPP

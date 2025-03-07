/**
* @author Tim Luchterhand
* @date 09.07.24
* @brief
*/

#ifndef TEMPO_HEURISTIC_FACTORIES_HPP
#define TEMPO_HEURISTIC_FACTORIES_HPP

#include <nlohmann/json.hpp>
#include <string>
#include <iostream>

#include "LRB.hpp"
#include "RandomBinaryValue.hpp"
#include "SolutionGuided.hpp"
#include "TightestValue.hpp"
#include "heuristic_interface.hpp"

#include "util/factory_pattern.hpp"
#include "util/Options.hpp"
#include "util/traits.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {

    namespace detail {

    auto getVarHName(Options::ChoicePointHeuristics heuristic) -> std::string;
    auto getValHName(Options::PolarityHeuristic heuristic) -> std::string;

    template <typename... Heuristics>
    struct VariantHeuristicWrapper : std::variant<Heuristics...> {
      using std::variant<Heuristics...>::variant;
      using std::variant<Heuristics...>::emplace;

      template <concepts::scalar T>
        requires(variable_heuristic<Heuristics, T> && ...)
      DYNAMIC_DISPATCH(nextVariable, const Solver<T> &solver, &, solver, EMPTY)

      template <concepts::scalar T>
        requires(value_heuristic<Heuristics, Solver<T>> && ...)
      DYNAMIC_DISPATCH(valueDecision,
                       ESCAPE(VariableSelection x, const Solver<T> &solver), &,
                       ESCAPE(x, solver), EMPTY)

      template <concepts::scalar T>
        requires(heuristic<Heuristics, T> && ...)
      DYNAMIC_DISPATCH(branch, const Solver<T> &solver, &, solver, EMPTY)
    };
    }


    // --- Value Heuristics ---

    // Add further heuristic types here

    using TightestSolutionGuided = SolutionGuided<TightestValue>;
    using RandomSolutionGuided = SolutionGuided<RandomBinaryValue>;
    using ValueHeuristic = detail::VariantHeuristicWrapper<TightestValue, RandomBinaryValue, TightestSolutionGuided, RandomSolutionGuided>;

    // Define heuristic factory types here

    MAKE_TEMPLATE_FACTORY(TightestValue, concepts::scalar T, Solver<T> &solver) {
            return TightestValue(solver);
        }
    };

//MAKE_FACTORY(TightestValue, const Options &options) {
//        return TightestValue(options.polarity_epsilon);
//    }
//};

    MAKE_DEFAULT_TEMPLATE_FACTORY(RandomBinaryValue, concepts::scalar T, Solver<T> &)

MAKE_TEMPLATE_FACTORY(TightestSolutionGuided, concepts::scalar T, Solver<T> &solver) {
            return TightestSolutionGuided(solver);
        }
    };

MAKE_TEMPLATE_FACTORY(RandomSolutionGuided, concepts::scalar T, Solver<T> &solver) {
            return RandomSolutionGuided(solver);
        }
    };

    
    MAKE_FACTORY_PATTERN(ValueHeuristic, TightestValue, RandomBinaryValue, TightestSolutionGuided, RandomSolutionGuided)


    // --- Use these factory methods ---

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
                std::cout << "-- no learning, cannot use " << detail::getVarHName(hType) << std::endl;
            }

            hType = Options::ChoicePointHeuristics::WeightedDegree;
        }
        if (options.verbosity >= Options::QUIET) {
            const auto name = detail::getVarHName(hType);
            std::cout << "-- using variable selection strategy '" << name << "'" << std::endl;
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
    auto make_value_heuristic(Solver<T> &solver) {
        const auto &options = solver.getOptions();
        const auto name = detail::getValHName(options.polarity_heuristic);
        if (options.verbosity >= Options::QUIET) {
            std::cout << "-- using value selection strategy '" << name << "'" << std::endl;
        }

        return ValueHeuristicFactory::getInstance().create(name, solver);
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

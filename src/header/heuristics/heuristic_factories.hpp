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

#include "heuristic_interface.hpp"
#include "Tightest.hpp"
#include "VSIDS.hpp"
#include "RandomVariableSelection.hpp"
#include "WeightedDegree.hpp"
#include "RandomBinaryValue.hpp"
#include "TightestValue.hpp"
#include "SolutionGuided.hpp"

#include "util/factory_pattern.hpp"
#include "util/Options.hpp"
#include "util/traits.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {

    /// --- Variable Heuristics ---

    namespace detail {
        using VSIDS_M = MovableHeuristic<VSIDS>;
        using WeightedDegree_M = MovableHeuristic<WeightedDegree>;

        auto getVarHName(Options::ChoicePointHeuristics heuristic) -> std::string;
        auto getValHName(Options::PolarityHeuristic heuristic) -> std::string;

        template<typename ...Heuristics>
        struct VariantHeuristicWrapper : std::variant<Heuristics...> {
            using std::variant<Heuristics...>::variant;
            using std::variant<Heuristics...>::emplace;

            template<concepts::scalar T> requires(variable_heuristic<Heuristics, T> && ...)
            DYNAMIC_DISPATCH(nextVariable, const Solver<T> &solver, &, solver, EMPTY)

            template<concepts::scalar T> requires(value_heuristic<Heuristics, Solver<T>> && ...)
            DYNAMIC_DISPATCH(valueDecision, ESCAPE(VariableSelection x, const Solver<T> &solver), &, ESCAPE(x, solver), EMPTY)

            template<concepts::scalar T> requires(heuristic<Heuristics, T> && ...)
            DYNAMIC_DISPATCH(branch, const Solver<T> &solver, &, solver, EMPTY)
        };


    }

    using VariableHeuristic = detail::VariantHeuristicWrapper<Tightest, detail::VSIDS_M, detail::WeightedDegree_M,
                                                              RandomVariableSelection>;

    /// Define heuristic factory types here

    MAKE_TEMPLATE_FACTORY(Tightest, concepts::scalar T, const Solver<T> &) {
            return Tightest{};
        }
    };

    MAKE_TEMPLATE_P_FACTORY(VSIDS, VariableHeuristic, concepts::scalar T, Solver<T> &solver) {
            if(solver.getOptions().learning) {
                return detail::VSIDS_M(solver);
            } else {
                return detail::WeightedDegree_M(solver);
            }
        }
    };

    MAKE_TEMPLATE_FACTORY(WeightedDegree, concepts::scalar T, Solver<T> &solver) {
            return detail::WeightedDegree_M(solver);
        }
    };

    MAKE_TEMPLATE_FACTORY(RandomVariableSelection, concepts::scalar T, Solver<T> &) {
            return RandomVariableSelection{};
        }
    };

    MAKE_T_FACTORY_PATTERN(VariableHeuristic, template<concepts::scalar T>, Solver<T> &, Tightest, VSIDS,
                           WeightedDegree, RandomVariableSelection)


    /// --- Value Heuristics ---

    /// Add further heuristic types here

    using TightestSolutionGuided = SolutionGuided<TightestValue>;
    using ValueHeuristic = detail::VariantHeuristicWrapper<TightestValue, RandomBinaryValue, TightestSolutionGuided>;

    /// Define heuristic factory types here

    MAKE_FACTORY(TightestValue, const Options &options) {
            return TightestValue(options.polarity_epsilon);
        }
    };

    MAKE_DEFAULT_FACTORY(RandomBinaryValue, const Options&)

    MAKE_FACTORY(TightestSolutionGuided, const Options &options) {
            return TightestSolutionGuided(options.polarity_epsilon, 0);
        }
    };

    MAKE_FACTORY_PATTERN(ValueHeuristic, const Options&, TightestValue, RandomBinaryValue, TightestSolutionGuided)


    /// --- Use these factory methods ---

    /**
     * Infers variable selection heuristic from solver options
     * @tparam T timing type
     * @param solver solver for which to create variable selection heuristic
     * @return variable selection heuristic inferred from solver options
     */
    template<concepts::scalar T>
    auto make_variable_heuristic(Solver<T> &solver) {
        const auto &options = solver.getOptions();
        const auto name = detail::getVarHName(options.choice_point_heuristics);
        if (options.verbosity >= Options::QUIET) {
            std::cout << "-- using variable selection strategy '" << name << "'" << std::endl;
        }

        return VariableHeuristicFactory::getInstance().create(name, solver);
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

        return ValueHeuristicFactory::getInstance().create(name, options);
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

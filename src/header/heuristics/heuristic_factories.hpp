/**
* @author Tim Luchterhand
* @date 09.07.24
* @brief
*/

#ifndef TEMPO_HEURISTIC_FACTORIES_HPP
#define TEMPO_HEURISTIC_FACTORIES_HPP

#include <nlohmann/json.hpp>
#include <string>

#include "heuristic_interface.hpp"
#include "Tightest.hpp"
#include "VSIDS.hpp"
#include "WeightedDegree.hpp"
#include "RandomBinaryValue.hpp"
#include "TightestValue.hpp"

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
    }

    using VariableHeuristic = PolymorphicHeuristic<Tightest, detail::VSIDS_M, detail::WeightedDegree_M>;

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

    MAKE_T_FACTORY_PATTERN(VariableHeuristic, template<concepts::scalar T>, Solver<T>&, Tightest, VSIDS, WeightedDegree)

    /// --- Value Heuristics ---

    /// Add further heuristic types here

    using ValueHeuristic = PolymorphicHeuristic<TightestValue, RandomBinaryValue>;

    /// Define heuristic factory types here

    MAKE_FACTORY(TightestValue, const Options &options) {
            return TightestValue(options.polarity_epsilon);
        }
    };

    MAKE_DEFAULT_FACTORY(RandomBinaryValue, const Options&)

    MAKE_FACTORY_PATTERN(ValueHeuristic, const Options&, TightestValue, RandomBinaryValue)

    using BranchingHeuristic = CompoundHeuristic<VariableHeuristic, ValueHeuristic>;

    /// --- Use this factory method ---

    template<concepts::scalar T>
    auto make_heuristic(Solver<T> &solver) -> BranchingHeuristic {
        const auto &options = solver.getOptions();
        const auto &varHF = VariableHeuristicFactory::getInstance();
        const auto &valHF = ValueHeuristicFactory::getInstance();
        return {varHF.create(detail::getVarHName(options.choice_point_heuristics), solver),
                valHF.create(detail::getValHName(options.polarity_heuristic), options)};
    }
}

#endif //TEMPO_HEURISTIC_FACTORIES_HPP

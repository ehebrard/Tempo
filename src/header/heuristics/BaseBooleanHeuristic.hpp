/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief Base value selection heuristic
 */

#ifndef TEMPO_BASEBOOLEANHEURISTIC_HPP
#define TEMPO_BASEBOOLEANHEURISTIC_HPP

#include <concepts>
#include <stdexcept>

#include "Global.hpp"
#include "util/crtp.hpp"
#include "Literal.hpp"

namespace tempo::heuristics {

/**
 * @brief Contains all information necessary for instantiating value selection
 * heuristics
 */
struct ValueHeuristicConfig {
    double epsilon;
};

/**
 * Interface for value selection heuristic implementations that derive from
 * BaseBooleanHeuristic
 * @tparam H heuristic type
 * @tparam Solver information provider (usually the Solver)
 */
template <typename H, typename Solver>
concept binary_heuristic_implementation = requires(H heuristic, const Solver &solver, Literal<int> l) {
    { heuristic.choose(l, solver) } -> std::convertible_to<bool>;
};

/**
 * Interface for value selection heuristics
 * @tparam H heuristic type
 * @tparam Solver information provider (usually the scheduler)
 */
template <typename H, typename Solver>
concept value_heuristic = requires(H heuristic, VariableSelection x, const Solver &solver) {
    { heuristic.valueDecision(x, solver) } -> concepts::same_template<Literal>;
};

template<typename Solver>
concept boolean_info_provider = requires(const Solver s, var_t x) {
    { s.boolean.getLiteral(true, x) } -> concepts::same_template<Literal>;
};

/**
 * @brief CRTP base class for value selection heuristics.
 * @details @copybrief
 * Implements epsilon-greedy mechanism and implements the value selection
 * heuristic interface
 * @tparam Impl
 */
template <typename Impl>
class BaseBooleanHeuristic : public crtp<Impl, BaseBooleanHeuristic> {
    static constexpr auto EpsScale = 10000ul;
    unsigned long epsilon;

public:
    /**
     * Ctor
     * @param epsilon epsilon value for epsilon greedy selection (0:
     * deterministic, 1: purely random)
     */
    explicit constexpr BaseBooleanHeuristic(double epsilon)
            : epsilon(static_cast<unsigned long>(epsilon * EpsScale)) {
        if (epsilon > 1 or epsilon < 0) {
            throw std::runtime_error("invalid epsilon value");
        }
    }

    /**
     * Value selection heuristic interface
     * @tparam Solver
     * @param selection
     * @param solver solver instance
     * @return
     */
    template<boolean_info_provider Solver>
    requires(binary_heuristic_implementation <Impl, Solver>)
    constexpr auto valueDecision(const VariableSelection &selection, const Solver &solver) {
        assert(selection.second == VariableType::Boolean);
        auto lit = solver.boolean.getLiteral(true, selection.first);
        auto rval = tempo::random();
        if (rval % EpsScale < epsilon) {
            return rval % 2 == 0 ? lit : ~lit;
        }

        return this->getImpl().choose(lit, solver) ? lit : ~lit;
    }
};
} // namespace tempo::heuristics

#endif // TEMPO_BASEBOOLEANHEURISTIC_HPP

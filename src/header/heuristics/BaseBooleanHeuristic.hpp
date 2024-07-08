/************************************************
 * Tempo BaseBooleanHeuristic.hpp
 *
 * Copyright 2024 Tim Luchterhand
 *
 * Tempo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 * Tempo is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tempo.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************/

#ifndef TEMPO_BASEBOOLEANHEURISTIC_HPP
#define TEMPO_BASEBOOLEANHEURISTIC_HPP

#include <concepts>
#include <stdexcept>
#include <cassert>

#include "Global.hpp"
#include "Literal.hpp"
#include "util/crtp.hpp"
#include "util/traits.hpp"

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
concept binary_heuristic_implementation = requires(H heuristic,
                                                   const Solver &solver,
                                                   var_t x) {
  { heuristic.choose(x, solver) } -> concepts::same_template<Literal>;
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
    { s.boolean.hasSemantic(x) } -> std::convertible_to<bool>;
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
      auto rval = tempo::random();
      if ((rval % EpsScale) < epsilon) {
        auto lit = solver.boolean.getLiteral((rval % 2) == 0, selection.first);
        assert(lit.isBoolean());
        return lit;
      }

      return this->getImpl().choose(selection.first, solver);
    }
};
} // namespace tempo::heuristics

#endif // TEMPO_BASEBOOLEANHEURISTIC_HPP

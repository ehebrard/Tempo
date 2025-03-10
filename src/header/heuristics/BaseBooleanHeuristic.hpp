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
#include "Solver.hpp"
#include "Literal.hpp"
#include "util/crtp.hpp"
#include "util/traits.hpp"
#include "util/random.hpp"
#include "heuristics/heuristic_interface.hpp"

namespace tempo::heuristics {

template<typename S>
concept boolean_info_provider = requires(const S s, var_t x) {
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
template <typename Impl, concepts::scalar T>
class BaseBooleanHeuristic : public crtp<Impl, BaseBooleanHeuristic, T>, public BaseValueHeuristic<T> {
    static constexpr auto EpsScale = 10000ul;
    unsigned long epsilon;
protected:

    template<boolean_info_provider S>
    auto valueDecisionImpl(const VariableSelection &selection, const S &solver) {
        assert(selection.second == VariableType::Boolean);
        auto rval = tempo::random();
        if ((rval % EpsScale) < epsilon) {
            auto lit = solver.boolean.getLiteral((rval % 2) == 0, selection.first);
            assert(lit.isBoolean());
            return lit;
        }

        return this->getImpl().choose(selection.first, solver);
    }
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
    auto valueDecision(const VariableSelection &selection, const Solver<T> &solver) -> Literal<T> override {
        return valueDecisionImpl(selection, solver);
    }
};
} // namespace tempo::heuristics

#endif // TEMPO_BASEBOOLEANHEURISTIC_HPP

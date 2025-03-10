/************************************************
 * Tempo RandomBinaryValue.hpp
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

#ifndef TEMPO_RANDOMBINARYVALUE_HPP
#define TEMPO_RANDOMBINARYVALUE_HPP

#include <cassert>

#include "heuristic_interface.hpp"
#include "util/random.hpp"
#include "Solver.hpp"

namespace tempo::heuristics {

/**
 * @brief Random value selection heuristic.
 * @details @copybrief
 * Chooses choice point polarity always randomly
 */
    template <concepts::scalar T>
    class RandomBinaryValue : public BaseValueHeuristic<T> {
    public:
        auto valueDecision(const VariableSelection &selection,
                           const Solver<T> &solver) noexcept -> Literal<T> override {
            assert(selection.second == VariableType::Boolean);
            return solver.boolean.getLiteral(tempo::random() % 2 == 0, selection.first);
        }
    };

    struct RandomBinaryValueFactory : MakeValueHeuristicFactory<RandomBinaryValueFactory> {
        RandomBinaryValueFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(const Solver<T> &) const -> ValueHeuristic<T> {
            return std::make_unique<RandomBinaryValue<T>>();
        }
    };

} // namespace tempo::heuristics

#endif // TEMPO_RANDOMBINARYVALUE_HPP

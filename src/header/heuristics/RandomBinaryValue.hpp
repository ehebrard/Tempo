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

#include "BaseBooleanHeuristic.hpp"
#include "Global.hpp"
#include "util/random.hpp"

namespace tempo::heuristics {

/**
 * @brief Random value selection heuristic.
 * @details @copybrief
 * Chooses choice point polarity always randomly
 */
    class RandomBinaryValue {
    public:
        
        explicit RandomBinaryValue() {}
        explicit RandomBinaryValue(double) {}
        
        template <concepts::scalar T>
        explicit RandomBinaryValue(Solver<T> &) {}
        
        /**
         * heuristic interface
         * @tparam Sched class that provides additional information for the actual
         * implementation
         * @param cp choice point
         * @param scheduler scheduler instance
         * @return either POS(cp) or NEG(cp)
         */
        template<boolean_info_provider Solver>
        static auto valueDecision(const VariableSelection &selection, const Solver &solver) noexcept {
            assert(selection.second == VariableType::Boolean);
            return solver.boolean.getLiteral(tempo::random() % 2 == 0, selection.first);
        }
    };

} // namespace tempo::heuristics

#endif // TEMPO_RANDOMBINARYVALUE_HPP

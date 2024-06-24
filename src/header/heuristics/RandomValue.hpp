/************************************************
 * Tempo RandomValue.hpp
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

#ifndef TEMPO_RANDOMVALUE_HPP
#define TEMPO_RANDOMVALUE_HPP
#include <cassert>

#include "BaseValueHeuristic.hpp"
#include "Global.hpp"
#include "util/factory_pattern.hpp"

namespace tempo::heuristics {

/**
 * @brief Random value selection heuristic.
 * @details @copybrief
 * Chooses choice point polarity always randomly
 */
class RandomValue {
public:
  /**
   * heuristic interface
   * @tparam Sched class that provides additional information for the actual
   * implementation
   * @param cp choice point
   * @param scheduler scheduler instance
   * @return either POS(cp) or NEG(cp)
   */
  template <typename Sched>
  static lit choosePolarity(var cp, const Sched &) noexcept {
    return tempo::random() % 2 == 0 ? POS(cp) : NEG(cp);
  }
};

MAKE_DEFAULT_FACTORY(RandomValue, const ValueHeuristicConfig &)

} // namespace tempo::heuristics

#endif // TEMPO_RANDOMVALUE_HPP

/************************************************
 * Tempo BaseValueHeuristic.hpp
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

#ifndef TEMPO_BASEVALUEHEURISTIC_HPP
#define TEMPO_BASEVALUEHEURISTIC_HPP

#include <concepts>
#include <stdexcept>

#include "Global.hpp"
#include "util/crtp.hpp"

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
 * BaseValueHeuristic
 * @tparam H heuristic type
 * @tparam Sched information provider (usually the scheduler)
 */
template <typename H, typename Sched>
concept value_heuristic_implementation = requires(H heuristic, var x,
                                                  const Sched &scheduler) {
  { heuristic.choose(x, scheduler) } -> std::same_as<lit>;
};

/**
 * Interface for value selection heuristics
 * @tparam H heuristic type
 * @tparam Sched information provider (usually the scheduler)
 */
template <typename H, typename Sched>
concept value_heuristic = requires(H heuristic, var x, const Sched &scheduler) {
  { heuristic.choosePolarity(x, scheduler) } -> std::same_as<lit>;
};

/**
 * @brief CRTP base class for value selection heuristics.
 * @details @copybrief
 * Implements epsilon-greedy mechanism and implements the value selection
 * heuristic interface
 * @tparam Impl
 */
template <typename Impl>
class BaseValueHeuristic : public crtp<Impl, BaseValueHeuristic> {
  static constexpr auto EpsScale = 10000ul;
  unsigned long epsilon;

public:
  /**
   * Ctor
   * @param epsilon epsilon value for epsilon greedy selection (0:
   * deterministic, 1: purely random)
   */
  explicit constexpr BaseValueHeuristic(double epsilon)
      : epsilon(static_cast<unsigned long>(epsilon * EpsScale)) {
    if (epsilon > 1 or epsilon < 0) {
      throw std::runtime_error("invalid epsilon value");
    }
  }

  /**
   * Value selection heuristic interface
   * @tparam Sched class that provides additional information for the actual
   * implementation
   * @param cp choice point
   * @param scheduler scheduler instance
   * @return either POS(cp) or NEG(cp)
   * @return polarity selected by the actual heuristic with probability
   * 1-epsilon, random polarity otherwise
   */
  template <typename Sched>
  requires(value_heuristic_implementation<Impl, Sched>) constexpr lit
      choosePolarity(var cp, const Sched &scheduler) {
    auto rval = tempo::random();
    if (rval % EpsScale < epsilon) {
      return rval % 2 == 0 ? POS(cp) : NEG(cp);
    }

    return this->getImpl().choose(cp, scheduler);
  }
};
} // namespace tempo::heuristics

#endif // TEMPO_BASEVALUEHEURISTIC_HPP

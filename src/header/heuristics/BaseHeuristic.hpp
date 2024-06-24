/************************************************
 * Tempo BaseHeuristic.hpp
 *
 * Copyright 2024 Tim Luchterhand and Emmanuel Hebrard
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

#ifndef TEMPO_BASEHEURISTIC_HPP
#define TEMPO_BASEHEURISTIC_HPP

//#define DEBUG_HEURISTICS

#include <concepts>

#include "util/crtp.hpp"

namespace tempo {

template<typename T>
class Solver;
}

namespace tempo::heuristics {
/**
 * @brief Requirement for a class that derives from BaseHeuristic
 * @details @copybrief
 * Requires a member function named getCost with valid signature
 * @tparam Impl
 */
template <typename Impl, typename T>
concept HeuristicImplementation = requires(Impl instance, var_t x, const Solver<T>& s) {
  { instance.getCost(x, s) } -> std::convertible_to<double>;
};

/**
 * @brief CRTP base class used to build heuristics that select a choice point by
 * minimizing a cost function
 * @tparam Impl type of concrete Implementation
 */
template <typename Impl>
class BaseHeuristic : public crtp<Impl, BaseHeuristic> {
public:
  /**
   * @brief Static polymorphic interface for the caller of the heuristic.
   * Selects a choice point from a scheduler by minimizing a cost function.
   * @details @copybrief
   * Also removes ground instances from the index sequence of the scheduler
   * @tparam T type of scheduler
   * @param solver solver for which to generate a choice point
   * @return selected choice point or DistanceConstraint::none if no further
   * choices can be made
   */
  template <typename T>
  requires(HeuristicImplementation<Impl, T>) auto nextChoicePoint(
      const Solver<T> &solver) {

    var_t best_var{Constant::NoVar};
    auto &indexSequence = solver.getBranch();
    double minCost = std::numeric_limits<double>::infinity();

    assert(not indexSequence.empty());

    for (auto x : indexSequence) {
      const auto cost = this->getImpl().getCost(x, solver);
      if (cost < minCost) {
        minCost = cost;
        best_var = x;
      }
    }

    assert(best_var != Constant::NoVar);
    return best_var;
  }

private:
  BaseHeuristic() = default;
  friend Impl;
};
}

#endif //TEMPO_BASEHEURISTIC_HPP
